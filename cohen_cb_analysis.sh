#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Использование: $0 <bam_dir> <markers_dir> <output_dir> [ref_dir]"
    exit 1
fi

BAM_DIR="$1"
MARKERS_DIR="$2"
OUT_DIR="$3"
REF_DIR="${4:-ref19}"
REF="${REF_DIR}/hg19.no_alt.fa"

MASTER_MARKERS_FILE="${MARKERS_DIR}/markers_master_table.tsv"

MIN_DEPTH=1
MIN_COVERAGE_RATIO=0.30
MIN_MATCH_RATIO=0.50

TMP_DIR=$(mktemp -d)
trap "rm -rf $TMP_DIR" EXIT
mkdir -p "$OUT_DIR"

echo "Cohen Branch Analysis Pipeline"
echo "BAM: $BAM_DIR"
echo "Маркеры: $MASTER_MARKERS_FILE"
echo "Вывод: $OUT_DIR"
echo "Референс: $REF"

echo "[0/4] Проверка входных данных..."
[[ ! -f "$REF" ]] && { echo "Ошибка: Референс не найден: $REF"; exit 1; }
[[ ! -f "${REF}.fai" ]] && samtools faidx "$REF"
[[ ! -f "$MASTER_MARKERS_FILE" ]] && { echo "Ошибка: Таблица не найдена: $MASTER_MARKERS_FILE"; exit 1; }

BAM_COUNT=$(find "$BAM_DIR" -maxdepth 1 -name "*.bam" -type f | wc -l)
[[ $BAM_COUNT -eq 0 ]] && { echo "Ошибка: BAM-файлы не найдены"; exit 1; }

echo "  BAM: $BAM_COUNT, Референс: OK, Маркеры: OK"

echo "[1/4] Создание BED-файла..."
BED_FILE="$TMP_DIR/markers.bed"
awk -F'\t' 'NR>1 && $3!="" && $3!="NA" && $3~/^[0-9]+$/ {print "chrY\t"$3-1"\t"$3"\t"$1"\t0\t+"}' "$MASTER_MARKERS_FILE" | sort -u -k1,1 -k2,2n > "$BED_FILE"
echo "  Позиций: $(wc -l < "$BED_FILE")"

echo "[2/4] Поиск BAM-файлов..."
BAM_LIST="$TMP_DIR/bam_list.txt"
SAMPLE_NAMES="$TMP_DIR/sample_names.txt"
find "$BAM_DIR" -maxdepth 1 -name "*.bam" -type f | sort > "$BAM_LIST"

> "$SAMPLE_NAMES"
while read -r bam; do
    sample_name=$(samtools view -H "$bam" 2>/dev/null | grep "^@RG" | head -1 | sed 's/.*SM:\([^ \t]*\).*/\1/' || echo "")
    if [[ -z "$sample_name" ]]; then
        sample_name=$(basename "$bam" .bam)
    fi
    echo -e "${sample_name}\t${bam}" >> "$SAMPLE_NAMES"
done < "$BAM_LIST"

echo "  Образцов: $(wc -l < "$BAM_LIST")"

echo "[3/4] Извлечение аллелей..."
RAW_PILEUP="$TMP_DIR/raw.pileup"
bam_args=$(cat "$BAM_LIST" | tr '\n' ' ')

samtools mpileup -f "$REF" -l "$BED_FILE" -q 0 -Q 0 -d 10000 $bam_args > "$RAW_PILEUP"
echo "  Позиций в pileup: $(wc -l < "$RAW_PILEUP")"

echo "[4/4] Генерация отчётов..."
FINAL_VCF="${OUT_DIR}/cohen_cb_genotypes.vcf.gz"
TSV_OUT="${OUT_DIR}/cohen_cb_summary.tsv"

python3 - "$MASTER_MARKERS_FILE" "$RAW_PILEUP" "$FINAL_VCF" "$TSV_OUT" "$SAMPLE_NAMES" "$MIN_DEPTH" "$MIN_COVERAGE_RATIO" "$MIN_MATCH_RATIO" "$REF" << 'PYTHON'
import pysam
from collections import defaultdict, Counter
import sys
import re

MASTER_MARKERS_FILE, RAW_PILEUP, FINAL_VCF, TSV_OUT, SAMPLE_NAMES = sys.argv[1:6]
MIN_DEPTH, MIN_COVERAGE_RATIO, MIN_MATCH_RATIO = int(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8])
REF_FA = sys.argv[9]

markers, markers_by_pos, cb_markers = {}, defaultdict(list), defaultdict(list)
with open(MASTER_MARKERS_FILE, 'r') as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            snp, cb, pos, ref, alt = parts[0], parts[1], int(parts[2]), parts[3], parts[4]
            markers[snp] = {'cb': cb, 'pos': pos, 'ref': ref.upper(), 'alt': alt.upper() if alt and alt != '.' else None}
            markers_by_pos[pos].append(snp)
            cb_markers[cb].append(snp)

print(f"  Загружено маркеров: {len(markers)}")

sample_names = [line.strip().split('\t')[0] for line in open(SAMPLE_NAMES) if line.strip()]
if not sample_names:
    sample_names = ['unknown']

pileup_by_pos = {}
with open(RAW_PILEUP, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 5:
            continue
        pos = int(parts[1])
        depth = int(parts[3])
        bases = parts[4]
        clean = re.sub(r'\^.', '', re.sub(r'$', '', re.sub(r'[0-9]+[+-][0-9]+[ACGTNacgtn]*', '', bases)))
        alt_bases = [b.upper() for b in re.findall(r'[ACGTNacgtn]', clean)]
        observed_alt = Counter(alt_bases).most_common(1)[0][0] if alt_bases else None
        pileup_by_pos[pos] = {'depth': depth, 'observed_alt': observed_alt, 'gt': 1 if observed_alt else 0}

print(f"  Покрыто позиций: {len(pileup_by_pos)}")

chrY_length = 59373566
try:
    for line in open(REF_FA + '.fai'):
        parts = line.split('\t')
        if parts[0] in ['chrY', 'Y']:
            chrY_length = int(parts[1])
            break
except:
    pass

header = pysam.VariantHeader()
header.add_line('##fileformat=VCFv4.2')
header.add_line(f'##contig=<ID=chrY,length={chrY_length}>')
header.add_line('##INFO=<ID=CB,Number=1,Type=String,Description="Cohen Branch">')
header.add_line('##INFO=<ID=EXPECTED_ALT,Number=1,Type=String,Description="Expected alternative allele">')
header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
header.add_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">')
header.add_line('##FORMAT=<ID=OBS,Number=1,Type=String,Description="Observed: REF, ., or A/C/G/T">')
header.add_line('##FORMAT=<ID=CBM,Number=1,Type=String,Description="CB match: YES/NO/NA">')
header.add_line('##FORMAT=<ID=CBS,Number=1,Type=String,Description="CB status">')

for s in sample_names:
    header.add_sample(s)

with pysam.VariantFile(FINAL_VCF, 'w', header=header) as vcf:
    no_cov_total = 0
    for snp, data in sorted(markers.items(), key=lambda x: x[1]['pos']):
        pos, cb, exp_ref, exp_alt = data['pos'], data['cb'], data['ref'], data['alt']
        rec = vcf.new_record()
        rec.chrom = 'chrY'
        rec.pos = pos
        rec.id = snp
        rec.ref = exp_ref
        rec.alts = (exp_alt,) if exp_alt and exp_alt != '.' else ()
        rec.filter.add('PASS')
        rec.info['CB'] = cb
        rec.info['EXPECTED_ALT'] = exp_alt if exp_alt else '.'
        for sample in sample_names:
            if pos in pileup_by_pos:
                p = pileup_by_pos[pos]
                obs_alt, gt, dp = p['observed_alt'], p['gt'], p['depth']
                if gt == 1 and obs_alt:
                    match = 'YES' if obs_alt == exp_alt else 'NO'
                    status = 'MATCH' if match == 'YES' else 'WRONG_ALLELE'
                    obs_value = obs_alt
                else:
                    match, status, obs_value = 'NA', 'REF_ONLY', 'REF'
            else:
                obs_alt, gt, dp = None, 0, 0
                match, status, obs_value = 'NA', 'NO_COVERAGE', '.'
                no_cov_total += 1
            rec.samples[sample]['GT'] = (gt,)
            rec.samples[sample]['DP'] = dp
            rec.samples[sample]['OBS'] = obs_value
            rec.samples[sample]['CBM'] = match
            rec.samples[sample]['CBS'] = status
        vcf.write(rec)
    print(f"  NO_COVERAGE: {no_cov_total}")

pysam.tabix_index(FINAL_VCF, preset='vcf', force=True)

results = {s: {cb: {'found': 0, 'covered': 0} for cb in cb_markers} for s in sample_names}
for snp, data in markers.items():
    cb, pos, exp_alt = data['cb'], data['pos'], data['alt']
    if pos in pileup_by_pos:
        p = pileup_by_pos[pos]
        if p['depth'] >= MIN_DEPTH:
            for sample in sample_names:
                results[sample][cb]['covered'] += 1
                if p['gt'] == 1 and p['observed_alt'] == exp_alt:
                    results[sample][cb]['found'] += 1

def classify(res, cb_m):
    if sum(d['covered'] for d in res.values()) == 0:
        return "NO_CB_COVERAGE"
    cands = []
    for cb, d in res.items():
        total = len(cb_m[cb])
        if total == 0:
            continue
        cov_r = d['covered'] / total
        match_r = d['found'] / d['covered'] if d['covered'] > 0 else 0
        if cov_r >= MIN_COVERAGE_RATIO and match_r >= MIN_MATCH_RATIO and d['found'] > 0:
            cands.append((cb, match_r, cov_r))
    if not cands:
        return "NO_CB_MATCH"
    cands.sort(key=lambda x: (-x[1], -x[2]))
    return cands[0][0] if len(cands) == 1 or cands[0][1] - cands[1][1] > 0.2 else ",".join(c[0] for c in cands[:3])

with open(TSV_OUT, 'w') as out:
    cbs = sorted(cb_markers.keys())
    out.write('\t'.join(['Sample'] + [f'{cb}_markers' for cb in cbs] + ['CB_Assignment']) + '\n')
    for s in sorted(results.keys()):
        row = [s] + [f"{results[s][cb]['found']}/{results[s][cb]['covered']}" for cb in cbs]
        row.append(classify(results[s], cb_markers))
        out.write('\t'.join(row) + '\n')

print(f"  Готово: VCF и TSV")
PYTHON

echo ""
echo "Анализ завершён"
echo "Выходные файлы:"
echo "  $FINAL_VCF"
echo "  $TSV_OUT"