# 3003

**cohen_cb_analysis.sh** analyzes sample BAM files for the presence of 160 Y-chromosome markers from a reference table, determining genotypes based on the most frequent alternative allele at each position. The output includes a VCF file with detailed genotypes for each sample and a TSV summary classifying samples into Cohen Branches (CB-01 — CB-09 from [doi:10.64898/2025.12.08.692646](https://doi.org/10.64898/2025.12.08.692646)) based on the proportion of matching markers.

---

## Input Data

| # | File | Description |
|---|------|-------------|
| 1 | `*.bam` + `*.bai` | Indexed sample BAM files |
| 2 | `markers_master_table.tsv` | Marker table (SNP_Name, CB_Branch, POS_hg19, REF, ALT) |
| 3 | `hg19.no_alt.fa` + `.fai` | Indexed hg19 reference sequence |

---

## Output Data

| File | Description |
|------|-------------|
| `cohen_cb_genotypes.vcf.gz` + `.tbi` | Indexed VCF documenting the allele observed at each marker position for each sample |
| `cohen_cb_summary.tsv` | Summary table counting detected markers for each haplogroup |

---

## Usage

```bash
./cohen_cb_analysis.sh bam CB output ref19
```

| Argument | Description |
|----------|-------------|
| `bam` | Directory containing sample BAM files |
| `CB` | Directory containing the marker table (`markers_master_table.tsv`) |
| `output` | Output directory for results |
| `ref19` | Directory containing the hg19 reference genome |

---

## Details: `cohen_cb_genotypes.vcf` Structure

### Header and Data Example

```vcf
##fileformat=VCFv4.2
##contig=<ID=chrY,length=59373566>
##INFO=<ID=CB,Number=1,Type=String,Description="Cohen Branch">
##INFO=<ID=EXPECTED_ALT,Number=1,Type=String,Description="Expected alternative allele">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=OBS,Number=1,Type=String,Description="Observed: REF, ., or A/C/G/T">
##FORMAT=<ID=CBM,Number=1,Type=String,Description="CB match: YES/NO/NA">
##FORMAT=<ID=CBS,Number=1,Type=String,Description="CB status">
#CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO                      FORMAT              sample1         sample2
chrY    2709209 ZS222   T    C    .     PASS    CB=CB-01;EXPECTED_ALT=C   GT:DP:OBS:CBM:CBS   0:5:REF:NA:REF  1:8:C:YES:MATCH
```

### VCF Fields

| Field | Description |
|-------|-------------|
| `ID` | Marker name |
| `ALT` | Expected marker allele |
| `INFO:CB` | Marker's Cohen Branch assignment |
| `INFO:EXPECTED_ALT` | Duplicate of the expected marker allele |

### FORMAT Fields (per sample)

| Field | Description | Values |
|-------|-------------|--------|
| `GT` | Genotype | `0` = reference, `1` = non-reference |
| `DP` | Read depth | Number of reads |
| `OBS` | Observed allele in sample | `REF` = reference allele, `.` = no coverage, `A/C/G/T` = alternative allele |
| `CBM` | Match to marker allele | `YES` / `NO` / `NA` |
| `CBS` | Verdict regarding marker allele | `MATCH` = marker allele detected, `WRONG_ALLELE` = alternative non-marker allele, `REF_ONLY` = reference allele only, `NO_COVERAGE` = no coverage |

---

## Details: `cohen_cb_summary.tsv` Structure

### Example

```tsv
Sample	CB-01_markers	CB-02_markers	CB-03_markers	CB-04_markers	CB-05_markers	CB-06_markers	CB-07_markers	CB-08_markers	CB-09_markers	CB_Assignment
sample1	18/22	0/4	2/74	0/6	0/3	0/3	0/5	0/3	0/27	CB-01
sample2	5/22	8/11	1/74	0/6	0/3	0/3	0/5	0/3	0/27	CB-02
sample3	0/22	0/4	0/74	0/6	0/3	0/3	0/5	0/3	0/27	NO_CB_MATCH
sample4	0/22	0/4	0/74	0/6	0/3	0/3	0/5	0/3	0/27	NO_CB_COVERAGE
```

### Columns

| Column | Description |
|--------|-------------|
| `Sample` | Sample name |
| `CB-XX_markers` | Branch markers in format `detected/read` |
| `CB_Assignment` | Analysis verdict |

### Marker Cell Format

```
detected markers / read positions
```

| Component | Description |
|-----------|-------------|
| `detected` | Number of markers where the expected allele was found (MATCH) |
| `read` | Number of markers with coverage ≥1 read |

### CB_Assignment Values

| Value | Description |
|-------|-------------|
| `NO_CB_COVERAGE` | Coverage threshold not met for assignment |
| `NO_CB_MATCH` | Match ratio threshold not met — no Cohen Branch assignment |
| `CB-01` — `CB-09` | One or more branches meet classification thresholds |

---

## Classification Thresholds

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MIN_DEPTH` | 1 | Minimum depth for a position to be considered |
| `MIN_COVERAGE_RATIO` | 0.30 | ≥30% of branch markers must be covered |
| `MIN_MATCH_RATIO` | 0.50 | ≥50% of covered markers must match the expected allele |

---

## Script Versions

| Version | File | Description |
|---------|------|-------------|
| **Basic** | `cohen_cb_analysis.sh` | Sequential BAM processing; regex patterns compiled in each loop iteration |
| **Optimized** | `cohen_cb_analysis2.sh` | **Recommended** — parallel mpileup + precompiled regex patterns |

### Key Differences in the Optimized Version:

1. **Parallel mpileup** — each BAM file is processed by a separate process simultaneously

**Before (sequential):**
```bash
samtools mpileup -f "$REF" -l "$BED_FILE" -q 0 -Q 0 -d 10000 $bam_args > "$RAW_PILEUP"
```

**After (parallel):**
```bash
PILEUP_DIR="$TMP_DIR/pileups"
mkdir -p "$PILEUP_DIR"
while IFS=$'\t' read -r sample bam; do
    samtools mpileup -f "$REF" -l "$BED_FILE" -q 0 -Q 0 -d 10000 "$bam" > "$PILEUP_DIR/${sample}.pileup" &
done < "$SAMPLE_NAMES"
wait
cat "$PILEUP_DIR"/*.pileup | sort -k1,1 -k2,2n > "$RAW_PILEUP"
```

2. **Precompiled regular expressions** — patterns are compiled once before the loop

**Before (compilation inside loop):**
```python
clean = re.sub(r'\^.', '', re.sub(r'$', '', re.sub(r'[0-9]+[+-][0-9]+[ACGTNacgtn]*', '', bases)))
alt_bases = [b.upper() for b in re.findall(r'[ACGTNacgtn]', clean)]
```

**After (compilation once):**
```python
# Before loop:
PILEUP_CLEAN_RE = re.compile(r'(?:\^.)|(?:\$)|(?:[0-9]+[+-][0-9]+[ACGTNacgtn]*)')
ALT_BASE_RE = re.compile(r'[ACGTNacgtn]')

# Inside loop:
clean = PILEUP_CLEAN_RE.sub('', bases)
alt_bases = [b.upper() for b in ALT_BASE_RE.findall(clean)]
```

### Runtime Benchmarks (3 samples, 160 markers):

| Version | real | user | sys | Speedup |
|---------|------------------|------------------|-----|---------|
| **Basic** (sequential) | 3m 54s | 3m 11s | 8.5s | 1× |
| **Optimized** (parallel) | **1m 20s** | 3m 16s | 7.2s | **~3×** |


## Academic Context

This repository contains supplementary materials for the Master's thesis:

**"Identification of Cohen Haplogroups and their Application to Modern Jewish Subpopulations"**

The repository includes:
- Analysis pipeline for Cohen branch identification (CB-01 to CB-09)
- Interactive visualizations of population genetic structure
- Marker tables and reference data

### Interactive Visualizations

The HTML files (`mtdna_combined.html` and `ydna_combined.html`) provide interactive exploration of:

**mtDNA Visualization:**
- Mitochondrial haplogroup distribution across five Jewish subpopulations (Ashkenazi, Mountain, Georgian, Bukharan, Kurdistani)
- Composition of prior cohort (n=132) and internal Genotek cohort (n=1,960)
- Maternal lineage diversity patterns supporting findings in Section 3.2

[Open Interactive Plot](https://mortew.github.io/3003/mtdna_combined.html)

**Y-DNA Visualization:**
- Y-chromosome haplogroup distribution across five Jewish subpopulations (Ashkenazi, Mountain, Georgian, Bukharan, Kurdistani)
- Composition of prior cohort (n=132) and internal Genotek cohort (n=1,960)
- Population-specific Y-chromosome signatures discussed in Section 3.3

[Open Interactive Plot](https://mortew.github.io/3003/ydna_combined.html)

These visualizations complement haplogroup distribution charts (Figs. 6-17) presented in the thesis.

> **Note:**  
> Clicking on the `.html` files directly in the repository will open them as code/text.  
> For interactive viewing with Plotly graphs, use the links above.
