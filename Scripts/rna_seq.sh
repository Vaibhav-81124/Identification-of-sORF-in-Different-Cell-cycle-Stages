#!/usr/bin/env bash
set -euo pipefail

# ─── CONFIGURATION ────────────────────────────────────────────────────────────
BASE="$HOME/multiomics_project"
DATA="${BASE}/data"
RESULTS="${BASE}/results"
REF="${DATA}/reference"

RNA_FQ="${DATA}/rna_seq/SRR3306578.fastq.gz"
GENOME_FA="${REF}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GENOME_GTF="${REF}/Homo_sapiens.GRCh38.115.gtf"
SAMPLE_NAME="HeLa_M_RNA_rep2"
THREADS=12
ADAPTER="AGATCGGAAGAGCACACGTCTGAAC"

STAR_INDEX="${REF}/star_index"

TRIMMED_FQ="${RESULTS}/02_trimmed/${SAMPLE_NAME}_trimmed.fastq.gz"
SORTED_BAM="${RESULTS}/03_aligned_rnaseq/${SAMPLE_NAME}_sorted.bam"

# ─── CREATE DIRECTORIES ───────────────────────────────────────────────────────
mkdir -p "${RESULTS}/01_qc/rnaseq_raw"
mkdir -p "${RESULTS}/01_qc/rnaseq_trimmed"
mkdir -p "${RESULTS}/02_trimmed"
mkdir -p "${RESULTS}/03_aligned_rnaseq"

echo "===================================================="
echo " RNA-seq Pipeline Starting: $(date)"
echo "===================================================="

# ─── STEP 1: FastQC RAW ───────────────────────────────────────────────────────
RAW_QC="${RESULTS}/01_qc/rnaseq_raw/$(basename ${RNA_FQ%.gz})_fastqc.zip"

if [ ! -f "$RAW_QC" ]; then
    echo "── STEP 1: FastQC (raw) ──"
    fastqc --outdir "${RESULTS}/01_qc/rnaseq_raw" \
           --threads ${THREADS} "${RNA_FQ}"
else
    echo "✓ STEP 1 skipped (FastQC raw exists)"
fi


# ─── STEP 2: CUTADAPT ─────────────────────────────────────────────────────────
if [ ! -f "$TRIMMED_FQ" ]; then
    echo "── STEP 2: Adapter trimming ──"
    cutadapt \
        -a "${ADAPTER}" \
        --quality-cutoff 20 \
        --minimum-length 30 \
        --cores ${THREADS} \
        -o "$TRIMMED_FQ" \
        "${RNA_FQ}" \
        2>&1 | tee "${RESULTS}/02_trimmed/${SAMPLE_NAME}_cutadapt.log"
else
    echo "✓ STEP 2 skipped (trimmed FASTQ exists)"
fi


# ─── STEP 2B: FastQC TRIMMED ───────────────────────────────────────────────────
TRIM_QC="${RESULTS}/01_qc/rnaseq_trimmed/${SAMPLE_NAME}_trimmed_fastqc.zip"

if [ ! -f "$TRIM_QC" ]; then
    echo "── STEP 2B: FastQC (trimmed) ──"
    fastqc --outdir "${RESULTS}/01_qc/rnaseq_trimmed" \
           --threads ${THREADS} "$TRIMMED_FQ"
else
    echo "✓ STEP 2B skipped (FastQC trimmed exists)"
fi


# ─── STEP 3: STAR INDEX ───────────────────────────────────────────────────────
if [ ! -d "${STAR_INDEX}" ] || [ -z "$(ls -A ${STAR_INDEX})" ]; then
    echo "── STEP 3: Building STAR index ──"
    mkdir -p "${STAR_INDEX}"

    STAR \
        --runMode genomeGenerate \
        --runThreadN ${THREADS} \
        --genomeDir "${STAR_INDEX}" \
        --genomeFastaFiles "${GENOME_FA}" \
        --sjdbGTFfile "${GENOME_GTF}" \
        --sjdbOverhang 50 \
        --genomeSAindexNbases 10 

    echo "✓ STAR index built"
else
    echo "✓ STEP 3 skipped (STAR index exists)"
fi


# ─── STEP 4: STAR ALIGNMENT ───────────────────────────────────────────────────
if [ ! -f "$SORTED_BAM" ]; then
    echo "── STEP 4: STAR alignment ──"

    STAR \
        --runThreadN ${THREADS} \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "$TRIMMED_FQ" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM MD \
        --outFileNamePrefix "${RESULTS}/03_aligned_rnaseq/${SAMPLE_NAME}_" \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000

    mv "${RESULTS}/03_aligned_rnaseq/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam" \
       "$SORTED_BAM"

    echo "✓ Alignment complete"

else
    echo "✓ STEP 4 skipped (BAM exists)"
fi


echo "===================================================="
echo " PIPELINE COMPLETE"
echo "===================================================="
