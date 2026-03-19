#!/usr/bin/env bash
set -euo pipefail

# ───────────────── CONFIGURATION ─────────────────
BASE="$HOME/multiomics_project"
DATA="${BASE}/data"
RESULTS="${BASE}/results"
REF="${DATA}/reference"

SORTMERNA_DB="${DATA}/rRNA_db"

RIBO_FQ="${DATA}/ribo_seq/SRR3306586.fastq.gz"
GENOME_GTF="${REF}/Homo_sapiens.GRCh38.115.gtf"
STAR_INDEX="${REF}/star_index"

SAMPLE="HeLa_M_RIBO_rep2"   # ← change per replicate

THREADS=$(nproc)
ADAPTER="TGGAATTCTCGGGTGCCAAGG"
RPF_MIN=26
RPF_MAX=34

mkdir -p ${RESULTS}/{01_qc,02_trimmed,03_rRNA_removed,04_aligned,05_psite}

echo "========== RIBO-SEQ PIPELINE START =========="
echo "Using ${THREADS} threads"

# ───────────────── SANITY CHECKS ─────────────────
[[ -f ${RIBO_FQ} ]] || { echo "Missing input FASTQ"; exit 1; }
[[ -f ${STAR_INDEX}/SA ]] || { echo "STAR index not found"; exit 1; }

echo "All required files found."

# ───────────────── STEP 1: TRIM + SIZE SELECT ─────────────────
TRIMMED="${RESULTS}/02_trimmed/${SAMPLE}_trimmed.fastq.gz"

cutadapt \
  -a ${ADAPTER} \
  --quality-cutoff 20 \
  --minimum-length ${RPF_MIN} \
  --maximum-length ${RPF_MAX} \
  --cores ${THREADS} \
  -o ${TRIMMED} \
  ${RIBO_FQ}

# ───────────────── STEP 2: rRNA REMOVAL (Rebuild Index Each Time) ─────────────────
NORRNA_PREFIX="${RESULTS}/03_rRNA_removed/${SAMPLE}_norRNA"
NORRNA="${NORRNA_PREFIX}.fq.gz"

WORKDIR="${RESULTS}/03_rRNA_removed/sortmerna_${SAMPLE}_workdir"
rm -rf ${WORKDIR}
mkdir -p ${WORKDIR}

sortmerna \
  --ref ${SORTMERNA_DB}/silva-euk-28s-id98.fasta \
  --ref ${SORTMERNA_DB}/silva-euk-18s-id95.fasta \
  --ref ${SORTMERNA_DB}/rfam-5s-database-id98.fasta \
  --ref ${SORTMERNA_DB}/rfam-5.8s-database-id98.fasta \
  --reads ${TRIMMED} \
  --workdir ${WORKDIR} \
  --other ${NORRNA_PREFIX} \
  --fastx \
  --threads ${THREADS}

# ───────────────── STEP 3: STAR ALIGNMENT ─────────────────
STAR_PREFIX="${RESULTS}/04_aligned/${SAMPLE}_"

STAR \
  --runThreadN ${THREADS} \
  --genomeDir ${STAR_INDEX} \
  --readFilesIn ${NORRNA} \
  --readFilesCommand zcat \
  --sjdbGTFfile ${GENOME_GTF} \
  --outFileNamePrefix ${STAR_PREFIX} \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS NM MD \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --outFilterMatchNmin 20 \
  --limitBAMsortRAM 20000000000

BAM="${STAR_PREFIX}Aligned.sortedByCoord.out.bam"
samtools index ${BAM}

# ───────────────── STEP 4: FILTER CANONICAL RPF LENGTHS ─────────────────
RPF_BAM="${RESULTS}/04_aligned/${SAMPLE}_RPFs.bam"

samtools view -h ${BAM} | \
awk -v min=${RPF_MIN} -v max=${RPF_MAX} '
BEGIN{OFS="\t"}
/^@/ {print; next}
{
  if (length($10) >= min && length($10) <= max)
    print
}' | samtools sort -@ ${THREADS} -o ${RPF_BAM}

samtools index ${RPF_BAM}

# ───────────────── STEP 5: P-SITE ASSIGNMENT ─────────────────
PSITE_BED="${RESULTS}/05_psite/${SAMPLE}_psites.bed"

samtools view ${RPF_BAM} | \
awk 'BEGIN{OFS="\t"}
{
  chrom=$3
  start=$4-1
  strand = and($2,16) ? "-" : "+"
  readlen=length($10)

  if (readlen==28 || readlen==29) offset=12
  else if (readlen==30 || readlen==31) offset=13
  else offset=12

  if (strand=="+")
    psite=start+offset
  else
    psite=start+readlen-offset-1

  print chrom, psite, psite+1, ".", ".", strand
}' > ${PSITE_BED}

echo "========== PIPELINE COMPLETE =========="
