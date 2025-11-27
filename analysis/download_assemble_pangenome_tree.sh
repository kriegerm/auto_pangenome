#!/usr/bin/env bash
set -euo pipefail


# RUN: bash scripts/project_setup.sh <PROJECT_ID>

# Grab the project ID from the first argument
PROJECT_ID="${1:-}"

if [[ -z "$PROJECT_ID" ]]; then
  echo "Usage: $0 <PROJECT_ID>" >&2
  exit 1
fi


echo "Running pipeline for project: ${PROJECT_ID}"

# Initialize conda for non-interactive shells
source "$(conda info --base)/etc/profile.d/conda.sh"

# This code activates and exports the conda environment to a YAML file if it doesn't already exist
activate_and_export() {
    local ENV_NAME="$1"
    local ENV_DIR="envs"
    local YAML_PATH="${ENV_DIR}/${ENV_NAME}.yml"

    mkdir -p "${ENV_DIR}"

    echo "Activating conda env: ${ENV_NAME}"
    conda activate "${ENV_NAME}"

    if [[ ! -f "${YAML_PATH}" ]]; then
        echo "No YAML found for ${ENV_NAME}. Exporting to ${YAML_PATH}..."
        conda env export --no-builds > "${YAML_PATH}"

        if [[ $? -eq 0 ]]; then
            echo "Environment exported to ${YAML_PATH}"
        else
            echo "Failed to export ${ENV_NAME}" >&2
            exit 1
        fi
    else
        echo "Found existing env file: ${YAML_PATH}"
    fi
}

# =============================================================================
# Load configuration
# =============================================================================

#Access config files
CONFIG_FILE="config/project_${PROJECT_ID}.yaml"

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "ERROR: Config file not found: ${CONFIG_FILE}" >&2
  exit 1
fi

export CONFIG_FILE

#Use yq to read variables from the config file 
PROJECT_NAME=$(yq -r '.project.name' "$CONFIG_FILE")


# Normalize possible "null" strings to empty
for v in REF_TAX REF_ASSEMBLY_LEVEL REFERENCE_ONLY; do
  [[ "${!v:-}" == "null" ]] && declare "$v"=""
done

export PROJECT_NAME
export REF_TAX
export REF_ASSEMBLY_LEVEL

# Grab the accession numbers from SRA and Reference Genomes as a list. 
ensure_list() {
  local YAML_PATH="$1"
  local FILE="$2"

  yq -r "
    ${YAML_PATH} as \$v |
    if \$v == null then
      empty
    elif (\$v | type) == \"!!seq\" then
      \$v[]
    else
      \$v
    end
  " "$FILE"
}

REF_GENOMES=()
while IFS= read -r line; do
  [[ -z "$line" ]] && continue
  REF_GENOMES+=("$line")
done < <(yq -r '.reference_genomes.accession_list[]?' "$CONFIG_FILE")


SRA_GENOMES=()
while IFS= read -r line; do
  [[ -z "$line" ]] && continue
  SRA_GENOMES+=("$line")
done < <(yq -r '.SRA_download.accession_list[]?' "$CONFIG_FILE")


#####################################
# Make project-specific directories
#####################################

activate_and_export "ncbi_datasets"

mkdir -p \
  "outputs/${PROJECT_NAME}/reference_genomes" \
  "outputs/${PROJECT_NAME}/logs" \
  "outputs/${PROJECT_NAME}/all_genomes" \
 
BASE_DIR="outputs/${PROJECT_NAME}"
LOG_DIR="outputs/${PROJECT_NAME}/logs"

# =================================================
# =================================================
# Download Reference Genomes
# =================================================
# =================================================


#Keep this outside the loop so we can access it later even if it's empty
REF_GENOMES_DIR="outputs/${PROJECT_NAME}/reference_genomes"
export REF_GENOMES_DIR

#### ----> Only run if there are actually reference genomes listed! 
if [[ ${#REF_GENOMES[@]} -gt 0 ]]; then

    mkdir -p \
    "outputs/${PROJECT_NAME}/reference_genomes" \

    REF_LOG_FILE="${LOG_DIR}/ref_genomes_log.txt"

    # Write header if file does not already exist
    if [[ ! -f "$REF_LOG_FILE" ]]; then
    echo "accession,organism,assembly_level,contig_count,total_length,gc_percent" > "$REF_LOG_FILE"
    fi

    ########################################
    # Download specific reference genomes


    for item in "${REF_GENOMES[@]}"; do
        [[ -z "$item" ]] && continue

        fasta_path="$REF_GENOMES_DIR/$item.fasta"

        # Skip if genome already exists
        if [[ -s "$fasta_path" ]]; then
            echo "[INFO] Genome for $item already exists at $fasta_path – skipping download."
            continue
        fi

        echo "[INFO] Downloading specific reference genome: $item"

        datasets download genome accession "$item" \
            --include genome \
            --filename "$REF_GENOMES_DIR/$item.zip"

        unzip -o "$REF_GENOMES_DIR/$item.zip" \
            -d "$REF_GENOMES_DIR/${item}_dir"

        mv "$REF_GENOMES_DIR/${item}_dir/ncbi_dataset/data/"*/*.fna \
            "$fasta_path"

        rm -rf "$REF_GENOMES_DIR/$item.zip" "$REF_GENOMES_DIR/${item}_dir"

        # Log metadata for JUST DOWNLOADED genome
        datasets summary genome accession "$item" \
        --as-json-lines | jq -r '
            [
            .accession,
            (.organism.organism_name // "NA"),
            (.assembly_info.assembly_level // "NA"),
            (.assembly_stats.number_of_contigs // "NA"),
            (.assembly_stats.total_sequence_length // "NA"),
            (.assembly_stats.gc_percent // "NA")
            ] | @csv
        ' >> "${REF_LOG_FILE}"

        echo "[INFO] Logged metadata for $item"
    done


    # Rename all the files because they have big ugly names
    for f in "${REF_GENOMES_DIR}"/*.fasta; do
    base=$(basename "$f")
    dir=$(dirname "$f")

    new=$(echo "$base" | sed -E '
        s/^(ERR[0-9]+).*/\1.fasta/;
        s/^((GC[AF]_[0-9]+\.[0-9]+)).*/\1.fasta/
    ')

    if [[ "$base" != "$new" ]]; then
        mv -v "$f" "$dir/$new"
    fi
    done


    # ########################################
    # Create summary list 
    # ########################################

    OUTFILE="$LOG_DIR/ref_genome_list.txt"
    echo "genome_accession description" > "$OUTFILE"

    shopt -s nullglob 
    for f in "$REF_GENOMES_DIR"/*.fasta; do
    [ -f "$f" ] || continue
    genome_acc=$(basename "$f" | cut -d_ -f1,2)
    header=$(grep -m1 '^>' "$f" | sed 's/^>//; s/;.*$//')
    description=$(printf '%s\n' "$header" | cut -d' ' -f2-)
    printf '%s %s\n' "$genome_acc" "$description" >> "$OUTFILE"
    done


fi #Ends check for empty REF_GENOMES


# =================================================
# =================================================
# Download SRA genomes and assemble
# =================================================
# =================================================


#Keep this outside the loop so we can access it later even if it's empty
SRA_DIR="outputs/${PROJECT_NAME}/sra_genomes"
GENOME_ASSEMBLY_BASE="outputs/${PROJECT_NAME}/sra_assemblies"
FINAL_ASSEMBLY="outputs/${PROJECT_NAME}/sra_assemblies/final_assemblies"
export SRA_DIR FINAL_ASSEMBLY GENOME_ASSEMBLY_BASE


#### ----> Only run if there are actually reference genomes listed! 
if [[ ${#SRA_GENOMES[@]} -gt 0 ]]; then

    #DOWNLOADING SRA GENOMES
    mkdir -p \
        "${SRA_DIR}" \
        "${FINAL_ASSEMBLY}"

    for ACC in "${SRA_GENOMES[@]}"; do
        export ACC

        [[ -z "$ACC" ]] && continue

        # Skip the download if it's already there
        if compgen -G "${SRA_DIR}/${ACC}_*.fastq.gz" > /dev/null; then
            echo "[INFO] FASTQs already exist for $ACC – skipping download."
        else
        # Otherwise, download it
            echo "[INFO] Downloading SRA for ${ACC}"

            prefetch "${ACC}" --output-directory "${SRA_DIR}"

            fasterq-dump "${SRA_DIR}/${ACC}/${ACC}.sra" \
            --split-files \
            --outdir "${SRA_DIR}"

            gzip -f "${SRA_DIR}/${ACC}"_*.fastq
            rm -rf "${SRA_DIR}/${ACC}"
            rm -f "${SRA_DIR}/${ACC}.sra"
        fi
    

        #ASEEMBLING SRA GENOMES
        # Only assemble if it hasn't already been done
        if compgen -G "$FINAL_ASSEMBLY/${ACC}.fasta" > /dev/null; then
            echo "[INFO] Assemblies already exist for $ACC – skipping assembly step."
        else
            conda deactivate # I don't understand why but I can't use wgs_spades for SRA download, oh well.
            activate_and_export "wgs_spades"

            RAW_R1="$SRA_DIR/${ACC}_1.fastq.gz"
            RAW_R2="$SRA_DIR/${ACC}_2.fastq.gz"

            # Check that FASTQs exist
            if [[ ! -f "$RAW_R1" || ! -f "$RAW_R2" ]]; then
                echo "[ERROR] FASTQ files not found: $RAW_R1 or $RAW_R2" >&2
                exit 1
            fi

            POLISH_DIR="$GENOME_ASSEMBLY_BASE/polishing"    
            SPADES_DIR="$GENOME_ASSEMBLY_BASE/spades/$ACC"
            POLISH_DIR_STATS="$GENOME_ASSEMBLY_BASE/polishing/stats"
            PILON_DIR="$POLISH_DIR/pilon1"

            mkdir -p \
            "$GENOME_ASSEMBLY_BASE/qc/raw" \
            "$GENOME_ASSEMBLY_BASE/qc/trimmed" \
            "$GENOME_ASSEMBLY_BASE/trimmed" \
            "$POLISH_DIR_STATS" \
            "$PILON_DIR" \
            "$GENOME_ASSEMBLY_BASE/logs" \
            "$GENOME_ASSEMBLY_BASE/spades/genomes"

            TRIM_R1="$GENOME_ASSEMBLY_BASE/trimmed/${ACC}_R1.trimmed.fastq.gz"
            TRIM_R2="$GENOME_ASSEMBLY_BASE/trimmed/${ACC}_R2.trimmed.fastq.gz"
            SPADES_CONTIGS="$GENOME_ASSEMBLY_BASE/spades/genomes/${ACC}_contigs.fasta"
            BAM="$POLISH_DIR/${ACC}.bam"


            # --------- QC RAW ---------
            fastqc "$RAW_R1" "$RAW_R2" \
            --outdir "$GENOME_ASSEMBLY_BASE/qc/raw" \
            --threads 4

            # --------- TRIMMING ---------
            fastp \
            -i "$RAW_R1" \
            -I "$RAW_R2" \
            -o "$TRIM_R1" \
            -O "$TRIM_R2" \
            --detect_adapter_for_pe \
            --cut_front \
            --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality 20 \
            --length_required 50 \
            --thread 8 \
            --html "$GENOME_ASSEMBLY_BASE/qc/trimmed/${ACC}_fastp.html" \
            --json "$GENOME_ASSEMBLY_BASE/qc/trimmed/${ACC}_fastp.json" \
            > "$LOG_DIR/${ACC}_fastp.log" 2>&1

            fastqc "$TRIM_R1" "$TRIM_R2" \
            --outdir "$GENOME_ASSEMBLY_BASE/qc/trimmed" \
            --threads 4

            # --------- SPADES ASSEMBLY ---------
            mkdir -p "$SPADES_DIR"

            spades.py \
            -1 "$TRIM_R1" \
            -2 "$TRIM_R2" \
            -o "$SPADES_DIR" \
            -k 21,33,55,77,99,127 \
            --careful \
            -t 8 \
            -m 24

            mv ${SPADES_DIR}/contigs.fasta ${SPADES_CONTIGS}
            mv "${SPADES_DIR}/spades.log" "${LOG_DIR}/${ACC}_spades.log"
            rm -rf "$SPADES_DIR"


            # --------- READ MAPPING ---------
            bwa index "$SPADES_CONTIGS"

            bwa mem -t 8 \
            "$SPADES_CONTIGS" \
            "$TRIM_R1" \
            "$TRIM_R2" \
            | samtools sort -o "$BAM"

            samtools index "$BAM"

            samtools flagstat "$BAM" > "$POLISH_DIR_STATS/${ACC}_flagstat.txt"

            samtools depth "$BAM" \
            | awk '{sum+=$3; n++} END {print "Mean coverage:", sum/n}' \
            > "$POLISH_DIR_STATS/${ACC}_mean_coverage.txt"

            samtools depth "$BAM" \
            | awk '{cov[$1]+=$3; len[$1]++} END{for (c in cov) print c, cov[c]/len[c]}' \
            | sort -k2,2n > "$POLISH_DIR_STATS/${ACC}_contig_coverage.txt"

            # --------- PILON POLISHING ---------
            pilon \
            --genome "$SPADES_CONTIGS" \
            --frags "$BAM" \
            --output "${ACC}_pilon1" \
            --outdir "$PILON_DIR" \
            --threads 8 \
            --fix all

            # --------- ASSEMBLY STATS ---------
            MASTER="${LOG_DIR}/all_samples_spades_pilon_assembly_stats.tsv"

            # If the master file doesn't exist yet, write header + rows
            if [[ ! -f "${MASTER}" ]]; then
            seqkit stats -T --all \
                "${SPADES_CONTIGS}" \
                "${PILON_DIR}/${ACC}_pilon1.fasta" \
            | awk -v acc="$ACC" '
                BEGIN { OFS = "\t" }
                NR==1 { print "ACC", $0; next }   # header row
                NR>1 { print acc, $0 }            # data rows
                ' > "$MASTER"

            # If it does exist, append just the data rows (skip header)
            else
            seqkit stats -T --all \
                "${SPADES_CONTIGS}" \
                "${PILON_DIR}/${ACC}_pilon1.fasta" \
            | awk -v acc="${ACC}" 'NR>1 { print acc, $0 }' >> "${MASTER}"
            fi

            echo "==> Finished ${ACC}"

            # Clean temporary bam + index and fastq files
            rm -f "$POLISH_DIR"/*.bam*
            rm -f "$SPADES_CONTIGS".*

            # Copy finished genome to final assemblies dir
            cp "$PILON_DIR/${ACC}_pilon1.fasta" "$FINAL_ASSEMBLY/${ACC}.fasta"

            conda deactivate

        fi #Ends check for existing final assembly
    done #Ends loop through SRA accessions

fi #Ends check for empty SRA_GENOMES

rm -rf "$GENOME_ASSEMBLY_BASE/trimmed"

conda deactivate

######################################################
# Collect all your genomes in one place
#######################################################

ALL_GENOMES_DIR="outputs/${PROJECT_NAME}/all_genomes"
mkdir -p "$ALL_GENOMES_DIR"

echo "[INFO] Collecting genomes into $ALL_GENOMES_DIR"
# From SRA_DIR
if [[ -d "$FINAL_ASSEMBLY" ]]; then
  echo "[INFO] Moving genomes from $FINAL_ASSEMBLY"
  find "$FINAL_ASSEMBLY" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) \
    -exec cp -n {} "$ALL_GENOMES_DIR"/ \;
else
  echo "[INFO] FINAL_ASSEMBLY does not exist — skipping"
fi

# From REF_GENOMES_DIR
if [[ -d "$REF_GENOMES_DIR" ]]; then
  echo "[INFO] Moving genomes from $REF_GENOMES_DIR"
  find "$REF_GENOMES_DIR" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) \
    -exec cp -n {} "$ALL_GENOMES_DIR"/ \;
else
  echo "[INFO] REF_GENOMES_DIR does not exist — skipping"
fi

# Rename files so that Anvi'o is happy in the next step
## It gets mad if there are dots
for f in "${ALL_GENOMES_DIR}"/*.fasta; do
  base=$(basename "$f")
  dir=$(dirname "$f")

  new=$(
    echo "$base" | \
      sed -E 's/^(ERR[0-9]+).*/\1.fasta/; s/^(GCF_[0-9]+).*/\1.fasta/'
  )

  if [[ "$base" != "$new" ]]; then
    mv -v "$f" "$dir/$new"
  fi
done


#====================================
#====================================
#Run Anvio phylogenomics 
#====================================
#====================================


activate_and_export "anvio-8"

FUN_HOM_INDEX_MAX=$(yq -r '.anvio_parameters.functional_homogeneity_index_threshold' "$CONFIG_FILE")


WORKING_DIR="${BASE_DIR}/anvio"
ANVIO_DB_DIR="${WORKING_DIR}/dbs"
ANVIO_PAN_DIR="${WORKING_DIR}/pangenome"
ANVIO_TREE_DIR="${WORKING_DIR}/tree"
CLEAN_GENOMES_DIR="${ANVIO_DB_DIR}/genomes"

mkdir -p \
  "${WORKING_DIR}/genomes" \
  "${ANVIO_DB_DIR}" \
  "${ANVIO_PAN_DIR}" \
  "${ANVIO_TREE_DIR}" \
  "${CLEAN_GENOMES_DIR}"

# Create genomes.txt without extension
ls "${ALL_GENOMES_DIR}"/*.fasta \
  | xargs -n1 basename \
  | sed 's/\.fasta$//' > "${WORKING_DIR}/genomes.txt"


###############################
# Reformat FASTA & build contigs DBs
###############################

while read -r g; do
  echo "Reformatting FASTA for $g ..."
  anvi-script-reformat-fasta \
    "${ALL_GENOMES_DIR}/${g}.fasta" \
    -o "${CLEAN_GENOMES_DIR}/${g}_contigs.fasta" \
    --simplify-names
done < "${WORKING_DIR}/genomes.txt"

# Make contigs databases
for g in $(cat "${WORKING_DIR}/genomes.txt"); do
  echo
  echo "Generating contigs database for $g ..."
  echo
  anvi-gen-contigs-database \
    -f "${CLEAN_GENOMES_DIR}/${g}_contigs.fasta" \
    -o "${ANVIO_DB_DIR}/${g}.db" \
    --num-threads 4 \
    -n "${g}"
done

# Get rid of clean contigs to save space
rm -rf "${CLEAN_GENOMES_DIR}"


###############################
# Annotate DBs (HMMs, tRNAs, SCG taxonomy)
###############################

# Read values from config
RUN_TRNAS=$(yq -r '.anvio_parameters.annotate_trnas' "$CONFIG_FILE")
RUN_COGS=$(yq -r '.anvio_parameters.annotate_ncbi_cogs' "$CONFIG_FILE")
RUN_KEGG=$(yq -r '.anvio_parameters.annotate_kegg_kofams' "$CONFIG_FILE")

for g in "${ANVIO_DB_DIR}"/*.db; do
  echo "Processing ${g}"

  # ALWAYS run HMMs and SCG taxonomy for trees
  echo "  ➜ Running HMMs"
  anvi-run-hmms -c "$g"
  anvi-run-scg-taxonomy -c "$g"

  if [[ "$RUN_TRNAS" == "YES" ]]; then
    echo "  ➜ Scanning tRNAs"
    anvi-scan-trnas -c "$g"
  fi

  if [[ "$RUN_COGS" == "YES" ]]; then
    echo "  ➜ Running NCBI COGs"
    anvi-run-ncbi-cogs -c "$g"
  fi

  if [[ "$RUN_KEGG" == "YES" ]]; then
    echo "  ➜ Running KEGG KOfams"
    anvi-run-kegg-kofams -c "$g"
  fi

done


###############################
# Genomes file + completeness
###############################

anvi-script-gen-genomes-file \
  --input-dir "${ANVIO_DB_DIR}" \
  -o "${WORKING_DIR}/final-genomes.txt"

anvi-estimate-genome-completeness \
  -e "${WORKING_DIR}/final-genomes.txt" \
  > "${WORKING_DIR}/genome_completeness.txt" 2>&1

rm -f "${WORKING_DIR}/genomes.txt"

###############################
# Pangenome
###############################

anvi-gen-genomes-storage \
  -e "${WORKING_DIR}/final-genomes.txt" \
  -o "${ANVIO_PAN_DIR}/GENOMES.db"

anvi-pan-genome \
  -g "${ANVIO_PAN_DIR}/GENOMES.db" \
  --project-name "PANGENOME" \
  --num-threads 4

mv PANGENOME/ "${ANVIO_PAN_DIR}/PANGENOME/"


###############################
# SCG alignment + tree
###############################

N=$(($(wc -l < "${WORKING_DIR}/final-genomes.txt") - 1))
echo "Total number of genomes: ${N}"

FILE_NAME="${ANVIO_TREE_DIR}/SCG_${FUN_HOM_INDEX_MAX}"

anvi-get-sequences-for-gene-clusters \
  -g "${ANVIO_PAN_DIR}/GENOMES.db" \
  -p "${ANVIO_PAN_DIR}/PANGENOME/PANGENOME-PAN.db" \
  -o "${FILE_NAME}.fa" \
  --concatenate-gene-clusters \
  --min-num-genomes-gene-cluster-occurs "${N}" \
  --max-num-genes-from-each-genome 1 \
  --max-functional-homogeneity-index "${FUN_HOM_INDEX_MAX}" \
  --align-with mafft

anvi-gen-phylogenomic-tree \
  -f "${FILE_NAME}.fa" \
  -o "${FILE_NAME}.nwk"


