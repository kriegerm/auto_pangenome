
##############################################################################
# ------------- Download a CSV of reference genomes summary ---------------- #
##############################################################################

# Conda env ncbi_datasets
PROJECT_ID="001"
REF_TAXA="Limosilactobacillus fermentum"
REF_ASSEMBLY_LEVEL="complete"

REF_TAXA_DIR="outputs/${PROJECT_ID}/taxa_acc_search/$(echo "$REF_TAXA" | tr ' ' '_')"
mkdir -p "${REF_TAXA_DIR}"

OUTFILE="${REF_TAXA_DIR}/all_genomes"

{
  # Header
  echo "accession,organism,assembly_level,contig_count,total_length,gc_percent"

  # Body
  datasets summary genome taxon "${REF_TAXA}" \
    --assembly-level "${REF_ASSEMBLY_LEVEL}" \
    --assembly-source refseq \
    --as-json-lines \
  | jq -r '
      [
        .accession,
        (.organism.organism_name // "NA"),
        (.assembly_info.assembly_level // "NA"),
        (.assembly_stats.number_of_contigs // "NA"),
        (.assembly_stats.total_sequence_length // "NA"),
        (.assembly_stats.gc_percent // "NA")
      ] | @csv
    '
} > "${OUTFILE}.csv"


##############################################################################
# Randomly subset CSV to get a list of accessions
##############################################################################  

# Get a random subset of n number of genomes from the CSV
number_of_genomes_to_select=40

{
  head -n 1 ${OUTFILE}.csv
  tail -n +2 ${OUTFILE}.csv | sort -R | head -n ${number_of_genomes_to_select}
} > "${REF_TAXA_DIR}/${number_of_genomes_to_select}_selected_genomes.csv"

# Create a text file of just the accessions for use in downloading genomes - copy/paste into the yaml
awk -F',' '
  NR>1 {
    gsub(/"/, "", $1);                 # remove all double quotes from field 1
    if ($1 ~ /^(GCF|GCA|ERR|SRR)/) {
      print "- \"" $1 "\""
    }
  }
' "${REF_TAXA_DIR}/${number_of_genomes_to_select}_selected_genomes.csv" \
  > "${REF_TAXA_DIR}/${number_of_genomes_to_select}_selected_genomes_for_yaml.txt"


