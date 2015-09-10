#!/bin/sh

#
# Create a GEMINI database by (a) decomposing and normalizing variants;
# (b) annotating variants with VEP or snpEff; (c) loading annotated variants
# into GEMINI; (d) adding custom annotations to the database.
#
# Dependencies:
#   samtools/bgzip/tabix
#   GNU parallel
#   bcftools
#   VEP
#   GEMINI
#

# Parameter checking.
if [ $# -lt "2" ]
then
  echo "Usage: `basename $0` <directory of VCFs or single VCF> <db_name> [genome_reference] [VEP/snpEff] [Annotator directory] [custom annos directory]>"
  exit -1
fi

# Set up home directory and default settings.
HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${HOME_DIR}/default_settings.sh

# Set parameters.
INPUT=$1
GEMINI_DB=$2
REFERENCE=${3:-$DEFAULT_REFERENCE}
ANNOTATOR=${4:-$DEFAULT_ANNOTATOR}   # VEP
ANNOTATOR_DIR=${5:-$DEFAULT_ANNOTATOR_DIR}
CUSTOM_ANNOS_DIR=${6:-$DEFAULT_CUSTOM_ANNOS_DIR}


# If there is a directory of VCFs, combine them into a VCF.
CUR_DIR=${PWD}
if [ -d ${INPUT} ]; then
    # Combine all files in directory into a single file.

    # Go to directory.
    pushd ${INPUT}

    # Compress and index VCFs.
    ${HOME_DIR}/bgzip_and_tabix.sh

    # Merge VCFs.
    bcftools merge -m none *.vcf.gz > ${CUR_DIR}/all.vcf

    # Input VCF is all VCFs in directory.
    INPUT_VCF=all.vcf

    popd
else
    # Single BCF.
    INPUT_VCF=${INPUT}
fi

# Annotate.
BASE=$(basename "${INPUT_VCF}" .vcf)
TEMP_VCF="${BASE}.anno.vcf.temp"
perl ${ANNOTATOR_DIR}/variant_effect_predictor.pl \
    -i $INPUT_VCF \
    -o $TEMP_VCF \
    --offline \
    --sift b \
    --polyphen b \
    --symbol \
    --numbers \
    --regulatory \
    --biotype \
    --total_length \
    --vcf \
    --fields Consequence,BIOTYPE,IMPACT,Codons,Amino_acids,Protein_position,EXON,PolyPhen,SIFT,Gene,Feature,SYMBOL

# Set up name for annotated VCF.
ANNO_VCF="${BASE}.anno.vcf"

# reorder VEP annotations; change 'AF' to 'GQ'
perl VEPtoGemini.pl $TEMP_VCF $ANNO_VCF -gq

# Load into GEMINI.
gemini load -v $ANNO_VCF -t ${ANNOTATOR} ${GEMINI_DB}

#
# HACK specfic to cancer amplicons: annotate with HP field from VCF.
#
bgzip $INPUT_VCF
tabix -p vcf ${INPUT_VCF}.gz
gemini annotate -f ${INPUT_VCF}.gz -t integer -a extract -c HP -e HP -o first ${GEMINI_DB}

#
# Add annotations to the database.
#
${HOME_DIR}/annotate_gemini_db.sh ${GEMINI_DB} ${CUSTOM_ANNOS_DIR}
