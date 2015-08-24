#!/bin/sh

#
# Create a GEMINI database by (a) decomposing and normalizing variants;
# (b) annotating variants with VEP or snpEff; (c) loading annotated variants
# into GEMINI; (d) adding custom annotations to the database.
#
# Dependencies:
#   samtools/bgzip/tabix
#   bcftools
#   vt
#   VEP and/or snpEff
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
ANNOTATOR=${4:-$DEFAULT_ANNOTATOR}
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
    bcftools merge *.vcf.gz > ${CUR_DIR}/all.vcf

    # Input VCF is all VCFs in directory.
    INPUT_VCF=all.vcf

    popd
else
    # Single BCF.
    INPUT_VCF=${INPUT}
fi

#
# HACK specfic to cancer amplicons: replace AF with GQ to get AF (allele frequency) into database as genotype quality.
#
sed -i.bak 's/AF/GQ/g' all.vcf

# Set up name for annotated VCF.
BASE=$(basename "${INPUT_VCF}" .vcf)
ANNO_VCF="${BASE}.anno.vcf"

# Normalize and split.
vt decompose -s ${INPUT_VCF} | vt normalize -r ${REFERENCE} -o ${BASE}_decnorm.vcf -

# Annotate.
if [ ${ANNOTATOR} = "VEP" ]; then
    # Can add --refseq to use refseq annotations.
    perl ${ANNOTATOR_DIR}/variant_effect_predictor.pl -i ${BASE}_decnorm.vcf \
    --cache \
    --offline \
    --assembly GRCh37 \
    --sift b \
    --polyphen b \
    --symbol \
    --numbers \
    --biotype \
    --total_length \
    -o ${ANNO_VCF} \
    --vcf \
    --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE
elif [ ${ANNOTATOR} = "snpEff" ]; then
	java -jar ${ANNOTATOR_DIR}/snpEff.jar -i vcf -o vcf GRCh37.75 ${BASE}_decnorm.vcf > ${ANNO_VCF}
fi

# Load into GEMINI.
gemini load -v ${ANNO_VCF} -t ${ANNOTATOR} ${GEMINI_DB}

#
# HACK specfic to cancer amplicons: annotate with HP field from VCF.
#
bgzip ${BASE}_decnorm.vcf && tabix -p vcf ${BASE}_decnorm.vcf.gz
gemini annotate -f ${BASE}_decnorm.vcf.gz -t integer -a extract -c HP -e HP -o first ${GEMINI_DB}

#
# Add annotations to the database.
#
${HOME_DIR}/annotate_gemini_db.sh ${GEMINI_DB} ${CUSTOM_ANNOS_DIR}
