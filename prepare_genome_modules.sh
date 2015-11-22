#!/bin/bash
if [ $# -lt 3 ] ; then
        echo ""
        echo "usage: prepare_genome_modules.sh <genome_name> <build_number> <reference.fasta> <gene_annotation.gff>"
        echo "build index for several aligners and writes a module file on condo"
        echo "does not overwrite any exisitng files"
        echo ""
        exit 0
fi

#load requried programs to create genome module files
module purge
module load perl
module load gmap-gsnap/2015-09-29
module load parallel
module load bowtie2
module load bwa
module load gatk
module load bedtools
module load samtools


#create local variables 
NAME="$1"
BUILD="$2"
REF="$3"
GFF="$4"
GSEQ="/data003/GIF/genomes/sequences"
GMOD="/data003/GIF/genomes/modules"
#intervals set to 100kb
WINDOW=100000

mkdir -p ${GSEQ}/${NAME}/${BUILD}
mkdir -p ${GMOD}/${NAME}
touch ${GMOD}/${NAME}/${BUILD}


# write a module file

cat <<MODULEFILE > ${GMOD}/${NAME}/${BUILD}
#%Module1.0#####################################################################
##
module-whatis   "${NAME}"

#unset the general variables to make sure we don't accidently load two unrelated genomes
unsetenv   GENOME
unsetenv   GMAPDB
unsetenv   GNAME

#set general variable names when working with only one genome
setenv  GENOMEDIR        ${GSEQ}/${NAME}/${BUILD}/
setenv  GENOMEFASTA	 ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.fasta
setenv	GENOMEINTERVALS	${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_100kb_coords.bed
setenv  GNAME         ${NAME}_${BUILD}
setenv  GMAPDB        ${GSEQ}/${NAME}/${BUILD}/$GNAME
setenv  modulefile	${GMOD}/${NAME}/${BUILD}
setenv  "${NAME}_${BUILD}_genome" ${GSEQ}/${NAME}/${BUILD}/
setenv  "${NAME}_${BUILD}_GMAPDB" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}
setenv  "${NAME}_${BUILD}_GNAME" ${NAME}_${BUILD}

setenv  "${NAME}_${BUILD}_intervals100k" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_100kb_coords.bed 

setenv  "${NAME}_${BUILD}_cdna" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.fasta
setenv  "${NAME}_${BUILD}_cdna" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.gff3
setenv  "${NAME}_${BUILD}_cdna" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.cdna.fasta
setenv  "${NAME}_${BUILD}_cds" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.cds.fasta
setenv  "${NAME}_${BUILD}_gene" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.gene.fasta
setenv  "${NAME}_${BUILD}_pep" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.pep.fasta
setenv  "${NAME}_${BUILD}_upstream3000" ${NAME}_${BUILD}.upstream3000.fasta

MODULEFILE



# Create Fasta files from GFF file
/data005/GIF2/severin/isugif/common_scripts/gff2fasta.pl ${REF} ${GFF} ${NAME}_${BUILD}
mv ${NAME}_${BUILD}* ${GSEQ}/${NAME}/${BUILD}/

# build index for GSNAP, Bowtie2, BWA and SAMTOOLS
module unload perl
parallel <<FIL

gmap_build -d ${NAME}_${BUILD} -D ${GSEQ}/${NAME}/${BUILD} ${REF}
bowtie2-build ${REF} ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}
samtools faidx ${REF}
bwa index -p ${NAME}_${BUILD} -a bwtsw ${REF}
java -Xmx100G -jar /data003/GIF/software/packages/picard_tools/1.130/picard.jar CreateSequenceDictionary \
  REFERENCE=${REF} \
  OUTPUT=${NAME}_${BUILD}.dict
FIL
# cleanup
mv ${REF}.fai ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.fai
mv ${NAME}_${BUILD}* ${GSEQ}/${NAME}/${BUILD}/
mv ${NAME}_${BUILD}.dict ${GSEQ}/${NAME}/${BUILD}/
ln -s ${NAME}_${BUILD}.fai ${NAME}_${BUILD}.fasta.fai

# build intervals and cleanup
fasta_length.py ${REF} > ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_length.txt
bedtools makewindows -w ${WINDOW} -g  ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_length.txt |  awk '{print $1"\t"$2+1"\t"$3}' >  ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_100kb_coords.bed
java -Xmx100G -jar $PICARD/picard.jar BedToIntervalList \
  INPUT=${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_100kb_coords.bed \
  SEQUENCE_DICTIONARY=${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_100kb_coords.dict \
  OUTPUT=${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_100kb_gatk_intervals.list

#move reference and GFF file to genome module locations

mv ${REF} ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.fasta
mv ${GFF} ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.gff3


