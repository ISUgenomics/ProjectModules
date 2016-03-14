#!/bin/bash

main() {
	if [ $# -lt 3 ] ; then 
	displayhelp 
	fi
#	createVarFun
	writeModuleFile
	convertRef
	refStats
	createFastaFilesFromGFF $@
	createGMAPDB
	createBowtie2DB
	createFaidx
	createBWADB
	createBLASTDBREF
	createBLASTDBPEP $@
	createREFdict
	createRefIntervals
	createCDBFasta $@
}


displayhelp () {
	echo ""
        echo "usage: prepare_genome_modules.sh <genome_name> <build_number> <reference.fasta> <gene_annotation.gff>"
        echo "build index for several aligners and writes a module file on condo"
        echo "does not overwrite any exisitng files"
        echo ""
        exit 0
}


#does command exist
checkCommand () {
        if command -v $1 > /dev/null 2>&1; then
        #  echo given-command is available
          echo TRUE 
        else
        #  echo given-command is not available
          echo FALSE
        fi
}

#create Variable names 
#createVarFun () {
module load `whoami`
export	NAME="$1"
export	BUILD="$2"
export	BUILD=$(echo ${BUILD//./p})
export	REFNAME="$3"
export	GFF="$4"
export	GSEQ="$MODBASE/genome/sequences"
export	GMOD="$MODBASE/genome/modules"
	#intervals set to 100kb
export	WINDOW=100000
	mkdir -p ${GSEQ}/${NAME}/${BUILD}
	mkdir -p ${GMOD}/${NAME}
	touch ${GMOD}/${NAME}/${BUILD}
#}
echo $BUILD
echo $GSEQ
echo $GMOD
# write basic module file
#Module file is written first in order to keep everything consistent.
#May decide to remove variables where files fail to be written
writeModuleFile () {
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
setenv  GNAME         ${NAME}_${BUILD}
setenv  modulefile	${GMOD}/${NAME}/${BUILD}
setenv	VERSION		${BUILD}
setenv  ${PROG}_${BUILD}_createdate    $(date '+%m/%d/%y_%H:%M:%S')
setenv  ${PROG}_${BUILD}_creator       $(whoami)
setenv  ${PROG}_${BUILD}_dir           $(pwd -P)
MODULEFILE
}


#convert Reference into a standard format
convertRef () {
	module load emboss
	echo "Format Reference"
	local commandcheck=`checkCommand seqret`
	if  [ $commandcheck = "TRUE" ]; then 
	seqret -sequence ${REFNAME} -outseq ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.fasta
	REF="${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.fasta"
cat <<MODULEFILE >> ${GMOD}/${NAME}/${BUILD}
setenv  "${NAME}_${BUILD}_genomefasta" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.fasta
setenv  "${NAME}_${BUILD}_GNAME" ${NAME}_${BUILD}
MODULEFILE
	else
	echo "seqret program not found"
	echo "emboss program may not be installed"
	echo "or module load emboss may have failed"
	exit 1;
	fi
}


#get statistics on the genome
echo "get reference stats"
refStats () {
#	module load perl
	module load $(whoami)
	perl $COMMON_SCRIPTS/new_Assemblathon.pl ${REF}> ${NAME}_${BUILD}_stat
	mv ${NAME}_${BUILD}_stat ${GSEQ}/${NAME}/${BUILD}/
cat << MODULEFILE >> ${GMOD}/${NAME}/${BUILD}
setenv "${NAME}_${BUILD}_stat" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_stat
MODULEFILE
}


# Create Fasta files from GFF file

createFastaFilesFromGFF () {
#	module load perl
	module load $(whoami)
	if [ $# -eq 4 ] ; then
	echo "Generate Fasta files from GFF file"
	perl $COMMON_SCRIPTS/gff2fasta.pl ${REF} ${GFF} ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}
cat <<MODULEFILE >> ${GMOD}/${NAME}/${BUILD}
setenv  GENOMEGFF3      ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.gff3
setenv  "${NAME}_${BUILD}_gff3" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.gff3
setenv  "${NAME}_${BUILD}_cdna" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.cdna.fasta
setenv  "${NAME}_${BUILD}_cds" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.cds.fasta
setenv  "${NAME}_${BUILD}_gene" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.gene.fasta
setenv  "${NAME}_${BUILD}_pep" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.pep.fasta
setenv  "${NAME}_${BUILD}_upstream3000" ${NAME}_${BUILD}.upstream3000.fasta
MODULEFILE
	fi
}


# build index for GSNAP, Bowtie2, BWA and SAMTOOLS

createGMAPDB () {
	module load gmap-gsnap
	local commandcheck=`checkCommand gmap_build`
	if  [ $commandcheck = "TRUE" ]; then
	gmap_build -d ${NAME}_${BUILD} -D ${GSEQ}/${NAME}/${BUILD} ${REF}
cat <<MODULEFILE >> ${GMOD}/${NAME}/${BUILD}
setenv  "${NAME}_${BUILD}_GMAPDB" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}
setenv  GMAPDB        ${GSEQ}/${NAME}/${BUILD}/$GNAME
MODULEFILE
	else
	echo "gmap_build script not found"
	echo "gsnap may not be installed"
	echo "or module load gmap-gsnap failed"
	fi
}


createBowtie2DB () {
	module load bowtie2
        local commandcheck=`checkCommand bowtie2-build`
        if  [ $commandcheck = "TRUE" ]; then
	bowtie2-build ${REF} ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}
cat <<MODULEFILE >> ${GMOD}/${NAME}/${BUILD}
setenv  "${NAME}_${BUILD}_BowtieDB" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}
MODULEFILE
	else
        echo "bowtie2-build script not found"
        echo "bowtie2 may not be installed"
        echo "or module load bowtie failed"
        fi
}

createFaidx () {
	module load samtools
        local commandcheck=`checkCommand samtools`
        if  [ $commandcheck = "TRUE" ]; then
	samtools faidx ${REF}
cat <<MODULEFILE >> ${GMOD}/${NAME}/${BUILD}
setenv  "${NAME}_${BUILD}_genomefasta_fai" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.fasta.fai
MODULEFILE
	else
	echo "samtools script not found"
        echo "samtools may not be installed"
        echo "or module load samtools failed"
        fi
}

createBWADB () {
	module load bwa
        local commandcheck=`checkCommand bwa`
        if  [ $commandcheck = "TRUE" ]; then
	bwa index -p ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD} -a bwtsw ${REF}
cat <<MODULEFILE >> ${GMOD}/${NAME}/${BUILD}
setenv  "${NAME}_${BUILD}_genomefasta_bwaDB" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_bwaDB
MODULEFILE
	else
        echo "bwa script not found"
        echo "bwa may not be installed"
        echo "or module load bwa failed"
        fi
}


createBLASTDBREF () {
	module load ncbi-blast
        local commandcheck=`checkCommand makeblastdb`
        if  [ $commandcheck = "TRUE" ]; then
	makeblastdb -in ${REF} -dbtype 'nucl' -out ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.blastdb
cat <<MODULEFILE >> ${GMOD}/${NAME}/${BUILD}
setenv  "${NAME}_${BUILD}_genomefasta_blastDB" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_blastDB
MODULEFILE
	else
        echo "makeblastdb script not found"
        echo "makeblastdb may not be installed"
        echo "or module load makeblastdb failed"
        fi
}

createBLASTDBPEP () {
        module load ncbi-blast
        local commandcheck=`checkCommand makeblastdb`
        if  [ $commandcheck = "TRUE" ]; then
	 if [ $# -eq 4 ] ; then
        makeblastdb -in ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.pep.fasta  -dbtype 'prot' -out ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.pep.blastdb
cat <<MODULEFILE >> ${GMOD}/${NAME}/${BUILD}
setenv  "${NAME}_${BUILD}_pep_blastDB" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.pep_blastDB
MODULEFILE
	 fi
	else
        echo "makeblastdb script not found"
        echo "makeblastdb may not be installed"
        echo "or module load makeblastdb failed"
        fi
}

createREFdict () {
	module load picard_tools
	java -Xmx100G -jar $PICARD/picard.jar CreateSequenceDictionary \
	  REFERENCE=${REF} \
	  OUTPUT=${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.dict
cat <<MODULEFILE >> ${GMOD}/${NAME}/${BUILD}
setenv  "${NAME}_${BUILD}_genomefasta_dict" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.dict
MODULEFILE
}




# build intervals and cleanup

createRefIntervals () {
	module load bedtools
        local commandcheck=`checkCommand bedtools`
        if  [ $commandcheck = "TRUE" ]; then
	cut -f 1,2 ${REF}.fai > ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_length.txt
	bedtools makewindows -w ${WINDOW} -g  ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_length.txt |  awk '{print $1"\t"$2+1"\t"$3}' >  ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_100kb_coords.bed
	java -Xmx100G -jar $PICARD/picard.jar BedToIntervalList \
	  INPUT=${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_100kb_coords.bed \
	  SEQUENCE_DICTIONARY=${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.dict \
	  OUTPUT=${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_100kb_gatk_intervals.list
cat <<MODULEFILE >> ${GMOD}/${NAME}/${BUILD}
setenv  GENOMEINTERVALS ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_100kb_gatk_coords.bed
setenv  "${NAME}_${BUILD}_genome_intervals" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}_100kb_gatk_coords.bed 
MODULEFILE
	else
        echo "bedtools script not found"
        echo "bedtools may not be installed"
        echo "or module load bedtools failed"	
	fi
}
createCDBFasta () {
	if [ $# -eq 4 ] ; then
		module load cdbfasta
		cdbfasta ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.cdna.fasta
		cdbfasta ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.cds.fasta
		cdbfasta ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.gene.fasta
		cdbfasta ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.pep.fasta
		cdbfasta ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.upstream3000.fasta 
		cat <<MODULEFILE >> ${GMOD}/${NAME}/${BUILD}
setenv  "${NAME}_${BUILD}_cdna_cidx" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.cdna.fasta.cidx
setenv  "${NAME}_${BUILD}_cds_cidx" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.cds.fasta.cidx
setenv  "${NAME}_${BUILD}_gene_cidx" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.gene.fasta.cidx
setenv  "${NAME}_${BUILD}_pep_cidx" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.pep.fasta.cidx
setenv  "${NAME}_${BUILD}_upstream3000_cidx" ${GSEQ}/${NAME}/${BUILD}/${NAME}_${BUILD}.upstream3000.fasta.cidx
MODULEFILE
	fi

}

cleanup () {
rm $REFNAME.index 
}


#start from a clean module list
module purge
module load `whoami`
main $@
