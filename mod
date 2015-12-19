#!/bin/bash


scriptName="mod"
function printUsage () {
    cat <<EOF

Synopsis

    $scriptName [-h | --help] COMMAND [ARGS]

The most commonly used mod commands are:
	init:		Initialize a software, project, user or genome module
        add:		Add an environmental variable to a module    
	rm:		Remove an environmental variable from a module
	prepend:	Add Prepend statement for a directory to a variable such as PATH, LIB, PERL5_LIB, LD_LIBRARY_PATH etc

Author

    Andrew Severin, Genome Informatics Facilty, Iowa State University
    severin@iastate.edu
    15 December, 2015


EOF
}


if [ $# -lt 1 ] ; then
        printUsage
        exit 0
fi

BASE=/data003/GIF/
#BASE=/data006a/GIF_2a/
#software
SOFTMODDIR=$BASE/software/modules/
#Project
PROJMODDIR=$BASE/project/modules
#GENOME
GENMODDIR=$BASE/genomes/modules/
#USER
USERMODDIR=$BASE/user/modules

if [ "$1" = "rm" ]; then 
	shift;
	modrm $@
fi

if [ "$1" = "add" ]; then
	shift;
        modadd $@
fi

if [ "$1" = "modprepend" ]; then
        shift;
        modprepend $@
fi

if [ "$1" = "init" ]; then
	shift;
        modinit $@
fi
