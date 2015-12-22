# ProjectModules

mod
```
Synopsis

    mod [-h | --help] COMMAND [ARGS]

The most commonly used mod commands are:
        init:           Initialize a software, project, user or genome module
        add:            Add an environmental variable to a module    
        rm:             Remove an environmental variable from a module
        prepend:        Add Prepend statement for a directory to a variable such as PATH, LIB, PERL5_LIB, LD_LIBRARY_PATH etc
```

mod init
```
Synopsis

    modinit [-h | --help] [-s -p]

Description:
        -s:     to start a software module      modinit -s
        -p:     to start a project module       modinit -p
        -g:     to start a genome module        modinit -g

Software modules:       Version control for Software. Standard use of unix modules made easy.
Project modules:        Version control for Projects. Create env variables for important files or directories.  Summarize project.
Genome Module           Version control for Genomes.  Creates commonly used databases and useful fasta files from genome and GFF.
```
mod add
```
Synopsis

    ma [-h | --help] <modulename> <VariableNAME> <Filename>

Description:
        modulename:     is the current module but does not include the version
        VariableNAME:   is the variable you want to add to the module file
        FileName:       is the name of the file or text you wish you put into a variable
```
mod prepend
```
Synopsis

    mod prepend [-h | --help] <modulename> <VariableNAME> <DirName>

Description:
        modulename:     is the current module but does not include the version
        VariableNAME:   is the variable you want to add to the module file
        DirName:       is the name of the directory or text you wish you put into a variable

```
mod rm
```
Synopsis

    modrm [-h | --help] <modulename> <VariableNAME> 

Description:
        modulename:     is the current module but does not include the version
        VariableNAME:   is the variable you want to remove from the module file


```
Author

    Andrew Severin, Genome Informatics Facilty, Iowa State University
    severin@iastate.edu
    15 December, 2015


