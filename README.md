
# ProjectModules
This program is a wrapper that makes module file creation fast, efficient and reproducible for a multi-member lab. This wrapper script also extends the concept of modules from software version control to project, user and genome version control. The goal of the project is to make it super easy to generate module files and add/remove environmental variables that point to files, folders and text to the module file all from the command line.  It works on all module files to which you have permission to write.  It requires the location of a directory to which you have permission to write (MODBASE).  See 

https://github.com/ISUgenomics/ProjectModules/blob/master/INSTALLATION.md and

http://gif.biotech.iastate.edu/Tutorial/doku.php?id=blog:introducing_mod_wrapper_script_for_unix_module_command 

for more information

The first time you run "mod init -u" it will create a user module and module folders in the specified directory ($MODBASE) as follows.
```
$MODBASE
`--user
    `-- modules
        `-- $(whoami)
`--software
    `--modules
    `--packages
`--genomes
    `--modules
`--project
```

#COMMANDS

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

    modinit [-h | --help] [-s -p -g -u]

Description:
        -s:     to start a software module      modinit -s
        -p:     to start a project module       modinit -p
        -g:     to start a genome module        modinit -g
        -u:     to start a user module          modinit -u

Software Modules:       Version control for Software. Standard use of unix modules made easy.
Project Modules:        Version control for Projects. Create env variables for important files or directories.  Summarize project.      
Genome Module:          Version control for Genomes.  Creates commonly used databases and useful fasta files from genome and GFF.
User Module:            Version control for Users. Sets MODBASE variable, offloads environmental variables to a loadable module instead of .bashrc

```
mod add
```
Synopsis

    mod add [-h | --help] <modulename> <VariableNAME> <Filename>

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

    mod rm [-h | --help] <modulename> <VariableNAME> 

Description:
        modulename:     is the current module but does not include the version
        VariableNAME:   is the variable you want to remove from the module file


```
Author

    Andrew Severin, Genome Informatics Facilty, Iowa State University
    severin@iastate.edu
    15 December, 2015


