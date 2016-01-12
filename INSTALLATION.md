# INSTALLATION for mod wrapper

Download the ProjectModules repository

this repository has dependencies of other isugif repositories which we will download and place in your home directory

it is important that they can be found in ~/isugif/

cd
mkdir isugif
cd isugif
git clone git@github.com:ISUgenomics/ProjectModules.git 
git clone git@github.com:ISUgenomics/common_scripts.git
git clone git@github.com:ISUgenomics/common_analyses.git


##Create a user module
cd ProjectModules
./mod init -u

This will ask you to set a module file directory
I put my module files in /data003/GIF on my supercomputer at ISU called condo.
```
The first time you run this command it will create the following folders
/data003/GIF/
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
##Add the following to your .bashrc

module load $(whoami)


