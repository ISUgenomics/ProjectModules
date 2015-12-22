# INSTALLATION for mod wrapper

## Download the ProjectModules repository
git clone https://github.com/trinityrnaseq/trinityrnaseq.git

##Create a user module
cd ProjectModules
./mod init -u

This will ask you to set a module file directory
I put my module files in /data003/GIF on my supercomputer at ISU called condo.

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

##Add the following to your .bashrc

module load $(whoami)


