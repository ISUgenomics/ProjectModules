# INSTALLATION for mod wrapper

Download the ProjectModules repository

this repository has dependencies of other isugif repositories which we will download and place in your home directory

it is important that they can be found in ~/isugif/

```
cd
mkdir isugif
cd isugif
git clone git@github.com:ISUgenomics/ProjectModules.git 
git clone git@github.com:ISUgenomics/common_scripts.git
git clone git@github.com:ISUgenomics/common_analyses.git
#The isugif fold can be located anywhere so long as a softlink to that location is found in your home directory.

#if this times out with the following error message 
ssh: connect to host github.com port 22: Connection timed out
then perform these commands

cd
mkdir isugif
cd isugif
git clone https://github.com/ISUgenomics/ProjectModules.git
git clone https://github.com/ISUgenomics/common_scripts.git
git clone https://github.com/ISUgenomics/common_analyses.git
#The isugif fold can be located anywhere so long as a softlink to that location is found in your home directory.```
```

##Create a user module
```
cd ProjectModules
./mod init -u
```


This will ask for a module file directory
Consider using /home/$(whoami)/privatemodules/
I put my module files in /data003/GIF on my supercomputer at ISU calld condo.
This is a local directory that I have permission to write to.

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
module use /home/$(whoami)/privatemodules/project/modules/
module use /home/$(whoami)/privatemodules/genomes/modules/
module use /home/$(whoami)/privatemodules/user/modules/
module use /home/$(whoami)/privatemodules/software/modules/
module load $(whoami)


