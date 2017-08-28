HardRod

Michael Stefferson

Description: Uses DDFT to solve active or passive continuum equations for colloid particles (disks, rods, spheres) with various short and long range interactions and external potentials.

Set-up:
> git pull https://github.com/mstefferson/HardRodML
> ./cpParams

Parameter File to Edit:
initParams.m: This should be the only file that needs to be edited. This acts as an input file to the code. If you do not see initParams in the working directory, run ./cpParams in the terminal (or matlab).

Executeables:
runHardRod.m: matlab executeable that solves for the time evolution of the system. Uses initParams.m as the file input. Outputs are saved in either ./runfiles (MakeOP = 0, MakeMovies = 0), ./runOPfiles (MakeOP = 1, MakeMovies = 0), or analyzedfiles (MakeOP = 1, MakeMovies = 1) depending on flags.
opHardRod.m: Makes the order parameters for files in ./runfiles/ and moves them to ./runOPfiles/
movieHardRod.m: Makes movies and some figures for files in ./runOPfiles/ and moves them to ./analyzedfiles/
continueHr.m: Continues an unfinished run.

Cluster jobs:
summitSlurmRun.sh: A slurm job submit file to run runHardRod. 
To run:
> sbatch summitSlurmRun.sh
summitSlurmCont.sh: A slurm job submit file to continue unfinished runs. 
To run: 
> sbatch summitSlurmRun.sh


Folders:
src/: Contains all of the source code
runfiles/, runOPfiles/, analyzedfiles/: output directories
