#!/bin/bash
#
# Laucher batch script file for TACC systems (like Frontera, Stampede2, etc.)
# Si Liu
# July 13, 2020
#
# Simple SLURM script for submitting multiple serial
# jobs (e.g. parametric studies) using a script wrapper
# to launch the jobs.
#
# To use, build the launcher executable and your
# serial application(s) and place them in your WORKDIR
# directory.  Then, edit the LAUNCHER_JOB_FILE to specify 
# each executable per process.
#-------------------------------------------------------
# 
#         <------ Setup Parameters ------>
#
#SBATCH -J NP_C8
#SBATCH -N 2 
#SBATCH -n 64
#SBATCH -p normal
#SBATCH -o Parametric.%j.out
#SBATCH -e Parametric.%j.err
#SBATCH --mail-user=rathbunc@dickinson.edu
#SBATCH --mail-type=all
#SBATCH -t 01:30:00
#SBATCH -A TG-BIO220002

#------------------------------------------------------

export ROSETTA_DATABASE=$WORK/rosetta_database/database
module load launcher
module load intel/18.0.2
module load impi/18.0.2
module load rosetta/3.9

export LAUNCHER_WORKDIR=$WORK/runs/NP_C8
export LAUNCHER_JOB_FILE='job.job'

$LAUNCHER_DIR/paramrun

