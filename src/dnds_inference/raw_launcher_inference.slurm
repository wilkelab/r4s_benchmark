#!/bin/bash
#
# Simple SLURM script for submitting multiple serial
# jobs (e.g. parametric studies) using a script wrapper
# to launch the jobs.
#
# To use, build the launcher executable and your
# serial application(s) and place them in your WORKDIR
# directory.  Then, edit the CONTROL_FILE to specify 
# each executable per process.
#-------------------------------------------------------
#-------------------------------------------------------
# 
#         <------ Setup Parameters ------>
#
#SBATCH -J inference          # Job name
#SBATCH -n 10       
#SBATCH -N 10           
#SBATCH -p normal             # Queue name
#SBATCH -o o.%j      # Name of stdout output file (%j expands to jobid)
#SBATCH -e e.%j      # Name of stdout output file (%j expands to jobid)
#SBATCH -t 09:00:00            # Run time (hh:mm:ss)
#SBATCH -A A-bio7
#------------------------------------------------------
module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=PLACEHOLDER
export LAUNCHER_WORKDIR=`pwd`
 
$LAUNCHER_DIR/paramrun

