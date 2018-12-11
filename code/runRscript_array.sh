#! /bin/bash

#$ -cwd
#$ -S /bin/bash
#$ -j y

args=("$@") # all arguments
unset args[0] # remove first argument (R script name)

len=${#args[@]} # number of arguments

# e.g. to export R library path:
export R_LIBS=$R_LIBS:/projects/resources/gactools/R_packages/library

R -q --vanilla --args ${args[@]} $SGE_TASK_ID < $1
