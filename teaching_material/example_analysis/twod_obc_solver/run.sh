#!/usr/bin/env zsh
#BSUB -J LEVELSETFeSiW
#BSUB -o LEVELSETFeSiW.%J
#BSUB -e LEVELSETFeSiW.e%J
#BSUB -M 2000000
#BSUB -W 72:00
#BSUB -u kuehbach@imm.rwth-aachen.de
#BSUB -N
#BSUB -n 128
#BSUB -a "bcs openmp" 
#BSUB -x
###BSUB -P jara0076
#BSUB -R "select[hpcwork_fast]"

###no_numa_balancing LevelSet_IMM FeSiWithSEECuboidsProduction.xml
twod_obc_solver FeSiWithSEECuboidsProduction.xml 1>FeSiWithSEECuboidsProduction.xml.STDOUT.txt 2>FeSiWithSEECuboidsProduction.xml.STDERR.txt

