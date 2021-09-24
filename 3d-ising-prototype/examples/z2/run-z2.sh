#!/bin/bash
#SBATCH -n 1                   # Total number of mpi tasks requested
#SBATCH -c 4
#SBATCH -t 30             # Run time (hh:mm:ss)

julia z2-magnetization-wolff-from-ones-compute.jl 