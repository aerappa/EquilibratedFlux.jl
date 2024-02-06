#!/bin/bash

#SBATCH --job-name=eqflux_scaling    # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task	
#SBATCH --time=02:00:00              # Time limit hrs:min:sec
#SBATCH --exclusive            
#SBATCH --mem=150G
#SBATCH --partition=cpu_homogen

# Your Julia script name
julia_script="test/ParallelTest.jl"

# Loop through the thread values
echo "Running Julia script with $nb_threads thread(s) and order $p"
julia --project --threads $nb_threads $julia_script --output-dir $out_dir -p $p
echo "Finished running with $nb_threads thread(s) and order $p"
