#!/bin/bash

# Total number of runs
NUM_RUNS=4

# Get number of available cores
NUM_CORES=4

for i in $(seq 1 $NUM_RUNS); do
    # Compute core number (cycling through available cores)
    CORE=$(( (i - 1) % NUM_CORES ))

    echo "Submitting job $i to core $CORE..."
    ts taskset -c $CORE /home/yazminbentancour/Documents/MATLAB/bin/matlab \
        -nodisplay -nosplash -nodesktop -singleCompThread \
        -r "parameter_analysis_fun($i); exit;"
done
