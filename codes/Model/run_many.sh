#!/bin/bash

# Total number of runs
NUM_RUNS=30

# Get number of available cores
NUM_CORES=30

for i in $(seq 1 $NUM_RUNS); do
    # Compute core number (cycling through available cores)
    CORE=$(( (i - 1) % NUM_CORES ))

    echo "Submitting job $i to core $CORE..."
    ts taskset -c $CORE /opt/matlab/R2023b/bin/matlab \
        -nodisplay -nosplash -nodesktop -singleCompThread \
        -r "parallel_fitting(5, 'pattern_search'); exit;"
done
