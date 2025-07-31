#!/bin/bash
NUM_RUNS5=29
NUM_RUNS1=29

#Number of cores
NUM_CORES=$((NUM_RUNS5 + NUM_RUNS))


NUM_RUNS=60

for i in $(seq 1 $NUM_RUNS5); do
    CORE=$((i - 1))
    echo "Submitting job $i (5) to core $CORE..."
    ts taskset -c $CORE /opt/matlab/R2023b/bin/matlab \
       -nodisplay -nosplash -nodesktop -singleCompThread \
       -r "parallel_fitting(5, 'pattern_search', 'last'); exit;"
done

for j in $(seq 1 $NUM_RUNS1); do
    CORE=$((NUM_RUNS5 + j -1)) 
    echo "Submitting job $IDX (1) to core $CORE..."
    ts taskset -c $CORE /opt/matlab/R2023b/bin/matlab \
        -nodisplay -nosplash -nodesktop -singleCompThread \
        -r "parallel_fitting(1, 'pattern_search'); exit;"
done
