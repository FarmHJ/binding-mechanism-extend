MODEL='Grandi'
LOG="./log/${MODEL}"
mkdir -p $LOG

python3 scripts/create_param_space.py ${MODEL}

IDS=(`grep -v -e '^#.*' ./simulation_data/parameter_space_exploration/SA_space/${MODEL}/missing_id.txt`)
ID_LEN=${#IDS[@]}
WORKERS=6

echo "## ${MODEL}" >> log/ran_id.log

MAX_REPS=$((${ID_LEN}/${WORKERS} + 1))
for ((reps=0; reps<${MAX_REPS}; reps++)); do
    current=$(date)
    echo "Current time: $current"
    SECONDS=0
    for ((x=0; x<${WORKERS}; x++));
    do
        n=$((reps*${WORKERS} + x))
        echo "${IDS[n]}"
        python3 scripts/SA_sequential.py ${MODEL} ${IDS[n]} > $LOG/${IDS[n]}.log 2>&1 &
        echo "# ${IDS[n]}" >> log/ran_id.log
        echo $! >> log/ran_id.log
    done
    wait
    echo "Elapsed Time: $(($SECONDS/60)) minutes"
done
wait

python3 scripts/combine_results.py ${MODEL}