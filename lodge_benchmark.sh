#!/bin/bash
# lodge_benchmark.sh
# Benchmarks LAMMPS MPI/OpenMP combinations for in.polymers_lifetime_eq
# Tests all combinations where np * omp_threads <= 12
# Output: benchmark_results.txt

OUTPUT_FILE="benchmark_results.txt"
A=130
CHAIN_LENGTH=31
N_BEADS=24800
PATCH_SPACING=5
RANDOM_SEED=12345

# MPI x OMP combinations where product <= 12
NP_LIST=(12 6 4 3 2 1)
OMP_LIST=(1 2 3 4 6 12)

echo "=======================================" > "$OUTPUT_FILE"
echo " LAMMPS Benchmark - $(date)" >> "$OUTPUT_FILE"
echo " System: $(hostname)" >> "$OUTPUT_FILE"
echo " CPUs: $(nproc)" >> "$OUTPUT_FILE"
echo " A=${A}, N_beads=${N_BEADS}, chain_length=${CHAIN_LENGTH}" >> "$OUTPUT_FILE"
echo "=======================================" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
printf "%-8s %-12s %-12s %-20s %-20s\n" "NP" "OMP_THREADS" "TOTAL_CORES" "WALL_TIME(s)" "TIME/STEP(ms)" >> "$OUTPUT_FILE"
echo "-----------------------------------------------------------------------" >> "$OUTPUT_FILE"

for i in "${!NP_LIST[@]}"; do
    NP=${NP_LIST[$i]}
    OMP=${OMP_LIST[$i]}
    TOTAL=$((NP * OMP))

    echo "Running: mpirun -np $NP, OMP threads=$OMP (total cores=$TOTAL)..."

    # Unique log file per run to avoid conflicts
    LOG_FILE="logs/benchmark_np${NP}_omp${OMP}.log"
    mkdir -p logs

    # Time the run
    START=$(date +%s%N)

    mpirun -np "$NP" lmp \
        -pk omp "$OMP" \
        -sf omp \
        -in "in.polymers_lifetime_eq" \
        -log "$LOG_FILE" \
        -var chain_length "$CHAIN_LENGTH" \
        -var N_beads "$N_BEADS" \
        -var patch_spacing "$PATCH_SPACING" \
        -var patch_gauss_A_strength "$A" \
        -var random_seed "$RANDOM_SEED" \
        2>/dev/null

    EXIT_CODE=$?
    END=$(date +%s%N)

    if [ $EXIT_CODE -ne 0 ]; then
        printf "%-8s %-12s %-12s %-20s %-20s\n" "$NP" "$OMP" "$TOTAL" "FAILED" "FAILED" >> "$OUTPUT_FILE"
        echo "  -> FAILED (exit code $EXIT_CODE)"
        continue
    fi

    # Wall time in seconds
    WALL_TIME=$(echo "scale=3; ($END - $START) / 1000000000" | bc)

    # Extract total steps and timing from LAMMPS log
    # LAMMPS reports "Loop time of X secs for Y steps"
    LOOP_TIME=$(grep "Loop time" "$LOG_FILE" | tail -1 | awk '{print $4}')
    TOTAL_STEPS=$(grep "Loop time" "$LOG_FILE" | tail -1 | awk '{print $8}')

    if [ -n "$LOOP_TIME" ] && [ -n "$TOTAL_STEPS" ] && [ "$TOTAL_STEPS" -gt 0 ]; then
        # Time per step in milliseconds
        TIME_PER_STEP=$(echo "scale=6; $LOOP_TIME / $TOTAL_STEPS * 1000" | bc)
    else
        TIME_PER_STEP="N/A"
    fi

    printf "%-8s %-12s %-12s %-20s %-20s\n" "$NP" "$OMP" "$TOTAL" "$WALL_TIME" "$TIME_PER_STEP" >> "$OUTPUT_FILE"
    echo "  -> Done: wall=${WALL_TIME}s, time/step=${TIME_PER_STEP}ms"
done

echo "" >> "$OUTPUT_FILE"
echo "=======================================" >> "$OUTPUT_FILE"
echo " Benchmark complete - $(date)" >> "$OUTPUT_FILE"
echo "=======================================" >> "$OUTPUT_FILE"

echo ""
echo "Done! Results saved to $OUTPUT_FILE"
cat "$OUTPUT_FILE"