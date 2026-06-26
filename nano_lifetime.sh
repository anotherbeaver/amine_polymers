A_list=(130)
# A_list=(90 100 110 120 130)
cutoff_list=(0.3 0.4 0.5 0.6 0.7)
cutoff_list=(0.1 0.2 0.3 0.4 0.5 0.6 0.7)

set -e

for A in "${A_list[@]}"; do
    echo "Simulating A=${A}"
    # mpirun -np 6 ../lammps/build/lmp \
    #     -in "in.polymers_lifetime_eq" \
    #     -var chain_length "31" \
    #     -var N_beads "24800" \
    #     -var patch_spacing "5" \
    #     -var patch_gauss_A_strength "$A" \
    #     -var random_seed "12345"

    # wait
    
    mpirun -np 6 ../lammps/build/lmp \
        -in "in.polymers_lifetime_prod" \
        -var chain_length "31" \
        -var N_beads "24800" \
        -var patch_spacing "5" \
        -var patch_gauss_A_strength "$A" \
        -var random_seed "12345"

    wait
    

    for cutoff in "${cutoff_list[@]}"; do
        echo "postprocessing A=${A}, cutoff=${cutoff}"

        ../lammps/tools/binary2txt \
            "data/patches/${cutoff}prodpatch_chainlength_31_patchspacing_5_N_beads_24800_gaussA_${A}_r0_0.25_gaussB_10_seed_12345.bin" &
    done

    wait

    for cutoff in "${cutoff_list[@]}"; do
        echo "removing binary file for A=${A}, cutoff=${cutoff}"
        rm "data/patches/${cutoff}prodpatch_chainlength_31_patchspacing_5_N_beads_24800_gaussA_${A}_r0_0.25_gaussB_10_seed_12345.bin"
    done

    wait

done