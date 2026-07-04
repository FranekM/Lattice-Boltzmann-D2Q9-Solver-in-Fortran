#!/bin/bash
#SBATCH --job-name=lbm-fast-scaling-test
#SBATCH --partition=cpu_il
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=256mb
#SBATCH --output=output.%j.out
#SBATCH --error=error.%j.err

echo "ID/NAME:    $SLURM_JOB_ID / $SLURM_JOB_NAME"
echo "USER:       $SLURM_JOB_USER"
echo "PARTITION:  $SLURM_JOB_PARTITION"
echo "TASKS:      $SLURM_TASKS_PER_NODE tasks/node x $SLURM_JOB_NUM_NODES nodes = $SLURM_NTASKS tasks"
echo "NODES:      $SLURM_JOB_NODELIST"

module purge
module load compiler/intel
module load mpi/impi

export FOR_COARRAY_CONFIG_FILE=./lbm.conf

EXE=./fast

echo "Compiling program ..."

make -f Makefile.linux clean
make -f Makefile.linux fast

rc=$?
if [ $rc -ne 0 ]; then
    echo "Compilation failed (exit code $rc)"
    exit $rc
fi

for i in {1..11}; do
    n=$((i*i))

    echo "Running experiment with $n images"

    echo "-genvall -genv I_MPI_FABRICS=shm:ofi -n $n $EXE" > "$FOR_COARRAY_CONFIG_FILE"

    $EXE

    rc=$?
    if [ $rc -ne 0 ]; then
        echo "Experiment with $n images failed (exit code $rc)"
        exit $rc
    fi
done

echo "All experiments completed successfully."