#!/bin/bash
#
#SBATCH --job-name=warpx
#SBATCH --output=WarpX_%j.txt
#SBATCH --error=WarpX_%j.txt
#
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=100
#SBATCH --partition=submit-gpu
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-gpu=1


EXE=/work/submit/$USER/warpx/build_gpu_a30/bin/warpx.3d
INPUTS=inputs

srun --cpu-bind=cores $EXE $INPUTS > output.txt
