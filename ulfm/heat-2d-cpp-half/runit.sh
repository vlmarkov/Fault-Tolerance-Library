#!/bin/bash

./makeit.sh && cd build && sbatch slurm.job && sbatch slurm-kill.job
