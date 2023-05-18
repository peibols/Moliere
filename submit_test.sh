#!/usr/bin/env bash
#SBATCH --job-name=table_test
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16096M
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --account=def-jeon

./table_test
