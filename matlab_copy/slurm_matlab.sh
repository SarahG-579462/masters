#!/bin/bash
#SBATCH --time=00-01:00
#SBATCH --mem=8G
#SBATCH --account=rrg-kirshbau-ac
#SBATCH --job-name=matlab_analysis
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=sarah.gammon@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-10%10
./wrapper.sh
