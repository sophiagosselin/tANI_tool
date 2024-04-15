#!/usr/bin/env bash
#SBATCH --job-name=auk_tani
#SBATCH --nodes=1
#SBATCH --qos=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48gb
#SBATCH --mail-type=END
#SBATCH --mail-user=sophia.gosselin@uconn.edu
#SBATCH -o tani_%j.out
#SBATCH -e tani_%j.err

echo hostname

#dependencies
module load blast/2.11.0
module load perl

perl tANI_tool.pl -id .7 -cv .7 -bt 100 -t 16
