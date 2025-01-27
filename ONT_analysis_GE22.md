# ONT data analysys for GE22
I followed the previous github DL1 analysis for basecalling. By default, dorado basecaller will attempt to detect any adapter or 
primer sequences at the beginning and ending of reads, and remove them from the output sequence.
```
#!/bin/bash
####### Reserve computing resources #############
#SBATCH --time=24:00:00
#SBATCH --mem=80G
#SBATCH --partition=bigmem
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

####### Set environment variables ###############
module load python/3.10.4
####### Run your script #########################
dorado basecaller --min-qscore 8 sup pod5/ > GE22_SodaLakes_basecalling.bam
```
