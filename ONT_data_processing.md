In this section, I present the analysis conducted on ONT sequencing data obtained from High Molecular Weight DNA extracted from microbiome culture 
samples collected from Alberta Soda Lake. The analysis was performed using a R10.4 MinION flow cell.

### Basecalling
The sequencing run was performed using the Ligation Sequencing Kit V14 (SQK-LSK114) and 1 µg of HMW DNA. 73 pod5 files were obtained (~12.5 Gb).
For basecalling, I run the sup model and set a quality limit of 8 using GPU partition `bigmem gpu:1`
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
dorado basecaller --min-qscore 8 sup pod5/ > DL1_SodaLakes_basecalling.bam
```
