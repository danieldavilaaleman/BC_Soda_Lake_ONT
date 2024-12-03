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

Basecalling output was:
- Simplex reads basecalled: 1,891,493
- Simplex reads filtered: 569,618

### Convert .bam output to .fastq
The output of `dorado basecaller` is *.bam file. To convert to fastq for next steps, I used `bedtools bamtofastq`

```
#!/bin/bash
####### Reserve computing resources #############
#SBATCH --time=24:00:00
#SBATCH --mem=40G
#SBATCH --partition=bigmem
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

####### Set environment variables ###############
module load biobuilds/2017.11
####### Run your script #########################
bedtools bamtofastq -i DL1_SodaLakes_basecalling.bam -fq DL1_SodaLakes_LongReads.fastq
```

### Long-reads QC using chopper
One recommendation from doi: 10.1371/journal.pcbi.1010905 suggest filter reads shorter than 1kb before performing assembly.

```
#! /bin/bash
# ======================================================================
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=0-02:00:00
#SBATCH --partition=cpu2023
# ======================================================================
source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate lr_assemblers

gunzip -c DL1_SodaLakes_LongReads.fastq.gz | chopper -q 10 -l 500 | gzip > Filtered_500_10_DL1_SodaLakes_LongReads.fastq.gz
```



