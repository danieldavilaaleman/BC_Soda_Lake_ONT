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

The output generated for dorado
```
[2025-01-27 14:26:20.374] [info]  - downloading dna_r10.4.1_e8.2_400bps_sup@v4.3.0 with curl
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  105M  100  105M    0     0  20.9M      0  0:00:05  0:00:05 --:--:-- 23.5M
[2025-01-27 14:26:27.845] [info] > Creating basecall pipeline
[2025-01-27 14:26:40.375] [info]  - set batch size for cuda:0 to 1664
[2025-01-27 18:44:56.531] [info] > Simplex reads basecalled: 2137384
[2025-01-27 18:44:56.533] [info] > Simplex reads filtered: 331000
[2025-01-27 18:44:56.533] [info] > Basecalled @ Samples/s: 1.487425e+07
[2025-01-27 18:44:57.424] [info] > Finished
```
-  Total of reads basecalled: 2,137,384
-  Total of reads filtered: 331,000

# Convert .bam to fastq
The next step is to convert [dorado](https://github.com/nanoporetech/dorado) .bam file to fastq using [bedtools bamtofastq](https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html)

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
#module load biobuilds/2017.11
####### Run your script #########################
bedtools bamtofastq -i GE22_SodaLakes_basecalling.bam -fq GE22_SodaLakes_LongReads.fastq
```

#### NOTE: module biobuilds/2017.11 were not available anymore after ARC software update. I dowloaded badtools.static.binary [version2.30.0](https://github.com/arq5x/bedtools2/releases) on my ```/bin``` directory 

### ARC Dmitri option:  Another option is to source bioconda - This comes with the update in ARC
```
source /global/software/bioconda/init-2024-10
```

# Long-read QC and length trimming
I used [chopper](https://github.com/wdecoster/chopper) for filtering <10 QC reads and <500 bp
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

cat GE22_SodaLakes_LongReads.fastq | chopper -q 10 -l 500 | gzip > Filtered_500_10_GE22_SodaLakes_LongReads.fastq.gz
```

The output of chopper was:

Kept 3,186,680 reads out of 4,307,968 reads ### NOTE: that this number is different from dorado basecalling Simplex reads basecalled. 
Why? ### Possible answer [here](https://github.com/nanoporetech/dorado/issues/992)

# Long-reads Assembly

Assembly was performed using [metaMDBG](https://github.com/GaetanBenoitDev/metaMDBG) tool 

```
#! /bin/bash
# ======================================================================
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40gb
#SBATCH --time=3-00:00:00
#SBATCH --partition=cpu2023
# ======================================================================
module load gcc/10.2.0 cmake/3.13.4 lib/zlib/1.2.11 openmpi/4.1.1-gnu
metaMDBG asm --out-dir metaMDBG_assembly_GE22 --in-ont Filtered_500_10_GE22_SodaLakes_LongReads.fastq.gz --threads 8
```


















