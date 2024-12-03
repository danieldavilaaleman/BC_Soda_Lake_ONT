## Soda Lake microbiome assembly
In this section, I present the analysis conducted on ONT sequencing data obtained from High Molecular Weight DNA extracted from microbiome culture 
samples collected from Alberta Soda Lake. The analysis was performed using a R10.4 MinION flow cell.

### Basecalling
The sequencing run was performed using the Ligation Sequencing Kit V14 (SQK-LSK114) and 1 µg of HMW DNA. 73 pod5 files were obtained (~12.5 Gb) an N50 aprrox 7.54Kb.
For basecalling, I run [dorado basecaller](https://github.com/nanoporetech/dorado) with the `sup` model and set a quality limit of 8 using GPU partition `bigmem gpu:1`
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
Following the recomendation on ["Assembling the perfect bacterial genome using Oxford Nanopore and Illumina sequencing"](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010905) suggest filter Nanopore shorter reads than <1kb and low quality score reads (<10) before performing assembly. For this purpose, I used [chopper](https://github.com/wdecoster/chopper).

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
The output of the [chopper](https://github.com/wdecoster/chopper)filter is `Kept 1589641 reads out of 1904350 reads`

### Long-read Assembly
The assembly was performed using [metaMDBG](https://github.com/GaetanBenoitDev/metaMDBG) software from [Benoit, G., Raguideau, S., James, R. et al. High-quality metagenome assembly from long accurate reads with metaMDBG. Nat Biotechnol 42, 1378–1383 (2024). https://doi.org/10.1038/s41587-023-01983-6](https://www.nature.com/articles/s41587-023-01983-6#Abs1) due to its version 1.0 can handle R10.4+ Nanopore data.

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
metaMDBG asm --out-dir metaMDBG_assembly_DL1 --in-ont Filtered_500_10_DL1_SodaLakes_LongReads.fastq.gz --threads 8
```

The output of [metaMDGB](https://github.com/GaetanBenoitDev/metaMDBG)got 11 circular contigs >1MB
        Run time:                   3h 32min 17sec
        Peak memory:                8.11041 GB
        Assembly length:            362362931
        Contigs N50:                172643
        Nb contigs:                 8393
        Nb Contigs (>1Mb):          38
        Nb circular contigs (>1Mb): 11

The next step is polishing using, [MEDAKA](https://github.com/nanoporetech/medaka), [Polypolish](https://github.com/rrwick/Polypolish), and [Pypolca](https://github.com/gbouras13/pypolca).

### Assembly polishing using Long-reads with MEDAKA

The input for [MEDAKA](https://github.com/nanoporetech/medaka) are the filtered Long-reads used for the assembly and the assembly directory. I added the flag `--bacteria` to allow the usage of a research model that improve consensus accuracy to metagenomic samples.

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
module load python/3.10.4
####### Run your script #########################
medaka_consensus -i Filtered_500_10_DL1_SodaLakes_LongReads.fastq.gz -d metaMDBG_assembly_DL1/contigs.fasta.gz \
-o medaka.DL1.assembly.out -t 6 --bacteria
```



