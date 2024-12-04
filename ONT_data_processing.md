## Soda Lake microbiome assembly
In this section, I present the analysis conducted on ONT sequencing data obtained from High Molecular Weight DNA extracted from microbiome culture 
samples collected from British Columbia, Soda Lakes. The analysis was performed using a R10.4 MinION flow cell.

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
The output of [dorado basecaller](https://github.com/nanoporetech/dorado) is *.bam file. To convert to fastq for next steps, I used `bedtools bamtofastq`

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
The output of the [chopper](https://github.com/wdecoster/chopper) filter is `Kept 1589641 reads out of 1904350 reads`

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

The output of [metaMDGB](https://github.com/GaetanBenoitDev/metaMDBG) got 11 circular contigs >1MB
- Run time:                   3h 32min 17sec
- Peak memory:                8.11041 GB
- Assembly length:            362362931
- Contigs N50:                172643
- Nb contigs:                 8393
- Nb Contigs (>1Mb):          38
- Nb circular contigs (>1Mb): 11

The next step is polishing using, [MEDAKA](https://github.com/nanoporetech/medaka), [Polypolish](https://github.com/rrwick/Polypolish), and [Pypolca](https://github.com/gbouras13/pypolca).

### Assembly polishing using Long-reads with MEDAKA

The input for [MEDAKA](https://github.com/nanoporetech/medaka) are the filtered Long-reads used for the assembly and the assembly directory. I added the flag `--bacteria` to allow the usage of a research model that improve consensus accuracy to metagenomic samples. **NOTE:** Aseembly file needs to be unzipped. This step is to solve structure erros (misassemblies), and the only errors remained will be single base pair substitutions, deletion or insertions.

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
gzip -d metaMDBG_assembly_DL1/contigs.fasta.gz
medaka_consensus -i Filtered_500_10_DL1_SodaLakes_LongReads.fastq.gz -d metaMDBG_assembly_DL1/contigs.fasta \
-o medaka.DL1.assembly.out -t 6 --bacteria
```

### Assembly polishing using short-reads

After Long-read polishing, I used two different softwares for polishing assembly using short-reads which its goal is to correct for those single bp errors left by the long-read polishing step,for example, long homopolymers tempt to be difficult to correct with Nanopore but not with Illumina sequencing data.

The first step is polishing using [Polypolish](https://github.com/rrwick/Polypolish)

```
#!/bin/bash
####### Reserve computing resources #############
#SBATCH --time=12:00:00
#SBATCH --mem=80G
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

####### Set environment variables ###############
module load biobuilds/2017.11
####### Run your script #########################
bwa index medaka.DL1.assembly.out/consensus.fasta
bwa mem -t 16 -a medaka.DL1.assembly.out/consensus.fasta ../../../RB_6/SR/Li50127-RS-DL-1-RT_S16_R1.fastq.gz > alignments_1.sam
bwa mem -t 16 -a medaka.DL1.assembly.out/consensus.fasta ../../../RB_6/SR/Li50127-RS-DL-1-RT_S16_R2.fastq.gz > alignments_2.sam

###### Polypolish insert size filter ############
polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
polypolish polish medaka.DL1.assembly.out/consensus.fasta filtered_1.sam filtered_2.sam > medaka.polypolish.DL1.assembly.fasta
```

The second step was using [Pypolca](https://github.com/gbouras13/pypolca) on top of the polypolish output.

```
#!/bin/bash
####### Reserve computing resources #############
#SBATCH --time=12:00:00
#SBATCH --mem=250G
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

####### Set environment variables ###############
source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate pypolca
####### Run your script #########################
pypolca run -a medaka.polypolish.DL1.assembly.fasta \
-1 ../../../RB_6/SR/Li50127-RS-DL-1-RT_S16_R1.fastq.gz -2 ../../../RB_6/SR/Li50127-RS-DL-1-RT_S16_R2.fastq.gz \
-t 12 -o medaka.polypolish.polca.DL1.assembly.fasta --careful
```

The report output of pypolca is the following:

Stats BEFORE polishing:
- Substitution Errors Found: 296,847
- Insertion/Deletion Errors Found: 113,225
- Assembly Size: 362,803,214
- Consensus Quality Before Polishing: 99.89
- Consensus QV Before Polishing: 29.47

### Binning and Refinement of the assembly

For binning and refinement, I used [metaWRAP](https://github.com/bxlab/metaWRAP) tool. I will keep all the bins twith completness higher than 50% and contamination lower than 10% with the following command:

```
#!/bin/bash
####### Reserve computing resources #############
#SBATCH --time=96:00:00
#SBATCH --mem=45G
#SBATCH --partition=cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

####### Set environment variables ###############
source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate metawrap-env
###### Run your script #########################
##GunZIp in case are in .gz
gunzip -c ../../../RB_6/SR/Li50127-RS-DL-1-RT_S16_R1.fastq.gz > DL1_SR_R1.fastq
gunzip -c ../../../RB_6/SR/Li50127-RS-DL-1-RT_S16_R2.fastq.gz > DL1_SR_R2.fastq

## BINNING ##
metawrap binning -o Binning_pypolca_DL1 -t 12 -a medaka.polypolish.polca.DL1.assembly.fasta/pypolca_corrected.fasta \
--metabat2 --maxbin2 --concoct -m 40 DL1_SR_R1.fastq DL1_SR_R2.fastq

## BIN REFINEMENT ##
metawrap bin_refinement -o Refinement_pypolca_DL1 -t 12 -A Binning_pypolca_DL1/metabat2_bins/ \
-B Binning_pypolca_DL1/maxbin2_bins/ -C Binning_pypolca_DL1/concoct_bins/ -c 50 -x 10 -m 40
```

metaWRAP Refinement module generates **62 "good" bins** with contaminations score < 10% and completeness score > 50%. Seven bins were classified as Cyanobacteria from [CheckM](https://github.com/Ecogenomics/CheckM) dependency in metaWRAP. To know the number of contigs per bin in the assembled metagenome, I used to following unix command:

```

```











