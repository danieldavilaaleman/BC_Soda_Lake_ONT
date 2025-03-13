## Soda Lake microbiome assembly
In this section, I present the analysis conducted on ONT sequencing data obtained from High Molecular Weight DNA extracted from microbiome culture 
samples collected from British Columbia, Soda Lakes sample **PL4**. The analysis was performed using a R10.4 MinION flow cell.

### Basecalling
The sequencing run was performed using the Ligation Sequencing Kit V14 (SQK-LSK114) and 1 µg of HMW DNA. 122 pod5 files were obtained.
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
dorado basecaller --min-qscore 8 sup pod5/ > PL4_SodaLakes_basecalling.bam
```

Basecalling output was:
- Simplex reads basecalled: 1,542,398
- Simplex reads filtered: 491,815

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
bedtools bamtofastq -i PL4_SodaLakes_basecalling.bam -fq PL4_SodaLakes_LongReads.fastq
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

gunzip -c PL4_SodaLakes_LongReads.fastq.gz | chopper -q 10 -l 500 | gzip > Filtered_500_10_PL4_SodaLakes_LongReads.fastq.gz
```
The output of the [chopper](https://github.com/wdecoster/chopper) filter is `Kept 2,634,870 reads out of 3,101,656 reads`

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
metaMDBG asm --out-dir metaMDBG_assembly_PL4 --in-ont Filtered_500_10_PL4_SodaLakes_LongReads.fastq.gz --threads 8
```

The output of [metaMDGB](https://github.com/GaetanBenoitDev/metaMDBG) got 11 circular contigs >1MB
|Run time:| 22h 39min 5sec |
|----------|-------------|
|Peak memory: | 10.1851 GB |
|Assembly length: | 300,003,546 |
|Contigs N50: | 95,734 |
|Nb contigs: | 8331 |
|Nb Contigs (>1Mb): | 34 |
|Nb circular contigs (>1Mb): | 2|

The next step is polishing using, [MEDAKA](https://github.com/nanoporetech/medaka), [Polypolish](https://github.com/rrwick/Polypolish), and [Pypolca](https://github.com/gbouras13/pypolca).

### Assembly polishing using Long-reads with MEDAKA

The input for [MEDAKA](https://github.com/nanoporetech/medaka) are the filtered Long-reads used for the assembly and the assembly directory. I added the flag `--bacteria` to allow the usage of a research model that improve consensus accuracy to metagenomic samples.  This step is to solve structure erros (misassemblies), and the only errors remained will be single base pair substitutions, deletion or insertions.

**NOTE:** Aseembly file needs to be unzipped before MEDAKA polishing!
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
gzip -d metaMDBG_assembly_PL4/contigs.fasta.gz
medaka_consensus -i Filtered_500_10_PL4_SodaLakes_LongReads.fastq.gz -d metaMDBG_assembly_PL4/contigs.fasta \
-o medaka.PL4.assembly.out -t 6 --bacteria
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
bwa index medaka.PL4.assembly.out/consensus.fasta
bwa mem -t 16 -a medaka.PL4.assembly.out/consensus.fasta ../../../RS/PL4/SR/Li49152-RS-PL4-30C_S2_R1.fastq.gz > alignments_1.sam
bwa mem -t 16 -a medaka.PL4.assembly.out/consensus.fasta ../../../RS/PL4/SR/Li49152-RS-PL4-30C_S2_R2.fastq.gz > alignments_2.sam

###### Polypolish insert size filter ############
polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
polypolish polish medaka.PL4.assembly.out/consensus.fasta filtered_1.sam filtered_2.sam > medaka.polypolish.PL4.assembly.fasta
```

The output of **polypolish**

Finished! (2025-03-12 15:01:15)    
Alignments before filtering: 420,424,874    
Alignments after filtering:  350,927,764    

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
pypolca run -a medaka.polypolish.PL4.assembly.fasta \
-1 ../../../RS/PL4/SR/Li49152-RS-PL4-30C_S2_R1.fastq.gz -2 ../../../RS/PL4/SR/Li49152-RS-PL4-30C_S2_R2.fastq.gz \
-t 12 -o medaka.polypolish.polca.PL4.assembly.fasta --careful
```

The report output of pypolca is the following:

|Stats BEFORE polishing| Value |
|--------|--------|
|Substitution Errors Found:| 182,504|
|Insertion/Deletion Errors Found:| 81,511|
|Assembly Size:| 300,315,466 |
|Consensus Quality Before Polishing:| 99.91|
|Consensus QV Before Polishing:| 30.56|

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
gunzip -c ../../../RS/PL4/SR/Li49152-RS-PL4-30C_S2_R1.fastq.gz > PL4_SR_R1.fastq
gunzip -c ../../../RS/PL4/SR/Li49152-RS-PL4-30C_S2_R2.fastq.gz> PL4_SR_R2.fastq

## BINNING ##
metawrap binning -o Binning_pypolca_PL4 -t 12 -a medaka.polypolish.polca.PL4.assembly.fasta/pypolca_corrected.fasta \
--metabat2 --maxbin2 --concoct -m 40 PL4_SR_R1.fastq PL4_SR_R2.fastq

## BIN REFINEMENT ##
metawrap bin_refinement -o Refinement_pypolca_PL4 -t 12 -A Binning_pypolca_PL4/metabat2_bins/ \
-B Binning_pypolca_PL4/maxbin2_bins/ -C Binning_pypolca_PL4/concoct_bins/ -c 50 -x 10 -m 40
```

metaWRAP Refinement module generates **62 "good" bins** with contaminations score < 10% and completeness score > 50%. Seven bins were classified as Cyanobacteria from [CheckM](https://github.com/Ecogenomics/CheckM) dependency in metaWRAP. To know the number of contigs per bin in the assembled metagenome, I used to following unix command:

```
for bin in *.fa;do grep ">" $bin | wc -l;done
```

This shows that eight bins contains one single circular contig.

### Taxonomic classification

For taxonomic classification, I used [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) v2.4.0 which the most recent reference data version r220 with the following command

```
#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=08:00:00
#SBATCH --mem=100G
#SBATCH --partition=bigmem

####### Set environment variables ###############
source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate gtdbtk
####### Run your script #########################
gtdbtk classify_wf --genome_dir Refinement_pypolca_PL4/metawrap_50_10_bins/ \
--out_dir gtdbtk_classify_wf --cpus 8 -x fa --skip_ani_screen
```




