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

#### ARC Dmitri option:  Another option is to source bioconda - This comes with the update in ARC
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

Kept 3,186,680 reads out of 4,307,968 reads 
#### NOTE: The input read number is different from dorado basecalling Simplex reads basecalled output (4,307,968 vs 2,137,384). Why? 
Possible answer [here](https://github.com/nanoporetech/dorado/issues/992) 

Solution output: 33938 - Does this make sense?

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

[metaMDBG](https://github.com/GaetanBenoitDev/metaMDBG) tool output:

- Total read bps:  28,579,858,182
- N50 read length: 12,846
  
| Run time: | 23h 44min 16sec|
|----------- | -------------- |
|Peak memory: |9.94197 GB |
|Assembly length:|449,096,913 |
|Contigs N50:|30,868 |
|Nb contigs:|19,080 |
|Nb Contigs (>1Mb):|16 |
|Nb circular contigs (>1Mb):| 6 |

The next step is polishing using, [MEDAKA](https://github.com/nanoporetech/medaka), [Polypolish](https://github.com/rrwick/Polypolish), and [Pypolca](https://github.com/gbouras13/pypolca).

### Assembly polishing using Long-reads with MEDAKA

The input for [MEDAKA](https://github.com/nanoporetech/medaka) are the filtered Long-reads used for the assembly and the assembly directory. I added the flag `--bacteria` to allow the usage of a research model that improve consensus accuracy to metagenomic samples. This step is to solve structure erros (misassemblies), and the only errors remained will be single base pair substitutions, deletion or insertions.  

**NOTE:** Aseembly file needs to be unzipped before MEDAKA polishing.

```
#!/bin/bash
####### Reserve computing resources #############
#SBATCH --time=24:00:00
#SBATCH --mem=40G
#SBATCH --partition=bigmem
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1

####### Set environment variables ###############
module load python/3.10.4
####### Run your script #########################
gzip -d metaMDBG_assembly_GE22/contigs.fasta.gz
medaka_consensus -i Filtered_500_10_GE22_SodaLakes_LongReads.fastq.gz -d metaMDBG_assembly_GE22/contigs.fasta \
-o medaka.GE22.assembly.out -t 10 --bacteria
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
source /global/software/bioconda/init-2024-10  #Notice the use of bioconda global software

####### Run your script #########################
bwa index medaka.GE22.assembly.out/consensus.fasta
bwa mem -t 16 -a medaka.GE22.assembly.out/consensus.fasta /work/ebg_lab/gm/RS/RB_3/SR/Li49155-RS-GE2022-RT_S5_R1.fastq.gz > alignments_1.sam
bwa mem -t 16 -a medaka.GE22.assembly.out/consensus.fasta /work/ebg_lab/gm/RS/RB_3/SR/Li49155-RS-GE2022-RT_S5_R2.fastq.gz > alignments_2.sam

###### Polypolish insert size filter ############
polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
polypolish polish medaka.GE22.assembly.out/consensus.fasta filtered_1.sam filtered_2.sam > medaka.polypolish.GE22.assembly.fasta
```

The output of Polypolish filter:

```
Finished! (2025-02-03 12:03:16)
Alignments before filtering: 433,159,301
Alignments after filtering:  368,148,351

Time to run: 1:13:09.315608
```

The output of Polypolish polish was too long and only got the running time:

  
  ```Time to run: 0:49:25.735863```

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
pypolca run -a medaka.polypolish.GE22.assembly.fasta \
-1 /work/ebg_lab/gm/RS/RB_3/SR/Li49155-RS-GE2022-RT_S5_R1.fastq.gz -2 /work/ebg_lab/gm/RS/RB_3/SR/Li49155-RS-GE2022-RT_S5_R2.fastq.gz \
-t 12 -o medaka.polypolish.polca.GE22.assembly.fasta --careful
```

The report output of pypolca is the following:

|Stats BEFORE polishing| Value |
|--------|--------|
|Substitution Errors Found:| 1,315,251|
|Insertion/Deletion Errors Found:| 273,951 |
|Assembly Size:| 449,890,817|
|Consensus Quality Before Polishing:| 99.65|
|Consensus QV Before Polishing:| 24.52|

# Binning
## metaWRAP

