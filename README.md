#Documentation for bsExpress v0.5cegx

# Introduction #

`bsExpress` is a program to perform analysis and quality control of bisulfite (BS-Seq) and oxidised bisulfite (oxBS-Seq) sequencing libraries. `bsExpress` is designed to perform quality control of oxBS-Seq libraries using _ad hoc_ control sequences where cytosine modification are known. However, the pipeline is not limited to control sequences but is also suitable for processing (ox)BS-Seq data from raw fastq files to genome-wide methylation calls.

# CEGX customisations from [bsExpress-0.4.1b](https://code.google.com/p/oxbs-sequencing-qc/wiki/bsExpressDoc) #
1. Updated bismark to version 14.2
2. Updated trim_galore to version 0.4.0
3. Modified `bam2methylation.py` to work with SE reads
4. Modified `oxbs_qc.py` with new bismark, trim_galore and SE options
5. Modified `oxbs_qc_func.py` with SE options
6. Modified `validate_args.py` to update `mpile2methylation.py` references to `bam2methylation.py`

# Requirements #

All of the components behind `bsExpress` are freely available. Some of the wrapper scripts and programs ([FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/), [bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/), [MarkDuplicates](http://picard.sourceforge.net/command-line-overview.shtml#MarkDuplicates)) are shipped with `bsExpress`. Others need to be installed by the user. Currently only Linux/Unix systems are suitable, Windows equipped with a Unix emulator like [Cygwin](http://www.cygwin.com/) should be able to run `bsExpress` although it has not been tested.

  * *perl*, *python (2.7)*, *Java 1.6*: Most Linux\Unix systems should have all these components already installed.
  * [R](https://www.r-project.org/)
  * [bedtools](http://bedtools.readthedocs.org/en/latest/index.html)
  * [cutadapt](http://code.google.com/p/cutadapt/) 
  * [samtools](http://samtools.sourceforge.net/) (version 0.1.18+ but older than 1.0)
  * [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 
  * [clipOverlap](http://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap)  _Note_: Only required for paired-end reads

# Installation #

## Via Ubuntu 14.04 tar.gz ##

A tar.gz package is provided in the Downloads section with the cegx-bsexpress code and dependencies that will work on a Linux Ubuntu 14.04 environment.

### Install the dependencies below on a terminal window:

```
sudo apt-get install wget
sudo apt-get install make
sudo apt-get install build-essential
sudo apt-get install zlib1g
sudo apt-get install zlib1g-dev
sudo apt-get install libncurses5-dev
sudo apt-get install unzip
sudo apt-get install git
sudo apt-get install openjdk-6-jre
sudo apt-get install bedtools
sudo apt-get install r-base
sudo apt-get install python-pip
sudo pip install cutadapt==1.11
```

### Install with the instructions below on a terminal window:

```
cd $HOME
cegx_bsexpress_version="cegx-bsexpress-ubuntu-14.04.4"
tar xzf ${cegx_bsexpress_version}.tar.gz
```

```
cd $HOME/${cegx_bsexpress_version}/cegx_bsexpress
sudo python setup.py install --install-scripts $HOME/${cegx_bsexpress_version}/
cd $HOME
```

```
$HOME/${cegx_bsexpress_version}/bsExpress --help
```

### Place your _R1_ fastq.gz reads in this directory:

```
mkdir -p $HOME/in/reads_fastqgz/
cp -v /original/location/of/*_R1_*.fastq.gz $HOME/in/reads_fastqgz/
```

### Download a copy of run_cegx_bsexpress.pl from cegxtools and make a copy of it in your $HOME/ directory:

http://bitbucket.org/cegx-bfx/cegxtools/downloads/

### Now run using the run_cegx_bsexpress.pl script:

```
suffix=.fastq.gz
perl $HOME/run_cegx_bsexpress.pl -d $HOME/in/reads_fastqgz/ -suffix $suffix \
 -bsexpress $HOME/${cegx_bsexpress_version}/bsExpress \
 -SQ_controls $HOME/${cegx_bsexpress_version}/cegx_bsexpress/control_reference/oxBS_controls-v1.0.fa \
 -DC_controls $HOME/${cegx_bsexpress_version}/cegx_bsexpress/digestion_reference/DC_controls-v1.0.fa
```


## Via Docker container ##

(In progress)

## Via BitBucket code download ##
`bsExpress` is a python package installable via [setup.py](http://docs.python.org/2/install/#) . Download `bsExpress` from [BitBucket](https://bitbucket.org/russellshamilton/cegx_bsexpress) source forge]. Assuming the archive is in `~/Downloads`, unpack the archive, move inside it and install:

```
cd ~/Downloads
tar zxvf bsExpress-x.x.x.tar.gz
cd bsExpress-x.x.x
python setup.py install --user --install-scripts ~/bin/
```

The actual installation is done by the last command. In this example, the `--user` option makes `bsExpress` to be installed only for the current user, therefore no administrator rights are required. Remove the `--user` option for a system-wide installation. 

The parameter `--install-scripts` determines where the front-end program `bsExpress` is going to be copied to. In this example it goes to the bin directory in the user's home (`~/bin`). You can choose any directory of your convenience (e.g. one on your PATH or `/usr/local/bin` for system-wide installation)  

# `bsExpress`: Workflow overview #

`bsExpress` is effectively a wrapper around several independent tools and analysis steps. The workflow is designed to be flexible and modular so that several entry points are allowed. Below is the the full pipeline:

| Step | Description                         | Underlying program                               | Option to skip or enable step  
| :--- | :---------------------------------- | :----------------------------------------------- | :----------------------------
| 1    | Raw sequence QC                     | `FastQC`                                         | `--skip_fastqc/-Sf` 
| 2    | Trim adapters and low quality ends  | `cutadapt`/`trim_galore`                         | `--skip_trim/-St` 
| 2.1  | Trim reads from MspI/RRBS libraries | `trim_galore`                                    | `--rrbs` 
| 3    | Shorten reads<sup>1</sup>           | `ShortenFastq.jar`                               | `--skip_shorten/-Ss` 
| 4    | Align reads to reference            | `bowtie2/bismark`                                | `--skip_aln/-Sa` 
| 5    | Convert SAM to BAM, sort and index  | `samtools`                                       | `--skip_sam2bam/-Sb`
| 6    | Mark duplicate reads<sup>2</sup>             | `MarkDuplicates.jar`                             | `--mark_duplicates/-Md` 
| 7    | Clip overlapping PE reads<sup>3</sup>        | `clipOverlap`                                    |  `--skip_clip/-Sc `
| 8    | Call methylation<sup>4</sup>                 | `samtools mpileup &#124; methylation2mpileup.py` | `--skip_mcall/-Sm`
| 9    | Report QC on control sequences<sup>5</sup>   | `oxbs_report.R`                                  | `--skip_report/-Sr` 


*Notes*: <sup>1</sup> Only use this option when reads are longer than the shortest (control) reference minus 2bp. When the reference sequence(s) are considerably longer than the reads disable this step. <sup>2</sup> Skip MarkDuplicates when aligning against the control reference since most if not all of the reads would be marked as duplicate (See note 4 also). <sup>4</sup> Only relevant to paired-end reads. <sup>4</sup> _NB:_ The `mpileup` engine is hard-coded (at least up to version 0.1.18) to skip reads marked as duplicate (see also note 3 and this post on [SeqAnswers](http://seqanswers.com/forums/showthread.php?t=22752 SEQanswers)). <sup>5</sup> Skip this report when not aligning to control reference sequences. 

# Tutorial #

## Control sequences ##

We want to align paired-end reads to the reference control sequences: 

```
bsExpress -i fastq/oxBS_R1.fq.gz fastq/oxBS_R2.fq.gz -r bsExpress-x.x.x/control_reference/oxBS_controls-v1.0.fa -p runqc
```

The `bsExpress` archive contains the directory `control_reference` with thr reference fasta and all the other files to perform the QC. If you move this directory outside the `bsExpress` archove or you want to use a different control change the argument to `-r/-ref` accordingly.

The analysis output will be in directory `oxbs_qc` and output file names will have prefix `runqc`.

## Genome analysis ##

The analysis of "real" genomic reads is virtually the same as for the control sequences. We only need to skip the last reporting step (`--skip_report/-Sr`) and enable the marking of duplicate reads (`--mark_duplicates/-Md`, not necessary but recommended). Remember that reads marked as duplicate are not taken into account for computing methylation counts. 

We also pass to a bed file of the positions of CpG islands. Only cytosines overlapping CpG islands will be produced in output. 

```
bsExpress -i fastq/oxBS_R1.fq.gz fastq/oxBS_R2.fq.gz -r genome/mm9.fa -p oxbs_qc \
    -l cpgi.bed -Md -Sr
```

Again, directory `genome` is expected to contain the indexes for `bowtie2` created by `bismark_genome_preparation`.

If reads have been trimmed in a previous run of the pipeline (e.g. for the control sequences), we can skip the quality and adapter trimming by setting the `-St` option and passing to the `-i` options the trimmed fastq files. Also the FastQC report can be skipped (`-Sf`)
 
```
bsExpress -Sf -St -i oxbs_qc/oxbs_qc.R1.trim.fq.gz oxbs_qc/oxbs_qc.R2.trim.fq.gz ...
```

# Output files #

The type of output can be recognized by the extension and most of the output files will begin with the prefix assigned with the `--prefix` option.

The order if the description below reflects the order in which the `bsExpress` produces them, assuming all the workflow options are enabled.

#### `*fastqc.zip` ####

Basic quality control of the raw reads as produced by [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Open `fastqc_report.html` inside the unzipped folder for visual inspection. 

#### `*.trim.fq.gz` & `*.trimming_report.txt` ####

Respectively the trimmed fastq file(s) and the trimming report produced by `cutadapt/trim_galore`. 

#### `*.short.fq.gz` ####

Fastq reads shortened to the length set by `--maxlen`.

#### `<prefix>.bam/.bam.bai` & `Bismark_mapping_report.txt` ####

`.bam` is the aligned file produced by `bowtie2/bismark` _sorted_ and _indexed_ (index is `<prefix>.bam.bai`). `Bismark_mapping_report.txt` gives overall statistics about alignment and methylation.

#### `*.mdup.bam/.bai` & `*.mdup.metrics.txt` ####

Bam file, sorted and indexed, produced by bismark with duplicate reads marked by `MarkDuplicates`. Duplication metrics are in `*.mdup.metrics.txt`. See also [FAQs on MarkDuplicates](https://sourceforge.net/apps/mediawiki/picard/index.php?title=Main_Page#Q:_How_does_MarkDuplicates_work.3F)

#### `*.clip.bam/.bam.bai` & `*.clipStats` ####

BAM files from previous step with overalapping paired reads clipped to avoid double counting coverage and methylation. Clipping statistics are in `*.clip.Stats`. See also [clipOverlap](http://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap clipOverlap).

#### `.mcall.bdg.gz` ####

Tab delimited file without header produced by Step 8 above. The format is similar to bed or bedgraph. It reports for each cytosine the methylation percentage and count. Columns are:

  * Chromosome or control sequence name
  * Cytosine position 0-based. In bed format this would be the "start" position
  * Cytosine position 1-based. In bed format this would be the "end" position.
  * Percentage unconverted C (% methylated) expressed as 0 to 100. I.e. [column 5 div column 6] x 100
  * Number of unconverted C (count methylated)
  * Total number of converted and unconverted C. Note that bases that are neither C or T (+ strand) or G or A (- strand) are not counted 
  * Strand: + for C on forward strand, - for C on reverse strand

This file will always have positions sorted *within* chromosomes. Currently (version 0.1.4a) it is also sorted by chromosome but chromosome sorting might change in future versions if sorting large files takes too long.

By default, the positions included in this file are *all* the cytosines in the reference file on forward and reverse strand. This means 1) That cytosines with 0 coverage are not included in this file. Therefore the mcall files produced by different input files with the same reference might not align side by side. 2) For a mammalian genome (mouse, human) this file might contain, with default settings, hundreds of millions of rows. To reduce the size consider using the option `--listpos/-l <bed-file>'` where `bed file` is a bed file of regions to scan (e.g. only CpG sites, promoters, CpG islands).

_Example output:_

```
chr1     424     425     0.0     0       54      +
chr1     425     426     1.8519  1       54      +
chr1     427     428     0.0     0       22      -
chr1     428     429     0.0     0       22      -
```

## For control sequences only ##

####  `*.coverage.pdf` & `*.conversion.pdf` ####

Graphical representation of the read coverage and percent unconverted C to T at each cytosine in the control sequences. 

#### `*.oxqc_summary.txt` ####

Summary of conversions summarized by control sequence (chromosome) and modification. This is a file in tabular format with columns: 

| Field     | Description 
| :-------- | :----------------------------------------------------------- |
| chrom     | Control sequence name (_SQ1hmC, SQ3hmC, etc._)               |
| mod       | Cytosine modification (_5mC, 5hmC, 5fC, or C_)               |
| pct.met   | Percentage unconverted cytosines (as _cnt.met / tot_reads_)  |
| cnt.met   | Count of unconverted cytosines                               |
| tot_reads | Total number of reads                                        |

#### `*.oxqc.txt` ####

File similar in format to `.mcall.bdg.gz` with additional annotation of each position in the control sequences:

| Field             | Description 
| :---------------- | :----------------------------------------------------------- 
| chrom             | Control sequence name (_SQ1hmC, SQ3hmC, etc._) 
| pos               | Cytosine position 
| pct.met           | Percentage of unconverted cytosines 
| cnt.met           | Count of unconverted cytosines 
| tot_reads         | Total number of reads 
| strand            | + or - for cytosine on forward strand (+) or reverse (-) strand 
| base_iupac        | C or G 
| short_description | Description of modification: _5mC+, 5mC-, 5hmC+, 5hmC-, 5fC+, 5fC-, C+, C-_ 

# Additional notes #

  * Passing options to underlying programs like `samtools mpileup` is done by passing a string to the appropriate parameter (e.g. `--mpileup_opts`). This string must start with a blank space after the quotes. E.g.:

`--mpileup_opts ' -l cpgi.bed'`  *OK*

`--mpileup_opts '-l cpgi.bed'`   *WRONG*

(This is because python would interpret -l as an additional parameter)

  * The scripts `trim_galore` and `bismark` shipped with `bsExpress` have been edited to replace `zcat` with `gunzip -c` since some MacOS systems do not have `zcat`.