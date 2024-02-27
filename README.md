## ChIP-seq raw data analysis 

## Alignment and Peak Calling

This repository contains workflows and tools for performing alignment and peak calling on paired-end ChIP-seq data. 

The workflows can be run on paired-end (PE) or single-end (SE) next-generation sequcning (NGS) files in `.fastq` format.

The workflows are written to run as batch shell scripts using `SLURM` on a high performance computing (HPC) cluster. 

The first `ChIP_align.sh` (or step1) workflow uses `bowtie2` (v2.3.1), and the second `ChIP_peakcall.sh` (or step2) uses `MACS2` (v2.1.0), based on the GRCh38/hg38 assembly of the human genome:
  
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/

**NOTE:** `ChIP_align.sh` step1 and `ChIP_peakcall.sh` step2 workflows have to be run sequentially, as step2 input files are step1 output files.

Reference genome files in `<.fa>` and `<.fasta>` formats for `bowtie2` and `MACS2` are required for generating the reference indexes. 

Example scripts for generating `bowtie2` and `MACS2` indexes from UCSC references genome files are provided for the hg38 assembly of the human genome.

https://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/

For documentation on a ChIP-seq pipeline from the ENCODE consortium see:

https://github.com/ENCODE-DCC/chip-seq-pipeline2

In addition to `bowtie2` and `MACS2`, the modules in the workflows will require: 

`python`, `java`, `cutadapt`, `bedGraphToBigWig`, `align2rawsignal`, `MCR`, `ucsc_tools`, `picard-tools`, and `samtools` programs. 

## align2rawsignal documentation:

Below are system and installation requirements for `align2rawsignal` taken directly from: https://github.com/akundaje/align2rawsignal?tab=readme-ov-file

align2rawsignal (aka. WIGGLER .. because it generates wiggle files) reads in a set of tagAlign/BAM files, filters out multi-mapping tags and creates a consolidated genome-wide signal coverage file using various tag-shift and smoothing parameters as well as various normalization schemes

The method accounts for the following:

* depth of sequencing
* the mappabilty of the genome (based on read length and ambiguous bases)
* differentiates between positions that shown 0 signal simply because they are unmappable vs positions that are mappable by have no reads. The former are not represented in the output wiggle or bedgraph files while the latter are represented as 0s.
* different tag shifts for the different datasets being combined
* Several types of normalization are implemented. (See usage below)

Wiggler does NOT implement input-DNA/control corrections which can be important for some applications and genomes. 

This tool is primarily used with the following kinds of functional sequencing data
* TF and histone ChIP-seq
* DNase and FAIRE-seq
* MNase-seq for nucleosome positioning

## Compatibility
* Linux 64 bit system with an installation of the free Matlab Compiler Runtime (MCR)
* Supported input formats are tagAlign and BAM files (Single end reads ONLY)
* Supported output formats are bedgraph, wiggle and .mat (matlab)

## Directory structure
```
align2rawsignal/
        /src/                                               --- MATLAB .m source files                  
                *.m
        /bin/                                               --- binary files
                align2rawsignal

        /umap/<version>                                     --- Optional directory containing uniqueness maps
                globalmap_k<mink>tok<maxk>/chr*.unique

        /seq/<version>                                      --- Optional directory containing chromosome sequences
                chr*.fa
```

## Required Files
1. The align2rawsignal executable
2. A directory containing all chromosome sequences in fasta format e.g. align2rawsignal/seq/hg18Female. There should be one file per chromosome in this directory chr*.fa. You can download them from
  * hg18: http://hgdownload.cse.ucsc.edu/goldenPath/hg18/bigZips/chromFa.zip
  * encodeHg19Female: http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/referenceSequences/femaleByChrom/ (*.gz files)
  * encodeHg19Male: http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/referenceSequences/maleByChrom/ (*.gz files)
  * **NOTE:** Remove the _random*.fa and other contig files
  * **NOTE:** If you are working with a female genome, delete or move the chrY.fa file out of the directory pointed to by --seq-dir
3. A directory containing the global uniqueness/mappability map for a range of kmer lengths e.g. umap/encodeHg19Female/globalmap_k20tok54/
  * There should be one file per chromosome in this directory chr*.uint8.unique.
  * NOTE: The mappability file name prefix for each chromosome MUST correspond to the prefix of the sequence file for chromosome e.g. chr1.fa <=> chr1.uint8.unique
  * The mappability tracks can be downloaded from http://www.broadinstitute.org/~anshul/projects/umap/
  * You can unzip these using tar -xvzf globalmap_k`*`tok`*`.tgz
  * The mappability tracks currently support reads upto 54 bp for most organisms. I will be upgrading these to support upto 300 bp reads by Dec 2013.

## Installation Instructions
**NOTE:** These are installation/running instructions for 64-bit LINUX distributions. If you need executables for other platforms please contact Anshul (anshul_at_kundaje_dot_net)

### 1. MCR Installation
In order to run the align2rawsignal code and/or any MATLAB compiled code, you will need the MATLAB runtime library. Please only use the MCR version referenced in this README. This version of the executable was compiled using MCR V7.14 which is equivalent to R2010b release. You can download the MCR here http://mitra.stanford.edu/kundaje/software/MCR2010b.bin

If you haven't installed the MCR, you MUST do that using this command

`./MCR2010b.bin -console`

If you need to specify a specific temp directory then also use the option `-is:tempdir <tempdirname>`

The installer will prompt you to select the directory (<MCR_ROOT>) you want to install the MCR into. e.g. `/home/user/MCR/v714`

**NOTE:** Make sure your installation directory has write permissions and has atleast 500 MB of disk space.

The installation should go smoothly with the above command. However, if you are interested in other installation options you can consult http://www.mathworks.com/access/helpdesk/help/toolbox/compiler/bru23df-1.html

**NOTE:** You need to install the MCR ONLY once on the machine/cluster you plan to run MATLAB compiled code.

### 2. Setting paths
You need to set the following environment variables for the compiled MATLAB code to run correctly. These environment variables MUST be set before calling the align2rawsignal executable or any other MATLAB compiled code.

You can add the following lines to your `.bashrc` or `.cshrc` file if you want to avoid settings these variables everytime you want to run the code

If you are using the bash shell or modifying .bashrc then use

```
MCRROOT=<MCR_ROOT>/v714
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/runtime/glnxa64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64
MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}
XAPPLRESDIR=${MCRROOT}/X11/app-defaults
export LD_LIBRARY_PATH
export XAPPLRESDIR
```

If you are using the csh shell or modifying .cshrc then use
```
setenv MCRROOT <MCR_ROOT>/v714
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${MCRROOT}/runtime/glnxa64
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64 
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64/native_threads
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64/server
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64
setenv XAPPLRESDIR ${MCRROOT}/X11/app-defaults
```

### 3. Samtools
If you are working with BAM input files then the samTools executable is required to be in your system PATH. You can get samtools from here http://samtools.sourceforge.net/

You will need to add samtools to your path so that align2rawsignal can call samtools.

`export PATH="<directory_containing_samtools_binary>:${PATH}"`

If you are running align2rawsignal on a cluster make sure the cluster nodes have this directory in their $PATH. You can set it in your submit script.

### 4. Now you can run the align2rawsignal executable
For options/help simply type

`./align2rawsignal`

OR

`./align2rawsignal --help`

OR
```
export $PATH=<directory_containing_align2rawsignal>/bin:$PATH # add the /bin directory to your path`
align2rawsignal # call align2rawsignal
```

### 5. Running instructions on a cluster

  * Make sure the shell you use in your submit script is bash
  `#!/bin/bash`
  * Make sure all the MCR environment variables have been explicitly set in your submit script (See above.)
  * You can prevent bottlenecking by also defining the following environment variable
    * `export MCR_CACHE_ROOT=$TMPDIR   # If $TMPDIR is defined` OR
    * `export MCR_CACHE_ROOT=<ValidTmpDir> # where <ValidTmpDir> can be for example /tmp or some other temporary directory on the local cluster node`
  * Your submit script must have the path to samtools in your $PATH variable if you are reading in BAM files
  * I also recommend that you set the `$TMP` environment variable to a valid temporary directory with sufficient disk space (incase it isn't already set).
`export TMP=<ValidTmpDir>`

To generate mappability tracks using `gem-mappability` in `MACS2_ref.sh`, you would need to install GEM (Genome Multi-tool) from: 
https://github.com/smarco/gem3-mapper

Please see below for information on other modules used in this workflow.

## Documentation and References:

https://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads

https://hgdownload.soe.ucsc.edu/downloads.html#source_downloads

http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/

https://www.htslib.org/download/

https://broadinstitute.github.io/picard/
  
https://bowtie-bio.sourceforge.net/index.shtml

https://bowtie-bio.sourceforge.net/bowtie2/index.shtml

https://pypi.org/project/MACS2/

