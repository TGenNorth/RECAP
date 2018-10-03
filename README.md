# snpDensity
`snpDensityMatrix` is a tool used to analize **SNP pipeline output** and build graphical output for **SNP density**, **read depth**, and **PHI statistic**. 

---
## Getting Started
### Prerequisites:
* [Python](https://www.python.org/downloads/): Version 3.4 or greater.
* [SAMTools](https://github.com/samtools/samtools): Version 1.4.1 or greater.
* [VCFTools](https://github.com/vcftools/vcftools): Version 0.1.15 or greater.

### Basic Installation:
`snpDensityMatrix` primarily runs by calling the `snpDensityMatrix.py` script and can be cloned basically anywhere without the need to compile anything. All of the accompanying **HTML**, **CSS**, and **JavaScript** files are a barebones offline web application that will be copied to the working directory upon use and therefore must be kept in their current locations after cloning the repo.

**The following is a strict step by step guide that must be followed in the proper order:**  
1. Clone this repo to the directory of your choosing.
2. Done.

---
## Usage
### Command Line:
`snpDensityMatrix` is called with the following command line prompt:  
`$ snpDensityMatrix.py [-h, --help] [--window WINDOW] [--step STEP] [--nasp | --cfsan | --lyve | --snvphyl] DIR`

### Arguments:
* `[-h, --help]`  
  - Displays the help menu for snpDensityMatrix before exiting.  
* `[--window WINDOW]`  
  - Optional argument with a default value of **10000**. Set size of capture window to measure snp density.  
* `[--step STEP]`  
  - Optional argument with a default value of **5000**. Set step size to move window down snp matrix.  
* `[--nasp | --cfsan | --lyve | --snvphyl]`  
  - Required argument must be chosen and must match the intended snp pipeline format. Flag used for snp pipeline format.  
* `DIR`  
  - Required Positional argument. Root directory for snp pipeline output.

### Examples:
Running `snpDensityMatrix` to get the snp density visualization for snp pipeline output located at **'./example/snpOut/'**, with a window size of **1000** and a step size of **500**:
```
$ /installDirectory/snpDensityMatrix.py --window 1000 --step 500 --nasp ./example/snpOut/
```  
Would result in the command line output:
```
====================( Job Starting: snpDensityMatrix.py )
(1/4) gathering coverage from BAM/VCF files
    creating file "./snpDensityOut/snpPositions.txt"
    creating file "./snpDensityOut/Sample1-Depth.txt" (1 of 4)
    creating file "./snpDensityOut/Sample2-Depth.txt" (2 of 4)
    creating file "./snpDensityOut/Sample3-Depth.txt" (3 of 4)
    creating file "./snpDensityOut/Sample4-Depth.txt" (4 of 4)
    waiting for threads to return ▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉
(2/4) writing export and excise data files
    creating EXCISE file ...
    creating EXPORT file ▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉
(3/4) processing hash table
    processing 1 of 1 contigs ▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉
    processing 2 of 2 contigs ▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉
(4/4) writing results csv file
    creating file "./snpDensityOut/snpDensityMatrix.csv"
====================( Total Time: 0:0:0:30 )
====================( Job Complete: snpDensityMatrix.py )
```

---
## Technologies
### Visualization:
* [Dygraph](http://dygraphs.com): Used for graphical visualization of **SNP**, **depth**, and **Phi Statistic** data.
* [PhiPack](https://www.maths.otago.ac.nz/~dbryant/software.html): Used to generate **Phi Statistic** from sample alignment data.

### Supported SNP Pipelines:
* [NASP](https://github.com/TGenNorth/NASP): Northern Arizona SNP Pipeline.
* [Cfsan](https://github.com/CFSAN-Biostatistics/snp-pipeline): Center for Food Safety and Applied Nutrition SNP Pipeline.
* [Lyve-SET](https://github.com/lskatz/lyve-SET): Listeria, Yersinia, Vibrio, and Enterobacteriaceae SNP Extraction Tool.
<!-- * [SNVPhyl](https://snvphyl.readthedocs.io/en/latest/): Single Nucleotide Variant Phylogenomics Pipeline -->
