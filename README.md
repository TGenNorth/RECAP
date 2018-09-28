# snpDensity
## Usage:
 The `snpDensityMatrix` script is used to read files from a snp pipeline output directory and output an interactive web page based visualization of snp density, depth of reads, and possible recombination represented with a phi statistic. This script takes the optional arguments window size, and step size, an argument for the snp pipeline format and a positional argument for the root directory of the snp pipeline output. The **window** size is the range in which a snp density will be captured per window, the step **size** is the length of spaces the window will be moved forward for the next density capture window, and the **DIR** is the location of the snp pipeline output that will be measured.

### Command Line:
`snpDensityMatrix` is called with the following command line prompt:
```
$ snpDensityMatrix.py [-h, --help] [--window WINDOW] [--step STEP] [--nasp | --cfsan | --lyve | --snvphyl] DIR
```

### Arguments:
The following is a list of arguments for *snpDensityMatrix*:
`[-h, --help]` Displays the help menu for snpDensityMatrix before exiting.
`[--window WINDOW]` Set size of capture window to measure snp density.
`[--step STEP]` Set step size to move window down snp matrix.
`[--nasp | --cfsan | --lyve | --snvphyl]` Flag used for snp pipeline format.
`DIR`Root directory for snp pipeline output.

### Example:
Running `snpDensityMatrix` to get the snp density visualization for snp pipeline output located at **'./example/snpOut/'**, with a window size of **1000** and a step size of **500**:
```
$ snpDensityMatrix.py --window 1000 --step 500 --nasp ./example/snpOut/
```
Would result in the command line output:</p>
```
====================( Job Starting: snpDensityMatrix.py )
(1/3) gathering coverage from BAM/VCF files
    creating file "./snpDensityOut/snpPositions.txt"
    creating file "./snpDensityOut/Sample1-Depth.txt" (1 of 10)
    creating file "./snpDensityOut/Sample2-Depth.txt" (2 of 10)
    creating file "./snpDensityOut/Sample3-Depth.txt" (3 of 10)
    creating file "./snpDensityOut/Sample4-Depth.txt" (4 of 10)
    creating file "./snpDensityOut/Sample5-Depth.txt" (5 of 10)
    creating file "./snpDensityOut/Sample6-Depth.txt" (6 of 10)
    creating file "./snpDensityOut/Sample7-Depth.txt" (7 of 10)
    creating file "./snpDensityOut/Sample8-Depth.txt" (8 of 10)
    creating file "./snpDensityOut/Sample9-Depth.txt" (9 of 10)
    creating file "./snpDensityOut/Sample10-Depth.txt" (10 of 10)
(2/3) processing hash table
    processing 1 of 1 contigs |▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉| (100%)
(3/3) writing results csv file
    creating file "./snpDensityOut/snpDensityMatrix.csv"
====================( Total Time: 0:0:1:11 )
====================( Job Complete: snpDensityMatrix.py )</div>
```
