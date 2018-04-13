#!/usr/bin/env python3.4
#==================================================================================================== ( description )
# snpPositions.txt format:
#     contig1 position1
#     contig1 position2
#     contig2 position1
#
# sample1-Depth.txt
#     contig1 position1 depth
#     contig1 position2 depth
#     contig2 position1 depth
#
# hash table format:
#     {contig1: {position1: [['A', 'A', 'A', 'A', 'G'], ['42', '42', '42', '42']], position2: [['T', 'T', 'T', 'T', 'C'], ['42', '42', '42', '42']]},
#      contig2: {position1: [['A', 'A', 'A', 'A', 'G'], ['42', '42', '42', '42']], position2: [['T', 'T', 'T', 'T', 'C'], ['42', '42', '42', '42']]}}
#
# example rolling window:
#     contig,from,to,phi,aggregate,sample1,sample2,sample3,sample4
#     contig1,1,1000,--,0,0,0,0,1
#     contig1,1001,2001,--,0,0,0,0,1
#
#==================================================================================================== ( import )
### import
# for parsing options and arguments
import sys
import argparse
import math

# for python string substitution
import re

# for getting time taken for jobs
import time

# for running terminal commands
import subprocess
import os

#==================================================================================================== ( classes )
#============================================================( getArgs )
### parse and return all arguments to NAMESPACE
# all arguments are collected from command line
#
# 'renameSamples -h' should run in any event that the arguments are invalid

class ArgParse:
  def __init__(self):

    # initialize argument parser
    self.parser=argparse.ArgumentParser(usage='''%(prog)s [-h, --help] [ --nasp | --cfsan | --lyve | --snvphyl ] [ --window ] [ --step ] FILE''', formatter_class=argparse.RawDescriptionHelpFormatter, description='''''',epilog='''''')

    # set input file type
    self.input = self.parser.add_mutually_exclusive_group(required=True)
    self.input.add_argument('--nasp', action='store_true', help='Northern Arizona SNP Pipeline')
    self.input.add_argument('--cfsan', action='store_true', help='CFSAN SNP Pipeline')
    self.input.add_argument('--lyve', action='store_true', help='LYVE SNP Extraction Tool')
    self.input.add_argument('--snvphyl', action='store_true', help='SNV Phylogenomics Pipeline')

    # set window size and step size
    self.parser.add_argument('--window', default=1000, type=int, help='set size of capture window to measure snp density')
    self.parser.add_argument('--step', default=500, type=int, help='set step size to move window down snp matrix')
    self.parser.add_argument('DIR', nargs='*', help='SNP pipeline output directory')
    self.args = self.parser.parse_args()
    self.args.window = abs(self.args.window)-1
    self.args.step = abs(self.args.step)

  def getArgs(self):
    return self.args

  def getHelp(self):
    self.parser.print_help()
    exit()

#============================================================( jobTimer )
### allows for tracking the length of time for running a certain job
# sets object creation time, returns difference between creation and current time
#
# for job length of 12 minutes and 24 seconds
# self.getTime() --> 0:0:12:24

class JobTimer:
  def __init__(self):

    # set start time and initializes all time partitions
    self.startTime = time.time()
    self.days = 0
    self.hours = 0
    self.minutes = 0
    self.seconds = 0

  def updateTime(self):

    # updates all time partitions based on difference between start and current time
    totalTime = time.time()-self.startTime
    self.days = int(totalTime / 86400)
    totalTime -= self.days*86400
    self.hours = int(totalTime / 3600)
    totalTime -= self.hours*3600
    self.minutes = int(totalTime / 60)
    totalTime -= self.minutes*60
    self.seconds = int(totalTime)

  def getTime(self):

    # update and print time
    self.updateTime()
    print("====================( Total Time: " + str(self.days) + ":" + str(self.hours) + ":" + str(self.minutes) + ":" + str(self.seconds) + " )")

#==================================================================================================== ( helper functions )
#============================================================( readFile )
def readFile(file):

  inFile = []

  # if target file exists, cat and return file, else exit
  if os.path.isfile(file):
    f = open(file, 'r')
    file = f.read()
    f.close()
  else:
    print("====================( \'" + str(file) + "\' does not exist )")
    exit(1)

  return file

#============================================================( naspCoverage )
def naspHash(DIR):

  # open "matrices/bestsnp.tsv"
  file = readFile(DIR + "matrices/bestsnp.tsv")
  file = file.split("\n")[:-1]

  # get sample names
  samples = file[0].split("#SNPcall")[0].split("\t")[1:-1]
  for i in range(0, len(samples)):
    samples[i] = samples[i].split("::")[0]

  # write samples and positions to file
  print("    creating file \"./snpDensityOut/snpPositions.txt\"") 
  f = open("snpDensityOut/snpPositions.txt", 'w')
  snpPositions = subprocess.Popen("awk \'NR > 1 {print $1}\' " + DIR + "matrices/bestsnp.tsv", universal_newlines=True, shell=True, stdout=subprocess.PIPE)
  snpPositions = snpPositions.stdout.read().split("\n")[:-1]
  for i in range(0, len(snpPositions)):
    snpPositions[i] = "::".join(snpPositions[i].split("::")[:-1]) + '\t' + snpPositions[i].split("::")[-1]
    f.write(snpPositions[i] + "\n")
  f.close()

  # build position hash table for each contig
  hash = {}
  for line in snpPositions:
    line = line.split("\t")
    if line[0] not in hash:
      hash[line[0]] = {}
    if int(line[1]) not in hash[line[0]]:
      hash[line[0]][int(line[1])] = [[],[]]

  # write SNP data to hash table
  for i in range(1, len(file)):
    line = file[i].split("\t")
    contig = "::".join(line[0].split("::")[:-1])
    position = int(line[0].split("::")[-1])
    for j in range(1, len(samples)+1):
      hash[contig][position][0].append(line[j])

  # write coverage to file
  for i in range(1, len(samples)):
    if not os.path.exists("snpDensityOut/" + "".join(samples[i].split(" ")) + "-Depth.txt"):
      print("    creating file \"./snpDensityOut/" + "".join(samples[i].split(" ")) + "-Depth.txt\" (" + str(i) + " of " + str(len(samples) - 1) + ")")
      coverage = subprocess.Popen("module load samtools; samtools depth -b snpDensityOut/snpPositions.txt " + DIR + "bwamem/" + samples[i] + "-bwamem.bam", universal_newlines=True, shell=True, stdout=subprocess.PIPE)
      coverage = coverage.stdout.read()
      f = open("snpDensityOut/" + "".join(samples[i].split(" ")) + "-Depth.txt", 'w')
      f.write(coverage)
      f.close()
    else:
      coverage = readFile("snpDensityOut/" + "".join(samples[i].split(" ")) + "-Depth.txt")
    coverage = coverage.split("\n")[:-1]

    # assemble hash
    for line in coverage:
      line = line.split("\t")
      hash[line[0]][int(line[1])][1].append(int(line[2]))

  return [samples, hash]

#============================================================( cfsanHash )
def cfsanHash(DIR):

  # get sample names
  samples = sorted(readFile(DIR + "sampleDirectories.txt").split("\n")[:-1])
  for i in range(0, len(samples)):
    samples[i] = samples[i][2:]
  samples.insert(0, "Reference")

  # write samples and positions to file
  print("    creating file \"./snpDensityOut/snpPositions.txt\"") 
  snpPositions = subprocess.Popen("awk \'{print $1 \"\t\" $2}\' " + DIR + "snplist.txt", universal_newlines=True, shell=True, stdout=subprocess.PIPE)
  snpPositions = snpPositions.stdout.read().split("\n")[:-1]
  f = open("snpDensityOut/snpPositions.txt", 'w')
  for i in range(0, len(snpPositions)):
    f.write(snpPositions[i] + "\n")
  f.close()

  # build position hash table for each contig
  hash = {}
  for line in snpPositions:
    line = line.split("\t")
    if line[0] not in hash:
      hash[line[0]] = {}
    if int(line[1]) not in hash[line[0]]:
      hash[line[0]][int(line[1])] = [[],[]]

  # write reference SNP bases to hash table
  ref = readFile(DIR + "referenceSNP.fasta").split(">")[1:]
  for i in range (0,len(ref)):
    contig = ref[i].split("\n")[0]
    positions = sorted(list(hash[contig].keys()))
    bases = "".join(ref[i].split("\n")[1:])
    for j in range(0, len(positions)):
      hash[contig][positions[j]][0].append(bases[j])

  # for each sample
  for i in range(1, len(samples)):

    # get SNP data
    vcf = readFile(DIR + samples[i] + "/consensus.vcf")
    vcf = vcf.split(samples[i])[-1].split("\n")[1:-1]
    SNPs = {}
    for line in vcf:
      line = line.split("\t")
      if line[0] not in SNPs:
        SNPs[line[0]] = {}
      if int(line[1]) not in SNPs[line[0]]:
        if line[4] != ".":
          SNPs[line[0]][int(line[1])] = line[4]
        else:
          SNPs[line[0]][int(line[1])] = line[3]

    # write coverage to file
    if not os.path.exists("snpDensityOut/" + "".join(samples[i].split(" ")) + "-Depth.txt"):
      print("    creating file \"./snpDensityOut/" + "".join(samples[i].split(" ")) + "-Depth.txt\" (" + str(i) + " of " + str(len(samples) - 1) + ")")
      coverage = {}
      read = subprocess.Popen("module load samtools; samtools depth -b snpDensityOut/snpPositions.txt " + DIR + "/" + samples[i] + "/reads.sorted.bam", universal_newlines=True, shell=True, stdout=subprocess.PIPE)
      read = read.stdout.read().split("\n")[:-1]
      for line in read:
        line = line.split("\t")
        if line[0] not in coverage:
          coverage[line[0]] = {}
        if int(line[1]) not in coverage[line[0]]:
          coverage[line[0]][int(line[1])] = 0
        coverage[line[0]][int(line[1])] += int(line[2])
      coverageOut = ""
      for contig in sorted(list((coverage.keys()))):
        for position in sorted(list((coverage[contig].keys()))):
          coverageOut += "\t".join([contig,str(position),str(coverage[contig][position])]) + "\n"
      f = open("snpDensityOut/" + "".join(samples[i].split(" ")) + "-Depth.txt", 'w')
      f.write(coverageOut)
      f.close()
    else:

      # read coverage from file and build coverage hash table
      read = readFile("snpDensityOut/" + "".join(samples[i].split(" ")) + "-Depth.txt")
      read = read.split("\n")[:-1]
      coverage = {}
      for line in read:
        line = line.split("\t")
        if line[0] not in coverage:
          coverage[line[0]] = {}
        if int(line[1]) not in coverage[line[0]]:
          coverage[line[0]][int(line[1])] = 0
        coverage[line[0]][int(line[1])] = int(line[2])

    # add SNP data to hash table
    for contig in sorted(list(hash.keys())):
      for position in sorted(list(hash[contig].keys())):
        if int(position) not in SNPs[contig].keys():
          hash[contig][position][0].append("-")
        else:
          hash[contig][position][0].append(SNPs[contig][int(position)])
        if int(position) not in coverage[contig].keys():
          hash[contig][position][1].append(0)
        else:
          hash[contig][position][1].append(coverage[contig][int(position)])

  return [samples, hash]

#============================================================( lyveHash )
def lyveHash(DIR):

  # open "matrices/bestsnp.tsv"
  file = readFile(DIR + "msa/out.snpmatrix.tsv")
  file = file.split("\n")[:-1]

  # get sample names
  samples = file[0].split("REF")[-1].split("\t")[1:-1:2]
  for i in range(0, len(samples)):
    samples[i] = "]".join(samples[i].split("]")[1:])
    samples[i] = samples[i].split("_R1_")[0]
  samples.insert(0, "Reference")

  # write samples and positions to file
  print("    creating file \"./snpDensityOut/snpPositions.txt\"") 
  f = open("snpDensityOut/snpPositions.txt", 'w')
  snpPositions = subprocess.Popen("awk \'NR > 1 {print $1 \"\t\" $2}\' " + DIR + "msa/out.snpmatrix.tsv", universal_newlines=True, shell=True, stdout=subprocess.PIPE)
  snpPositions = snpPositions.stdout.read().split("\n")[:-1]
  for i in range(0, len(snpPositions)):
    f.write(snpPositions[i] + "\n")
  f.close()

  # build position hash table for each contig
  hash = {}
  for line in snpPositions:
    line = line.split("\t")
    if line[0] not in hash:
      hash[line[0]] = {}
    if int(line[1]) not in hash[line[0]]:
      hash[line[0]][int(line[1])] = [[],[]]

  # write SNP data to hash table
  for i in range(1, len(file)):
    line = file[i].split("\t")
    contig = str(line[0])
    position = int(line[1])
    hash[contig][position][0].append(line[2])
    for j in range(3, (len(samples))*2, 2):
      if line[j] == line[j+1] and line[j] == ".":
        hash[contig][position][0].append("-")
      elif line[j] == ".":
        hash[contig][position][0].append(line[j+1])
      else:
        hash[contig][position][0].append(line[j])

  # write coverage to file
  for i in range(1, len(samples)):
    if not os.path.exists("snpDensityOut/" + "".join(samples[i].split(" ")) + "-Depth.txt"):
      coverage = {}
      print("    creating file \"./snpDensityOut/" + "".join(samples[i].split(" ")) + "-Depth.txt\" (" + str(i) + " of " + str(len(samples) - 1) + ")")

      # get coverage from both forward and backwards bam files
      bam = subprocess.Popen("ls " + DIR + "bam/" + samples[i] + "*.bam", universal_newlines=True, shell=True, stdout=subprocess.PIPE)
      bam = bam.stdout.read().split("\n")[:-1]
      coverage1 = subprocess.Popen("module load samtools; samtools depth -b snpDensityOut/snpPositions.txt " + bam[0], universal_newlines=True, shell=True, stdout=subprocess.PIPE)
      coverage1 = coverage1.stdout.read().split("\n")[:-1]
      coverage2 = subprocess.Popen("module load samtools; samtools depth -b snpDensityOut/snpPositions.txt " + bam[1], universal_newlines=True, shell=True, stdout=subprocess.PIPE)
      coverage2 = coverage2.stdout.read().split("\n")[:-1]

      # build coverage hash table and save coverage to file
      for read in [coverage1, coverage2]:
        for line in read:
          line = line.split("\t")
          if line[0] not in coverage:
            coverage[line[0]] = {}
          if int(line[1]) not in coverage[line[0]]:
            coverage[line[0]][int(line[1])] = 0
          coverage[line[0]][int(line[1])] += int(line[2])
      coverageOut = ""
      for contig in sorted(list((coverage.keys()))):
        for position in sorted(list((coverage[contig].keys()))):
          coverageOut += "\t".join([contig,str(position),str(coverage[contig][position])]) + "\n"
      f = open("snpDensityOut/" + "".join(samples[i].split(" ")) + "-Depth.txt", 'w')
      f.write(coverageOut)
      f.close()
    else:

      # read coverage from file and build coverage hash table
      read = readFile("snpDensityOut/" + "".join(samples[i].split(" ")) + "-Depth.txt")
      read = read.split("\n")[:-1]
      coverage = {}
      for line in read:
        line = line.split("\t")
        if line[0] not in coverage:
          coverage[line[0]] = {}
        if int(line[1]) not in coverage[line[0]]:
          coverage[line[0]][int(line[1])] = 0
        coverage[line[0]][int(line[1])] = int(line[2])

    # assemble output hash table
    for contig in sorted(list((hash.keys()))):
      for position in sorted(list((hash[contig].keys()))):
        try:
          hash[contig][int(position)][1].append(int(coverage[contig][position]))
        except:
          hash[contig][int(position)][1].append(0)

  return [samples, hash]
 
#============================================================( snpPHYL )
  ### TODO ###
  # 1. support snpPHYL output
  exit()

#============================================================( processHash )
def processHash(samples, hash, windowSize, stepSize):

  # process hash table
  print("(2/3) processing hash table")

  contigs = sorted(list(hash.keys()))

  # build header
  results = "Contig,FromPos,ToPos,Phi,AggregateDist"
  for sampleType in ["::frequency", "::coverage"]:
    for sample in samples[1:]:
      results += "," + sample + sampleType

  for i in range(0, len(contigs)):
    # set start and end values for hash table
    start = 1
    indexes = sorted(list(hash[contigs[i]].keys()))
    max = indexes[-1]
    # sort all of the hash keys into window ranges
    rangeDict = {}
    for index in indexes:
      windowStart = ((index//stepSize)*stepSize)+start

      # move window back to capture all instances of index
      while index <= windowStart+windowSize-stepSize:
        windowStart -= stepSize

      # move window forward and add all instances of index to hash table
      while index >= windowStart and index <= windowStart+windowSize:
        if str(windowStart)+"::"+str(windowStart+windowSize) not in rangeDict:
          rangeDict[str(windowStart)+"::"+str(windowStart+windowSize)]=[]
        rangeDict[str(windowStart)+"::"+str(windowStart+windowSize)].append(index)
        windowStart += stepSize

    # create empty fasta data for phi
    fasta = {}
    blankSample = ''*(windowSize+1)
    for sample in samples:
      fasta[sample] = ['>'+sample,str(blankSample)]

    # create counter array to count snps
    snpEmpty = []
    for j in range(0, len(samples)):
      snpEmpty.append(0)

    # create array for coverage
    coverageEmpty = []
    for j in range(1, len(samples)):
      coverageEmpty.append(0)

    maxCount = ((max-windowSize-1)//stepSize)+1
    if ((max-windowSize-1)%stepSize) != 0:
      maxCount += 1
    if max < windowSize+1:
      maxCount = 1

    # while start of window i < maximum possible hash table size
    for j in range(1, maxCount+1):

      # reset values for results
      if start+windowSize+1 > max:
        string = str(contigs[i]) + "," + str(start) + "," + str(max)
      else:
        string = str(contigs[i]) + "," + str(start) + "," + str(start+windowSize)
      phi = "--"
      snpCount = list(snpEmpty)
      coverageCount = list(coverageEmpty)

      if str(start)+"::"+str(start+windowSize) in rangeDict:
        for index in rangeDict[str(start)+"::"+str(start+windowSize)]:
          for k in range(0, len(snpCount)):
            fasta[samples[k]][1] += hash[contigs[i]][index][0][k]
            if hash[contigs[i]][index][0][k] != hash[contigs[i]][index][0][0]:
              snpCount[k] += 1
          if sum(snpCount[1:]) > 1:
            snpCount[0] += 1
          for k in range(0, len(coverageCount)):
            coverageCount[k] += hash[contigs[i]][index][1][k]
        for k in range(0, len(coverageCount)):
          coverageCount[k] //= len(rangeDict[str(start)+"::"+str(start+windowSize)])

        ### TODO ###
        # 1. viciously improve run time by solving PHI
        # 2. variable window size for phi statistic to increase capture rate of informative sites
        # 3.
        #tempFasta = ""
        #for sample in samples:
        #  tempFasta += fasta[sample][0]+'\n'+fasta[sample][1]+'\n'
        #f = open("tempFasta.fasta", 'w')
        #f.write(tempFasta)
        #f.close()
        #phi = subprocess.Popen("module load phipack/phipack; Phi -w 1 -f tempFasta.fasta", universal_newlines=True, shell=True, stdout=subprocess.PIPE)
        #phi = str(phi.stdout.read().split(":")[-1].strip())

      for sample in samples:
        fasta[sample][1] = blankSample

      # add all counts and coverage to results
      string += ","+phi
      for k in range(0, len(snpCount)): 
        string += ","+str(snpCount[k])
      for k in range(0, len(coverageCount)): 
        string += ","+str(coverageCount[k])
      results += "\n"+string

      if j%((maxCount//50)+1) == 0 or j == maxCount:
        echo = str("processing " + str(i+1) + " of " + str(len(contigs)) + " contigs \|" + "\u2589"*int(j/maxCount*25) + "-"*(25-int((j)/maxCount*25)) + "\| \(" + str(int((j)/maxCount*100)) + "%\)\ \ ")
        subprocess.call("echo -ne \r\'    \'" + echo, universal_newlines=True, shell=True)

      start += stepSize
  print()
  return results

#============================================================( printResults )
def printResults(results):
  
  f = open("snpDensityOut/snpDensityMatrix.csv", 'w')
  f.write(results)
  f.close()

  results = results.split("\n")
  for i in range(0, len(results)):
    results[i] = "\"" + results[i] + "\","

  results = "buildGraphs([" + "\n".join(results) + "])"

  f = open("snpDensityOut/snpDensityMatrix.jsonp", 'w')
  f.write(results)
  f.close()

  return

#==================================================================================================== ( main )
#============================================================( main )
###### parse arguments and run program

def main():

  print("\n====================( Job Starting: "+sys.argv[0]+" )")
  # time length of script
  timer = JobTimer()

  # parse all arguments on initial call
  args = ArgParse().getArgs()
  if args.DIR[0][-1] != "/":
    args.DIR[0] += "/"

  # create "snpDensityOut" DIR in the working DIR if one does not exist
  if not os.path.exists("snpDensityOut"):
    os.mkdir("snpDensityOut")
    for file in os.listdir("/".join(sys.argv[0].split("/")[:-1]) + "/html"):
      subprocess.call("cp -r " + "/".join(sys.argv[0].split("/")[:-1]) + "/html/" + file + " snpDensityOut/" + file, universal_newlines=True, shell=True)
  ### build hash table ###
  # samplePositions.txt format:
  #     contig1 position1
  #     contig1 position2
  #     contig2 position1
  # 
  # example hash table:
  #     "{contig1: {position1: [['A', 'A', 'A', 'A', 'G'], ['42', '42', '42', '42']], position2: [['T', 'T', 'T', 'T', 'C'], ['42', '42', '42', '42']]},
  #       contig2: {position1: [['A', 'A', 'A', 'A', 'G'], ['42', '42', '42', '42']], position2: [['T', 'T', 'T', 'T', 'C'], ['42', '42', '42', '42']]}}"

  ### TODO ###
  # 1. test cfsan and lyveset with multiple contig data
  print("(1/3) gathering coverage from bam files")
  if args.nasp == True:
    results = naspHash(args.DIR[0])
  elif args.cfsan == True:
    results = cfsanHash(args.DIR[0])
  elif args.lyve == True:
    results = lyveHash(args.DIR[0])

  ### rolling window ###
  # example rolling window:
  #     contig,from,to,phi,aggregate,sample1,sample2,sample3,sample4
  #     contig1,1,1000,--,0,0,0,0,1
  #     contig1,1001,2001,--,0,0,0,0,1
  results = processHash(results[0], results[1], args.window, args.step)

  # print results to file
  print("(3/3) writing results csv file")
  print("    creating file \"./snpDensityOut/snpDensityMatrix.csv\"")
  printResults(results)

  timer.getTime()
  print("====================( Job Complete: "+sys.argv[0]+" )\n")

#============================================================( run main )
###### runMain
if __name__ == "__main__":
  main()

