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
import threading

# for python string substitution
import re
from collections import Counter

# for getting time taken for jobs
import time

# for running terminal commands
import subprocess
import os

# for getting the phi statistic
import phiStat

#==================================================================================================== ( classes )
#------------------------------------------------------------( getArgs )
### parse and return all arguments to NAMESPACE
# all arguments are collected from command line
#
# 'renameSamples -h' should run in any event that the arguments are invalid

class ArgParse:
  def __init__(self):

    # initialize argument parser
    self.parser=argparse.ArgumentParser(usage='''%(prog)s [-h, --help] [--nasp | --cfsan | --lyve | --snvphyl] [--window WINDOW] [--step STEP] DIR''', formatter_class=argparse.RawDescriptionHelpFormatter, description='''''',epilog='''''')

    # set input file type
    self.input = self.parser.add_mutually_exclusive_group(required=True)
    self.input.add_argument('--nasp', action='store_true', help='Northern Arizona SNP Pipeline')
    self.input.add_argument('--cfsan', action='store_true', help='CFSAN SNP Pipeline')
    self.input.add_argument('--lyve', action='store_true', help='LYVE SNP Extraction Tool')
    self.input.add_argument('--snvphyl', action='store_true', help='SNV Phylogenomics Pipeline')

    # set window size and step size
    self.parser.add_argument('-p', default=False, action='store_true', help='set size of capture window to measure snp density')
    self.parser.add_argument('-d', default=False, action='store_true', help='set size of capture window to measure snp density')
    self.parser.add_argument('--window', default=10000, type=int, help='set size of capture window to measure snp density')
    self.parser.add_argument('--step', default=5000, type=int, help='set step size to move window down snp matrix')
    self.parser.add_argument('DIR', nargs='*', help='SNP pipeline output directory')
    self.args = self.parser.parse_args()
    self.args.window = abs(self.args.window)-1
    self.args.step = abs(self.args.step)

  def getArgs(self):
    return self.args

  def getHelp(self):
    self.parser.print_help()
    exit()

#------------------------------------------------------------( jobTimer )
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
#------------------------------------------------------------( readFile )
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

def bamThread(DIR, hash, snpPositions, sample):
  if not os.path.exists("snpDensityOut/" + "".join(sample.split(" ")) + "-Depth.txt"):
    #print("    creating file \"./snpDensityOut/" + "".join(sample.split(" ")) + "-Depth.txt\" (" + str(i) + " of " + str(len(samples) - 1) + ")")

    # get depth from BAM file (faster than VCF)
    if os.path.exists(DIR + "bwamem/" + sample + "-bwamem.bam"):
      coverage = subprocess.Popen("module load samtools; samtools depth -b snpDensityOut/snpPositions.txt " + DIR + "bwamem/" + sample + "-bwamem.bam", universal_newlines=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      coverage, coverageERR = coverage.communicate()

    # if no BAM file, get depth from VCF file
    elif os.path.exists(DIR + "gatk/" + sample + "-bwamem-gatk.vcf"):
      coverage = subprocess.Popen("module load vcftools; vcftools --vcf " + DIR + "gatk/" + sample + "-bwamem-gatk.vcf --positions snpDensityOut/snpPositions.txt --site-depth --stdout | awk 'NR > 1 {print $1 \"\t\" $2 \"\t\" $3}'", universal_newlines=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      coverage, coverageERR = coverage.communicate()

    # if no BAM or VCF file, then write in coverage as 0
    else:
      coverage = list(snpPositions)
      for j in range(0, len(coverage)):
        coverage[j] += "\t0"
      coverage = "\n".join(coverage)
    f = open("snpDensityOut/" + "".join(sample.split(" ")) + "-Depth.txt", 'w')
    f.write(coverage)
    f.close()
  else:
    #print("    reading file \"./snpDensityOut/" + "".join(sample.split(" ")) + "-Depth.txt\" (" + str(i) + " of " + str(len(samples) - 1) + ")")
    coverage = readFile("snpDensityOut/" + "".join(sample.split(" ")) + "-Depth.txt")
  if coverage.split("\n")[-1] == "":
    coverage = coverage.split("\n")[:-1]
  else:
    coverage = coverage.split("\n")

  # assemble hash
  for line in coverage:
    line = line.split("\t")
    hash[line[0]][int(line[1])][1].append(int(line[2]))


#------------------------------------------------------------( naspCoverage )
def naspHash(args):

  DIR = args.DIR[0]
 

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
  f.write("\n".join(snpPositions))
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
  print("    building hash table") 
  for i in range(1, len(file)):
    line = file[i].split("\t")[:len(samples)+1]
    ID = line[0].split("::")
    contig = "::".join(ID[:-1])
    position = int(ID[1])
    hash[contig][position][0] = line[1:]
    if args.d:
      hash[contig][position][1] = (len(samples)-1) * [0,]

  # write SNP coverage to files
  if args.d:
    return [samples, hash]

  for i in range(1, len(samples)):
    print("    creating file \"./snpDensityOut/" + "".join(samples[i].split(" ")) + "-Depth.txt\" (" + str(i) + " of " + str(len(samples) - 1) + ")")
    while threading.activeCount() > 10:
      time.sleep(1)
    t = threading.Thread(target=bamThread, args = [DIR, hash, snpPositions, samples[i]])
    t.daemon = True
    t.start()

  while threading.activeCount() > 1:
    time.sleep(1)



  return [samples, hash]


#------------------------------------------------------------( cfsanHash )
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

      # get depth from BAM file
      if os.path.exists(DIR + "/" + samples[i] + "/reads.sorted.bam"):
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

      # if no BAM file, then write in coverage as 0
      else:
        coverageOut = list(snpPositions)
        for j in range(0, len(coverageOut)):
          coverageOut[j] += "\t0"
          line = coverageOut[j].split("\t")
          if line[0] not in coverage:
            coverage[line[0]] = {}
          if int(line[1]) not in coverage[line[0]]:
            coverage[line[0]][int(line[1])] = 0
          coverage[line[0]][int(line[1])] += int(line[2])
        coverageOut = "\n".join(coverageOut)

      f = open("snpDensityOut/" + "".join(samples[i].split(" ")) + "-Depth.txt", 'w')
      f.write(coverageOut)
      f.close()

    # read coverage from file and build coverage hash table
    else:
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
    # if SNP information is missing, add 0 to hash
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

#------------------------------------------------------------( lyveHash )
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

      if bam != []:
        coverage1 = subprocess.Popen("module load samtools; samtools depth -b snpDensityOut/snpPositions.txt " + bam[0], universal_newlines=True, shell=True, stdout=subprocess.PIPE)
        coverage1 = coverage1.stdout.read().split("\n")[:-1]
        coverage2 = subprocess.Popen("module load samtools; samtools depth -b snpDensityOut/snpPositions.txt " + bam[1], universal_newlines=True, shell=True, stdout=subprocess.PIPE)
        coverage2 = coverage2.stdout.read().split("\n")[:-1]

      else:
        vcf = subprocess.Popen("ls " + DIR + "vcf/" + samples[i] + "*.vcf.gz", universal_newlines=True, shell=True, stdout=subprocess.PIPE)
        vcf = vcf.stdout.read().split("\n")[:-1]

        if vcf != []:
          coverage1 = subprocess.Popen("module load vcftools; vcftools --gzvcf " + vcf[0] + " --positions snpDensityOut/snpPositions.txt --site-depth --stdout | awk 'NR > 1 {print $1 \"\t\" $2 \"\t\" $3}'", universal_newlines=True, shell=True, stdout=subprocess.PIPE)
          coverage1 = coverage1.stdout.read().split("\n")[:-1]
          coverage2 = subprocess.Popen("module load vcftools; vcftools --gzvcf " + vcf[1] + " --positions snpDensityOut/snpPositions.txt --site-depth --stdout | awk 'NR > 1 {print $1 \"\t\" $2 \"\t\" $3}'", universal_newlines=True, shell=True, stdout=subprocess.PIPE)
          coverage2 = coverage2.stdout.read().split("\n")[:-1]

        else:
          coverage1 = list(snpPositions) 
          for j in range(0, len(coverage1)):
            coverage1[j] += "\t0"
          coverage2 = list(coverage1)

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
 
#------------------------------------------------------------( snpPHYL )
  ### TODO ###
  # 1. support snpPHYL output
  exit()

#------------------------------------------------------------( processHash )
def processHash(samples, hash, windowSize, stepSize):

  # process hash table
  print("(3/4) processing hash table")

  contigs = sorted(list(hash.keys()))

  # build header
  results = ["{\"samples\":","\"".join(str(samples[1:]).split("\'")),",","\"contigs\":[","]}"]
    
  contigList= []
  for i in range(0, len(contigs)):
    contig = ["{\"name\":\"" + str(contigs[i]) + "\",\"data\":[","]}"]

    # set start and end values for hash table
    indexes = sorted(list(hash[contigs[i]].keys()))
    start = 1
    while start + stepSize < indexes[0]:
      start += stepSize
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

    dataList = []

    # while start of window i < maximum possible hash table size
    for j in range(((start-windowSize-1)//stepSize)+1, maxCount+1):
      # reset values for results
      if start+windowSize+1 > max:
        position = [start,max]
      else:
        position = [start,start+windowSize]
      phi = 1.0
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

        if len(samples) >= 4 and len(fasta[samples[0]][1]) >= 5:
          tempFasta = ""
          for sample in samples:
            tempFasta += fasta[sample][0]+'\n'+fasta[sample][1]+'\n'
          phi = phiStat.fasta(tempFasta)
          if phi < 0.0000000001:
            phi = 0.0000000001

      for sample in samples:
        fasta[sample][1] = blankSample

      dataList.append("{\"position\":"+"\"".join(str(position).split("\'"))+",\"phi\":\""+str(phi)+"\",\"aggregate\":"+str(snpCount[0])+",\"SNPs\":"+str(snpCount[1:])+",\"depth\":"+str(coverageCount)+"}")

      if j%((maxCount//50)+1) == 0 and j > 0 or j == maxCount and j > 0:
        echo = str("processing " + str(i+1) + " of " + str(len(contigs)) + " contigs \|" + "\u2589"*int(j/maxCount*25) + "-"*(25-int((j)/maxCount*25)) + "\| \(" + str(int((j)/maxCount*100)) + "%\)\ \ ")
        subprocess.call("echo -ne \r\'    \'" + echo, universal_newlines=True, shell=True)

      start += stepSize
    contig.insert(-1, ",".join(dataList))
    contigList.append("".join(contig))

  results.insert(-1, ",".join(contigList))
  results = "".join(results)

  print()
  return results

#------------------------------------------------------------( printResults )
def printSnps(results):
  out = "exciseContig({"
  contigOut = ""
  for contig in results:
    contigOut += "\"" + contig + "\":{"
    dataOut = ""
    for index in sorted(results[contig]):
      dataOut += str(index) + ":\"" + "".join(results[contig][index][0]) + "\","
    contigOut += dataOut[:-1] + "},"
  out += contigOut[:-1] + "})"

  f = open("snpDensityOut/excise.jsonp", 'w')
  f.write(out)
  f.close()

  return

#------------------------------------------------------------( printResults )
def printResults(results):
  
  results = "buildGraphs(" + results + ")"

  f = open("snpDensityOut/snpDensityMatrix.jsonp", 'w')
  f.write(results)
  f.close()

  return

#==================================================================================================== ( main )
#------------------------------------------------------------( main )
###### parse arguments and run program
def main():

  # time length of script
  print("\n====================( Job Starting: "+sys.argv[0]+" )")
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

  print("(1/4) gathering coverage from BAM/VCF files")
  if args.nasp == True:
    results = naspHash(args)
  elif args.cfsan == True:
    results = cfsanHash(args.DIR[0])
  elif args.lyve == True:
    results = lyveHash(args.DIR[0])

  print("(2/4) writing export and excise data files")
  printSnps(results[1])
  quit()

  results = processHash(results[0], results[1], args.window, args.step)

  ### rolling window ###
  # example rolling window:
  # {"contigs":[{"name":"gi|514064966|ref|NC_021554.1|",
  #              "data":[{"position": [1,500],
  #                       "phi": "--",
  #                       "depth": [1,2,3,4,5],
  #                       "coverage": [1,2,3,4,5]}]}]
  # }

  # print results to file
  print("(4/4) writing results csv file")
  print("    creating file \"./snpDensityOut/snpDensityMatrix.csv\"")
  printResults(results)

  timer.getTime()
  print("====================( Job Complete: "+sys.argv[0]+" )\n")

#------------------------------------------------------------( run main )
###### runMain
if __name__ == "__main__":
  main()

