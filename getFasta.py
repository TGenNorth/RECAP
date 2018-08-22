#!/usr/bin/env python3
import sys
import math
import subprocess

def main():
  print("getting fasta and duplicates")

  # get DIR
  DIR = sys.argv[1]

  # get reference file
  f = open(DIR + "reference/reference.fasta", 'r')
  ref = f.read().split(">")[1:]
  f.close()
  contigs = {}
  for i in range(0, len(ref)):
    # sequences begin with space to match bestsnp count
    contigs[ref[i].split("\n")[0]] = "".join(ref[i].split("\n")[1:])

  # range file
  RANGE = sys.argv[2]
  f = open(RANGE, 'r')
  newRange = f.read().split("\n")[0].split("\t")
  f.close()
  contig = newRange[0]
  newRange = newRange[1:]
  if int(newRange[0]) < 0:
    newRange[0] = 0
  if int(newRange[1]) > len(contigs[contig]):
    newRange[1] = len(contigs[contig])

  # get snps
  header = subprocess.Popen("head -n 1 " + DIR + "matrices/bestsnp.tsv", universal_newlines=True, shell=True, stdout=subprocess.PIPE)
  header = header.stdout.read().strip()
  snps = subprocess.Popen("cat " + DIR + "matrices/bestsnp.tsv | grep \""+ contig +"\"", universal_newlines=True, shell=True, stdout=subprocess.PIPE)
  snps = snps.stdout.read().strip().split("\n")[:-1]

  #f = open(DIR + "matrices/bestsnp.tsv", 'r')
  #snps = f.read().split("\n")[:-1]
  #f.close()
  samples = header.split("\t#SNPcall")[0].split("\t")[2:]
  for i in range(0, len(samples)):
    samples[i] = samples[i].split("::")[0]
  #snps = snps[1:]
  newSamples = {}
  for sample in samples:
    newSamples[sample+"::"+contig] = list(contigs[contig])

  # build samples
  for i in range(0, len(snps)):
    pos = snps[i].split("\t")[0].split("::")
    line = snps[i].split("\t")[2:len(samples)+2]
    if pos[0] == contig:
      for j in range(0, len(line)):
        newSamples[samples[j]+"::"+pos[0]][int(pos[1])-1]=line[j]

  # build fasta output
  fasta = ""
  for sample in samples:
    fasta += ">"+sample+"::"+contig+"\n"
    newSamples[sample+"::"+contig] = newSamples[sample+"::"+contig][int(newRange[0]):int(newRange[1])]
    for i in range(0,len(newSamples[sample+"::"+contig]),80):
      fasta += "".join(newSamples[sample+"::"+contig][i:i+80]) + "\n"

  f = open("out.fasta", 'w')
  f.write(fasta)
  f.close()

  # write duplicates
  f = open(DIR + "reference/duplicates.txt", 'r')                                                                                                             
  dup = f.read().strip()
  f.close()

  dup = dup.split(">")[1:]
  dupContigs = {}
  for line in dup:
    dupContigs[line.split("\n")[0]] = "".join(line.split("\n")[1:])

  f = open(sys.argv[2], 'r')
  rangeFile = f.read().strip().split("\t")
  f.close()

  # open range file
  f = open("outDup.txt", 'w')
  for index in sorted(dupContigs.keys()):
    if index == rangeFile[0]:
      min = int(rangeFile[1])
      if min < 0:
        min = 0 
      max = int(rangeFile[2])
      string = ("1"*min)+dupContigs[index][min:max]+("1"*(len(dupContigs[index])-max))
      if max <= min:
        string = dupContigs[index]
    else:
      string = dupContigs[index]
    newString=">"+index+"\n"
    for i in range(0,len(string),80):
      newString+=string[i:i+80]+"\n"
    #print(len(newString))
    f.write(newString)
  f.close()

  subprocess.call("module load nasp; nasp matrix --dto-file " + DIR + "matrix_dto.xml --reference-dups outDup.txt", universal_newlines=True, shell=True)

#============================================================( run main )
###### runMain
if __name__ == "__main__":
  main()
