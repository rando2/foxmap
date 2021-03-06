#!/usr/bin/env python2
import csv
import gzip

"""This script reassigns positions in a vcf file called against vv2.2 (Kukekova et al., 2018)
Adjusted positions are on the chromosome fragments of vv2.4 (Rando et al., 2018) 
The output is a new vcf file. Only the first two columns are modified
The vv2.4 assembly does not include all scaffolds, so any SNV on scaffolds without positions on the chromosomes are skipped
"""
#change this to your file (it assumes gzip, but can be modified for .vcf by changing syntax on line 29 to open(infile, 'rb')
infile='fox5k.quantfilt.recode.vcf.gz'

#outfile is not set up to be compressed, you can modify using gzip (https://docs.python.org/3/library/gzip.html) 
outfile='fox5k.quantfilt.vv3.vcf'

scaff_pos = {}
with open('scripts/vv2.2-vv2.4.csv', 'rb') as posfile:
  #The format of the csv file is: newchrom, newlow, newhigh, oldchrom, oldlow, oldhigh, dir
    refdata = csv.reader(posfile, delimiter=',')
    for line in refdata:
        scaff = line[3]
        if scaff == "GAPS" or scaff == "INTERGAP" or scaff=="Scaffold" or scaff=='': #won't have SNV in gaps
            continue
        scdata = scaff_pos.get(scaff, [])
        scdata.append([int(line[4]),int(line[5]),line[6], line[0], int(line[1]), int(line[2])])
        scaff_pos[scaff]=scdata

header = []
with gzip.open(infile) as vcffile, open(outfile, 'wb') as outfile:
    vcf = csv.reader(vcffile, delimiter='\t')
    writerbot = csv.writer(outfile, delimiter='\t')
    lastscaff= ""
    scaffinfo = []
    for line in vcf:
        if line[0][0]=="#": #all the headerlines
            writerbot.writerow(line)
            continue
        scaff, pos = line[0:2]
        pos = int(pos)

        if scaff not in scaff_pos.keys():
            print "no conversion for", scaff
            continue
        elif scaff != lastscaff or pos > match[1] : #if it's a new scaffold or new syntenic block, udpate the information
            if scaff != lastscaff: #if we're on a new scaffold, update that information
                print scaff
                lastscaff = scaff
                scaffinfo = scaff_pos[scaff]
            match = [] #identifiy new syntenic block
            for synt in scaffinfo: #for each syntenic block associated with this scaffold
                if synt[0] <= pos and synt[1] >= pos: #check the position of this position relative to the sb bounds
                    match = synt

            if len(match) ==0: #if it didn't find an overlapping match
                print "no match", scaff, pos, scaffinfo
                continue

        chrom = match[3] #now, the chromosome can be identified from the syntenic block
        if match[2] == "+": #based on direction, calculate the position given the sb position
            newpos = match[4] + (pos - match[0]) # if it's position 1 in vcf, it should be good (+0)
        elif match[2] == "-":
            newpos = match[5] - (pos - match[0] - 1) #if it's position 1 in vcf, it should be last - 0
        writerbot.writerow([chrom, newpos] + line[2:])
