import csv
import re
from essential_mods import *

"""Note: this script is designed specifically for the red fox project and requires information that may
not be available for other assembly projects (e.g., meiotic linkage map). You are welcome to try to adapt
it for your purposes, but it will not run. It is here to show the logic we used for assembly."""

# Setup parameters
scaff_dir = '/home/lab/fox/vv2/indivs' #where are the scaffold files stored?
raca_frag_dir = '/home/rando2/genome/RACA/raca_frags' #where RACA fasta should be written
rec_chrs_file = '/home/rando2/genome/RACA/rec_chrs.refined_40kbp.txt' #from RACA

def loadsc(scaff):
    # Open the scaffold specified. Assumes multifasta has been split into individual files for each scaffold.
    # Accepts: basename of scaffold file (e.g., "scaffold1" if the file is "scaffold1.fa"
    # Returns: contents of "scaffold1.fa"
    fin = open(scaff_dir + "/" + scaff+'.fa','r')
    data = fin.read()
    data = data.replace('>' + scaff ,'')
    data = data.replace('\n','')
    return data

def initfrag(folder, frag):
    # Initialize fasta file for RACA fragment
    # Accepts: name of directory that these fragments should be build in within the raca_frag_dir and
    #          name of fragment 
    # Returns: Name of file
    fname=raca_frag_dir + '/' + folder + '/' + frag + '.fa'
    outfile = open(fname,'w')
    outfile.write('>' + frag + '\n')
    return fname

def grabseq(scaff, slow, shigh, sdir):
    # Returns sequence based on start and end positions and direction
    # Accepts: scaffold name (format: "scaffold1"), start position (low) as integer, end position (high) as integer
    #          direction of sequence desired ("+" or "-")
    # Returns: specified sequence, reverse and comped if negative direction
    seq = loadsc(scaff)[slow:shigh] #raca is using zero indexed values as they run 0 to i where i the length of the scaff
    if sdir == "-":
        seq = rev_comp(seq)
    return seq

def openfrag(folder, frag):
    # Open an existing RACA fragment
    # Accepts: Name of folder containing RACA fragments within raca_frag_dir
    # Returns: Sequence of RACA fragment
    fopen = open(raca_frag_dir+ "/" +folder+'/' + frag + '.fa','r')
    seq = fopen.read()
    seq = seq.replace('>' + frag, '')
    seq = seq.replace('\n', '')
    return seq

with open(rec_chrs_file, 'rb') as csvfile:
    data = csv.reader(csvfile, delimiter='\t')
    last = "" #raca
    outfile = ""
    for line in data:
        if line[3]=="GAPS": #if it's a gap, don't try to unpack, just write and move on
            outfile.write("N"*100)
            continue
        #Here we make individual fragments corresponding to each RACA Fragment
        raca, rlow, rhigh, scaff, slow, shigh, sdir = line  #format of data used: 12_33   0       12472085        scaffold55      0       12472085        -
        rlow, rhigh, slow, shigh = [int(i) for i in rlow, rhigh, slow, shigh]
        if len(last)==0 : #if it's the first RACA frag
            outfile = initfrag('raca_orig',raca)
        elif last != raca: #if we're starting a new RACA frag
            outfile.write('\n\n')
            outfile.close()
            outfile = initfrag('raca_orig',raca)
        last = raca
        seq = grabseq(scaff, slow, shigh, sdir)
        outfile.write(seq)

"""If you don't have a meiotic linkage map, don't try to run the code below
        
#Check whether each raca fragment goes + or - in dog
racadog = {}
with open('/home/rando2/genome/RACA/rec_chrs.canFam3.segments.formatted.csv','rb') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    for line in data:
        raca, rlow, rhigh, rdir, dog, dlow, dhigh, ddir = line #1a,1949,582677,+,chr1,24994866,25534824,+
    
        if raca not in racadog.keys():
            racadog[raca] = [[int(dlow), int(dhigh), ddir]]
        else:
            racadog[raca].append([int(dlow), int(dhigh), ddir])

final_dir = {}
for rd,ls in racadog.items():
    if len(ls) == 1:
        final_dir[rd] = ls[0][2]
    else:
        maxcvg = 0
        maxdir = ""
        for block in ls:
            if block[1]-block[0] > maxcvg:
                maxcvg = block[1]-block[0]
                maxdir = block[2]
        final_dir[rd] = maxdir

with open('/home/rando2/genome/RACA/racafrag_foxfrag.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    last = [] #raca, fox #15a,10b
    #outfile = ""
    head = 0

    for line in data:        
        if head == 0:
            head = 1
            continue
        #Here we make individual fragments corresponding to each Dog Fragment
        raca, fox = line  #26a, 10a
        print line
        if len(last)==0 : #if it's the first RACA frag
            fake = 2 #just for annotation
            #outfile = initfrag('dog_syntblocks',fox)
        elif last[1] != fox: #if we're starting a new RACA frag
            #outfile.write('\n\n')
            #outfile.close()
            #outfile = initfrag('dog_syntblocks',fox)
            fake=2 #just for writing annotation right now
        else: # if the same raca chr, write buffer
            prevdog = [int(i[0]) for i in racadog[last[0]]] #what are the dog positions of the last fragment?
            prevdog.extend([int(i[1]) for i in racadog[last[0]]])
            currentdog = [int(i[0]) for i in racadog[raca]] #what are the dog positions of the current fragment?
            currentdog.extend([int(i[1]) for i in racadog[raca]])
            if min(currentdog) > max(prevdog):
                #outfile.write("N"*(min(currentdog) - max(prevdog)))
                print "intergap", last, "to", raca, fox, ":", min(currentdog) - max(prevdog)
            elif max(currentdog) < min(prevdog):
                #outfile.write("N"*(min(prevdog) - max(currentdog)))
                print "intergap", last, "to", raca, fox, ":", min(prevdog) - max(currentdog)

        last = [raca, fox]
        seq = openfrag('raca_orig',raca)

        dogsynt = {'1a':'-','1b':'+','1c':'+','2a':'+','2b':'-','2c':'+','3a':'+','3b':'-','3c':'+','4a':'+','4b':'-','4c':'+','5a':'+','5b':'-','5c':'-','5d':'+','5e':'-','6a':'+','6b':'+','7a':'+','7b':'+','8a':'+','8b':'-','9a':'-','9b':'+','10a':'+','10b':'+','11a':'-','11b':'+','12a':'-','12b':'-','12c':'+','13a':'-','13b':'-','13c':'-','14a':'-','14b':'+','15a':'-','15b':'+','15c':'+','16a':'-','16b':'+','X':'+'}
        if (final_dir[raca]=="+" and dogsynt[fox] == "-") or (final_dir[raca]=="-" and dogsynt[fox] == "+"):
            print "reversing", raca, "to use on", fox
            seq = rev_comp(seq)
        else:
            print raca, "runs same direction on", fox
"""
