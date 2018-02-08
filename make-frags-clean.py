import csv
import re
from essential_mods import *

def loadsc(scaff):
    fin = open('/home/lab/fox/vv2/indivs/'+scaff+'.fa','r')
    data = fin.read()
    data = data.replace('>' + str(scaff) ,'')
    data = data.replace('\n','')
    return data

def initfrag(folder, frag):
    outfile = open('/home/rando2/genome/RACA/raca_frags/'+folder+'/' + frag + '.fa','w')
    outfile.write('>' + frag + '\n')
    return "Null" #outfile

def grabseq(scaff, slow, shigh, sdir):
    seq = loadsc(scaff)[slow:shigh] #raca is using zero indexed values as they run 0 to i where i the length of the scaff
    if sdir == "-":
        seq = rev_comp(seq)
    return seq

def openfrag(folder, frag):
    fopen = open('/home/rando2/genome/RACA/raca_frags/'+folder+'/' + frag + '.fa','r')
    seq = fopen.read()
    seq = seq.replace('>' + frag, '')
    seq = seq.replace('\n', '')
    return seq

with open('/home/rando2/genome/RACA/rec_chrs.refined_40kbp.txt', 'rb') as csvfile:
    data = csv.reader(csvfile, delimiter='\t')
    last = "" #raca
    outfile = ""
    for line in data:
        if line[3]=="GAPS": #if it's a gap, don't try to unpack, just write and move on
            outfile.write("N"*100)
            continue
        #Here we make individual fragments corresponding to each RACA Fragment
        raca, rlow, rhigh, scaff, slow, shigh, sdir = line  #12_33   0       12472085        scaffold55      0       12472085        -
        rlow, rhigh, slow, shigh = [int(i) for i in rlow, rhigh, slow, shigh]
        if len(last)==0 : #if it's the first RACA frag
            outfile = initfrag('raca_orig',raca)
            fake=2 #just need to write annotation so skip this
        elif last != raca: #if we're starting a new RACA frag
            outfile.write('\n\n')
            outfile.close()
            outfile = initfrag('raca_orig',raca)
            fake = 2 #just trying to get annotation rn
        #else: # if the same raca chr
        last = raca
        seq = grabseq(scaff, slow, shigh, sdir)
        outfile.write(seq)

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
        #elif + and + don't 
        #elif - and - don't
        #outfile.write(seq)


#dogsynt = {'1a':'-','1b':'+','1c':'+','2a':'+','2b':'-','2c':'+','3a':'+','3b':'-','3c':'+','4a':'+','4b':'-','4c':'+','5a':'+','5b':'-','5c':'-','5d':'+','5e':'-','6a':'+','6b':'+','7a':'+','7b':'+','8a':'+','8b':'-','9a':'-','9b':'+','10a':'+','10b':'+','11a':'-','11b':'+','12a':'-','12b':'-','12c':'+','13a':'-','13b':'-','13c':'-','14a':'-','14b':'+','15a':'-','15b':'+','15c':'+','16a':'-','16b':'+','X':'+'}

#chrom = ["chr" + x for x in range(1,39)]
#chrom.append("X")

#    open('./dog_syntblocks/' + raca + '.fa','w')
