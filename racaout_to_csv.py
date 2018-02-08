import csv
"""This program converts RACA's output files to a csv-readable format"""

#Change these parameters
RACA_directory = "./" #where is the RACA output located?
target = 'vv2' #name of target genome as provided to RACA
reference = 'canFam3' #name of reference genome as provided to RACA
outgroup = 'felCat5' #name of outgroup genome as provided to RACA

#Conversion below
genomes = {"target":target, "reference":reference, "outgroup":outgroup}

for gtype, gname in genomes.items(): 
    fileroot = 'rec_chrs.'+gname+'.segments.refined'
    filein = open(fileroot + '.txt', 'r')
    fileread = filein.read()
    
    with open(fileroot + '.csv', 'wb') as fileout:
        writerbot = csv.writer(fileout, delimiter = ',')

        header = ["blocknum"]
        blocks = fileread.split('>')
        for block in blocks[1:]:
            lines = block.split('\n')
            blocknum = lines[0] #RACA assigns each block a number
            outlist = [blocknum]
            
            for line in lines[1:]:
                if line != "":
                    chrominfo, direction = line.split(' ')
                    if len(header) >0:
                        header.append(chrominfo.split('.')[0])
                        header.extend(['low','high','direction'])
                    fragname = chrominfo.split('.')[1].split(':')[0]
                    low, high = (chrominfo.split(':')[1]).split('-')
                    outlist.extend([fragname, low, high, direction])
            if len(header) > 0:
                writerbot.writerow(header)
                header = []
            writerbot.writerow(outlist)
