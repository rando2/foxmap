def rev_comp(seq):
    subs = {"A":"T","C":"G","G":"C","T":"A", "N":"N", 'n':'n', 'a':'t', 'c':'g', 'g':'c', 't':'a'}
    new_seq = ""
    for i in range(1,len(seq)+1):
        new_seq += subs[seq[-i]]
    return new_seq

def clean_seq(seq, name):
    seq = seq.replace('>' + name,'')
    seq = seq.replace('\n','')
    return seq

def get_scaff(scaff, start, end, direct):
    fin = open('/home/lab/fox/vv2/indivs/' + scaff + '.fa','r')
    seq = fin.read()
    seq = clean_seq(seq, scaff)
    target = seq[int(start):int(end)]
    if direct == "-":
        target = rev_comp(target)
    return target
