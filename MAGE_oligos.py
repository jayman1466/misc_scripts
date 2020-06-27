from RNAfold_wrapper import RNAfold

#output reverse complement
def reversecomplement(inpt):
    codonr=inpt[::-1]
    basecomplement = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'a':'t', 'c':'g', 't':'a', 'g':'c'}
    codonrc=''
    for base in codonr:
        codonrc += basecomplement[base]
    return codonrc

#minimize dG of the oligo
def optmage(oligo):
    i=0
    dG=-100
    optoligo=""
    while i<60:
        this_oligo = oligo[i:i+90]
        this_dG = RNAfold(this_oligo).folding[1]
        if this_dG > dG:
            optoligo = this_oligo
            dG = this_dG
        i = i+1
    return optoligo,dG

f=open('essentials.txt')
ess=f.read().split()
f.close

#open list of genes for which to design TGA -> TAA codon MAGE oligos
f=open('20190710.txt')
relevant=f.read().split()
f.close

f=open('overlapslist.txt')
overlapping=f.read().split()
f.close

f=open('exoverlapslist.txt')
exoverlapping=f.read().split()
f.close

seq=open('mg1655_seq.txt').read().replace('\n','').upper()


import csv
#import gene coordinate file
with open('gene_coords.txt') as csvfile:
    coords=csv.reader(csvfile, delimiter='\t')
    for row in coords:
        row[1].replace("'","")
        if row[1] in relevant:
            if row[2]=="Clockwise":
                start=int(row[3])
                stop=int(row[4])
                essential = 0
                ovrlp= 0
                exovrlp=0
                if seq[stop-3:stop]=='TGA':
                    if row[1] in overlapping:
                        left=stop-3
                        right=stop
                        mutation='taatga'
                        ovrlp=1
                        if row[1] in exoverlapping:
                            exovrlp=1
                    else:
                        left=stop-3
                        right=stop
                        mutation='taa'
                    
                    #create a +/- 75bp window around the mutation site to scan for minimal oligo structure
                    oligo=seq[stop-76:stop-3]
                    oligo+=mutation
                    oligo+=seq[stop:stop+74]
                    if stop not in range(1597981,3932974):
                        oligo=reversecomplement(oligo)
                    optimize = optmage(oligo)
                    oligo = optimize[0]
                    oligo = oligo[0] + "*" + oligo[1] + "*" + oligo[2:]
                    dG = optimize[1]
                    
                    replicore = 1
                    if start in range(1597981,3923998):
                        replicore = 2

                    if row[1] in ess:
                        essential=1
                    
                    print(row[1],'\t',oligo,'\t',dG,'\t','+','\t',replicore,'\t',left,'\t',right,'\t', 'M', '\t', mutation, '\t', essential, '\t', ovrlp, '\t',exovrlp)
            else:
                start=int(row[3])
                stop=int(row[4])
                essential = 0
                ovrlp= 0
                exovrlp=0
                if seq[start-1:start+2]=='TCA':
                    if row[1] in overlapping:
                        left=start-1
                        right=start+2
                        mutation='taatga'
                        ovrlp=1
                        if row[1] in exoverlapping:
                            exovrlp=1
                    else:
                        left=start-1
                        right=start+2
                        mutation='taa'

                    oligo=seq[start-74:start-1]
                    oligo+=reversecomplement(mutation)
                    oligo+=seq[start+2:start+75]
                    if start not in range(1597981,3932974):
                        oligo=reversecomplement(oligo)
                    optimize = optmage(oligo)
                    oligo = optimize[0]
                    oligo = oligo[0] + "*" + oligo[1] + "*" + oligo[2:]
                    dG = optimize[1]

                    replicore = 1
                    if start in range(1597981,3923998):
                        replicore = 2

                    if row[1] in ess:
                        essential = 1

                    print(row[1],'\t',oligo,'\t',dG,'\t','-','\t',replicore,'\t',left,'\t',right,'\t', 'M', '\t', mutation, '\t', essential, '\t', ovrlp, '\t',exovrlp)




