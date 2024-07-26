import os,sys,pymol
from numpy import mean

#TODO: This program needs to be run in the folder of the clusterstructures and not in the main library folder than the other programs... at some point you wanna fix this to adopt.
clusterlist=[item for item in os.listdir('.') if '.pdb' in item and not 'all' in item]
clusterlist.sort()
#Making a dictionary giving all my clusters a running number as key. Need that to store the key as bfactor, so I can retrieve the cluster later.
numberdict={}
for i in range(0,len(clusterlist)):
    numberdict[clusterlist[i]]=i
# Ok, conserve that dictionary by printing it out as textfile.
outfile=open('dictionary_clusterID_numbers.txt','w')
for item in numberdict:
    outfile.write(str(numberdict[item])+':'+str(item)+'\n')
outfile.close()



calist=[]
clist=[]
nlist=[]
felist=[]

for item in clusterlist:
    f=open(item,'r')
    fe=[]
    for line in f:
        if 'ATOM' in line[0:4]:
            if 'CA  CYS' in line[13:20]: #This should find all Calpha of cys residues.
                calist.append(line)
                newline=calist[-1][0:60]+' '+str(numberdict[item]).zfill(5)+calist[-1][66:] # Here was a problem with the correct row counting leading to weird 6 digit clusternumbers when the original b factor at this position was in the hundreds.
                calist[-1]=newline
            if 'C   CYS' in line[13:20]: #This should find all C of cys residues.
                clist.append(line)
                newline=clist[-1][0:60]+' '+str(numberdict[item]).zfill(5)+clist[-1][66:]
                clist[-1]=newline
            if 'N   CYS' in line[13:20]: #This should find all n of cys residues.
                nlist.append(line)
                newline=nlist[-1][0:60]+' '+str(numberdict[item]).zfill(5)+nlist[-1][66:]
                nlist[-1]=newline
        if 'HETATM' in line[:6] and 'SF4' in line[17:20] and 'FE' in line[76:78].upper(): # This should grab all Fe atoms within SF4 clusters; note the 'upper', to avoid trouble with lower/upper case writing of element.
            fe.append(line)
    atoms=[]
    for atom in fe:
        atoms.append([atom[31:38],atom[38:46],atom[46:54]])
    #print item # debugging only
    #print [number[0] for number in atoms]
    averagefe=[mean([float(number[0]) for number in atoms]),mean([float(number[1]) for number in atoms]),mean([float(number[2]) for number in atoms])]
    try:
        newline1=fe[-1][0:30]+'%8.3f' %(averagefe[0])+'%8.3f' %(averagefe[1])+'%8.3f' %(averagefe[2])+fe[-1][54:]
        newline2=newline1[0:60]+' '+str(numberdict[item]).zfill(5)+newline1[66:]
        felist.append(newline2)
    except:
        pass
    f.close()

outfile=open('allca.pdb','w')
for item in calist:
    outfile.write(item+'\n')

outfile.close()

outfile=open('allc.pdb','w')
for item in clist:
    outfile.write(item+'\n')

outfile.close()

outfile=open('alln.pdb','w')
for item in nlist:
    outfile.write(item+'\n')

outfile.close()

outfile=open('allfe.pdb','w')
for item in felist:
    outfile.write(item+'\n')

outfile.close()
