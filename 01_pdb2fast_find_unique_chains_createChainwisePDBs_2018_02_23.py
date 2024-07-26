
import sys, os

# Latest version. This one should take a whole folder, look at every pdb file,
# and write a .fasta file with the sequence for every chain in every pdb file in a subfolder sequences.

# First get all pdb files in the folder into an array.
files=os.listdir('.')
pdbfiles=[]
for item in files:
    if '.pdb' in item:
        pdbfiles.append(item)
# Now make a subfolder sequences, but don't fail and stop the program if that already exists.
try:
    os.mkdir('sequences')
except:
    pass

#This dictionary defines which three letter code from the pdbfile is translated to which one letter code in our output fasta sequence.
letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
           'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
           'TYR':'Y','VAL':'V'}
#input_file = open(sys.argv[1])

#And now some functions we will need.

#Following function should take one pdb file and write one fasta file per chain.
def pdb2fasta_per_chain(input_pdb):
    input_file=open(input_pdb,'r')
    pdbarray=[]
    for line in input_file:
        if line[0:4] == 'ATOM': pdbarray.append(line) # Make an array, which contains only the lines with ATOM identifier in PDB file (we are looking for the sequence, so not interested in header or waters or hetatm)
    input_file.close()
    chains=[]
    for item in pdbarray: chains.append(item[21:22]) #These rows contain the chain information in a pdb file; writes an array with all chain info from all files.
    uniquechains=list(set(chains)) # set ... will give the unique members of a list. Unfortunately it is not returned as a list you can work with; therefore need to convert it to a list.
    chainnumber=len(uniquechains)
    for number in range(0,chainnumber):
        #print "Working on chain number "+str(number)+' on structure '+str(input_pdb)
        output_file=open('sequences/'+input_pdb[:-4]+'_'+str(uniquechains[number])+'.fasta','w')
        output_file.write('>'+input_pdb[:-4]+'_'+str(uniquechains[number])+'\n') #This is needed to indicate the start of a sequence in fasta file format. Alignement programs will not understand it otherwise.
        #print '>',sys.argv[1]
        prev = '-1'
        subarray=[]
        for item in pdbarray: # Goes through every line in the pdbarray (which are only the atom lines as defined earlier).
            if item[21:22] == uniquechains[number]: # In case the chain identifier is the one the routine is working on at the moment (number), then write that line into the current subarray.
                subarray.append(item)
        for line in subarray: # Now working on the subarray, which only contains the lines with ATOM from the chain we are interested in at the moment.
            if line[0:4] != 'ATOM': continue #Not sure that is needed since the pdbarray where these lines are taken from, should only contain lines with ATOM anyways as defined above.
            if int(line[23:26]) != prev: # These rows contain the residue number; since subarray contains lines for all atoms, this line checks if the program worked on an atom of this residue already or not.
                try:
                    output_file.write('%c' % letters[line[17:20]]) # In case the program has not yet worked on that residue, take the residue name (rows 17:20) and look up one letter code in my letters dict and write it into the file.
                except: # This is for cases, in which e.g. a DNA structure is present with atoms, which are not HETATM but still not part of recognizable amino acid. Just writes an X instead of residue one letter code.                    output_file.write('X')
                    pass
                prev = int(line[23:26]) # This sets the "prev" variable to the current residue number to make sure one residue is only visited once even though it has many atoms. This is a bit more precise than just looking for CA atoms in each residue since there might be two of those.
        output_file.write('\n')
        output_file.close()

#Following function takes 2 fasta files and determines, if these have roughly the same sequence. This is to check wether two chains in a PDB are actually two different proteins, or just two copies of the same.
from Bio import pairwise2 # Needs biopython installed and uses pairwise2 module of that. (sudo apt-get install python-biopyton)
def check_homodimers_ASU(fasta1,fasta2):
    score=0
    f1=open(fasta1,'r')
    f2=open(fasta2,'r')
    array1=[]
    array2=[]
    for line in f1: array1.append(line)
    for line in f2: array2.append(line)
    f1.close()
    f2.close()
    # Up to here just loaded the two files into two arrays to work with.
    if 0.7<(float(len(array1[1])/len(array2[1])))<1.42:  #This is a quick check if the sequences are roughly the same length within about 70%. Only if they are, there is a chance that they might be homodimers. Otherwise don't even align, because it takes time.
        normalscore=pairwise2.align.globalxx(array1[1],array2[1],score_only=True) #This gives the number of matching AA in the alignment.
        avglength=(len(array1[1])+len(array2[1]))/2 #Average length between the two input protein sequences.
        score=normalscore/avglength # And finally a score, which gives something like a percentage of aligned AA per length.
    else:
        pass
    if score>0.85: #Kinda arbitrarily set this to 85% alignment match. Probably could be much closer to 100% since I only want to make sure that sequences, which have some missing AA, but are otherwise identical are still seen as identical.
        match=True
    else: match=False
    return match #Just returns a boolean value for match or no match.


############################ Start with actual program here ########################
# Lets make fasta files for each chain of all pdb files in the folder.
print 'Generating sequence files for each chain of each of the ' +str(len(pdbfiles))+' input PDB files.'
for input_pdb in pdbfiles:
    pdb2fasta_per_chain(input_pdb)

# Get a list of all the files we just created.
sequencefiles=[item for item in os.listdir('sequences/') if '.fasta' in item]
sequencefiles.sort()
print 'Created '+str(len(sequencefiles))+ ' fasta files.'

# And now check if some of the chains within one PDB file are actually just homo oligomers or copies of the same protein in the ASU.
# This assumes that only pdb files have been selected from the PDB, which contain non-identical sequences; thus we don't have to check anymore if sequences between the different PDB files are identical or not.
#Let's make groups for all the sequences with the same PDB ID. (Probably could be done more elegantly when combining with the programs above, but these were developed independently, so we have to read in fasta files from harddrive instead taking from RAM).
seqnames=[item[:4] for item in sequencefiles] # This deliberately only takes the PDB ID into account and not the 'bundle number' if there is any. I need to compare all chains over the different bundles within one of those files.
uniqueseqnames=list(set(seqnames))
uniqueseqnames.sort
group=[]
for name in uniqueseqnames:
    group.append([name])
    for fasta in sequencefiles:
        if fasta[:4] == name:
            group[-1].append(fasta)
for item in group:
    del item[0]
print 'Found '+ str(len(uniqueseqnames))+' unique structures.'
#TODO: Problem is still that there are structures for which the chains are spread over multiple bundle files.

# Next step: within each group, check each sequence against each other.
unique=[]
duplicates=[]
from copy import deepcopy
workgroup=deepcopy(group) # Just making sure we work on an actual copy of groups and not just a link (pointer) of group.
for g in range(0,len(workgroup)): # for every group
#for g in range(0,len(workgroup)): # for every group
    print 'Now acting on group #' + str(g) + ' of '+ str(len(workgroup)) + ' groups.'
    for h in range(len(workgroup[g])-1,-1,-1): #For every item within this group; note the very odd loop running in reverse.
        print 'Now checking item '+str(h)+' in that group against other items in that group'
        go=True
        for f in range(h-1,-1,-1): # against every other item within the same group, again running this loop from reverse.
            if go==True:
                #print 'Checking against item # '+str(f)+' in that group'
                if not check_homodimers_ASU('sequences/'+str(workgroup[g][h]),'sequences/'+str(workgroup[g][f])):
                    #print 'nohit'
                    pass
                else:
                    duplicates.append(workgroup[g].pop(h))
                    #print 'hit'
                    go=False
            else:
                break
# Now I have to make a new folder with one PDB file for each chain that is in the workgroup i.e. that is a unique chain.
# Probably not necessary: from shutil import copy2
#alluniquechains=list(itertools.chain(*workgroup))
try:
    os.mkdir('uniqueChains')
except:
    pass
#Build an array, which for every member of workgroup (which is also every PDB in my original PDB list), it just contains the PDB ID and the identifiers for the unique chains.
chains=[]
for g in range(0,len(workgroup)):
    if len(workgroup[g][0]) == 12: #This is the case if it is a normal PDB file with PDB ID _ chain . fasta.
        chains.append([workgroup[g][i][-7] for i in range(0,len(workgroup[g]))]) # Takes the 7th last character from each member of a workgroup and saves it. This is the chain identifier.
        chains[-1].insert(0,workgroup[g][0][0:-8]) # Adds the file name (w/o '.pdb') in front of each workgroup in order to keep name and chain IDs connected.
    else:
        chains.append([workgroup[g][i][-9:-6] for i in range(0,len(workgroup[g]))]) # This is the case if it comes from a bundle file. In that case I need to keep the number of the bundle file together with the chain ID.
        chains[-1].insert(0,workgroup[g][0][0:-9]) # as above, but this time without the .pdb but also without the bundle number directly left to .pdb.
#Next lines just to produce some output file containing this information since the above step might take a good deal of time and I don't want to repeat that all the time while troubleshooting.
outfile=open('PDBswithuniquechains.txt','w')
for line in chains:
    outfile.write(str(line)+'\n')
outfile.close()


#Need to run this going through chains and not through the pdbfiles, because the pdb-bundle files might contain identical chains.
for i in range(0,len(chains)): # Looping over all structures (might be spread over multiple pdb bundle files)
    if len(chains[i][1]) == 1: # This is the case for single-file structures with normal chain ID.
        input_pdb=open(chains[i][0]+'.pdb','r') # Loading the corresponding input_pdb file into an array.
        pdbarray=[]
        for line in input_pdb:
            if line[0:4] == 'ATOM' or line[0:6] == 'HETATM': pdbarray.append(line) # I only want the ATOM or HETATM lines.
        input_pdb.close()
        for k in range(1,len(chains[i])): # Loops over all unique chains in this structure.
            outfile=open('uniqueChains/'+chains[i][0]+'_'+chains[i][k]+'.pdb','w') # Make a file for each uniqiue chain.
            for line in pdbarray:
                if line[21] == chains[i][k]: # Go through full ATOM and HETATM array of this pdb file and find those with matching chain ID.
                    outfile.write(line) # Write every matching line into the new file.
            outfile.close()
    else: # This is the case if the structure comes with bundle files, over which all chains are spread with potentially duplicated chain ID's.
        # Here we need to loop the file opening and array formation over every unique chain already, because it also contains the bundle file number.
        for k in range(1,len(chains[i])): # Loops over all unique chains in this structure.
            input_pdb=open(chains[i][0]+chains[i][k][0]+'.pdb','r') # Loading the corresponding input_pdb file into an array.
            pdbarray=[]
            for line in input_pdb:
                if line[0:4] == 'ATOM' or line[0:6] == 'HETATM': pdbarray.append(line) # I only want the ATOM or HETATM lines.
            input_pdb.close()
            outfile=open('uniqueChains/'+chains[i][0]+'_'+chains[i][k]+'.pdb','w') # Make a file for each uniqiue chain.
            for line in pdbarray:
                if line[21] == chains[i][k][-1]: # Go through full ATOM and HETATM array of this pdb file and find those with matching chain ID.
                    outfile.write(line) # Write every matching line into the new file.
            outfile.close()


#TODO: Oho... this might go wildly wrong, if the clusters have chain IDs different from the chains surrounding them ... or what about interfacial clusters?!?!?! I'm actually really interested in those, but they will not come through bc of the two different chains around it.


"""
Old code... maybe useful at some point.
# Here I should combine all bundle files into one file and move the single bundle files to a different subfolder.
bundlefiles=[item for item in pdbfiles if 'bundle' in item]
uniquebundles=list(set([item[:4] for item in bundlefiles]))
for name in uniquebundles:
    bundlegroup=[]
    for pdbbundle in bundlefiles:
        if pdbbundle[:4] == name:
            bundlegroup.append(pdbbundle)
    #print bundlegroup
    outfile=open(name+'.pdb','w')
    for ffile in bundlegroup:
        inputfile=open(ffile,'r')
        for line in inputfile:
            outfile.write(line)
        inputfile.close()
    outfile.close()
#TODO:Wait... now we a file, which might contain several chains with the same name, because they are concatenated from different bundle files.
try:
    os.mkdir('single_bundlefiles')
except:
    pass
"""
"""from shutil import move
for item in bundlefiles:
    move(item,'single_bundlefiles/'+item)
"""
