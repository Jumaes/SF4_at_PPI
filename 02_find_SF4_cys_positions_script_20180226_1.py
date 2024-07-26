import sys, os
import pymol

#This script should take a bunch of pdb structures with SF4 clusters as input (e.g. all with that cluster from the PDB) and generate as output the positions of the other three cys coordinating that cluster relative to the first one.

# Step one: Find clusters coordinated by four Cys and create a new pdb file for each cluster with the immediate surrounding of the cluster.
# Load pdb into pymol, find out how many clusters. Check which of them are 4 cys.
# From pymolwiki copied the following start sequence for pymol from python:

# pymol environment
moddir='/usr/lib/python2.7/dist-packages/pymol'
sys.path.insert(0, moddir)
os.environ['PYMOL_PATH'] = os.path.join(moddir, 'pymol/pymol_path')

#making subfolders to keep everything neat and clean.
try:
    os.mkdir('clusterstructures')
except:
    pass
# Unfotunately starting several pymol instances after each other from one python instance is not really supported and always messes up.
# Means this has to get a bit messy, opening pymol and then looping within that pymol instance through all tasks.

# pymol launching: quiet (-q), without GUI (-c) and with arguments from command line
pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
pymol.finish_launching()
cmd = pymol.cmd
from pymol import stored
#This is really made more complicated by the fact that we have python modules working on variables, which are filled with values from pymol. E.g. chainlist or resilist contain values spit out by pymol, while then we do set() or min() on them. Therefore the weird "stored" module is needed to make them accesible to both python and pymol.
#DONE: We actually should kick out any chains, which are just NCS chains of the same protein. So somehow I would need to differntiate between a heterooligomer in a structure and a homooligomer. Maybe using their sequence or so?
# Quite a bit of work, but the 01_pdb2fasta.... program does that.
#DONE: Along the same lines, we should try to get a database of structures, which is devoid of duplications to avoid bias by those proteins, which have been crystallized and deposited in the PDB multiple times.
# New database is devoid of duplicates.

#This is the function that makes single pdb files of every SF4 with 4 cys around from an input pdb file.
def make_SF4_pdbs(inputPDB): #For some bizar reason, inputPDB has to be given as "6fd1.pdb" and not '6fd1.pdb'?!?!?!?!
    cmd.load('uniqueChains/'+inputPDB,'full') # a bit hairy here that I have to move the folder path into the function... but otherwise have output problems.
    # Since there can be multiple chains, which might or might not contain clusters at the same resi, we have to make an outer loop over the chanis while keeping the counter counting.
    # UPDATE: When running the 01_pdb2fast.... program first, this program here should run on pdb files, which already contain only one chain each. Thus we probably could make this easier... might save time. But keep it for now.
    failedclusters=[]
    counter=0 # Counter used to count the number of the sf4 of THIS PDB we are working on at the moment.
    stored.chainlist=[]
    cmd.iterate('resn sf4','stored.chainlist.append(chain)') #for every atom in an sf4, put the chain ID in a list. set(chainlist) will condense that list to unique values.
    if not stored.chainlist:
        failedclusters.append('PDB '+str(inputPDB)+' has no SF4.')
    for chain in set(stored.chainlist):
        stored.resilist=[]
        cmd.iterate('resn sf4 and chain '+str(chain),'stored.resilist.append(resi)') #Puts the resi of all atoms in SF4s in a list.
        number_of_sf4=len(set(stored.resilist)) # Counts the unique entries in that list.
        for i in set(stored.resilist):
            cmd.select('coordinating','byres (resn cys and chain '+str(chain)+') within 3 of resi '+ str(i)) #Make selection of cys within 2.6angstrom of sf4. Hope that should be a good distance.
            stored.cyslist=[]
            cmd.iterate('coordinating and name sg','stored.cyslist.append(resn + " " + resi)') # Lists the resn of all sg atoms (of cysteines) in the vicinity of SF4.
            stored.coordinatinglist=[]
            cmd.iterate('(byres chain '+str(chain)+' within 3 of resi '+ str(i)+') and name ca' ,'stored.coordinatinglist.append(resn)')
            if not len(set(stored.cyslist)) == 4: # Only if there are 4 Cys nearby is this sf4 interesting for us. If not, print these error messages. Otherwise go to next and save that substructure. Since I now print RESN RESI, I can look for unique entries with set(). Somehow he was always counting some cys twice.
                print 'PDB '+str(inputPDB)+' does not have four cys coordinating the cluster.'
                failedclusters.append('PDB '+str(inputPDB)+' coordinated by '+ str(stored.coordinatinglist)+' having only cysteins '+str(stored.cyslist))
                #DONE: I want a file containing these not conform cluster ideally with their actual coordinating residues.
                # These are now stored in failedclusters list within this function. Need to return this list for every new inputPDB and append to a general list.
                pass
            else:
                print str(inputPDB) + ': ok'
                counter=counter+1 # Need this counter to make a seperate structure for each cluster within one pdb file.
                cmd.create(inputPDB[:-4]+'_'+str(counter), 'byres all within 8 of (resi '+ str(i) + ' and chain '+str(chain)+' )' ) # Create an object of the cluster with 8A environment.
                #Does make things down the road a lot easier, if I  delete all cys from here, which are NOT coordinating the cluster as defined in the coordinating selection.
                #Make a selection first, which contains all cys in the newly created object, which are not within 3 of the sf4 cluster in that chain. Line after removes those atoms.
                cmd.select('toremove',inputPDB[:-4]+'_'+str(counter)+' and resn cys and (! (byres (resn cys and chain '+str(chain)+') within 3 of resi '+ str(i)+'))')
                cmd.remove('toremove')
                cmd.save('clusterstructures/'+inputPDB[:-4]+'_'+str(counter)+'.pdb', inputPDB[:-4]+'_'+str(counter)) # Save that object under PBDID_clusternumber.pdb
                failedclusters.append(inputPDB[:-4]+'_'+str(counter)+'.pdb created.')
            cmd.delete('coordinating') # Should delete this selection before the next cluster is processed. Hope this solves the problem with him finding 5 amino acids in cyslist with repetition of one AA.
    cmd.delete('*') #Deletes all atoms' and objects from pymol. Needed, because reinitialize doesn't do it and memory usage piles up until pymol runs out of RAM.
    #DONE: Delete doesn't do the trick. Still RAM usage piles up.
    # Workaround: Just process a batch of 400 structures or so at a time. Change the script and process the next 400...
    cmd.reinitialize() #instead of quitting pymol, we reinitialize and use that same thread of pymol for other tasks.
    #DONE: There is sth wrong with 7 structures, that leads to much larger and pymol unfriendly clusterfiles, which cannot be processed by the align function.
    # These seem to be large cryo EM structures. There are a decent amount of them. Also some solution NMR. Should I filter them out before? EXPDTA is the keyword.
    # Now added a function to filter out the xray structures only to avoid these problems.
    return failedclusters # Returns the clusters without four cys and their actual coordinating residues.

# Step two: Pick a structure and cystein #0, #1,#2 or #3 and define CA as origin of coordinate system with C and N aligned in some way. Apply this coordinate shift to the rest of the structure (maybe align actually).
# UPDATED: Now aligns not only the first cys to the origin and saves that structure, but actually loops over all 4 cys and makes a structure after each alignement to get also the other vectors in the end.
def align_cysbase(structure): #Also here for troubleshooting or so, always run with "6fd1.pdb" and not '6fd1.pdb'.
    print 'Now aligning '+str(structure)
    cmd.load(str(clusterfolder)+str(structure))
    stored.clusterres=[]
    cmd.iterate('resn cys and name sg within 3 of resn sf4','stored.clusterres.append(resi)')
    stored.clusterres.sort()
    #Now create a pseudoatom cys backbone that is nicely aligned to origin and the basis vectors.
    cmd.pseudoatom('oricys',resi='1',name='ca',pos='[0,0,0]')
    cmd.pseudoatom('oricys',resi='1',name='n',pos='[-1.45,0,0]')
    cmd.pseudoatom('oricys',resi='1',name='cb',pos='[0.34,1.476,0]')
    for cysposition in range(0,len(stored.clusterres)): # This is the new loop to ge the vectors between all Cys and not only cys #1 towards the other 3.
        cmd.select('thecys','resi '+ str(stored.clusterres[cysposition]))
        #Using pair_fit for the fit, because the more complicated fitting routines don't like this.
        cmd.pair_fit('thecys and name ca','oricys and name ca', 'thecys and name n', 'oricys and name n', 'thecys and name cb', 'oricys and name cb')
        cmd.save('alignedclusterstructures/'+str(structure[:-4])+'_'+str(cysposition)+'.pdb')
    cmd.delete('*') #Deletes all atoms' and objects from pymol. Needed, because reinitialize doesn't do it and memory usage piles up until pymol runs out of RAM.
    #DOONE: Delete doesn't do the trick. Still RAM usage piles up.
    # Workaround: Just process a batch of 400 structures or so at a time. Change the script and process the next 400...
    cmd.reinitialize()

#Following function makes lists of pdb files in a folder sorted by their method of structure solution given in the EXPDTA line in the file.
# Returns all these lists, so it should be used like this: xraylist, emlist, nmrlist, otherlist = sort_pdb_by_struturesolution_method('.')
def sort_pdb_by_struturesolution_method(folder):
    fulllist=os.listdir(folder)
    pdblist=[]
    for item in fulllist:
        if '.pdb' in item:
            pdblist.append(item)
    xraylist=[]
    emlist=[]
    nmrlist=[]
    otherlist=[]
    for item in pdblist:
        f=open(item,'r')
        for line in f:
            if 'EXPDTA' in line:
                if 'X-RAY DIFFRACTION' in line:
                    xraylist.append(item)
                elif 'ELECTRON MICROSCOPY' in line:
                    emlist.append(item)
                elif 'SOLUTION NMR' in line:
                    nmrlist.append(item)
                else:
                    otherlist.append(item)
        f.close()
    xraylist.sort()
    emlist.sort()
    nmrlist.sort()
    otherlist.sort()
    return xraylist, emlist, nmrlist, otherlist

# Getting list of structures solved by respective method so I can later loop only over the ones solved by Xray.
#xraylist, emlist, nmrlist, otherlist = sort_pdb_by_struturesolution_method('.')
#Above function not needed anymore, because new databse only contains those solved by xray anyways.

#Because of the piling up RAM usage, can't run all structures in one run. Instead need to manually run this script several times with limited range.
#Run 1 for i in range(0,400):
#Run 2 for i in range(400,len(xraylist)):
#UPDATE: Probably not necessary anymore since new database devoid of duplicates and weeding out of NCS copies should decrease the number of files and clusters significantly.
pdblist=[item for item in os.listdir('uniqueChains/') if '.pdb' in item]

failedclusters=[]
for i in range(0,len(pdblist)):
    print 'Making clusterstructures for unique chain ' + str(i) + ' of ' + str(len(pdblist)) + ' chains.'
    failedclusters.append(make_SF4_pdbs(pdblist[i])) # Weird way to write it, but it does write all the pdb files anyways and just gives as output stream any cluster that didn't have 4 cys.
f=open('cluster_without_4cys.txt','w')
for item in failedclusters:
    f.write(str(item) + '\n')
f.close()


clusterfolder='clusterstructures/'
clusterlist=os.listdir(clusterfolder)

try:
    os.mkdir('alignedclusterstructures')
except:
    pass
for i in range(0,len(clusterlist)):
    align_cysbase(clusterlist[i])



#DONE IN NEXT PROGRAM STEP: Write a file with CA atoms only and then another one with C and another one with N atoms only. Make a list of clusterfile name and running number and use that running number as Bfactor for the atoms to keep the information where they come from.
#DONE IN NEXT PROGRAM STEP: Use those files to map xyz for CA, N and C. To define clusters.


#Use lenght of CA-C and CA-N difference vectors to scale those vectors to 1 and then somehow calculate rotation function from that.
# Use numpy for that. Defining a vector as VECTOR=numpy.array()[...]). The numpy.linalg.norm(VECTOR) gives you the length. VECTOR/THATLENGTH normalizes it.

# Step three: Spit out C, CA and N positions of the other Cys in a 11-value tuple (PDBIDclusterID,clusterNumber, 3xC,3xCA,3xN)
# As output, write all those positions into a pdb file with same resi for all atoms within each cluster and resi being the clusterID number (which I guess I have to create for this).
# Display in pymol, find spots of higher occupancy etc...
