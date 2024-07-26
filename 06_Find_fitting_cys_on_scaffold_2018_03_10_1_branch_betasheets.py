import os,sys,pymol
from pymol import stored
moddir='/usr/lib/python2.7/dist-packages/pymol'
sys.path.insert(0, moddir)
os.environ['PYMOL_PATH'] = os.path.join(moddir)
pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
pymol.finish_launching()
cmd = pymol.cmd

#Parameter for the run:
minhits=0 # Only those boxes of CA atoms count as hits, which have MORE than 'minhits' SF4 template structure CA's in it. I.e. minhits=2 selects boxes with 3 or more hits.
edge=0.5 # This is the edge-length of the caboxes. Necessary, because only the coordinates of the lower corner are stored with the box. So you need to calculate the higher corner with the edge value.
capadding=1.6 # Padding, which is allowed for CA atom in question to miss the CA box in any direction. Note that multiple boxes might be found for each CA atom.
fepadding=1 # Padding around the "average" Fe atom of the SF4 used in the search for clashes of the new SF4 with the backbone of the scaffold structure.
#Since we are looking at the positions of the average Fe atom of the SF4 (kinda the center of the cluster), we need to allow for a padding to accomodate the actual cluster. Rather want to keep some bad solutions in, than kick good ones out. So keep it a bit too small.

# Following function used to load files with box information produced as arrays in program before.
def load_boxfile(infile):
    #Following to load in an array file produced before.
    f=open(infile,'r')
    loaded=[]
    for line in f:
        loaded.append(line.strip('\n[]').split(','))
    for i in range(0,len(loaded)):
        for k in range(0,len(loaded[0])):
            loaded[i][k]=float(loaded[i][k])
    f.close()
    return loaded

# Since complete insanity just to read in the list of lists that is the caarray.txt file and convert it into what it was before.
def load_caboxfile(infile):
    #Following to load in an array file produced before.
    f=open(infile,'r')
    loaded=[]
    for line in f:
        loaded.append(line.replace('[','').replace(']','').strip('\n').split(',')) # Need to completely flatten that line by removing all [] and then splitting into single items by ,.
    f.close()
    loaded2=[[int(item[0]),[float(item[1]),float(item[2]),float(item[3])],int(item[4]),item[5:]] for item in loaded] # This makes a list of lists out of the subitems 1-3, which are the x,y,z coordinates.
    for item in loaded2: # Now much fun with the list of tuples, which are the bfactor and resi of the atoms in that box.
        newlist=[]
        length=len(item[3]) # Gets number of subitems.
        if length>1: # If there are none, don't do the stuff below, because it will fail.
            for i in range(0,length,2):
                newlist.append([int(item[3][i]),int(item[3][i+1])]) # Every two subitems in item[3] are forming one tuple.
            item[3]=newlist
        else:
            pass
    return loaded2

#Following function to load the dictionary that can translate between the cluster ID in the bafactors and the actual name of that cluster pdb file.
def load_dictionary(infile):
    f=open(infile,'r')
    temp=[]
    for line in f:
        temp.append(line.rstrip('\n').split(':'))
    dictio={}
    for item in temp:
        dictio[int(item[0])]=item[1]
    return dictio

# Loading the search boxes first from one txt file for each atom type and loading the dictionary.
folder_with_boxfiles='/media/tezcanlab/Data/Julian/SF4project/newSF4_PDB_library/alignedclusterstructures/'
caarray=load_caboxfile(str(folder_with_boxfiles)+'caarray.txt')
carray=load_boxfile(str(folder_with_boxfiles)+'carray.txt')
narray=load_boxfile(str(folder_with_boxfiles)+'narray.txt')
fearray=load_boxfile(str(folder_with_boxfiles)+'fearray.txt')
dictionary=load_dictionary(str(folder_with_boxfiles)+'dictionary_clusterID_numbers.txt')
#outfilefolder='/media/tezcanlab/Data/Julian/SF4project/Azurinproject/Cupredoxin_PDB_Structures/Hits/'
#outfilefolder='Hits/' # Now dynamically defined for each scaffold in the loop further down.
structurefilefolder='/media/tezcanlab/Data/Julian/SF4project/newSF4_PDB_library/alignedclusterstructures/'
scaffoldlibfolder='/media/tezcanlab/Data/Julian/SF4project/Azurinproject/Cupredoxin_PDB_Structures/'


#Defining subarray of caarray: This will define how many CA atoms must be in a CA search box to count as a hit. Making this subarray will drastically speed up the search below.
casubarray=[]
for box in caarray:
    if box[2]>minhits: # Picking all searchboxes, which contain more than 'minhits' atoms with minhits defined in Parameter section on top.
        casubarray.append(box)

#This program should go through all CA atoms of a scaffold protein (or a subset of those) and
#check if there are any other CA atoms, which fall into one of the predfined 11 boxes for typical
# CA atoms within SF4 clusters.
# If so, it should check of the corresponding C and N atom of that AA falls also into the respective
# N or C box. I.e. if the backbone orientation is also correct.
# It should then spit out the information for each amino acid:
# Are there any CA atoms in boxes and if so which boxes and which CA
# Which of those other AA have also matching C and N atoms.



# Now using pymol for this and define a sphere of maybe 15 A around each scaffold CA first,
# and screen atoms only within this sphere to speed things up. This is because I first need to
# translate and especially rotate the scaffold atom, for which the search is done to the
# origin and orientation of the origin cys used to make the boxes.
#scaffold='1dur.pdb'
#cmd.load('../../newSF4_PDB_library/1dur.pdb')


#Define a function to take a scaffold structure and check all the CA atoms in it for potential second CA atoms fitting SF4.
# At this point pymol needs to be loaded. Also the casubarray needs to be constructed already.
def checking_scaffold_protein(scaffold):
    print '%s%s' %(scaffoldlibfolder, scaffold)
    cmd.load('%s%s' %(scaffoldlibfolder, scaffold))
    scaffold_name=scaffold[:-4]
    try: # make subfolder where all output structure files for this scaffold should be stored.
        os.mkdir(scaffold_name)
    except:
        pass
    outfilefolder=('%s/' %(scaffold_name))
    try: # make subfolder in subfolder to store the actually SF4 clashing structures to keep everything a bit more tidy.
        os.mkdir('%s/SF4clashes' %(scaffold_name))
    except:
        pass
    cmd.select('start_ca','name ca and (resi 1-6 or resi 29-37 or resi 61-69)')
    #cmd.select('start_ca','name ca') # Define the selection of atoms, which should be taken into consideration as first cys.
    #TODO: For some reason with 3sbq he is looking into amino acids 1-57 or so, even though there is not residue of that number and thus there should be no CA atoms...
    # Make list of of all CA atoms in scaffold protein with their chain/resi/, which should be taken into consideration for first cys.
    stored.start_ca_list=[]
    cmd.iterate('start_ca and name ca','stored.start_ca_list.append(chain+"/"+resi+"/")')

    # Make a dummy Cys backbone at origin and with the same orientation used to make the search boxes from the SF4 database.
    cmd.pseudoatom('oricys',resi='1',name='ca',pos='[0,0,0]')
    cmd.pseudoatom('oricys',resi='1',name='n',pos='[-1.45,0,0]')
    cmd.pseudoatom('oricys',resi='1',name='cb',pos='[0.34,1.476,0]')
    cmd.pseudoatom('oricys',resi='1',name='c',pos='[0.69,-0.70, 1.13]') # Introduced to fit gly's. Hope the fit is as good as with the older cb fit.

    reportfile=open('detailed_report_'+str(scaffold_name)+'.log','w') # All report files go into the folder in which this script is run.
    reportfile.write('Parameters of this search: \n')
    reportfile.write('Scaffold: %s. Minimum hits for a positive CA box: %d. Edgelength: %f Angstrom. Padding for CA box: %f angstrom. Padding for average Fe: %f angstrom. \n' %(scaffold, minhits, edge, capadding, fepadding))
    outfilecounter=0 #This will be used to number the pse files for hits. Just because there might be multiple hits for each resi A/ resi B pair depending on the padding.
    cahits=0 # Next lines are setting a number of reporting variables to 0.
    nhits=0
    chits=0
    clusterclashes=0
    hits=0
    for starting_ca in stored.start_ca_list:
        reportfile.write('#'*40+'\n Now looking to find partners for CA of '+str(starting_ca)+'. \n')
        cmd.select('starting_cys',scaffold_name+' and '+str(starting_ca))
        #Now align our starting cys to the dummy cys backbone at the origin.
        print 'Attempting fit to origin for amino acid '+str(starting_ca)
        # Originally used ca,n and cb. Doesn't work for gly. Therefore newer version below. cmd.pair_fit('starting_cys and name ca','oricys and name ca', 'starting_cys and name n', 'oricys and name n', 'starting_cys and name cb', 'oricys and name cb')
        cmd.pair_fit('starting_cys and name ca','oricys and name ca', 'starting_cys and name n', 'oricys and name n', 'starting_cys and name c', 'oricys and name c')
        #Make a selection of 15 A around the starting cys so that only those atoms have to be searched.
        cmd.select('15sphere','(name ca) within 15 of '+str(starting_ca))
        stored.target_list=[]
        cmd.iterate('15sphere','stored.target_list.append(chain+"/"+resi+"/")')
        for query_resi in stored.target_list: # Going through all CA within 15 A of start Cys CA
            stored.query_resi_coor=[]
            cmd.iterate_state('1',str(query_resi)+' and name ca','stored.query_resi_coor.append([x,y,z])') # Getting x,y,z for this CA atom.
            for box in casubarray: # Let's see if it is in any of our CA boxes, with enough CA atoms as defined above.
                if box[1][0]-capadding<stored.query_resi_coor[0][0]<box[1][0]+edge+capadding and box[1][1]-capadding<stored.query_resi_coor[0][1]<box[1][1]+edge+capadding and box[1][2]-capadding<stored.query_resi_coor[0][2]<box[1][2]+edge+capadding:
                    # The plus 0.5 value is needed, because the coordinates stored with the box are the ones of the lower value corner and the edgelength is 0.5.
                    #print 'Query atom %s sits in CA box %5d . %3d CA atoms residing there.' %(query_resi,box[0],box[2])
                    reportfile.write('Match for CA: For ' +str(starting_ca)+' CA of '+str(query_resi)+'fits into cabox '+str(box[0])+'.\n')
                    cahits=cahits+1
                    for nbox in narray:
                        if nbox[0]==box[0]: #Find the nbox, which is associated with the CA box, where we had a hit.
                            stored.query_resi_Ncoor=[]
                            cmd.iterate_state('1',str(query_resi)+' and name n','stored.query_resi_Ncoor.append([x,y,z])')
                            # Now compare the position of the query_resi N atom with that box.
                            if nbox[4]<stored.query_resi_Ncoor[0][0]<nbox[3] and nbox[6]<stored.query_resi_Ncoor[0][1]<nbox[5] and nbox[8]<stored.query_resi_Ncoor[0][2]<nbox[7]: #TODO: When running on 3sbq, this line throws a list index out of range error... dunno what is wrong there.
                                    #TODO: Integrate a possibility for N box padding in line above.
                                    #print 'Query residue %s N atom sits in N box %5d . %3d N atoms residing there.' %(query_resi,nbox[0],nbox[2])
                                    reportfile.write('Match for N: For ' +str(starting_ca)+' N of '+str(query_resi)+' fits into nbox '+str(nbox[0])+'.\n')
                                    nhits=nhits+1
                                    for cbox in carray:
                                        if cbox[0]==nbox[0]:
                                            stored.query_resi_Ccoor=[]
                                            cmd.iterate_state('1',str(query_resi)+' and name c','stored.query_resi_Ccoor.append([x,y,z])')
                                            if cbox[4]<stored.query_resi_Ccoor[0][0]<cbox[3] and cbox[6]<stored.query_resi_Ccoor[0][1]<cbox[5] and cbox[8]<stored.query_resi_Ccoor[0][2]<cbox[7]:
                                                #TODO: Integrate a possibility for Cbox padding in line above.
                                                reportfile.write('Match for C:For ' +str(starting_ca)+' C of '+str(query_resi)+' fits into cbox '+str(cbox[0])+'.\n')
                                                chits=chits+1
                                                structures=[dictionary[item[0]] for item in box[3]]
                                                reportfile.write('%2d atoms in this cabox from clusters %s \n' %(box[2],structures))
                                                for febox in fearray: # Now lets look at the average Fe atom of the SF4's, which are in that cabox.
                                                    if febox[0]==box[0]: # Looking for the febox that corresponds to the cabox.
                                                        xmin=febox[4]-fepadding # Calculate the min and max values from the febox values with the padding.
                                                        xmax=febox[3]+fepadding
                                                        ymin=febox[6]-fepadding
                                                        ymax=febox[5]+fepadding
                                                        zmin=febox[8]-fepadding
                                                        zmax=febox[7]+fepadding
                                                        #print febox
                                                        #print xmin, xmax, ymin, ymax, zmin, zmax
                                                        #print scaffold_name
                                                        cmd.select('clasher','%s and (name ca or name c or name n) and x>%s and x<%s and y>%s and y<%s and z>%s and z<%s' %(scaffold_name, xmin, xmax, ymin, ymax, zmin, zmax))
                                                        # Above selects all backbone residues of the scaffold within that min/max +- padding box of the fe.
                                                        if cmd.count_atoms('clasher')>0: # Only raise alarm if we find a backbone atom in there. Maybe wanna raise the value to 2 or 3 atoms... just 1 might be accomodated somehow.
                                                            stored.clasherres=[]
                                                            cmd.iterate('clasher','stored.clasherres.append(resi)')
                                                            reportfile.write('Found a clash of prospective SF4 position %.1f %.1f %.1f %.1f %.1f %.1f with backbone atoms from residues %s. \n' %(xmin, xmax, ymin, ymax, zmin, zmax, stored.clasherres))
                                                            clusterclashes=clusterclashes+1
                                                            # Don't completely trust the clash finder routine. Partly, because the fe-boxes are sometimes very ill defined.
                                                            # Therefore output a pse file also for the clashing ones, just with a tag in the name that it clashes.
                                                            for pdbfile in structures:
                                                                cmd.load(str(structurefilefolder)+str(pdbfile))
                                                            cmd.show('cartoon',str(scaffold_name))
                                                            cmd.label('(starting_cys or %s) and name ca' %(query_resi),'resi')
                                                            cmd.set('label_size','25')
                                                            cmd.set('label_position','1,1,3')
                                                            cmd.show('sticks','starting_cys or %s' %(query_resi))
                                                            cmd.hide('everything','!(oricys or %s) and !(resn cys or resn SF4)' %(scaffold_name))
                                                            cmd.color('red','!(oricys or %s) and !(resn SF4)'%(scaffold_name))
                                                            cmd.save('%sSF4clashes/%s_%s_%s_%d_clash.pse'%(outfilefolder,scaffold_name,starting_ca.replace('/',''),query_resi.replace('/',''), outfilecounter)) # Outfilecounter just to avoid overwriting files, if there are multiple hits for a resi A/B pair.
                                                            for pdbfile in structures:
                                                                cmd.delete(pdbfile[:-4])
                                                            cmd.hide('sticks','all')
                                                            cmd.label()
                                                            outfilecounter=outfilecounter+1
                                                        else:
                                                            reportfile.write('HIT! Didn t find clashes of backbone atoms with SF4 position %.1f %.1f %.1f %.1f %.1f %.1f. \n' %(xmin, xmax, ymin, ymax, zmin, zmax))
                                                            hits=hits+1
                                                            for pdbfile in structures:
                                                                cmd.load(str(structurefilefolder)+str(pdbfile))
                                                            cmd.show('cartoon',str(scaffold_name))
                                                            cmd.label('(starting_cys or %s) and name ca' %(query_resi),'resi')
                                                            cmd.set('label_size','25')
                                                            cmd.set('label_position','1,1,3')
                                                            cmd.show('sticks','starting_cys or %s' %(query_resi))
                                                            cmd.hide('everything','!(oricys or %s) and !(resn cys or resn SF4)' %(scaffold_name))
                                                            cmd.color('red','!(oricys or %s) and !(resn SF4)'%(scaffold_name))
                                                            cmd.save('%s%s_%s_%s_%d.pse'%(outfilefolder,scaffold_name,starting_ca.replace('/',''),query_resi.replace('/',''), outfilecounter)) # Outfilecounter just to avoid overwriting files, if there are multiple hits for a resi A/B pair.
                                                            for pdbfile in structures:
                                                                cmd.delete(pdbfile[:-4])
                                                            cmd.hide('sticks','all')
                                                            cmd.label()
                                                            outfilecounter=outfilecounter+1
                            '''else: # Mostly for debugging: Spits out a structure for all non-hits, which do NOT fit an Nbox. Just to check  manually if the filter works reasonably.
                                structures=[dictionary[item[0]] for item in box[3]]
                                for pdbfile in structures:
                                    cmd.load(str(structurefilefolder)+str(pdbfile))
                                cmd.show('cartoon',str(scaffold_name))
                                cmd.label('(starting_cys or %s) and name ca' %(query_resi),'resi')
                                cmd.set('label_size','25')
                                cmd.set('label_position','1,1,3')
                                cmd.show('sticks','starting_cys or %s' %(query_resi))
                                cmd.hide('everything','!(oricys or %s) and !(resn cys or resn SF4)' %(scaffold_name))
                                cmd.color('red','!(oricys or %s) and !(resn SF4)'%(scaffold_name))
                                cmd.save('%sNclash_%s_%s_%s_%d.pse'%(outfilefolder,scaffold_name,starting_ca.replace('/',''),query_resi.replace('/',''), outfilecounter)) # Outfilecounter just to avoid overwriting files, if there are multiple hits for a resi A/B pair.
                                for pdbfile in structures:
                                    cmd.delete(pdbfile[:-4])
                                cmd.hide('sticks','all')
                                cmd.label()
                                outfilecounter=outfilecounter+1'''
    reportfile.close()
    cmd.delete('*')
    return [cahits, nhits, chits, clusterclashes, hits]
                #else: # For debugging... spams screen otherwise.
            #    print 'Query atom %s not in any CA box with enough CA atoms to count as hit.' %(query_resi)

# Wanna loop over the different scaffold files all the while using that function from above.
#scaffold=[item for item in os.listdir(scaffoldlibfolder) if '.pdb' in item]
scaffold=['8paz.pdb']
for item in scaffold:
    repsumfile=open('reportsummary.log','a')
    hitarray=checking_scaffold_protein(item)
    repsumfile.write('Stats for scaffold: %s \n' %(item))
    repsumfile.write('Hits for CA: %d \n' %(hitarray[0]))
    repsumfile.write('Hits for N: %d \n' %(hitarray[1]))
    repsumfile.write('Hits for C: %d \n' %(hitarray[2]))
    repsumfile.write('Of these clashes of potential SF4 with backbone atoms: %d \n' %(hitarray[3]))
    repsumfile.write('Real hits: %d \n' %(hitarray[4]))
    repsumfile.write('#'*40)
    repsumfile.close()



"""for pdbfile in stuctures:
    cmd.load(str(structurefilefolder)+str(pdbfile))
cmd.label('(starting_cys or query_resi) and name ca'.'resi')
cmd.set('label_size','25')
cmd.show('sticks','starting_cys or query_resi')
cmd.hide('everything','!(oricys or %s) and !(resn cys or resn SF4)' %(scaffold_name))
cmd.color('red','!(oricys or %s) and !(resn SF4)'%(scaffold_name))
cmd.save('%s%s_%s_%s_%d.pse'%(outfilefolder,scaffold_name,starting_cys.replace('/',''),query_resi.replace('/',''), outfilecounter)) # Outfilecounter just to avoid overwriting files, if there are multiple hits for a resi A/B pair.
cmd.delete('!(oricys or %s)' %(scaffold_name))
outfilecounter=outfilecounter+1"""


"""
In case of a HIT! I want to see a result .pse file.
That means, take the scaffold_name in the orientation it is in (aligned to ori with startcys), label the start_cys and the query atom with their resi at CA, show sticks for these two AA, load the clusterstructures from the cabox, hide everything from them except the lines of their cys and SF4, color them red.
Save as a .pse in a subfolder.
Remove labels, remove all non-scaffold_name or oricys entries from pymol (hopefully this is not crowding the RAM too much...)
"""

"""for febox in fearray:
    fepadding=1
    if febox[0]==box[0]:
        xmin=febox[4]-fepadding
        xmax=febox[3]+fepadding
        ymin=febox[6]-fepadding
        ymax=febox[5]+fepadding
        zmin=febox[8]-fepadding
        zmax=febox[7]+fepadding
        cmd.select('clasher','%s and (name ca or name c or name n) and x>%s and x<%s and y>%s and y<%s and z>%s and z<%s' %(scaffold, xmin, xmax, ymin, ymax, zmin, zmax))
        if cmd.count_atoms('clasher')>0:
            stored.clasherres=[]
            cmd.iterate('clasher','stored.clasherres.append(resi)')
            reportfile.write('Found a clash of prospective SF4 position with backbone atoms from residues %s.' %(stored.clasherres))

xmin=febox[4]-fepadding
xmax=febox[3]+fepadding
ymin=febox[6]-fepadding
ymax=febox[5]+fepadding
zmin=febox[8]-fepadding
zmax=febox[7]+fepadding
"""

"""

        for caboxnumber in range(0,len(cabox)): # Going through all targetboxes for CA's
             if boxquery(cabox[caboxnumber],str(query_resi)+' and name ca'):
                 reportfile.write('Match for CA: For ' +str(starting_ca)+' CA of '+str(query_resi)+'fits into cabox '+str(caboxnumber)+'.\n')
                 if boxquery(nbox[caboxnumber],str(query_resi)+' and name n'):
                     reportfile.write('Match for N: For ' +str(starting_ca)+' N of '+str(query_resi)+' fits into nbox '+str(caboxnumber)+'.\n')
                     if boxquery(cbox[caboxnumber],str(query_resi)+' and name c'):
                         reportfile.write('Match for C:For ' +str(starting_ca)+' C of '+str(query_resi)+' fits into cbox '+str(caboxnumber)+'.\n')
                     else: reportfile.write('C of ' +str(query_resi)+'does not fit into Cbox '+str(caboxnumber)+'.\n')
                 else: reportfile.write('N of ' +str(query_resi)+'does not fit into Nbox '+str(caboxnumber)+'.\n')
             else: reportfile.write('CA of ' +str(query_resi)+'does not fit into CAbox '+str(caboxnumber)+'.\n')
    reportfile.write('#'*40+'\n \n')

reportfile.close()

cmd.delete('*')"""
# Need to make some output from this. And need to test it!






# Load scaffold into pymol
# Create reference origin cys
# Make list of all CA atoms on scaffold, which should be searched
# Can I premake box definition searches, somehow?

#Loop over each CA atom in scaffold, which should be searched;
# Align to origin cys
# Make a subselection of 15 A around it with only CA, C and N atoms
# Make list of CA atoms in there.
# For every one by AA, check if it falls into one of the boxes
# If so, take the N atom of that AA and check if it falls into the same box.
# Then take C atom of that AA and check if that falls into the same box.

#Write report line by line on the go into a report file? Or how do I eventually like that to be reported?
# Maybe good also to have a digital list of hits, so I can do sth with them later.
