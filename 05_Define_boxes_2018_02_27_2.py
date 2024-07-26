import os,sys
import pymol
from numpy import mean
from math import ceil
moddir='/usr/lib/python2.7/dist-packages/pymol'
sys.path.insert(0, moddir)
os.environ['PYMOL_PATH'] = os.path.join(moddir, 'pymol/pymol_path')
from pymol import stored
pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
pymol.finish_launching()
cmd = pymol.cmd

#This program will use a allca.pdb and later alln.pdb and allc.pdb and allfe.pdb.
# It will spand a playing field of 20 x 20 x 20 A over around the middle and make boxes of 0.5 each edge.
# Every box will get a running number and a triple of its lower value corner.
# For every box it will count the number of CA atoms in there and give the sum as additional value.
# It will also note a multuple of unique tuple's for each atom consisting of the b factor (which now equals the cluster number) and the resi (because every cluster will have 4 cys).
# Finally those numbers will be used to look up boxes for the n atoms and c atoms and Fe atoms for those structures.


cmd.load('allca.pdb')
#Need to read in the file and define the box from there.. otherwise it becomes to big zif I apply ample padding.
stored.posx=[]
stored.posy=[]
stored.posz=[]
cmd.iterate_state('1','all','stored.posx.append(x)')
cmd.iterate_state('1','all','stored.posy.append(y)')
cmd.iterate_state('1','all','stored.posz.append(z)')
orix=ceil(mean([max(stored.posx),min(stored.posx)]))
oriy=ceil(mean([max(stored.posy),min(stored.posy)]))
oriz=ceil(mean([max(stored.posz),min(stored.posz)]))
lengthx=ceil(max(stored.posx)-min(stored.posx))+1
lengthy=ceil(max(stored.posy)-min(stored.posy))+1
lengthz=ceil(max(stored.posz)-min(stored.posz))+1
# length=20 # defines the length of the playing field edges in A. Since now I have dynamic field length in all dimensions to better accomodate the actual cluster, this is obsolete.
edge=0.5 # length of one edge of one box.
allboxes=int((lengthx*lengthy*lengthz)/(edge**3)) #calculate the number of boxes we will have.

caarray=['a']*allboxes # create the whole array of boxes with dummy.
startx=orix-(lengthx/2) #The playing field will span half of length in each direction from the origin of that direction.
starty=oriy-(lengthy/2)
startz=oriz-(lengthz/2)
counter=0 #internal counter for each box.
print 'Field: Length in X: %3d with origin at %3d. Length in Y: %3d with origin in %3d. Length in Z: %3d with origin in %3d. Boxedge %3.2f. Total of %6d boxes.' %(lengthx,orix,lengthy,oriy,lengthz,oriz,edge,allboxes)
for x in range(0,int(lengthx/edge)):
    for y in range(0,int(lengthy/edge)):
        for z in range(0,int(lengthz/edge)):
            print ' box number ' + str(counter) + ' of ' + str(allboxes)
            lowx=startx+edge*x #for each box define the lower x,y,z and the higher x,y,z
            highx=startx+edge*(x+1)
            lowy=starty+edge*y
            highy=starty+edge*(y+1)
            lowz=startz+edge*z
            highz=startz+edge*(z+1)
            stored.cas=[] # Need that stored ... function of pymol again to read/write values between pymol and python proper.
            # print stored.cas
            cmd.select('box','x>'+str(lowx)+' and x<'+str(highx)+' and y>'+str(lowy)+' and y<'+str(highy)+' and z>'+str(lowz)+' and z<'+str(highz)) # Select all atoms within the box boundaries.
            cmd.iterate('box','stored.cas.append([int(b),int(resi)])') # Write out the b factor of those atoms (contains the cluster ID) and the residue ID. Always 4 atoms will have the same cluster ID, because there are 4 cys.
            caarray[counter]=[counter, [lowx,lowy,lowz],len(stored.cas),stored.cas] # Write into the respective field of that box: Box number, lower corner of box, number of atoms in box, Bfactors and resi of each atom.
            cmd.delete('box') # delete selection to make sure nothing piles up.
            counter=counter+1 # count up and go to next box.

#Following code just to show as a pdb file kind of a heat map. Makes a large 64 000 atoms pdb file! A bit hard to load into pymol.
#TODO: Maybe better to load if those would not all be in one amino acid? Maybe unbonded types like H2O in different residues would be good so that he doesn't try to draw connections?
out=open('mapfile.pdb','w')
counter=0
for item in caarray:
    counter=counter+1
    #newline='ATOM  %5d  CA  GLY A 001   %8.3f%8.3f%8.3f  1.00%5.2f           C' %(counter,item[1][0],item[1][1],item[1][2],item[2]) # first attempt. Was writing every atom as one amino acid. Long time to load; maybe bc trying to draw bonds?
    newline='HETATM%5d  O   HOH A 001   %8.3f%8.3f%8.3f  1.00%5.2f           C' %(counter,item[1][0],item[1][1],item[1][2],item[2])
    out.write(newline+'\n')

out.close()

#This gives the actual .txt file we are going to use as input later.
outfile=open('caarray.txt','w')
for line in caarray:
    outfile.write(str(line)+'\n')

outfile.close()


cmd.remove('allca')
cmd.load('alln.pdb')
def find_nbox(box):
        #print 'Now working on box %4d' %(box[0])
        padding=0.5
        atoms=box[3] # This contains a tuple of b-factor(equals cluster ID) and resi for each CA atom in that box.
        cmd.select('nbox','none')
        for atom in atoms:
            cmd.select('nbox','nbox or (resi %4d and b=%5d)' %(atom[1],atom[0]))
        stored.xpos=[]
        stored.ypos=[]
        stored.zpos=[]
        cmd.iterate_state('1','nbox','stored.xpos.append(x)')
        cmd.iterate_state('1','nbox','stored.ypos.append(y)')
        cmd.iterate_state('1','nbox','stored.zpos.append(z)')
        xav=mean(stored.xpos)
        yav=mean(stored.ypos)
        zav=mean(stored.zpos)
        counter=0
        for i in range(len(stored.xpos)-1,-1,-1):
            if not ((xav-1)<stored.xpos[i]<(xav+1)):
                del stored.xpos[i]
                del stored.ypos[i]
                del stored.zpos[i]
                counter=counter+1
        #print '%10d outliers so far' %(counter)
        for i in range(len(stored.ypos)-1,-1,-1):
            if not ((yav-1)<stored.ypos[i]<(yav+1)):
                del stored.xpos[i]
                del stored.ypos[i]
                del stored.zpos[i]
                counter=counter+1
        #print '%10d outliers so far' %(counter)
        for i in range(len(stored.zpos)-1,-1,-1):
            if not ((zav-1)<stored.zpos[i]<(zav+1)):
                del stored.xpos[i]
                del stored.ypos[i]
                del stored.zpos[i]
                counter=counter+1
        #print '%10d outliers so far' %(counter)
        if not (counter == 0 or counter>box[2]):
            if float(counter)/(box[2]/float(100))>30:
                print 'Box: %4d: PROBLEM: More than 30 percent outliers rejected.' %(box[0])
                print '          %3d outliers of total of %3d atoms detected.' %(counter, box[2])
                xmax,xmin,ymax,ymin,zmax,zmin = [0,0,0,0,0,0]
                #TODO: If there are too many outliers, this should probably be reflected in the return output... right now that just keeps going.
            else:
                print 'Box: %4d: %3d outliers of total of %3d atoms detected.' %(box[0], counter, box[2])
        #else:
        #    print 'Box: %4d: No outliers of %3d atoms.' %(box[0],box[2])
        if not len(stored.xpos) == 0:
            xmax=max(stored.xpos)+padding
            xmin=min(stored.xpos)-padding
            ymax=max(stored.ypos)+padding
            ymin=min(stored.ypos)-padding
            zmax=max(stored.zpos)+padding
            zmin=min(stored.zpos)-padding
        else:
            xmax,xmin,ymax,ymin,zmax,zmin = [0,0,0,0,0,0]
        return [box[0],box[2],len(stored.xpos),xmax,xmin,ymax,ymin,zmax,zmin]

#Following also only works using state 0 in the iterate state command....
def find_febox(box): # Only needed because selection by b factor AND resi does not work for fe atoms, since the SF4 has a different resi ID.
        #print 'Now working on box %4d' %(box[0])
        padding=0.5
        atoms=box[3] # This contains a tuple of b-factor(equals cluster ID) and resi for each CA atom in that box.
        cmd.select('nbox','none')
        for atom in atoms:
            cmd.select('nbox','nbox or b=%5d' %(atom[0])) # Only difference really. Since there is only ONE SF4 fe average atom per aligned cluster, no need for resi and that would anyways be not the same as the one of the cys.
        stored.xpos=[]
        stored.ypos=[]
        stored.zpos=[]
        cmd.iterate_state('0','nbox','stored.xpos.append(x)') # For some reason this only works with state 0!!!!!
        cmd.iterate_state('0','nbox','stored.ypos.append(y)')
        cmd.iterate_state('0','nbox','stored.zpos.append(z)')
        #print 'Found %3d atoms.' %(len(stored.xpos))
        xav=mean(stored.xpos)
        yav=mean(stored.ypos)
        zav=mean(stored.zpos)
        counter=0
        for i in range(len(stored.xpos)-1,-1,-1):
            if not ((xav-1)<stored.xpos[i]<(xav+1)):
                del stored.xpos[i]
                del stored.ypos[i]
                del stored.zpos[i]
                counter=counter+1
        #print '%10d outliers so far' %(counter)
        for i in range(len(stored.ypos)-1,-1,-1):
            if not ((yav-1)<stored.ypos[i]<(yav+1)):
                del stored.xpos[i]
                del stored.ypos[i]
                del stored.zpos[i]
                counter=counter+1
        #print '%10d outliers so far' %(counter)
        for i in range(len(stored.zpos)-1,-1,-1):
            if not ((zav-1)<stored.zpos[i]<(zav+1)):
                del stored.xpos[i]
                del stored.ypos[i]
                del stored.zpos[i]
                counter=counter+1
        #print '%10d outliers so far' %(counter)
        if not (counter == 0 or counter>box[2]):
            if float(counter)/(box[2]/float(100))>30:
                print 'Box: %4d: PROBLEM: More than 30 percent outliers rejected.' %(box[0])
                print '          %3d outliers of total of %3d atoms detected.' %(counter, box[2])
                xmax,xmin,ymax,ymin,zmax,zmin = [0,0,0,0,0,0]
                #TODO: If there are too many outliers, this should probably be reflected in the return output... right now that just keeps going.
            else:
                print 'Box: %4d: %3d outliers of total of %3d atoms detected.' %(box[0], counter, box[2])
        #else:
        #    print 'Box: %4d: No outliers of %3d atoms.' %(box[0],box[2])
        if not len(stored.xpos) == 0:
            xmax=max(stored.xpos)+padding
            xmin=min(stored.xpos)-padding
            ymax=max(stored.ypos)+padding
            ymin=min(stored.ypos)-padding
            zmax=max(stored.zpos)+padding
            zmin=min(stored.zpos)-padding
        else:
            xmax,xmin,ymax,ymin,zmax,zmin = [0,0,0,0,0,0]
        return [box[0],box[2],len(stored.xpos),xmax,xmin,ymax,ymin,zmax,zmin]


narray=[]
for cabox in caarray:
    if not cabox[2]==0:
        narray.append(find_nbox(cabox))

cmd.remove('alln')

cmd.load('allc.pdb')
carray=[]
for cabox in caarray:
    if not cabox[2]==0:
        carray.append(find_nbox(cabox))
cmd.remove('allc')

cmd.load('allfe.pdb')
fearray=[]
for cabox in caarray:
    if not cabox[2]==0:
        fearray.append(find_febox(cabox))
cmd.remove('allfe')


outfile=open('carray.txt','w')
for line in carray:
    outfile.write(str(line)+'\n')

outfile.close()

outfile=open('narray.txt','w')
for line in narray:
    outfile.write(str(line)+'\n')

outfile.close()

outfile=open('fearray.txt','w')
for line in fearray:
    outfile.write(str(line)+'\n')

outfile.close()


"""
Following to load in an array file produced before.
nfile=open('narray.txt','r')
loaded=[]
for line in infile:
    loaded.append(line.strip('\n[]').split(','))
for i in range(0,len(loaded)):
    for k in range(0,len(loaded[0])):
        loaded[i][k]=float(loaded[i][k])
infile.close()
"""


"""
for item in caarray:

    if not item[2]==0:
            print item[0:3]

for item in caarray:
    if item[2]>10:
            print item[0:3]
"""
