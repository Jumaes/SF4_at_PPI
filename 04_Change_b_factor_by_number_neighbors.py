import pymol,os,sys
# This program will load a file allca.pdb in the same folder.
# Then calculate the number of close neighbors (0.5 A ) for all the atoms in there (when reasonably placed) and change b factors to number of neighbors.
moddir='/usr/lib/python2.7/dist-packages/pymol'
sys.path.insert(0, moddir)
os.environ['PYMOL_PATH'] = os.path.join(moddir, 'pymol/pymol_path')

pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
pymol.finish_launching()
cmd = pymol.cmd
from pymol import stored

cmd.load('allca.pdb')
stored.indices=[]
cmd.select('outlier','(x>20 or x<-20) or (y>20 or y<-20) or (z>20 or z<-20)') # Defines outliers, which are so far away from the origin, that they are most likely other cys not binding to the cluster.
cmd.select('core','(x<0.2 and x>-0.2) and (y<0.2 and y>-0.2) and (z<0.2 and z>-0.2)')
cmd.iterate('name ca and (! core ) and (! outlier)','stored.indices.append(index)') # Only those, which are not the origin cys (and thus super close to 0) and not outliers.
for item in stored.indices:
    stored.neighbors=[]
    cmd.iterate('all within 0.5 of index ' + str(item),'stored.neighbors.append(resi)') # For all those selected above, put all neighbors within 0.5 A in a list.
    cmd.alter('index '+ str(item),'b='+str(len(stored.neighbors))) # Check the length of list of neighbors and change the b factor of this atom to the number of its neighbors.
cmd.remove('core')
cmd.remove('outlier')
cmd.save('allca_neighbornumber_b.pdb') # Make a new pdb file, so that I can still look at the old one, in which the B factor represents the cluster number, which lets me look up the PDB ID it was coming from.
