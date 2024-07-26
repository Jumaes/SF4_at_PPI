import os,sys,pymol
from pymol import stored
moddir='/usr/lib/python2.7/dist-packages/pymol'
sys.path.insert(0, moddir)
os.environ['PYMOL_PATH'] = os.path.join(moddir)
pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
pymol.finish_launching()
cmd = pymol.cmd

PSEFILE='1JOI_A96_A94_588.pse'
radius=2 #Radius, within which other CA, N or C atoms of the original atoms are searched.
clashcutoff=7 #Number of backbone atoms of the two structures, which are maximally allowed to be within 1.5 A of each other, before structure is counted as massive clash, which does not validate manual inspection.
try:
    os.mkdir('ChainBs')
except:
    pass
subfolder='ChainBs'

cmd.load(PSEFILE)
clusterobjects=[item for item in cmd.get_object_list('all') if '_' in item]
originalstruc=[item for item in cmd.get_object_list('all') if len(item)==4]
cmd.create('copystruc',originalstruc[0])
cmd.select('othercys','name ca and resn cys and %s and !( all within 2 of rep sticks)' %(clusterobjects[0]))
stored.othercysID=[]
cmd.iterate('othercys','stored.othercysID.append(resi)')
for firstcystein in stored.othercysID:
    print 'Looking at cystein %s in clusterstructure %s now.' %(firstcystein,clusterobjects[0])
    secondcystein=[item for item in stored.othercysID if not item==firstcystein]
    print 'Second cystein to be searched for is cystein %s in clusterstructure %s' %(secondcystein, clusterobjects[0])
    cmd.select('queryresidue', '%s and resi %s' %(clusterobjects[0],firstcystein))
    stored.othercyspos=[]
    atoms=['ca','n','c']
    for atom in atoms:
        cmd.iterate_state('1','queryresidue and name %s' %(atom),'stored.othercyspos.append([x,y,z])')
    cmd.pseudoatom('ori2cys',resi='1',name='ca',pos='%s'%(stored.othercyspos[0]))
    cmd.pseudoatom('ori2cys',resi='1',name='n',pos='%s'%(stored.othercyspos[1]))
    cmd.pseudoatom('ori2cys',resi='1',name='c',pos='%s'%(stored.othercyspos[2]))
    print 'Constructed alternative origin at position of cystein %s' %(firstcystein)
    stored.queryresiduelist=[]
    cmd.iterate('copystruc and name ca','stored.queryresiduelist.append(resi)')
    queryresidues=['copystruc//A/%s/' %(item) for item in stored.queryresiduelist]
    for residue in queryresidues:
        print 'Looking for match for residue %s.' %(residue)
        cmd.pair_fit(str(residue)+' and name ca','ori2cys and name ca', str(residue)+' and name n', 'ori2cys and name n', str(residue)+' and name c', 'ori2cys and name c')
        cmd.select('CAs','copystruc and name ca within %s of (%s and resi %s and name ca)' %(radius, clusterobjects[0], secondcystein[0]))
        if cmd.count_atoms('CAs')>0:
            print 'Found a fitting CA.'
            cmd.select('Ns','copystruc and name n within %s of (%s and resi %s and name n)' %(radius, clusterobjects[0], secondcystein[0]))
            if cmd.count_atoms('Ns')>0:
                print 'Found a fitting N.'
                cmd.select('Cs','copystruc and name c within %s of (%s and resi %s and name c)' %(radius,clusterobjects[0], secondcystein[0]))
                if cmd.count_atoms('Cs')>0:
                    stored.match=[]
                    cmd.iterate_state('1','Cs','stored.match.append(resi)')
                    print 'HEUREKA! Found also a fitting C atom. Residue %s has a partner in residue %s!' %(residue,stored.match[0])
                    cmd.select('clashes','(%s and backbone) within 1.5 of (copystruc and backbone)'%(originalstruc[0]))
                    clashscore=cmd.count_atoms('clashes')
                    if clashscore>clashcutoff:
                        print 'Unfortunately a massive clash between the two chains having %s backbone atoms within 1.5 A of each other.' %(clashscore)
                    else:
                        print 'LOOKS GOOD! Seems the two chains will not have a massive clash!'
                        cmd.color('marine','copystruc and name c*')
                        cmd.show('sticks','queryresidue or copystruc//A/%s/' %(secondcystein[0]))
                        cmd.save('%s/%s_%s_%s.pse'%(subfolder,PSEFILE[:-4],residue.split('/')[-2],stored.match[0]))
    cmd.remove('ori2cys')




'''
clusterobjects=[item for item in cmd.get_object_list('all') if '_' in item]
for clusobject in clusterobjects:
    cmd.select('othercys','name ca and resn cys and %s and !( all within 2 of rep sticks)' %(clusobject))
    stored.othercysID=[]
    cmd.iterate('othercys','stored.othercysID.append(resi)')
    for cystein in stored.othercysID:
        cmd.select('queryresidue', '%s and resi %s' %(clusobject,cystein))
        cmd.pair_fit('queryresidue and name ca','ori2cys and name ca', 'queryresidue and name n', 'ori2cys and name n', 'queryresidue and name c', 'ori2cys and name c')'''
