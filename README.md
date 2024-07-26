# SF4_at_PPI
Program to extract relative coordinates of Cys holding 4Fe4S clusters in protein structures and using that data to find AA positions in other proteins, where potentially an SF4 could fit at a protein-protein interface. Written as part of my research as a postdoc in the Tezcan lab at UCSD around spring of 2018.
It uses pymol via commandline to do some of the tasks (could likely be ported to pure python, but ...)

## General architecture
(Written in retrospect in 2024)
- Whole program is split into 6 (+1 now obsolete) steps, which are individual python files, which need to be run separately and in order.
- Understand steps 1-5 as "training" the model on the empirical structures with SF4 files and step 6 and 7 as using that knowledge to find positions on a new structure of interest. I.e. you need to run 1-5 only once and then repeat 6 and 7 for every structure, where you want to find a position to engineer an SF4 onto.
- Beforehand, one should have searched for and downloaded all pdb files containing the ligand "SF4", which is the code for 4Fe4S clusters. Probably you want to take those from a non-redundant subset of the PDB, potentially apply some filters for resolution or such. This is not part of the script!
- Steps:
    1. Looks through all pdb files in a folder, reads out the protein sequence for every chain, keeps only one copy of almost identical chains (homodimers) and writes a single pdb file for each chain.
    2. Goes through list of single-chain pdb files, selects the SF4 clusters and cys around them and creates new structures containing just those. For each cluster, generates 4 pdb files, where in each a different Cys is aligned to 0,0,0 witht he C-alpha and two defined positions with the N and C-beta atoms. 
    3. Goes through all pdb-files, which only contain SF4's and their Cys, and collects all lines from all structures for each atom type in a different pdb file. So you have coordinates of all Cys C-alpha atoms, all Cys N-atoms etc..
    Note: B-factors represent a unique "Clusternumber" so one can backtrack to the original PDB Structure.
    4. ( SEEMS MADE OBSOLETE NOW AND REPLACED BY 5) 
    Looks at every atom (not too far away and not close to origin) in those C-alpha and C-beta and N-pdb-files, counts neighboring atoms within 0.5 A distance and changes the B-factor to that number. Spits out new pdb-file with those b-factors.
    5. Divides the space around the ori in boxes of 0.5 A edgelength, then goes through the pdb-file with all C-alpha and for each box, counts the number of atoms in there. Then for each box of C-alpha atoms, checks the respective N-atoms, rejects outliers of more than 1 A away from their average and counts how many there are and where. Returns that as N-box for the C-alpha box. Same then for C-atoms and the Fe atoms. Spits out one pdb-file with atoms in the box positions and b-factor with the count of atoms in that box for visualization. But also a .txt file containing a running number for the box, coordinates of the all-lower corner of the box, number of atoms, etc.... this file is going to be used in the next step(s).
    6. This uses the files from (5) and an input pdb-file, where you want to find positions to put an SF4 in, let's you select a subset of residues you are interested in (if you want the SF4 to be in a specific area) and then loops over all those residues:
    - Align whole structure to have that residue at the origin
    - Make a 15 A sphere around that residue to limit secondary residues, which could be of potential interest
    - Loop over all seconday residues in questions:
        - Go through all CA boxes from (5) and check for each, if the CA atom of the seconday residue is in that box +- considerable padding to allow for non-perfect hits. If it's box has more than X hits:
            - Find the box of N atoms, which are associated with the box of CA atoms, we just hit and see if the secondary residue's N atom falls into that. If so... do the same with the C atom... 
            - If still matching, use the associated average Fe atoms positions and check in our query structure, if any of those is within X A (tyically 1 A) of any other atoms in that structure (Clash!)
            - if not, report success, print out a couple of messages, structures and stats....
    7. Uses hits from 6 and tries to find a homodimer, which would allow placement of a SF4 at the interface and will not lead to massive clashes (smaller clashes allowed here, going to manual inspection)