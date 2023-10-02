from glob import glob
import Bio.PDB
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy
import numpy as np
z1=0
z2=0
zc1=0
zc2=0
box_cent=184.0
chains=['A','B','C','D','E','F','G','H','I','J','K']
def stitchy(dimer,trimer,N):
    start_id = 82
    end_id   = 102
    atoms_to_be_aligned = range(start_id, end_id + 1)
    atoms_to_be_aligned3 = range(start_id+237, end_id+237 + 1)
    atoms_to_be_aligned4 = range(start_id+237+237, end_id+237+237 + 1)
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
# Get the structures
    ref_structure = pdb_parser.get_structure("reference",trimer)
    ref_model    = ref_structure[0]

# Make a list of the atoms (in the structures) you wish to align.
# In this case we use CA atoms whose index is in the specified range
    ref_atoms = []
    ref_atoms3 = []
    ref_atoms5 = []
    for ref_chain in ref_model:
      atoms_to_be_aligned = range(start_id, end_id + 1)
  # Iterate of all residues in each model in order to find proper atoms
      for ref_res in ref_chain:
    # Check if residue number ( .get_id() ) is in the list
        if ref_res.get_id()[1] in atoms_to_be_aligned:
      # Append CA atom to list
            ref_atoms.append(ref_res['CA'])
        if ref_res.get_id()[1] in atoms_to_be_aligned3:
      # Append CA atom to list
            ref_atoms3.append(ref_res['CA'])
        if ref_res.get_id()[1] in atoms_to_be_aligned4:
      # Append CA atom to list
            ref_atoms5.append(ref_res['CA'])

    atoms_to_be_aligned = range(start_id, end_id + 1)
    atoms_to_be_aligned2 = range(start_id+237, end_id+237 + 1)
    #print(pdb)
# Start the parser
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
# Get the structures
    sample_structure = pdb_parser.get_structure("sample",dimer)

# Use the first model in the pdb-files for alignment
# Change the number 0 if you want to align to another structure
    sample_model = sample_structure[0]

# Make a list of the atoms (in the structures) you wish to align.
# In this case we use CA atoms whose index is in the specified range
    sample_atoms = []
    sample_atoms2=[]
# Iterate of all chains in the model in order to find all residues
# Do the same for the sample structure
    for sample_chain in sample_model:
      atoms_to_be_aligned = range(start_id, end_id + 1)
      for sample_res in sample_chain:
        if sample_res.get_id()[1] in atoms_to_be_aligned:
          sample_atoms.append(sample_res['CA'])
        if sample_res.get_id()[1] in atoms_to_be_aligned2:
          sample_atoms2.append(sample_res['CA'])

    l=len(sample_atoms)
    fixed_coord=numpy.zeros((l, 3))
    fixed_coord2=numpy.zeros((l, 3))
    moving_coord=numpy.zeros((l, 3))
    moving_coord2=numpy.zeros((l, 3))
    moving_coord3=numpy.zeros((l, 3))
    moving_coord4=numpy.zeros((l, 3))
    moving_coord5=numpy.zeros((l, 3))
    moving_coord6=numpy.zeros((l, 3))
    for i in range(0, len(sample_atoms)):
            fixed_coord[i]=sample_atoms[i].get_coord()
            fixed_coord2[i]=sample_atoms2[i].get_coord()
            moving_coord[i]=ref_atoms[i].get_coord()
            moving_coord3[i]=ref_atoms3[i].get_coord()
            moving_coord5[i]=ref_atoms5[i].get_coord()
# Now we initiate the superimposer:
    super_imposer = SVDSuperimposer()
    super_imposer.set(fixed_coord, moving_coord)
    super_imposer.run()
    #super_imposer.apply(sample_model.get_atoms())
    rms1=super_imposer.get_init_rms()
    if 10.0>super_imposer.get_init_rms():
        mergy(dimer,trimer,1,1,N)
        return
    super_imposer.set(fixed_coord, moving_coord3)
    super_imposer.run()
    if 10.0>super_imposer.get_init_rms():
        mergy(dimer,trimer,1,2,N)
        return
    super_imposer.set(fixed_coord, moving_coord5)
    super_imposer.run()
    if 10.0>super_imposer.get_init_rms():
        mergy(dimer,trimer,1,3,N)
        return
    #super_imposer.apply(sample_model.get_atoms())
    super_imposer = SVDSuperimposer()
    super_imposer.set(fixed_coord2, moving_coord)
    super_imposer.run()
    #super_imposer.apply(sample_model.get_atoms())
    if 10.0>super_imposer.get_init_rms():
        mergy(dimer,trimer,2,1,N)
        return
    super_imposer.set(fixed_coord2, moving_coord3)
    super_imposer.run()
    if 10.0>super_imposer.get_init_rms():
        mergy(dimer,trimer,2,2,N)
        return
    super_imposer.set(fixed_coord2, moving_coord5)
    super_imposer.run()
    if 10.0>super_imposer.get_init_rms():
        mergy(dimer,trimer,2,3,N)
        return

def mergy(pdbD,pdbT,dC,tC,N):
    cL=237
    xxx=np.zeros((1,3))
    chain=chains[N-1]
    fin=open(pdbD,'r')
    lline=fin.readline()
    fout=open('asu_%d.pdb'%N,'w')
    while len(lline)>0:
        if len(lline)>4 and lline[:4]=="ATOM":
            if (int(lline[22:26])<(93+(dC-1)*237)) and (int(lline[22:26])>((dC-1)*237)):
                fout.write(lline[:21]+chain+lline[22:])
        lline=fin.readline()
    fin.close()
    fin=open(pdbT,'r')
    lline=fin.readline()
    while len(lline)>0:
        if len(lline)>4 and lline[:4]=="ATOM":
            if (int(lline[22:26])>= (93+(tC-1)*237)) and (int(lline[22:26])<=(tC*237)):
                fout.write(lline[:21]+chain+lline[22:])
        lline=fin.readline()
    fin.close()
    fout.write('TER\n')
    fout.close()


stitchy('../D1.pdb','../T1.pdb',1)
stitchy('../D2.pdb','../T1.pdb',2)
stitchy('../D1.pdb','../T1b.pdb',3)
