from glob import glob
import numpy
import Bio.PDB
from Bio.SVDSuperimposer import SVDSuperimposer
dimers=glob('../dimer_1/selects/snap*.pdb')
trimers=glob('../trimer/selects/snap*.pdb')
fout=open('logger1.log','w')
for trimer in trimers:
  print(trimer)
  flag=False
  start_id = 82
  end_id   = 102
  atoms_to_be_aligned = range(start_id, end_id + 1)
  atoms_to_be_aligned3 = range(start_id+237, end_id+237 + 1)
  atoms_to_be_aligned4 = range(start_id+237+237, end_id+237+237 + 1)
  pdb_parser = Bio.PDB.PDBParser(QUIET = True)
  ref_structure = pdb_parser.get_structure("reference",trimer)
  ref_model    = ref_structure[0]

  ref_atoms = []
  ref_atoms3 = []
  ref_atoms5 = []
  for ref_chain in ref_model:
      atoms_to_be_aligned = range(start_id, end_id + 1)
      for ref_res in ref_chain:
        if ref_res.get_id()[1] in atoms_to_be_aligned:
            ref_atoms.append(ref_res['CA'])
        if ref_res.get_id()[1] in atoms_to_be_aligned3:
            ref_atoms3.append(ref_res['CA'])
        if ref_res.get_id()[1] in atoms_to_be_aligned4:
            ref_atoms5.append(ref_res['CA'])
  for dimer in dimers: 
    start_id = 82
    end_id   = 102
    atoms_to_be_aligned = range(start_id, end_id + 1)
    atoms_to_be_aligned2 = range(start_id+237, end_id+237 + 1)
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    sample_structure = pdb_parser.get_structure("sample",dimer)

    sample_model = sample_structure[0]

    sample_atoms = []
    sample_atoms2=[]
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
    super_imposer = SVDSuperimposer()
    super_imposer.set(fixed_coord, moving_coord)
    super_imposer.run()
    rms1=super_imposer.get_init_rms()
    super_imposer.set(fixed_coord, moving_coord3)
    super_imposer.run()
    if rms1>super_imposer.get_init_rms():
        rms1=super_imposer.get_init_rms()
    super_imposer.set(fixed_coord, moving_coord5)
    super_imposer.run()
    if rms1>super_imposer.get_init_rms():
        rms1=super_imposer.get_init_rms()
    super_imposer = SVDSuperimposer()
    super_imposer.set(fixed_coord2, moving_coord)
    super_imposer.run()
    if rms1>super_imposer.get_init_rms():
        rms1=super_imposer.get_init_rms()
    super_imposer.set(fixed_coord2, moving_coord3)
    super_imposer.run()
    if rms1>super_imposer.get_init_rms():
        rms1=super_imposer.get_init_rms()
    super_imposer.set(fixed_coord2, moving_coord5)
    super_imposer.run()
    if rms1>super_imposer.get_init_rms():
        rms1=super_imposer.get_init_rms()
    if rms1<1.5:
      print(rms1,trimer,dimer)
      fout.write('%3.2f '%(rms1)+trimer+' '+dimer+'\n')
      flag=True
fout.close()

