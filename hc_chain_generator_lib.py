import copy
import numpy as np
from numpy.random import choice, seed

seed(2020)

Print_genlog = True

Particles = ['CH4','CH3','CH2','CH1','CH0']

def vec_to_str(v):
	s = ''
	for i in v:
		s += str(i)+','
	return s[:-1]	

def find_key(arr):
	type = arr[0]
	pos = arr[1]
	return vec_to_str(np.append(pos,type))

def try_to_add_bond(atoms,ak,j):
	directions = np.array([	[	[ 0, 0, 0, 4], [-1,-1, 0, 5], [-1, 0,-1, 6], [ 0,-1,-1, 7]	],
							[	[ 0, 0, 0, 4], [ 0, 0, 0, 5], [ 0, 0,-1, 6], [ 0, 0,-1, 7]	],
							[	[ 0, 0, 0, 4], [ 0,-1, 0, 5], [ 0, 0, 0, 6], [ 0,-1, 0, 7]	],
							[	[ 0, 0, 0, 4], [-1, 0, 0, 5], [-1, 0, 0, 6], [ 0, 0, 0, 7]	],
							[	[ 0, 0, 0, 0], [ 0, 0, 0, 1], [ 0, 0, 0, 2], [ 0, 0, 0, 3]	],
							[	[ 1, 1, 0, 0], [ 0, 0, 0, 1], [ 0, 1, 0, 2], [ 1, 0, 0, 3]	],
							[	[ 1, 0, 1, 0], [ 0, 0, 1, 1], [ 0, 0, 0, 2], [ 1, 0, 0, 3]	],
							[	[ 0, 1, 1, 0], [ 0, 0, 1, 1], [ 0, 1, 0, 2], [ 0, 0, 0, 3]	]	])				
	atom = copy.deepcopy(atoms[ak])
	type = atom[0]
	pos = atom[1]
	near_atoms_presence = atom[2]
	bonds_presence = atom[3]
	neigb = directions[atom[0]][j]
	new_type = neigb[3]
	new_pos = neigb[0:3]+atom[1]
	new_ak = find_key([new_type,new_pos])
	if new_ak in atoms.keys():
		new_atom = copy.deepcopy(atoms[new_ak])
		near_atoms_presence_for_new = new_atom[2]
		bonds_presence_for_new = new_atom[3]
		bonds_presence_for_new[type%4] = True
		near_atoms_presence[type%4] = True
		bonds_presence[new_type%4] = True
		new_atom = [new_type, new_pos, near_atoms_presence_for_new, bonds_presence_for_new]
		atom = [type, pos, near_atoms_presence, bonds_presence]
		new_bond = [ak, atom, new_ak, new_atom]
	else:
		near_atoms_presence_for_new = [False,False,False,False]
		bonds_presence_for_new = [False,False,False,False]
		bonds_presence_for_new[type%4] = True
		for i, d in enumerate(directions[new_type]):
			tak = find_key([d[3],d[0:3]+new_pos])
			if tak in atoms.keys():
				near_atoms_presence_for_new[i] = True
		near_atoms_presence[new_type%4] = True
		bonds_presence[new_type%4] = True
		new_atom = [new_type, new_pos, near_atoms_presence_for_new, bonds_presence_for_new]
		atom = [type, pos, near_atoms_presence, bonds_presence]
		new_bond = [ak, atom, new_ak, new_atom]
	return new_bond

def add_bond(atoms,new_bond):
	ak = new_bond[0]
	new_ak = new_bond[2]
	atom = new_bond[1]
	new_atom = new_bond[3]
	adds_atom = not new_ak in atoms.keys()
	atoms[ak] = atom
	atoms[new_ak] = new_atom
	return (ak, new_ak, adds_atom)

def convert_atoms_to_xyz(atoms):
	cell_size = 0.357
	type_shifts = [	[0.00, 0.00, 0.00],[0.50, 0.50, 0.00],[0.50, 0.00, 0.50],[0.00, 0.50, 0.50],
					[0.25, 0.25, 0.25],[0.75, 0.75, 0.25],[0.75, 0.25, 0.75],[0.25, 0.75, 0.75]	]			
	lines = []
	for ak in atoms.keys():
		type = atoms[ak][0]
		pos = atoms[ak][1]
		x = pos[0]+type_shifts[type][0]
		y = pos[1]+type_shifts[type][1]
		z = pos[2]+type_shifts[type][2]
		bonds_num = atoms[ak][3].count(True)
		particle_type = Particles[bonds_num]
		l = particle_type+'\t'+'%f\t'*3 % (x*cell_size,y*cell_size,z*cell_size)
		lines.append(l)
	return lines

def output(atoms, bonds_list, fout, foutbonds):
	lines = convert_atoms_to_xyz(atoms)
	fout.write(str(len(lines))+'\n')
	fout.write('\n')
	for l in lines:
		fout.write(l+'\n')
	for b in bonds_list:
		ii = list(atoms.keys()).index(b[0])
		jj = list(atoms.keys()).index(b[1])
		foutbonds.write(str(ii+1)+'\t'+str(jj+1)+'\n')
	
def calculate_energy(pb,temperature,bond_energies):
	beta = 1./0.00831446/temperature
	bonds_num = pb[1][3].count(True)
	bonds_num_new = pb[3][3].count(True)				
	if bonds_num>=bonds_num_new:
		bond_type = Particles[bonds_num]+'-'+Particles[bonds_num_new]
	else:
		bond_type = Particles[bonds_num_new]+'-'+Particles[bonds_num]
	try:
		u = bond_energies[bond_type]
	except KeyError:
		print('ERROR')
		print(pb)
		fout = open('error.xyz', 'w')
		output(atoms, fout)
		fout.close()
		exit()
	return u*beta, bond_type

def calculate_weights(possible_bonds,temperature,bond_energies):
	w = np.zeros(len(possible_bonds))
	bt = ['' for i in range(len(w))]
	for i, pb in enumerate(possible_bonds):
		v, bond_type = calculate_energy(pb,temperature,bond_energies)
		if v>10**6:
			w[i] = 0.
		else:
			w[i] = np.exp(-v)
		bt[i] = bond_type
	s = np.sum(w)
	w = w/s
	return w, bt

########
########
########

def generate_chain(outname, max_atom_num, temperature, bond_energies):
	initial_atom_type = 0
	initial_atom_position = np.array([0,0,0])
	atoms = {vec_to_str([0,0,0,0]) : [initial_atom_type, initial_atom_position, [False,False,False,False], [False,False,False,False]]}

	if Print_genlog:
		genlog = open(outname+'_genlog.txt', 'w')
	bonds_list = []
	while len(atoms.keys())<max_atom_num:
		possible_bonds = []
		for ak in atoms.keys():
			for j in range(4):
				if atoms[ak][3][j]==False:
					new_bond = try_to_add_bond(atoms,ak,j)
					possible_bonds.append(new_bond)
		weights, bond_types = calculate_weights(possible_bonds, temperature, bond_energies)
		bond_choice = choice(range(len(possible_bonds)), p=weights)
		b = add_bond(atoms,possible_bonds[bond_choice])
		bonds_list.append([possible_bonds[bond_choice][0],possible_bonds[bond_choice][2]])
		if Print_genlog:
			genlog.write(str(len(bonds_list))+'\n')
			genlog.write(str(weights)+'\n')
			genlog.write(str(bond_types)+'\n')
			genlog.write(str(b[0])+'\t'+str(b[1])+'\t'+str(b[2])+'\n')
			genlog.write(str(bond_types[bond_choice])+'\t'+str(weights[bond_choice])+'\t'+str(1./len(weights))+'\n')
			genlog.write('\n')
	if Print_genlog:	
		genlog.close()
	
	fout = open(outname+'.xyz', 'w')
	foutbonds = open(outname+'_bonds.txt', 'w')
	output(atoms, bonds_list, fout, foutbonds)
	fout.close()
	foutbonds.close()
	
	return(len(bonds_list))