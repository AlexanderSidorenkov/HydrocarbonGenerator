from hc_chain_generator_lib import *
from pathlib import Path

bond_energies_equal = {	'CH0-CH3': 0.,
						'CH0-CH2': 0.,
						'CH0-CH1': 0.,
						'CH0-CH0': 0.,
						'CH1-CH3': 0.,
						'CH1-CH2': 0.,
						'CH1-CH1': 0.,
						'CH2-CH3': 0.,
						'CH2-CH2': 0.,
						'CH3-CH3': 0.	}

bond_energies_no_CH0 = {	'CH0-CH3': 1.e10,
							'CH0-CH2': 1.e10,
							'CH0-CH1': 1.e10,
							'CH0-CH0': 1.e10,
							'CH1-CH3': 0.,
							'CH1-CH2': 0.,
							'CH1-CH1': 0.,
							'CH2-CH3': 0.,
							'CH2-CH2': 0.,
							'CH3-CH3': 0.	}

bond_energies_dendrit = {	'CH0-CH3': 1.e10,
							'CH0-CH2': 1.e10,
							'CH0-CH1': 1.e10,
							'CH0-CH0': 1.e10,
							'CH1-CH3': 10.,
							'CH1-CH2': 1.e10,
							'CH1-CH1': 1.e10,
							'CH2-CH3': 0.,
							'CH2-CH2': 1.e10,
							'CH3-CH3': 0.	}

bond_energies_linear = {	'CH0-CH3': 1.e10,
							'CH0-CH2': 1.e10,
							'CH0-CH1': 1.e10,
							'CH0-CH0': 1.e10,
							'CH1-CH3': 1.e10,
							'CH1-CH2': 1.e10,
							'CH1-CH1': 1.e10,
							'CH2-CH3': 0.,
							'CH2-CH2': 1.e10,
							'CH3-CH3': 0.	}

bond_energies_archipelago = {	'CH0-CH3': 1.e10,
								'CH0-CH2': 1.e10,
								'CH0-CH1': 1.e10,
								'CH0-CH0': 1.e10,
								'CH1-CH3': 12.,
								'CH1-CH2': 3.,
								'CH1-CH1': 3.,
								'CH2-CH3': 0.,
								'CH2-CH2': 7.,
								'CH3-CH3': 0.	}


molecules_num = 1
temperature = 373.15
max_atom_num = 1000
twoD = False
outpath = f'Dendrit3D_Generated_C{max_atom_num:d}/'
Path(outpath).mkdir(parents=True, exist_ok=True)
outnames = [outpath+f'{i:04d}' for i in range(molecules_num)]

for outname in outnames:
	bonds_num = generate_chain(outname, max_atom_num, temperature, bond_energies_dendrit, twoD=twoD)
	print(outname+f'\tBonds number = {bonds_num:d}')
