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

bond_energies_no_CH0 = {	'CH0-CH3': 9999999.,
							'CH0-CH2': 9999999.,
							'CH0-CH1': 9999999.,
							'CH0-CH0': 9999999.,
							'CH1-CH3': 0.,
							'CH1-CH2': 0.,
							'CH1-CH1': 0.,
							'CH2-CH3': 0.,
							'CH2-CH2': 0.,
							'CH3-CH3': 0.	}

molecules_num = 10
temperature = 373.15
max_atom_num = 24
outpath = f'Generated_C{max_atom_num:d}/'
Path(outpath).mkdir(parents=True, exist_ok=True)
outnames = [outpath+f'{i:04d}' for i in range(molecules_num)]

for outname in outnames:
	bonds_num = generate_chain(outname, max_atom_num, temperature, bond_energies_equal)
	print(outname+f'\tBonds number = {bonds_num:d}')
