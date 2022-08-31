from pathlib import Path
from glob import glob
import pandas as pd
import numpy as np
from sparrow import Protein
import os

dtype = 'hydrophobic'


all_mean_data_list = [] 
all_cov_data_list = [] 

num_residues = 24
mean_singular_value_thresh = 8
cov_singular_value_thresh = 8


if dtype == 'kappa':

	print('getting kappa variants')
	home_dir = '/home/ishan/Lasker'
	data_dir = sorted(glob('%s/mpipi/ood_constructs/kappa_variants/output_data/*' % home_dir, recursive = True))

	for d in data_dir:
		
		#print(d)
		if os.path.exists(d):
			
			mean_ij_singular_values_df = pd.read_csv('%s/mean_ij_singular_value_data.csv' % d)
			cov_ij_singular_values_df = pd.read_csv('%s/cov_ij_singular_value_data.csv' % d)

			mean_ij_singular_values_df = mean_ij_singular_values_df[mean_ij_singular_values_df['num_singular_values'] <= mean_singular_value_thresh]
			mean_ij_singular_values_df = mean_ij_singular_values_df[['num_singular_values', 'singular_value', 'seq']]

			mean_ij_singular_values_df = mean_ij_singular_values_df.pivot(['seq'], columns='num_singular_values', values='singular_value').reset_index()
			mean_ij_singular_values_df = mean_ij_singular_values_df.rename_axis(None, axis=1)
			#mean_ij_singular_values_df['source'] = 'kappa_variants'
				
			all_mean_data_list.append(mean_ij_singular_values_df)

			cov_ij_singular_values_df = cov_ij_singular_values_df[cov_ij_singular_values_df['num_singular_values'] <= cov_singular_value_thresh]

			cov_ij_singular_values_df = cov_ij_singular_values_df[['num_singular_values', 'singular_value', 'seq']]

			cov_ij_singular_values_df = cov_ij_singular_values_df.pivot(['seq'], columns='num_singular_values', values='singular_value').reset_index()
			cov_ij_singular_values_df = cov_ij_singular_values_df.rename_axis(None, axis=1)
			#cov_ij_singular_values_df['source'] = 'kappa_variants'
				
			all_cov_data_list.append(cov_ij_singular_values_df)

		
	all_mean_data = pd.concat(all_mean_data_list)
	all_cov_data = pd.concat(all_cov_data_list)


elif dtype == 'hydrophobic':

	print('getting hydrophobic variants')
	home_dir = '/home/ishan/Lasker'
	data_dir = sorted(glob('%s/mpipi/ood_constructs/hydrophobic_variants/output_data/*' % home_dir, recursive = True))

	for d in data_dir:
		
		#print(d)
		if os.path.exists(d):

			mean_ij_singular_values_df = pd.read_csv('%s/mean_ij_singular_value_data.csv' % d)
			cov_ij_singular_values_df = pd.read_csv('%s/cov_ij_singular_value_data.csv' % d)

			mean_ij_singular_values_df = mean_ij_singular_values_df[mean_ij_singular_values_df['num_singular_values'] <= mean_singular_value_thresh]
			mean_ij_singular_values_df = mean_ij_singular_values_df[['num_singular_values', 'singular_value', 'seq']]

			mean_ij_singular_values_df = mean_ij_singular_values_df.pivot(['seq'], columns='num_singular_values', values='singular_value').reset_index()
			mean_ij_singular_values_df = mean_ij_singular_values_df.rename_axis(None, axis=1)
			#mean_ij_singular_values_df['source'] = 'hydrophobic_variants'
				
			all_mean_data_list.append(mean_ij_singular_values_df)

			cov_ij_singular_values_df = cov_ij_singular_values_df[cov_ij_singular_values_df['num_singular_values'] <= cov_singular_value_thresh]

			cov_ij_singular_values_df = cov_ij_singular_values_df[['num_singular_values', 'singular_value', 'seq']]

			cov_ij_singular_values_df = cov_ij_singular_values_df.pivot(['seq'], columns='num_singular_values', values='singular_value').reset_index()
			cov_ij_singular_values_df = cov_ij_singular_values_df.rename_axis(None, axis=1)
			#cov_ij_singular_values_df['source'] = 'hydrophobic_variants'
			 
			all_cov_data_list.append(cov_ij_singular_values_df)

		
	all_mean_data = pd.concat(all_mean_data_list)
	all_cov_data = pd.concat(all_cov_data_list)





#####

print('getting seq metadata')

ncpr_list_all = []
hydrophobicity_list_all = []
kappa_list_all = []

ncpr_list_per_seq = []
hydrophobicity_list_per_seq = []
kappa_list_per_seq = []

all_seq = list(all_mean_data['seq'])

seq_dict = {}

for i,seq in enumerate(all_seq):
  

	if i % 1000 == 0:
		print(i)
			
	path_to_seq_metadata = '/media/ishan/UNTITLED/mpipi/ood_constructs/%s_variants/%s/seq_metadata.csv' % (dtype, seq)
	#print(path_to_seq_metadata)
	seq_metadata = pd.read_csv(path_to_seq_metadata)

	full_seq = str(seq_metadata['seq'])

	if full_seq in seq_dict:
		
		ncpr = seq_dict[full_seq]['ncpr']
		hydrophobicity = seq_dict[full_seq]['hydrophobicity']
		kappa = seq_dict[full_seq]['kappa']
	else:
		
		ncpr = float(seq_metadata['ncpr'])
		hydrophobicity = float(seq_metadata['hydrophobicity'])
		kappa = float(seq_metadata['kappa'])
		
		ncpr_list_per_seq.append(ncpr)
		hydrophobicity_list_per_seq.append(hydrophobicity)
		kappa_list_per_seq.append(kappa)
		
		seq_dict[full_seq] = {}
		seq_dict[full_seq]['ncpr'] = ncpr
		seq_dict[full_seq]['hydrophobicity'] = hydrophobicity
		seq_dict[full_seq]['kappa'] = kappa
		
	ncpr_list_all.append(ncpr)
	hydrophobicity_list_all.append(hydrophobicity)
	kappa_list_all.append(kappa)



ncpr_list_per_seq = np.array(ncpr_list_per_seq)
kappa_list_per_seq = np.array(kappa_list_per_seq)
hydrophobicity_list_per_seq = np.array(hydrophobicity_list_per_seq)

ncpr_list_all= np.array(ncpr_list_all)
kappa_list_all = np.array(kappa_list_all)
hydrophobicity_list_all = np.array(hydrophobicity_list_all)

kappa_list_all[kappa_list_all == -1] = 0 #ok for the model as it will ignore this 
kappa_list_per_seq[kappa_list_per_seq == -1] = 0 

seq_metadata_df = pd.DataFrame({'ncpr': ncpr_list_all, 'kappa': kappa_list_all})
per_seq_metadata_df = pd.DataFrame({'ncpr': ncpr_list_per_seq, 'kappa': kappa_list_per_seq})




###

print('saving data')

save_path = './model_input_data_%s' % dtype
path = Path(save_path)
path.mkdir(parents=True, exist_ok=True)


all_mean_data.to_csv('%s/all_mean_data.csv' % save_path, index=False)
all_cov_data.to_csv('%s/all_cov_data.csv' % save_path, index=False)
seq_metadata_df.to_csv('%s/seq_metadata_df.csv' % save_path, index=False)
per_seq_metadata_df.to_csv('%s/per_seq_metadata_df.csv' % save_path, index=False)




