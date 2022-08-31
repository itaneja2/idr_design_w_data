from pathlib import Path
from glob import glob
import pandas as pd
import numpy as np
from sparrow import Protein
import os

print('getting polyampholyte')

home_dir = '/home/ishan/Lasker'
data_dir = sorted(glob('%s/mpipi/polyampholyte/num_charge=8/output_data/bootstrap/*' % home_dir, recursive = True))

all_mean_bootstrap_data_list = [] 
all_cov_bootstrap_data_list = [] 

num_residues = 24
mean_singular_value_thresh = 8
cov_singular_value_thresh = 8
num_bootstrap_samples = 1000 

for d in data_dir:
  
  #print(d)
	if os.path.exists(d):

		mean_ij_singular_values_df = pd.read_csv('%s/mean_ij_singular_value_data.csv' % d)
		cov_ij_singular_values_df = pd.read_csv('%s/cov_ij_singular_value_data.csv' % d)

		mean_ij_singular_values_df = mean_ij_singular_values_df[mean_ij_singular_values_df['num_bootstrap_samples'] == num_bootstrap_samples]
		mean_ij_singular_values_df = mean_ij_singular_values_df[mean_ij_singular_values_df['num_singular_values'] <= mean_singular_value_thresh]

		mean_ij_singular_values_df = mean_ij_singular_values_df[['num_singular_values', 'singular_value', 'bootstrap_rep', 'num_bootstrap_samples', 'seq']]

		mean_ij_singular_values_df = mean_ij_singular_values_df.pivot(['bootstrap_rep', 'num_bootstrap_samples', 'seq'], columns='num_singular_values', values='singular_value').reset_index()
		mean_ij_singular_values_df = mean_ij_singular_values_df.rename_axis(None, axis=1)
			
		all_mean_bootstrap_data_list.append(mean_ij_singular_values_df)

		cov_ij_singular_values_df = cov_ij_singular_values_df[cov_ij_singular_values_df['num_bootstrap_samples'] == num_bootstrap_samples]
		cov_ij_singular_values_df = cov_ij_singular_values_df[cov_ij_singular_values_df['num_singular_values'] <= cov_singular_value_thresh]

		cov_ij_singular_values_df = cov_ij_singular_values_df[['num_singular_values', 'singular_value', 'bootstrap_rep', 'num_bootstrap_samples', 'seq']]

		cov_ij_singular_values_df = cov_ij_singular_values_df.pivot(['bootstrap_rep', 'num_bootstrap_samples', 'seq'], columns='num_singular_values', values='singular_value').reset_index()
		cov_ij_singular_values_df = cov_ij_singular_values_df.rename_axis(None, axis=1)
			
		all_cov_bootstrap_data_list.append(cov_ij_singular_values_df)

  
all_mean_bootstrap_data = pd.concat(all_mean_bootstrap_data_list)
all_cov_bootstrap_data = pd.concat(all_cov_bootstrap_data_list)


#####

print('getting polyampholyte_kappa_variants')


home_dir = '/home/ishan/Lasker'
data_dir = sorted(glob('%s/mpipi/polyampholyte_kappa_variants/num_charge=8/output_data/bootstrap/*' % home_dir, recursive = True))
print(data_dir)


for d in data_dir:
  
  #print(d)
	if os.path.exists(d):

		mean_ij_singular_values_df = pd.read_csv('%s/mean_ij_singular_value_data.csv' % d)
		cov_ij_singular_values_df = pd.read_csv('%s/cov_ij_singular_value_data.csv' % d)


		mean_ij_singular_values_df['num_singular_values'] = np.tile(np.arange(1,num_residues+1), int(len(mean_ij_singular_values_df)/num_residues))
		cov_ij_singular_values_df['num_singular_values'] = np.tile(np.arange(1,(.5*num_residues*(num_residues-1))+1), int(len(cov_ij_singular_values_df)/(.5*num_residues*(num_residues-1))))


		mean_ij_singular_values_df = mean_ij_singular_values_df[mean_ij_singular_values_df['num_bootstrap_samples'] == num_bootstrap_samples]
		mean_ij_singular_values_df = mean_ij_singular_values_df[mean_ij_singular_values_df['num_singular_values'] <= mean_singular_value_thresh]

		mean_ij_singular_values_df = mean_ij_singular_values_df[['num_singular_values', 'singular_value', 'bootstrap_rep', 'num_bootstrap_samples', 'seq']]

		mean_ij_singular_values_df = mean_ij_singular_values_df.pivot(['bootstrap_rep', 'num_bootstrap_samples', 'seq'], columns='num_singular_values', values='singular_value').reset_index()
		mean_ij_singular_values_df = mean_ij_singular_values_df.rename_axis(None, axis=1)
			
		all_mean_bootstrap_data_list.append(mean_ij_singular_values_df)

		cov_ij_singular_values_df = cov_ij_singular_values_df[cov_ij_singular_values_df['num_bootstrap_samples'] == num_bootstrap_samples]
		cov_ij_singular_values_df = cov_ij_singular_values_df[cov_ij_singular_values_df['num_singular_values'] <= cov_singular_value_thresh]

		cov_ij_singular_values_df = cov_ij_singular_values_df[['num_singular_values', 'singular_value', 'bootstrap_rep', 'num_bootstrap_samples', 'seq']]

		cov_ij_singular_values_df = cov_ij_singular_values_df.pivot(['bootstrap_rep', 'num_bootstrap_samples', 'seq'], columns='num_singular_values', values='singular_value').reset_index()
		cov_ij_singular_values_df = cov_ij_singular_values_df.rename_axis(None, axis=1)

		all_cov_bootstrap_data_list.append(cov_ij_singular_values_df)

  
all_mean_bootstrap_data = pd.concat(all_mean_bootstrap_data_list)
all_cov_bootstrap_data = pd.concat(all_cov_bootstrap_data_list)



#####

print('getting seq metadata')


ncpr_list_all = []
hydrophobicity_list_all = []
kappa_list_all = []

ncpr_list_per_seq = []
hydrophobicity_list_per_seq = []
kappa_list_per_seq = []

all_seq = list(all_mean_bootstrap_data['seq'])

seq_dict = {}

for i,seq in enumerate(all_seq):
  
  if i % 1000 == 0:
    print(i)
    
  if '_' not in seq:
    
  
    full_seq = ''
    for s in seq:
        full_seq += 'GS%s' % s

    if len(full_seq) == 24:

      if full_seq in seq_dict:
        ncpr = seq_dict[full_seq]['ncpr']
        hydrophobicity = seq_dict[full_seq]['hydrophobicity']
        kappa = seq_dict[full_seq]['kappa']
      else:
        ncpr = round(Protein(full_seq).NCPR,2)
        hydrophobicity = round(Protein(full_seq).hydrophobicity,2)
        kappa = round(Protein(full_seq).kappa,2)

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
  else:
    
    seq_abbrev, kappa = seq.split('_')
    kappa_str = 'kappa=%s' % kappa
    path_to_seq_metadata = '/home/ishan/Lasker/mpipi/polyampholyte_kappa_variants/num_charge=8/%s/%s/seq_metadata.csv' % (seq_abbrev, kappa_str)
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

save_path = './model_input_data'
path = Path(save_path)
path.mkdir(parents=True, exist_ok=True)


all_mean_bootstrap_data.to_csv('%s/all_mean_bootstrap_data.csv' % save_path, index=False)
all_cov_bootstrap_data.to_csv('%s/all_cov_bootstrap_data.csv' % save_path, index=False)
seq_metadata_df.to_csv('%s/seq_metadata_df.csv' % save_path, index=False)
per_seq_metadata_df.to_csv('%s/per_seq_metadata_df.csv' % save_path, index=False)




