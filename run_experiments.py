from made_gauss_conditional_model import * 


import numpy as np
import pandas as pd

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

from pathlib import Path
from glob import glob
from matplotlib import pyplot as plt

from sparrow import Protein

import statsmodels.api as sm
from sklearn.preprocessing import PolynomialFeatures

num_residues = 24
mean_singular_value_thresh = 8
cov_singular_value_thresh = 8
num_bootstrap_samples = 1000 


def get_train_test_data(partition_by_feature, dataset_info_dict):

	conditional_str = dataset_info_dict['conditional_str']
	lower_val = dataset_info_dict['lower_val']
	upper_val = dataset_info_dict['upper_val']
	dataset_type = dataset_info_dict['dataset_type']

	if partition_by_feature:
		if conditional_str == 'ncpr':
			train_rows = ((seq_metadata_df['ncpr'] < lower_val)  | (seq_metadata_df['ncpr'] > upper_val))
			test_rows = ~train_rows
			train_rows_idx = seq_metadata_df.index[train_rows].tolist()
			test_rows_idx = seq_metadata_df.index[test_rows].tolist()
		elif conditional_str == 'kappa':
			train_rows = ((seq_metadata_df['kappa'] < lower_val)  | (seq_metadata_df['kappa'] > upper_val)) 
			test_rows = (~train_rows) 
			train_rows_idx = seq_metadata_df.index[train_rows].tolist()
			test_rows_idx = seq_metadata_df.index[test_rows].tolist()

	else:
	  num_seq = len(np.unique(np.array(all_mean_bootstrap_data.iloc[:, 2])))
	  num_bootstrap_rep = len(np.unique(np.array(all_mean_bootstrap_data.iloc[:, 0])))
	  train_rows_idx = np.arange(int(num_seq*.8)*num_bootstrap_rep)
	  test_rows_idx = np.arange(int(num_seq*.8)*num_bootstrap_rep,(num_seq*num_bootstrap_rep))
  
	
	num_seq = len(np.unique(np.array(all_mean_bootstrap_data.iloc[:, 2])))
	num_bootstrap_rep = len(np.unique(np.array(all_mean_bootstrap_data.iloc[:, 0])))


	if dataset_type == 'mean':
	  train_data = np.array(all_mean_bootstrap_data.iloc[train_rows_idx, 3:])
	  test_data = np.array(all_mean_bootstrap_data.iloc[test_rows_idx, 3:])
	  cov_singular_value_thresh = 0 
	elif dataset_type == 'cov':
	  train_data = np.array(all_cov_bootstrap_data.iloc[train_rows_idx, 3:])
	  test_data = np.array(all_cov_bootstrap_data.iloc[test_rows_idx, 3:])
	  mean_singular_value_thresh = 0
	elif dataset_type == 'mean_cov':
	  train_data = np.concatenate((np.array(all_mean_bootstrap_data.iloc[train_rows_idx, 3:]), np.array(all_cov_bootstrap_data.iloc[train_rows_idx, 3:])), axis=1)
	  test_data = np.concatenate((np.array(all_mean_bootstrap_data.iloc[test_rows_idx, 3:]), np.array(all_cov_bootstrap_data.iloc[test_rows_idx, 3:])), axis=1)

	conditional_train_data = np.array(seq_metadata_df.iloc[train_rows_idx, :])
	conditional_test_data = np.array(seq_metadata_df.iloc[test_rows_idx, :])

	labels = np.array(all_mean_bootstrap_data.iloc[:, 2])

	mean_train_data = np.mean(train_data, axis=0)
	std_train_data = np.std(train_data, axis=0)

	train_data = (train_data - mean_train_data)/std_train_data
	test_data = (test_data - mean_train_data)/std_train_data

	train_data = np.concatenate((train_data, conditional_train_data), axis=1)
	test_data = np.concatenate((test_data, conditional_test_data), axis=1)

	polynomial_features= PolynomialFeatures(degree=2)
	conditional_train_data_poly = polynomial_features.fit_transform(conditional_train_data)
	ols_model = sm.OLS(train_data[:,0], conditional_train_data_poly).fit() #predict s1 using ncpr/kappa

	return train_data, test_data, mean_train_data, std_train_data, ols_model




def eval_model(partition_by_feature, dataset_info_dict, use_conditional):

	print('partition by feature = %d' % partition_by_feature)
	print('use conditional = %d' % use_conditional)
	print(dataset_info_dict)
	
	conditional_str = dataset_info_dict['conditional_str']
	lower_val = dataset_info_dict['lower_val']
	upper_val = dataset_info_dict['upper_val']
	dataset_type = dataset_info_dict['dataset_type']


	train_data, test_data, mean_train_data, std_train_data, ols_model = get_train_test_data(partition_by_feature, dataset_info_dict)

	nin = train_data.shape[1]
	hiddens = [100, 100]
	sampling_method = 'gauss'
	num_conditional = 2

	if sampling_method == 'gauss':
	  num_components = None 
	  nout = 2*(nin-num_conditional)
	  
	train_loader = torch.utils.data.DataLoader(train_data, batch_size=100, shuffle=True)
	test_loader = torch.utils.data.DataLoader(test_data, batch_size=100)

	num_rep = 30
	num_epochs = 20
	test_losses_final_epoch = [] 

	for i in range(0,num_rep):
		model = MADE(nin, hiddens, nout, num_conditional, use_conditional, sampling_method, num_components)
		train_losses, test_losses = train_epochs(model, train_loader, test_loader, dict(epochs=num_epochs, lr=1e-3), quiet=False)
		test_losses_final_epoch.append(test_losses[-1])

	std_train_data_w_ones = np.append(std_train_data, np.ones(num_conditional))
	mean_train_data_w_zeros = np.append(mean_train_data, np.zeros(num_conditional))

	num_samples = 10000
	samples = model.sample_unconditional(num_samples, train_data).detach().numpy()
	samples_reconstructed_unconditional = samples*std_train_data_w_ones + mean_train_data_w_zeros

	
	ncpr_vals = ncpr_list_per_seq[np.random.choice(len(ncpr_list_per_seq), size=num_samples, replace=True)]
	kappa_vals = kappa_list_per_seq[np.random.choice(len(kappa_list_per_seq), size=num_samples, replace=True)]

	ncpr_vals = ncpr_vals[:,None]
	kappa_vals = kappa_vals[:,None]
	conditional_samples = np.concatenate((ncpr_vals, kappa_vals), axis=1)

	conditional_samples = torch.from_numpy(conditional_samples)
	samples = model.sample_conditional_ols(num_samples, train_data, conditional_samples, ols_model).detach().numpy()
	samples_reconstructed_conditional_all_data = samples*std_train_data_w_ones + mean_train_data_w_zeros



	if conditional_str == 'kappa':
	  ncpr_vals = ncpr_list_per_seq[np.random.choice(len(ncpr_list_per_seq), size=num_samples, replace=True)]
	  kappa_list_per_seq_subset = kappa_list_per_seq[(kappa_list_per_seq >= lower_val) & (kappa_list_per_seq <= upper_val)]
	  kappa_vals = kappa_list_per_seq_subset[np.random.choice(len(kappa_list_per_seq_subset), size=num_samples, replace=True)]
	  ncpr_vals = ncpr_vals[:,None]
	  kappa_vals = kappa_vals[:,None]
	elif conditional_str == 'ncpr':
	  ncpr_list_per_seq_subset = ncpr_list_per_seq[(ncpr_list_per_seq >= lower_val) & (ncpr_list_per_seq <= upper_val)]
	  ncpr_vals = ncpr_list_per_seq_subset[np.random.choice(len(ncpr_list_per_seq_subset), size=num_samples, replace=True)]
	  kappa_vals = kappa_list_per_seq[np.random.choice(len(kappa_list_per_seq), size=num_samples, replace=True)]
	  ncpr_vals = ncpr_vals[:,None]
	  kappa_vals = kappa_vals[:,None]


	conditional_samples = np.concatenate((ncpr_vals, kappa_vals), axis=1)
	conditional_samples = torch.from_numpy(conditional_samples)

	samples = model.sample_conditional_ols(num_samples, train_data, conditional_samples, ols_model).detach().numpy()
	samples_reconstructed_conditional_unseen_data = samples*std_train_data_w_ones + mean_train_data_w_zeros


	dataset_folder = '' 
	if partition_by_feature:
		dataset_folder = '%s_%.1f_%.1f' % (conditional_str, lower_val, upper_val)
	else:
		dataset_folder = '80_20_partition'
	
	if use_conditional:
		use_conditional_str = 'use_conditional=1'
	else:
		use_conditional_str = 'use_conditional=0'

	save_path = './model_output_data/%s/%s' % (dataset_folder, use_conditional_str)
	path = Path(save_path)
	path.mkdir(parents=True, exist_ok=True)
	print(save_path)


	colnames = [] 
	for i in range(1,mean_singular_value_thresh+1):
		colnames.append('S%d_mean' % i)
		
	if dataset_type == 'mean_cov':
		for i in range(1,cov_singular_value_thresh+1):
			colnames.append('S%d_cov' % i)
			
	colnames.extend(list(seq_metadata_df.columns))

	model_metadata = pd.DataFrame({'num_mean_singular_values': [mean_singular_value_thresh], 'num_cov_singular_values': [cov_singular_value_thresh], 'num_train_instances': [train_data.shape[0]], 'num_test_instances': [test_data.shape[0]], 'num_bootstrap_samples': [num_bootstrap_samples], 'dataset_info': [dataset_folder], 'use_conditional': [use_conditional]})
	model_metadata.to_csv('%s/made_metadata.csv' % save_path, index=False)

	
	test_losses_final_epoch_df = pd.DataFrame({'test_losses_final_epoch': test_losses_final_epoch})
	test_losses_final_epoch_df.to_csv('%s/test_losses_final_epoch.csv' % save_path, index=False)

	samples_reconstructed_df = pd.DataFrame(samples_reconstructed_unconditional, columns = colnames)
	samples_reconstructed_df.to_csv('%s/unconditional_samples.csv' % save_path, index=False)

	samples_reconstructed_df = pd.DataFrame(samples_reconstructed_conditional_all_data, columns = colnames)
	samples_reconstructed_df.to_csv('%s/conditional_all_data_samples.csv' % save_path, index=False)

	samples_reconstructed_df = pd.DataFrame(samples_reconstructed_conditional_unseen_data, columns = colnames)
	samples_reconstructed_df.to_csv('%s/conditional_unseen_data_samples.csv' % save_path, index=False)

	

model_input_data_path = './model_input_data'
all_mean_bootstrap_data = pd.read_csv('%s/all_mean_bootstrap_data.csv' % model_input_data_path)
all_cov_bootstrap_data = pd.read_csv('%s/all_cov_bootstrap_data.csv' % model_input_data_path)
seq_metadata_df = pd.read_csv('%s/seq_metadata_df.csv' % model_input_data_path)
per_seq_metadata_df = pd.read_csv('%s/per_seq_metadata_df.csv' % model_input_data_path)

ncpr_list_per_seq = np.array(per_seq_metadata_df['ncpr'])
kappa_list_per_seq = np.array(per_seq_metadata_df['kappa'])




print("A")
partition_by_feature = True

for use_conditional in [False, True]:

	for kappa_lower, kappa_upper in zip(list(np.arange(0,.8,.2)), list(np.arange(.2,1,.2))):

		dataset_info_dict = {}
		dataset_info_dict['conditional_str'] = 'kappa'
		dataset_info_dict['lower_val'] = round(kappa_lower,1)
		dataset_info_dict['upper_val'] = round(kappa_upper,1)
		dataset_info_dict['dataset_type'] = 'mean_cov'
		eval_model(partition_by_feature, dataset_info_dict, use_conditional)



print("B")
partition_by_feature = True

for use_conditional in [True, False]:

	for kappa_lower, kappa_upper in zip(list(np.arange(0,.9,.1)), list(np.arange(.1,1,.1))):

		dataset_info_dict = {}
		dataset_info_dict['conditional_str'] = 'kappa'
		dataset_info_dict['lower_val'] = round(kappa_lower,1)
		dataset_info_dict['upper_val'] = round(kappa_upper,1)
		dataset_info_dict['dataset_type'] = 'mean_cov'
		eval_model(partition_by_feature, dataset_info_dict, use_conditional)





print("C")
partition_by_feature = False

for use_conditional in [True, False]:
	dataset_info_dict = {}
	dataset_info_dict['conditional_str'] = 'NA'
	dataset_info_dict['lower_val'] = 'NA'
	dataset_info_dict['upper_val'] = 'NA'
	dataset_info_dict['dataset_type'] = 'mean_cov'
	eval_model(partition_by_feature, dataset_info_dict, use_conditional)


