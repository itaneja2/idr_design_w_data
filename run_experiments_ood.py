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
mean_singular_value_thresh = 4
cov_singular_value_thresh = 4
num_bootstrap_samples = 1000 


def get_train_test_data(dataset_type):
 
	if dataset_type == 'mean_cov':
		train_data = np.concatenate((np.array(all_mean_bootstrap_data_train.iloc[:, 3:(3+mean_singular_value_thresh)]), np.array(all_cov_bootstrap_data_train.iloc[:, 3:(3+mean_singular_value_thresh)])), axis=1)
		test_data = np.concatenate((np.array(all_mean_data_ood_test.iloc[:, 1:(1+mean_singular_value_thresh)]), np.array(all_cov_data_ood_test.iloc[:, 1:(1+mean_singular_value_thresh)])), axis=1)
	elif dataset_type == 'mean':
		train_data = np.array(all_mean_bootstrap_data_train.iloc[:, 3:(3+mean_singular_value_thresh)])
		test_data = np.array(all_mean_data_ood_test.iloc[:, 1:(1+mean_singular_value_thresh)])
	
	#test_data = np.concatenate((np.array(all_mean_data_ood_test.iloc[:, 3:]), np.array(all_cov_data_ood_test.iloc[:, 3:])), axis=1)

	conditional_train_data = np.array(seq_metadata_df_train)
	conditional_test_data = np.array(seq_metadata_df_ood_test)

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




def eval_model(use_conditional, dataset_type, dtype):

	print('use conditional = %d' % use_conditional)
	
	train_data, test_data, mean_train_data, std_train_data, ols_model = get_train_test_data(dataset_type)

	nin = train_data.shape[1]
	hiddens = [100, 100]
	sampling_method = 'gauss'
	num_conditional = 2

	if sampling_method == 'gauss':
	  num_components = None 
	  nout = 2*(nin-num_conditional)
	  
	train_loader = torch.utils.data.DataLoader(train_data, batch_size=100, shuffle=True)
	test_loader = torch.utils.data.DataLoader(test_data, batch_size=100)
		

	num_rep = 3
	num_epochs = 30
	test_losses_final_epoch = [] 

	for i in range(0,num_rep):
		model = MADE(nin, hiddens, nout, num_conditional, use_conditional, sampling_method, num_components)
		train_losses, test_losses = train_epochs(model, train_loader, test_loader, dict(epochs=num_epochs, lr=1e-3), quiet=False)
		test_losses_final_epoch.append(test_losses[-1])

	dataset_folder = 'unrestricted_%s_variants' % dtype 
	
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
	print(test_losses_final_epoch_df)
	

model_input_data_path = './model_input_data'
all_mean_bootstrap_data_train = pd.read_csv('%s/all_mean_bootstrap_data.csv' % model_input_data_path)
all_cov_bootstrap_data_train = pd.read_csv('%s/all_cov_bootstrap_data.csv' % model_input_data_path)
seq_metadata_df_train = pd.read_csv('%s/seq_metadata_df.csv' % model_input_data_path)

dtype = 'hydrophobic' 
model_input_data_ood_path = './model_input_data_%s' % dtype
all_mean_data_ood_test = pd.read_csv('%s/all_mean_data.csv' % model_input_data_ood_path)
all_cov_data_ood_test = pd.read_csv('%s/all_cov_data.csv' % model_input_data_ood_path)
seq_metadata_df_ood_test = pd.read_csv('%s/seq_metadata_df.csv' % model_input_data_ood_path)


#all_mean_data_ood_test = all_mean_bootstrap_data_train.iloc[100:150,]
#all_cov_data_ood_test = all_cov_bootstrap_data_train.iloc[100:150,]
#seq_metadata_df_ood_test = seq_metadata_df_train.iloc[100:150,]

dataset_type = 'mean_cov'
use_conditional = True
eval_model(use_conditional, dataset_type, dtype)

