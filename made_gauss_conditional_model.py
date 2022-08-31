import numpy as np 

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

import statsmodels.api as sm
from sklearn.preprocessing import PolynomialFeatures




class MaskedLinear(nn.Linear):
		""" same as Linear except has a configurable mask on the weights """

		def __init__(self, in_features, out_features, bias=True):
				super().__init__(in_features, out_features, bias)        
				self.register_buffer('mask', torch.ones(out_features, in_features))
				
		def set_mask(self, mask):        
				self.mask.data.copy_(torch.from_numpy(mask.astype(np.uint8).T))
				
		def forward(self, input):
				return F.linear(input, self.mask * self.weight, self.bias)

class MADE(nn.Module):
    
		def __init__(self, nin, hidden_sizes, nout, num_conditional, use_conditional, sampling_method, num_components, num_masks=1, natural_ordering=True):
				
				"""
				nin: integer; number of inputs
				hidden sizes: a list of integers; number of units in hidden layers
				nout: integer; number of outputs, which usually collectively parameterize some kind of 1D distribution
							note: if nout is e.g. 2x larger than nin (perhaps the mean and std), then the first nin
							will be all the means and the second nin will be stds. i.e. output dimensions depend on the
							same input dimensions in "chunks" and should be carefully decoded downstream appropriately.
							the output of running the tests for this file makes this a bit more clear with examples.
				num_masks: can be used to train ensemble over orderings/connections
				natural_ordering: force natural ordering of dimensions, don't use random permutations
				"""

				super().__init__()
				self.nin = nin
				self.num_conditional = num_conditional
				self.nout = nout
				self.hidden_sizes = hidden_sizes
				self.sampling_method = sampling_method 
				self.num_components = num_components 
				self.use_conditional = use_conditional
				assert self.nout % (self.nin-self.num_conditional) == 0, "nout must be integer multiple of nin"

				# define a simple MLP neural net
				self.net = []
				hs = [nin] + hidden_sizes + [nout]
				for h0,h1 in zip(hs, hs[1:]):
						self.net.extend([
										MaskedLinear(h0, h1),
										nn.ReLU()
								])
								
				self.net.pop() # pop ReLU

				self.net = nn.Sequential(*self.net)

				# seeds for orders/connectivities of the model ensemble
				self.natural_ordering = natural_ordering
				self.num_masks = num_masks
				self.seed = 0 # for cycling through num_masks orderings

				self.m = {}
				self.update_masks() # builds the initial self.m connectivity
				# note, we could also precompute the masks and cache them, but this
				# could get memory expensive for large number of masks.
					
		def update_masks(self):
      
				if self.m and self.num_masks == 1: return # only a single seed, skip for efficiency
				L = len(self.hidden_sizes)

				# fetch the next seed and construct a random stream
				rng = np.random.RandomState(self.seed)
				self.seed = (self.seed + 1) % self.num_masks

				# sample the order of the inputs and the connectivity of all neurons
				self.m[-1] = np.arange(self.nin) if self.natural_ordering else rng.permutation(self.nin)
				self.m['last'] = np.arange(self.nin-self.num_conditional) if self.natural_ordering else rng.permutation(self.nin)

				#if input is S1_m, S2_m, .. S1_c, S2_c, ..., NCPR, kappa
				#we only want to form connections from S1_m, S2_m, .. S1_c, .. S7_c (no point forming connections to S8_c as it can't be connected to anything in the final layer) 
				for l in range(L):
						self.m[l] = rng.randint(self.m[l-1].min(), self.nin-self.num_conditional-1, size=self.hidden_sizes[l]) 
							
				# construct the mask matrices
				# satisfy constraint that Di <= Dj
				masks = [self.m[l-1][:,None] <= self.m[l][None,:] for l in range(L)]

				#for conditional features, no restriction on their connections in the first layer 
				for i in range(1,(self.num_conditional+1)):
					masks[0][masks[0].shape[0]-i,:] = self.use_conditional
				 
				masks.append(self.m[L-1][:,None] < self.m['last'][None,:])

				# handle the case where nout = nin * k, for integer k > 1
				if self.nout >= self.nin:
						k = int(self.nout/(self.nin-self.num_conditional))
						# replicate the mask across the other outputs
						masks[-1] = np.concatenate([masks[-1]]*k, axis=1)

				# set the masks in all MaskedLinear layers
				layers = [l for l in self.net.modules() if isinstance(l, MaskedLinear)]
				for l,m in zip(layers, masks):
						l.set_mask(m)
		
		def forward(self, x):
			
				return self.net(x)
			
		def sample_unconditional(self, num_samples, train_data):
			
				samples = torch.zeros(num_samples, self.nin)
				order = self.m[-1]
				samples[:, order[0]] = torch.from_numpy(train_data[np.random.choice(train_data.shape[0], num_samples),order[0]])
				
				if self.sampling_method == 'gauss':
					for dim in order[1:(self.nin-self.num_conditional)]:
						out = self(samples)
						mu, log_std = torch.chunk(out, 2, dim=1)
						sample_output = torch.normal(mu[:,dim], torch.exp(log_std[:,dim]))
						samples[:, dim] = sample_output

				return(samples)
			
    
		#only based on kappa/ncpr
		def sample_conditional(self, num_samples, train_data, conditional_samples):
			
				samples = torch.zeros(num_samples, self.nin)
				order = self.m[-1]
				
				samples[:,(self.nin-self.num_conditional):] = conditional_samples

				if self.sampling_method == 'gauss':
					for dim in (order[0:(self.nin - self.num_conditional)]):
						out = self(samples)
						mu, log_std = torch.chunk(out, 2, dim=1)
						sample_output = torch.normal(mu[:,dim], torch.exp(log_std[:,dim]))
						samples[:, dim] = sample_output

				return(samples)

		#based on kappa/ncpr + s1 prediction
		def sample_conditional_ols(self, num_samples, train_data, conditional_samples, ols_model):
			
				samples = torch.zeros(num_samples, self.nin)
				order = self.m[-1]
							
				polynomial_features= PolynomialFeatures(degree=2)
				conditional_samples_poly = polynomial_features.fit_transform(conditional_samples)
				mean_s1_prediction = ols_model.get_prediction(conditional_samples_poly).summary_frame(alpha = .1) #predict s1 given kappa, ncpr
				
				samples[:, order[0]] = torch.from_numpy(np.random.uniform(mean_s1_prediction['obs_ci_lower'], mean_s1_prediction['obs_ci_upper'])) #sample using prediction interval      
				samples[:,(self.nin-self.num_conditional):] = conditional_samples

				if self.sampling_method == 'gauss':
					for dim in (order[1:(self.nin - self.num_conditional)]):
						out = self(samples)
						mu, log_std = torch.chunk(out, 2, dim=1)
						sample_output = torch.normal(mu[:,dim], torch.exp(log_std[:,dim]))
						samples[:, dim] = sample_output

				return(samples)
	




def calc_loss(model, x, out):
  
		x = x[:,0:(model.nin-model.num_conditional)]

		if model.sampling_method == 'gauss':
			loss = nn.GaussianNLLLoss(reduction='mean')
			mu, log_std = torch.chunk(out, 2, dim=1)
			nll_loss = loss(x, mu, torch.exp(log_std))

		return nll_loss
      
def train(model, train_loader, optimizer, epoch, grad_clip=None):
  
		model.train()
		train_losses = []
		for x in train_loader:
			out = model.forward(x.float())
			nll_loss = calc_loss(model, x, out)
			optimizer.zero_grad()
			nll_loss.backward()
			if grad_clip:
				torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
			optimizer.step()
			train_losses.append(nll_loss.item())
		return train_losses

def eval_loss(model, data_loader):

		model.eval()
		total_loss = 0
		with torch.no_grad():
			for x in data_loader:
				out = model.forward(x.float())
				nll_loss = calc_loss(model, x, out)
				total_loss += nll_loss * x.shape[0]
			avg_loss = total_loss / len(data_loader.dataset)

		return avg_loss.item()


def train_epochs(model, train_loader, test_loader, train_args, quiet):

		epochs, lr = train_args['epochs'], train_args['lr']
		grad_clip = train_args.get('grad_clip', None)
		optimizer = optim.Adam(model.parameters(), lr=lr)

		train_losses = []
		test_losses = [eval_loss(model, test_loader)]

		if not quiet:
				print('Initial test loss %.4f' % test_losses[0])
				
		for epoch in range(epochs):
			model.train()
			train_losses.extend(train(model, train_loader, optimizer, epoch, grad_clip))
			test_loss = eval_loss(model, test_loader)
			test_losses.append(test_loss)
			if not quiet:
				print(f'Epoch {epoch}, Test loss {test_loss:.4f}')

		return train_losses, test_losses
