library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS)
library(viridis)
library(reshape2)
library(data.table)
library(mvtnorm)

GS_ref_dir = paste('/home/ishan/Lasker/mpipi/GS_ref/num_rep=12',sep='/')
GS_pairwise_dist_output_save_dir = paste(GS_ref_dir, 'pairwise_dist_output', sep ='/')
GS_pairwise_dist_matrix = as.matrix(read.table(paste0(GS_pairwise_dist_output_save_dir, '/mean_pairwise_dist_matrix.csv')))
colnames(GS_pairwise_dist_matrix) = 1:nrow(GS_pairwise_dist_matrix)

pairwise_dist_cols = colnames(read.csv(paste(paste(GS_ref_dir,'rep_1', sep='/'), 'traj_analysis_data/pairwise_dist_matrix.csv', sep = '/'), check.names=FALSE))
pairwise_dist_cols = pairwise_dist_cols[2:length(pairwise_dist_cols)]

svd_output_save_dir = paste(GS_ref_dir, 'svd_output', sep ='/')
dir.create(svd_output_save_dir, recursive = TRUE, showWarnings = FALSE)

u_mean_pairwise_dist = as.matrix(read.table(paste0(svd_output_save_dir, '/mean_basis_u', '.csv')))
d_mean_pairwise_dist = as.matrix(read.table(paste0(svd_output_save_dir, '/mean_basis_d', '.csv')))
v_mean_pairwise_dist = as.matrix(read.table(paste0(svd_output_save_dir, '/mean_basis_v', '.csv')))

u_cov_pairwise_dist = as.matrix(read.table(paste0(svd_output_save_dir, '/cov_basis_u', '.csv')))
d_cov_pairwise_dist = as.matrix(read.table(paste0(svd_output_save_dir, '/cov_basis_d', '.csv')))
v_cov_pairwise_dist = as.matrix(read.table(paste0(svd_output_save_dir, '/cov_basis_v', '.csv')))





mean_singular_value_df_gsek = read.csv('/home/ishan/Lasker/mpipi/polyampholyte/num_charge=8/output_data/all_samples/mean_ij_singular_value_data.csv')
cov_singular_value_df_gsek = read.csv('/home/ishan/Lasker/mpipi/polyampholyte/num_charge=8/output_data/all_samples/cov_ij_singular_value_data.csv')

mean_singular_value_df_kappa_variants = read.csv('/home/ishan/Lasker/mpipi/polyampholyte_kappa_variants/num_charge=8/output_data/all_samples/mean_ij_singular_value_data.csv')
cov_singular_value_df_kappa_variants = read.csv('/home/ishan/Lasker/mpipi/polyampholyte_kappa_variants/num_charge=8/output_data/all_samples/cov_ij_singular_value_data.csv')

mean_singular_value_df_gsek$seq_id = mean_singular_value_df_gsek$full_seq
cov_singular_value_df_gsek$seq_id = cov_singular_value_df_gsek$full_seq

mean_singular_value_df_gsek = mean_singular_value_df_gsek[,c(1,2,3,4,7,8,9)]
cov_singular_value_df_gsek = cov_singular_value_df_gsek[,c(1,2,3,4,7,8,9)]

mean_singular_value_df_kappa_variants$seq_id = paste(mean_singular_value_df_kappa_variants$full_seq, mean_singular_value_df_kappa_variants$seq, sep='_')
cov_singular_value_df_kappa_variants$seq_id = paste(cov_singular_value_df_kappa_variants$full_seq, cov_singular_value_df_kappa_variants$seq, sep='_')

mean_singular_value_df_kappa_variants = mean_singular_value_df_kappa_variants[,c(1,2,3,4,8,9,10)]
cov_singular_value_df_kappa_variants = cov_singular_value_df_kappa_variants[,c(1,2,3,4,8,9,10)]

mean_singular_value_df_merged = bind_rows(mean_singular_value_df_gsek, mean_singular_value_df_kappa_variants)
cov_singular_value_df_merged = bind_rows(cov_singular_value_df_gsek, cov_singular_value_df_kappa_variants)

num_sequences = length(unique(mean_singular_value_df_merged$seq_id))


#convert to wide
mean_singular_value_df_merged$num_singular_values = paste0('S', mean_singular_value_df_merged$num_singular_values, '_mean')
mean_singular_value_df_merged = mean_singular_value_df_merged[,c('seq_id', 'ncpr', 'kappa', 'singular_value', 'num_singular_values')]
mean_singular_value_df_merged_wide = spread(mean_singular_value_df_merged, key = num_singular_values, value = singular_value) %>% as.data.frame()
mean_singular_value_df_merged_wide[mean_singular_value_df_merged_wide$kappa == -1, 'kappa'] = 0
mean_singular_value_df_merged_wide = mean_singular_value_df_merged_wide[,c('seq_id', 'ncpr', 'kappa', mean_singular_value_df_merged$num_singular_values[1:24])]

cov_singular_value_df_merged$num_singular_values = paste0('S', cov_singular_value_df_merged$num_singular_values, '_cov')
cov_singular_value_df_merged = cov_singular_value_df_merged[,c('seq_id', 'ncpr', 'kappa', 'singular_value', 'num_singular_values')]
cov_singular_value_df_merged_wide = spread(cov_singular_value_df_merged[,c('seq_id', 'ncpr', 'kappa', 'singular_value', 'num_singular_values')], key = num_singular_values, value = singular_value) %>% as.data.frame()
cov_singular_value_df_merged_wide[cov_singular_value_df_merged_wide$kappa == -1, 'kappa'] = 0
cov_singular_value_df_merged_wide = cov_singular_value_df_merged_wide[,c('seq_id', 'ncpr', 'kappa', cov_singular_value_df_merged$num_singular_values[1:276])]

tmp1 = mean_singular_value_df_merged_wide[,1:7]
tmp2 = cov_singular_value_df_merged_wide[,c(1,4:7)]
mean_cov_singular_value_df_merged_wide = left_join(tmp1, tmp2, by = 'seq_id') %>% as.data.frame()


mean_singular_value_df_merged_wide$dataset_info = 'ground_truth'
cov_singular_value_df_merged_wide$dataset_info = 'ground_truth'
mean_cov_singular_value_df_merged_wide$dataset_info = 'ground_truth'

mean_singular_value_df_merged_wide = subset(mean_singular_value_df_merged_wide, select = -seq_id)
cov_singular_value_df_merged_wide = subset(cov_singular_value_df_merged_wide, select = -seq_id)
mean_cov_singular_value_df_merged_wide = subset(mean_cov_singular_value_df_merged_wide, select = -seq_id)


mean_sd_ratio_ij_polyampholyte = read.csv('/home/ishan/Lasker/mpipi/polyampholyte/num_charge=8/output_data/all_samples/mean_sd_ratio_ij_data.csv')
mean_sd_ratio_ij_polyampholyte_kappa_variants = read.csv('/home/ishan/Lasker/mpipi/polyampholyte_kappa_variants/num_charge=8/output_data/all_samples/mean_sd_ratio_ij_data.csv')
mean_sd_ratio_ij_polyampholyte_kappa_variants$seq_w_kappa = paste(mean_sd_ratio_ij_polyampholyte_kappa_variants$seq, mean_sd_ratio_ij_polyampholyte_kappa_variants$kappa, sep='_')
mean_sd_ratio_ij_polyampholyte_kappa_variants = mean_sd_ratio_ij_polyampholyte_kappa_variants[,1:(ncol(mean_sd_ratio_ij_polyampholyte_kappa_variants))]


model_output_data_dir = '/home/ishan/Lasker/idr_design/model_output_data'
all_output_dirs = list.dirs(path = model_output_data_dir, full.names = TRUE, recursive = TRUE)

for (d in all_output_dirs)
{
  if ((length(str_split(d, '/')[[1]]) == 8) & (grepl('use_conditional=1', d)))
  {
    print(d)

    model_metadata = read.csv(paste(d, 'made_metadata.csv', sep='/'))
    conditional_all_data_samples = read.csv(paste(d, 'conditional_all_data_samples.csv', sep='/'))
    conditional_unseen_data_samples = read.csv(paste(d, 'conditional_unseen_data_samples.csv', sep='/'))
    
    dataset_info = model_metadata$dataset_info
    num_mean_singular_values = model_metadata$num_mean_singular_values
    num_cov_singular_values = model_metadata$num_cov_singular_values
    
    use_conditional =  model_metadata$use_conditional
    use_conditional_str = paste0('use_conditional=',use_conditional)
    
    if (grepl('kappa', dataset_info)) {
      kappa_lower_val = as.numeric(str_split(dataset_info, '_')[[1]][2])
      kappa_upper_val = as.numeric(str_split(dataset_info, '_')[[1]][3])
    } else {
      kappa_lower_val = 0
      kappa_upper_val = 1
    }
    
    #for plotting
    if (dataset_info == '80_20_partition') {
      dataset_info_str = 'rand_split'
    } else {
      dataset_info_str = dataset_info
    }
    
    conditional_all_data_samples$dataset_info = dataset_info_str
    conditional_unseen_data_samples$dataset_info = dataset_info_str

    mean_reconstruction_ratio_long_all_samples = data.frame()
    cov_reconstruction_long_all_samples = data.frame()
    
    for (i in 1:25)
    {
      mean_singular_values_pred = as.numeric(conditional_unseen_data_samples[i,grepl('mean', colnames(conditional_unseen_data_samples))])
      cov_singular_values_pred = as.numeric(conditional_unseen_data_samples[i,grepl('cov', colnames(conditional_unseen_data_samples))])
      
      mean_singular_values_combined = c(mean_singular_values_pred, d_mean_pairwise_dist[(num_mean_singular_values+1):length(d_mean_pairwise_dist)])
      cov_singular_values_combined = c(cov_singular_values_pred, d_cov_pairwise_dist[(num_cov_singular_values+1):length(d_cov_pairwise_dist)])
      
      mean_reconstruction = u_mean_pairwise_dist %*% diag(mean_singular_values_combined) %*% t(v_mean_pairwise_dist)
      cov_reconstruction = u_cov_pairwise_dist %*% diag(cov_singular_values_combined) %*% t(v_cov_pairwise_dist)
      
      mean_reconstruction_ratio_long = mean_reconstruction/GS_pairwise_dist_matrix
      mean_reconstruction_ratio_long[!is.finite(mean_reconstruction_ratio_long)] = 1
      
      mean_reconstruction_ratio_long = melt(mean_reconstruction_ratio_long) 
      colnames(mean_reconstruction_ratio_long) = c('atom_id1', 'atom_id2', 'mean')
      mean_reconstruction_ratio_long$sample_num = i 
      mean_reconstruction_ratio_long_all_samples = bind_rows(mean_reconstruction_ratio_long_all_samples, mean_reconstruction_ratio_long)
      
      colnames(cov_reconstruction) =  pairwise_dist_cols
      rownames(cov_reconstruction) = colnames(cov_reconstruction)
      
      cov_reconstruction_long = melt(cov_reconstruction) 
      colnames(cov_reconstruction_long) = c('atom_id1', 'atom_id2', 'cov')
      cov_reconstruction_long$sample_num = i 
      cov_reconstruction_long_all_samples = bind_rows(cov_reconstruction_long_all_samples, cov_reconstruction_long)
    }
    
    mean_reconstruction_ratio_long_all_samples$sample_num = factor(mean_reconstruction_ratio_long_all_samples$sample_num)
    mean_reconstruction_ratio_long_all_samples$atom_pair_id = paste(mean_reconstruction_ratio_long_all_samples$atom_id1, mean_reconstruction_ratio_long_all_samples$atom_id2, sep = '_')
    mean_reconstruction_ratio_long_all_samples$atom_id1 = as.numeric(mean_reconstruction_ratio_long_all_samples$atom_id1)
    mean_reconstruction_ratio_long_all_samples$atom_id2 = as.numeric(mean_reconstruction_ratio_long_all_samples$atom_id2)
    #upper triangular
    mean_reconstruction_ratio_long_all_samples = mean_reconstruction_ratio_long_all_samples %>% filter(atom_pair_id %in% unique(mean_sd_ratio_ij_polyampholyte_kappa_variants$atom_pair_id)) %>% as.data.frame()
    
    
    mean_sd_ratio_ij_polyampholyte_kappa_variants_subset = mean_sd_ratio_ij_polyampholyte_kappa_variants %>% filter(kappa >= kappa_lower_val) %>% filter(kappa <= kappa_upper_val) %>% as.data.frame()
    rel_seq = unique(mean_sd_ratio_ij_polyampholyte_kappa_variants_subset$seq_w_kappa)
    rel_seq_sample = rel_seq[sample(1:length(rel_seq),25)]
    mean_sd_ratio_ij_polyampholyte_kappa_variants_subset = mean_sd_ratio_ij_polyampholyte_kappa_variants_subset %>% filter(seq_w_kappa %in% rel_seq_sample) %>% as.data.frame()
    
    
    plot_save_dir = paste('/home/ishan/Lasker/idr_design/plots/generated_mean_ij/conditional_unseen_data', dataset_info, sep='/')
    dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)
    
    print(ggplot(mean_sd_ratio_ij_polyampholyte_kappa_variants_subset, aes(x = atom_id1, y = atom_id2, fill = mean)) + 
      geom_tile() + scale_fill_gradientn(colours = rainbow(5)) +
      facet_wrap(vars(seq_w_kappa)) + theme(text = element_text(size=20), axis.text.y = element_text(angle = 45), axis.text.x = element_text(angle = 45, hjust = 1)))
    ggsave(paste(plot_save_dir, paste0('mean_ij_actual_samples', '.png'), sep = '/'), width = 10, height = 10)
    
    print(ggplot(mean_reconstruction_ratio_long_all_samples, aes(x = atom_id1, y = atom_id2, fill = mean)) + 
      geom_tile() + scale_fill_gradientn(colours = rainbow(5), limits=c(min(mean_sd_ratio_ij_polyampholyte_kappa_variants_subset$mean), max(mean_sd_ratio_ij_polyampholyte_kappa_variants_subset$mean))) +
      facet_wrap(vars(sample_num)) + theme(text = element_text(size=20), axis.text.y = element_text(angle = 45), axis.text.x = element_text(angle = 45, hjust = 1)))
    ggsave(paste(plot_save_dir, paste0('mean_ij_gen_samples', '.png'), sep = '/'), width = 10, height = 10)
  }
}






    

#########simulated 2d pairwise distances#######


for (d in all_output_dirs)
{
  if ((length(str_split(d, '/')[[1]]) == 8) & (grepl('use_conditional=1', d)))
  {
    print(d)
    
    model_metadata = read.csv(paste(d, 'made_metadata.csv', sep='/'))
    conditional_all_data_samples = read.csv(paste(d, 'conditional_all_data_samples.csv', sep='/'))
    conditional_unseen_data_samples = read.csv(paste(d, 'conditional_unseen_data_samples.csv', sep='/'))
    
    dataset_info = model_metadata$dataset_info
    num_mean_singular_values = model_metadata$num_mean_singular_values
    num_cov_singular_values = model_metadata$num_cov_singular_values
    
    use_conditional =  model_metadata$use_conditional
    use_conditional_str = paste0('use_conditional=',use_conditional)
    
    if (grepl('kappa', dataset_info)) {
      kappa_lower_val = as.numeric(str_split(dataset_info, '_')[[1]][2])
      kappa_upper_val = as.numeric(str_split(dataset_info, '_')[[1]][3])
    } else {
      kappa_lower_val = 0
      kappa_upper_val = 1
    }
    
    #for plotting
    if (dataset_info == '80_20_partition') {
      dataset_info_str = 'rand_split'
    } else {
      dataset_info_str = dataset_info
    }
    
    conditional_all_data_samples$dataset_info = dataset_info_str
    conditional_unseen_data_samples$dataset_info = dataset_info_str
    
    pairwise_dist_cols_subset = pairwise_dist_cols[1:23]
    
 
    for (i in 1:10)
    {
      print(i)
      
      mean_singular_values_pred = as.numeric(conditional_unseen_data_samples[i,grepl('mean', colnames(conditional_unseen_data_samples))])
      cov_singular_values_pred = as.numeric(conditional_unseen_data_samples[i,grepl('cov', colnames(conditional_unseen_data_samples))])
      
      mean_singular_values_combined = c(mean_singular_values_pred, d_mean_pairwise_dist[(num_mean_singular_values+1):length(d_mean_pairwise_dist)])
      cov_singular_values_combined = c(cov_singular_values_pred, d_cov_pairwise_dist[(num_cov_singular_values+1):length(d_cov_pairwise_dist)])
      
      mean_reconstruction = u_mean_pairwise_dist %*% diag(mean_singular_values_combined) %*% t(v_mean_pairwise_dist)
      cov_reconstruction = u_cov_pairwise_dist %*% diag(cov_singular_values_combined) %*% t(v_cov_pairwise_dist)
      
      mean_reconstruction_long = melt(mean_reconstruction) 
      colnames(mean_reconstruction_long) = c('atom_id1', 'atom_id2', 'mean')
      mean_reconstruction_long$sample_num = i 
      mean_reconstruction_long$atom_pair_id = paste(mean_reconstruction_long$atom_id1, mean_reconstruction_long$atom_id2, sep='_')
      mean_reconstruction_long = mean_reconstruction_long %>% filter(atom_pair_id %in% unique(mean_sd_ratio_ij_polyampholyte_kappa_variants$atom_pair_id)) %>% as.data.frame()
      
      colnames(cov_reconstruction) =  pairwise_dist_cols
      rownames(cov_reconstruction) = colnames(cov_reconstruction)
      
      cov_reconstruction_long = melt(cov_reconstruction) 
      colnames(cov_reconstruction_long) = c('atom_id1', 'atom_id2', 'cov')
      cov_reconstruction_long$sample_num = i 
      
      simulated_2d_pairwise_distances = data.frame()
      
      for (j in 1:(length(pairwise_dist_cols_subset)-1))
      { 
        #print(paste('j=',j))
        for (k in (j+1):length(pairwise_dist_cols_subset))
        {
          atom_pair_id1 = as.character(pairwise_dist_cols[j])
          atom_pair_id2 = as.character(pairwise_dist_cols[k])
          
          mean_atom_id1 = as.numeric(mean_reconstruction_long[which(mean_reconstruction_long$atom_pair_id == atom_pair_id1),]['mean'])
          mean_atom_id2 = as.numeric(mean_reconstruction_long[which(mean_reconstruction_long$atom_pair_id == atom_pair_id2),]['mean'])
          
          cov_atom_id1 = as.numeric(cov_reconstruction_long[which((cov_reconstruction_long$atom_id1 == atom_pair_id1) & (cov_reconstruction_long$atom_id2 == atom_pair_id1)),]['cov'])
          cov_atom_id2 = as.numeric(cov_reconstruction_long[which((cov_reconstruction_long$atom_id1 == atom_pair_id2) & (cov_reconstruction_long$atom_id2 == atom_pair_id2)),]['cov'])
          cov_atom_id1_atom_id2 = as.numeric(cov_reconstruction_long[which((cov_reconstruction_long$atom_id1 == atom_pair_id1) & (cov_reconstruction_long$atom_id2 == atom_pair_id2)),]['cov'])
          
          cov_matrix = matrix(nrow=2,ncol=2)
          cov_matrix[1,1] =  cov_atom_id1      
          cov_matrix[2,2] =  cov_atom_id2         
          cov_matrix[1,2] =  cov_atom_id1_atom_id2         
          cov_matrix[2,1] =  cov_atom_id1_atom_id2         
          
          samples = rmvnorm(n=180, mean=c(mean_atom_id1,mean_atom_id2), sigma=cov_matrix)
          out = data.frame(dist1 = samples[,1], dist2 = samples[,2], atom_pair_id1 = atom_pair_id1, atom_pair_id2 = atom_pair_id2)
          simulated_2d_pairwise_distances = bind_rows(simulated_2d_pairwise_distances, out)
        }
      }
      
      simulated_2d_pairwise_distances$source = paste('gen_sample', i, sep='_')
      simulated_2d_pairwise_distances$atom_pair_merged = paste(simulated_2d_pairwise_distances$atom_pair_id1, simulated_2d_pairwise_distances$atom_pair_id2, sep='-')
      simulated_2d_pairwise_distances$atom_pair_id1 = factor(simulated_2d_pairwise_distances$atom_pair_id1, levels = pairwise_dist_cols_subset)
      simulated_2d_pairwise_distances$atom_pair_id2 = factor(simulated_2d_pairwise_distances$atom_pair_id2, levels = pairwise_dist_cols_subset)
      
      plot_save_dir = paste('/home/ishan/Lasker/idr_design/plots/simulated_2d_pairwise_distances/conditional_unseen_data', dataset_info, sep='/')
      dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)
      
      plot_name = paste('gen_sample', i, sep='_')
      print(ggplot(simulated_2d_pairwise_distances, aes(x = dist1, y = dist2)) + 
              geom_point(alpha=.05) +
              labs(title = paste('Generated Sample', i), x = 'Distance', y = 'Distance') + 
              scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) + 
              guides(colour = guide_legend(override.aes = list(alpha=1))) +
              facet_grid(atom_pair_id1 ~ atom_pair_id2) + theme(text = element_text(size=14), axis.text.x = element_text(angle = 45, hjust = 1)))
      ggsave(paste(plot_save_dir, paste0(plot_name, '.png'), sep = '/'), width=12, height=12)
      
    }  
  }
}


####actual sequence 2d pairwise distribution#######


mean_sd_ratio_ij_polyampholyte = read.csv('/home/ishan/Lasker/mpipi/polyampholyte/num_charge=8/output_data/all_samples/mean_sd_ratio_ij_data.csv')
mean_sd_ratio_ij_polyampholyte_kappa_variants = read.csv('/home/ishan/Lasker/mpipi/polyampholyte_kappa_variants/num_charge=8/output_data/all_samples/mean_sd_ratio_ij_data.csv')
mean_sd_ratio_ij_polyampholyte_kappa_variants$seq = paste(mean_sd_ratio_ij_polyampholyte_kappa_variants$seq, mean_sd_ratio_ij_polyampholyte_kappa_variants$kappa, sep='_')
mean_sd_ratio_ij_polyampholyte_all = bind_rows(mean_sd_ratio_ij_polyampholyte, mean_sd_ratio_ij_polyampholyte_kappa_variants)

pairwise_dist_cols_subset = pairwise_dist_cols[1:23]

atom_pair_merged_relevant = c()

for (j in 1:(length(pairwise_dist_cols_subset)-1))
{ 
  print(paste('j=',j))
  for (k in (j+1):length(pairwise_dist_cols_subset))
  {
    atom_pair_id1 = as.character(pairwise_dist_cols[j])
    atom_pair_id2 = as.character(pairwise_dist_cols[k])
    atom_pair_merged_relevant = c(atom_pair_merged_relevant, paste(atom_pair_id1, atom_pair_id2, sep='-'))
  }
} 

all_seq = unique(mean_sd_ratio_ij_polyampholyte_all$seq)

for (curr_seq in all_seq[seq(1, length(all_seq), 10)][27:196])
{
  print(curr_seq)
  
  mean_sd_ratio_ij_polyampholyte_all_subset = mean_sd_ratio_ij_polyampholyte_all %>% filter(seq == curr_seq) %>% as.data.frame()
  comparison_seq = mean_sd_ratio_ij_polyampholyte_all_subset[1,'seq']
  seq_str = comparison_seq
  
  if (grepl('_',comparison_seq)) {
    seq_wo_kappa = str_split(seq_str, '_')[[1]][1]
    kappa = str_split(seq_str, '_')[[1]][2]
    poly_type = 'polyampholyte_kappa_variant'
    comparison_seq_dir = paste('/media/ishan/UNTITLED/mpipi/polyampholyte_kappa_variants/num_charge=8', seq_wo_kappa, paste0('kappa=',kappa), sep='/')
  } else {
    poly_type = 'polyampholyte'
    comparison_seq_dir = paste('/media/ishan/UNTITLED/mpipi/polyampholyte/num_charge=8', seq_str, sep='/')
  }
  
  comparison_seq_pairwise_dist = data.frame()
  
  for (i in 1:3)
  {
    print(i)
    curr_dir = paste(comparison_seq_dir,paste0('rep_',i), sep='/')
    pairwise_dist = read.csv(paste(curr_dir, 'traj_analysis_data/pairwise_dist_matrix.csv', sep = '/'), check.names=FALSE)
    pairwise_dist$rep_num = i
    comparison_seq_pairwise_dist = bind_rows(comparison_seq_pairwise_dist, pairwise_dist)
  }
  
  
  pairwise_dist_long = comparison_seq_pairwise_dist %>% gather(atom_pair_id, dist, colnames(pairwise_dist)[2:(ncol(pairwise_dist)-1)]) %>% as.data.frame()
  pairwise_dist_long = pairwise_dist_long %>% filter(atom_pair_id %in% pairwise_dist_cols_subset) %>% as.data.frame()
  pairwise_dist_long = pairwise_dist_long %>% filter(frame_num %% 100 == 0) %>% as.data.frame()
  pairwise_dist_long = pairwise_dist_long %>% separate(col = atom_pair_id, into = c("atom_id1", "atom_id2"), sep = "_", remove=FALSE) %>% as.data.frame()
  
  pairwise_dist_long1 = pairwise_dist_long[,c(1,3,6)]
  pairwise_dist_long2 = pairwise_dist_long[,c(1,3,6)]
  
  colnames(pairwise_dist_long1) = c('frame_num', 'atom_pair_id1', 'dist1')
  colnames(pairwise_dist_long2) = c('frame_num', 'atom_pair_id2', 'dist2')
  
  comparison_seq_2d_pairwise_distances = pairwise_dist_long1 %>% left_join(pairwise_dist_long2, by = 'frame_num') %>% as.data.frame()
  comparison_seq_2d_pairwise_distances = comparison_seq_2d_pairwise_distances %>% filter(atom_pair_id1 != atom_pair_id2) %>% as.data.frame()
  comparison_seq_2d_pairwise_distances = comparison_seq_2d_pairwise_distances[,c('dist1', 'dist2', 'atom_pair_id1', 'atom_pair_id2')]
  
  comparison_seq_2d_pairwise_distances$seq = seq_str
  comparison_seq_2d_pairwise_distances$atom_pair_merged = paste(comparison_seq_2d_pairwise_distances$atom_pair_id1, comparison_seq_2d_pairwise_distances$atom_pair_id2, sep='-')
  comparison_seq_2d_pairwise_distances = comparison_seq_2d_pairwise_distances %>% filter(atom_pair_merged %in% atom_pair_merged_relevant) %>% as.data.frame()

  comparison_seq_2d_pairwise_distances$atom_pair_id1 = factor(comparison_seq_2d_pairwise_distances$atom_pair_id1, levels = pairwise_dist_cols_subset)
  comparison_seq_2d_pairwise_distances$atom_pair_id2 = factor(comparison_seq_2d_pairwise_distances$atom_pair_id2, levels = pairwise_dist_cols_subset)
  
  plot_save_dir = paste('/home/ishan/Lasker/idr_design/plots/actual_seq_2d_pairwise_distances', poly_type, sep='/')
  dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)
  
  plot_name = seq_str
  print(ggplot(comparison_seq_2d_pairwise_distances, aes(x = dist1, y = dist2)) + 
          geom_point(alpha=.01) +
          labs(title = seq_str, x = 'Distance', y = 'Distance') + 
          scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) + 
          guides(colour = guide_legend(override.aes = list(alpha=1))) +
          facet_grid(atom_pair_id1 ~ atom_pair_id2) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  ggsave(paste(plot_save_dir, paste0(plot_name, '.png'), sep = '/'), width=12, height=12)
  
}




