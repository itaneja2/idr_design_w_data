library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS)
library(viridis)
library(reshape2)
library(GGally)
library(FNN)

#grid.draw.gg <- function(x){
#  print(x)
#}

#read in singular values for original data
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




model_output_data_dir = '/home/ishan/Lasker/idr_design/model_output_data'
all_output_dirs = list.dirs(path = model_output_data_dir, full.names = TRUE, recursive = TRUE)


#conditional_all_data_samples_all_models = data.frame()
#conditional_unseen_data_samples_all_models = data.frame()
#unconditional_samples_all_models = data.frame()

test_losses_final_epoch_all_models = data.frame()
kl_divergence_all_models = data.frame()
cor_abs_diff_all_models = data.frame()

for (d in all_output_dirs)
{
  
  if (length(str_split(d, '/')[[1]]) == 8)
  {
    print(d)
    while (!is.null(dev.list()))  { dev.off() }
    
    model_metadata = read.csv(paste(d, 'made_metadata.csv', sep='/'))
    test_losses_final_epoch = read.csv(paste(d, 'test_losses_final_epoch.csv', sep='/'))
    conditional_all_data_samples = read.csv(paste(d, 'conditional_all_data_samples.csv', sep='/'))
    conditional_unseen_data_samples = read.csv(paste(d, 'conditional_unseen_data_samples.csv', sep='/'))
    unconditional_samples = read.csv(paste(d, 'unconditional_samples.csv', sep='/'))
    
    conditional_all_data_samples = conditional_all_data_samples[1:num_sequences,]
    unconditional_samples = unconditional_samples[1:num_sequences,]
    
    dataset_info = model_metadata$dataset_info
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
    } else{
      dataset_info_str = dataset_info
    }
   
    test_losses_final_epoch$dataset_info = dataset_info_str
    test_losses_final_epoch$model_conditional_str = use_conditional
    
    test_losses_final_epoch_all_models = bind_rows(test_losses_final_epoch_all_models, test_losses_final_epoch)
    
    conditional_all_data_samples$dataset_info = dataset_info_str
    conditional_unseen_data_samples$dataset_info = dataset_info_str
    unconditional_samples$dataset_info = dataset_info_str
    
    #conditional_all_data_samples_all_models = bind_rows(conditional_all_data_samples_all_models, conditional_all_data_samples)
    #conditional_unseen_data_samples_all_models = bind_rows(conditional_unseen_data_samples_all_models, conditional_all_data_samples)
    #unconditional_samples_all_models = bind_rows(unconditional_samples_all_models, conditional_all_data_samples)
    
    
    conditional_all_data_samples_mean = conditional_all_data_samples[,c('ncpr', 'kappa', 'S1_mean', 'S2_mean', 'S3_mean', 'S4_mean', 'S5_mean', 'S6_mean', 'S7_mean', 'S8_mean', 'dataset_info')]
    conditional_all_data_samples_cov = conditional_all_data_samples[,c('ncpr', 'kappa', 'S1_cov', 'S2_cov', 'S3_cov', 'S4_cov', 'S5_cov', 'S6_cov', 'S7_cov', 'S8_cov', 'dataset_info')]
    conditional_all_data_samples_mean_cov = conditional_all_data_samples[,c('ncpr', 'kappa', 'S1_mean', 'S2_mean', 'S3_mean', 'S4_mean', 'S1_cov', 'S2_cov', 'S3_cov', 'S4_cov', 'dataset_info')]
    
    conditional_all_data_samples_mean_w_gt = bind_rows(conditional_all_data_samples_mean, mean_singular_value_df_merged_wide[,colnames(conditional_all_data_samples_mean)])
    conditional_all_data_samples_cov_w_gt = bind_rows(conditional_all_data_samples_cov, cov_singular_value_df_merged_wide[,colnames(conditional_all_data_samples_cov)])
    conditional_all_data_samples_mean_cov_w_gt = bind_rows(conditional_all_data_samples_mean_cov, mean_cov_singular_value_df_merged_wide[,colnames(conditional_all_data_samples_mean_cov)])
    
    plot_save_dir = paste('/home/ishan/Lasker/idr_design/plots/singular_value_pairwise_distribution/conditional_all_data', dataset_info, use_conditional_str, sep='/')
    dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)

    
    plot_name = 'mean_ij_singular_values_distribution_pairwise'
    pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
    g = ggpairs(conditional_all_data_samples_mean_w_gt, aes(colour = dataset_info, alpha = 0.4), upper = list(continuous = "density")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(g, progress=F)
    while (!is.null(dev.list()))  { dev.off() }
    
    plot_name = 'cov_ij_singular_values_distribution_pairwise'
    pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
    g = (ggpairs(conditional_all_data_samples_cov_w_gt, aes(colour = dataset_info, alpha = 0.4), upper = list(continuous = "density")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    print(g, progress=F)
    while (!is.null(dev.list()))  { dev.off() }
    
    plot_name = 'mean_cov_ij_singular_values_distribution_pairwise'
    pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
    g = (ggpairs(conditional_all_data_samples_mean_cov_w_gt, aes(colour = dataset_info, alpha = 0.4), upper = list(continuous = "density")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    print(g, progress=F)
    while (!is.null(dev.list()))  { dev.off() }
    
    
    
    #####
    
    
    unconditional_samples_mean = unconditional_samples[,c('S1_mean', 'S2_mean', 'S3_mean', 'S4_mean', 'S5_mean', 'S6_mean', 'S7_mean', 'S8_mean', 'dataset_info')]
    unconditional_samples_cov = unconditional_samples[,c('S1_cov', 'S2_cov', 'S3_cov', 'S4_cov', 'S5_cov', 'S6_cov', 'S7_cov', 'S8_cov', 'dataset_info')]
    unconditional_samples_mean_cov = unconditional_samples[,c('S1_mean', 'S2_mean', 'S3_mean', 'S4_mean', 'S1_cov', 'S2_cov', 'S3_cov', 'S4_cov', 'dataset_info')]
    
    unconditional_samples_mean_w_gt = bind_rows(unconditional_samples_mean, mean_singular_value_df_merged_wide[,colnames(unconditional_samples_mean)])
    unconditional_samples_cov_w_gt = bind_rows(unconditional_samples_cov, cov_singular_value_df_merged_wide[,colnames(unconditional_samples_cov)])
    unconditional_samples_mean_cov_w_gt = bind_rows(unconditional_samples_mean_cov, mean_cov_singular_value_df_merged_wide[,colnames(unconditional_samples_mean_cov)])
    
    plot_save_dir = paste('/home/ishan/Lasker/idr_design/plots/singular_value_pairwise_distribution/unconditional', dataset_info, use_conditional_str, sep='/')
    dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)
    
                 
    plot_name = 'mean_ij_singular_values_distribution_pairwise'
    pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
    g = (ggpairs(unconditional_samples_mean_w_gt, aes(colour = factor(dataset_info), alpha = .4), upper = list(continuous = "density")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    print(g, progress=F)
    while (!is.null(dev.list()))  { dev.off() }
    
    plot_name = 'cov_ij_singular_values_distribution_pairwise'
    pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
    g = (ggpairs(unconditional_samples_cov_w_gt, aes(colour = dataset_info, alpha = 0.4), upper = list(continuous = "density")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    print(g, progress=F)
    while (!is.null(dev.list()))  { dev.off() }
    
    plot_name = 'mean_cov_ij_singular_values_distribution_pairwise'
    pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
    g = (ggpairs(unconditional_samples_mean_cov_w_gt, aes(colour = dataset_info, alpha = 0.4), upper = list(continuous = "density")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    print(g, progress=F)
    while (!is.null(dev.list()))  { dev.off() }
    
    
    #####
    
    mean_singular_value_df_merged_wide_kappa_subset = mean_singular_value_df_merged_wide %>% filter(kappa >= kappa_lower_val) %>% filter(kappa <= kappa_upper_val) %>% as.data.frame()
    cov_singular_value_df_merged_wide_kappa_subset = cov_singular_value_df_merged_wide %>% filter(kappa >= kappa_lower_val) %>% filter(kappa <= kappa_upper_val) %>% as.data.frame()
    mean_cov_singular_value_df_merged_wide_kappa_subset = mean_cov_singular_value_df_merged_wide %>% filter(kappa >= kappa_lower_val) %>% filter(kappa <= kappa_upper_val) %>% as.data.frame()
    
    
    #so sample sizes are equal between 
    conditional_unseen_data_samples = conditional_unseen_data_samples[1:nrow(mean_singular_value_df_merged_wide_kappa_subset),]
    conditional_unseen_data_samples_mean = conditional_unseen_data_samples[,c('ncpr', 'kappa', 'S1_mean', 'S2_mean', 'S3_mean', 'S4_mean', 'S5_mean', 'S6_mean', 'S7_mean', 'S8_mean', 'dataset_info')]
    conditional_unseen_data_samples_cov = conditional_unseen_data_samples[,c('ncpr', 'kappa', 'S1_cov', 'S2_cov', 'S3_cov', 'S4_cov', 'S5_cov', 'S6_cov', 'S7_cov', 'S8_cov', 'dataset_info')]
    conditional_unseen_data_samples_mean_cov = conditional_unseen_data_samples[,c('ncpr', 'kappa', 'S1_mean', 'S2_mean', 'S3_mean', 'S4_mean', 'S1_cov', 'S2_cov', 'S3_cov', 'S4_cov', 'dataset_info')]
    
    
    conditional_unseen_data_samples_mean_w_gt = bind_rows(conditional_unseen_data_samples_mean, mean_singular_value_df_merged_wide_kappa_subset[,colnames(conditional_unseen_data_samples_mean)])
    conditional_unseen_data_samples_cov_w_gt = bind_rows(conditional_unseen_data_samples_cov, cov_singular_value_df_merged_wide_kappa_subset[,colnames(conditional_unseen_data_samples_cov)])
    conditional_unseen_data_samples_mean_cov_w_gt = bind_rows(conditional_unseen_data_samples_mean_cov, mean_cov_singular_value_df_merged_wide_kappa_subset[,colnames(conditional_unseen_data_samples_mean_cov)])
    
    plot_save_dir = paste('/home/ishan/Lasker/idr_design/plots/singular_value_pairwise_distribution/conditional_unseen_data', dataset_info, use_conditional_str, sep='/')
    dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)
    
    if ((d != "/home/ishan/Lasker/idr_design/model_output_data/kappa_0.7_0.8/use_conditional=0") & (d != "/home/ishan/Lasker/idr_design/model_output_data/kappa_0.7_0.8/use_conditional=1")) {
      
      plot_name = 'mean_ij_singular_values_distribution_pairwise'
      pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
      g = (ggpairs(conditional_unseen_data_samples_mean_w_gt, aes(colour = dataset_info, alpha = 0.4), upper = list(continuous = "density")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
      print(g, progress=F)
      while (!is.null(dev.list()))  { dev.off() }
      
      plot_name = 'cov_ij_singular_values_distribution_pairwise'
      pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
      g = (ggpairs(conditional_unseen_data_samples_cov_w_gt, aes(colour = dataset_info, alpha = 0.4), upper = list(continuous = "density")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
      print(g, progress=F)
      while (!is.null(dev.list()))  { dev.off() }
      
      plot_name = 'mean_cov_ij_singular_values_distribution_pairwise'
      pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
      g = (ggpairs(conditional_unseen_data_samples_mean_cov_w_gt, aes(colour = dataset_info, alpha = 0.4), upper = list(continuous = "density")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
      print(g, progress=F)
      while (!is.null(dev.list()))  { dev.off() }
      
    } else {
      
      plot_name = 'mean_ij_singular_values_distribution_pairwise'
      pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
      g = (ggpairs(conditional_unseen_data_samples_mean_w_gt, aes(colour = dataset_info, alpha = 0.4), upper = list(continuous = wrap("cor", size = 2))) +
             theme(axis.text.x = element_text(angle = 45, hjust = 1)))
      print(g, progress=F)
      while (!is.null(dev.list()))  { dev.off() }
      
      plot_name = 'cov_ij_singular_values_distribution_pairwise'
      pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
      g = (ggpairs(conditional_unseen_data_samples_cov_w_gt, aes(colour = dataset_info, alpha = 0.4), upper = list(continuous = wrap("cor", size = 2))) +
             theme(axis.text.x = element_text(angle = 45, hjust = 1)))
      print(g, progress=F)
      while (!is.null(dev.list()))  { dev.off() }
      
      plot_name = 'mean_cov_ij_singular_values_distribution_pairwise'
      pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
      g = (ggpairs(conditional_unseen_data_samples_mean_cov_w_gt, aes(colour = dataset_info, alpha = 0.4), upper = list(continuous = wrap("cor", size = 2))) +
             theme(axis.text.x = element_text(angle = 45, hjust = 1)))
      print(g, progress=F)
      while (!is.null(dev.list()))  { dev.off() }
    }
    
    mean_cols = c('S1_mean', 'S2_mean', 'S3_mean', 'S4_mean')
    cov_cols = c('S1_cov', 'S2_cov', 'S3_cov', 'S4_cov')
    
    
    kl_divergence_all_models = bind_rows(kl_divergence_all_models, calc_kl_divergence(conditional_all_data_samples_mean, mean_singular_value_df_merged_wide, mean_cols, dataset_info_str, 'conditional_all_data', use_conditional_str))
    kl_divergence_all_models = bind_rows(kl_divergence_all_models, calc_kl_divergence(conditional_all_data_samples_cov, cov_singular_value_df_merged_wide, cov_cols, dataset_info_str, 'conditional_all_data', use_conditional_str))
    
    kl_divergence_all_models = bind_rows(kl_divergence_all_models, calc_kl_divergence(unconditional_samples_mean, mean_singular_value_df_merged_wide, mean_cols, dataset_info_str, 'unconditional', use_conditional_str))
    kl_divergence_all_models = bind_rows(kl_divergence_all_models, calc_kl_divergence(unconditional_samples_mean_cov, cov_singular_value_df_merged_wide, cov_cols, dataset_info_str, 'unconditional', use_conditional_str))
    
    kl_divergence_all_models = bind_rows(kl_divergence_all_models, calc_kl_divergence(conditional_unseen_data_samples_mean, mean_singular_value_df_merged_wide_kappa_subset, mean_cols, dataset_info_str, 'conditional_unseen_data', use_conditional_str))
    kl_divergence_all_models = bind_rows(kl_divergence_all_models, calc_kl_divergence(conditional_unseen_data_samples_cov, cov_singular_value_df_merged_wide_kappa_subset, cov_cols, dataset_info_str, 'conditional_unseen_data', use_conditional_str))
    
    ##
    
    cor_abs_diff_all_models = bind_rows(cor_abs_diff_all_models, calc_corr(conditional_all_data_samples_mean, mean_singular_value_df_merged_wide, mean_cols, dataset_info_str, 'conditional_all_data', use_conditional_str))
    cor_abs_diff_all_models = bind_rows(cor_abs_diff_all_models, calc_corr(conditional_all_data_samples_cov, cov_singular_value_df_merged_wide, cov_cols, dataset_info_str, 'conditional_all_data', use_conditional_str))
    
    cor_abs_diff_all_models = bind_rows(cor_abs_diff_all_models, calc_corr(unconditional_samples_mean, mean_singular_value_df_merged_wide, mean_cols, dataset_info_str, 'unconditional', use_conditional_str))
    cor_abs_diff_all_models = bind_rows(cor_abs_diff_all_models, calc_corr(unconditional_samples_mean_cov, cov_singular_value_df_merged_wide, cov_cols, dataset_info_str, 'unconditional', use_conditional_str))
    
    cor_abs_diff_all_models = bind_rows(cor_abs_diff_all_models, calc_corr(conditional_unseen_data_samples_mean, mean_singular_value_df_merged_wide_kappa_subset, mean_cols, dataset_info_str, 'conditional_unseen_data', use_conditional_str))
    cor_abs_diff_all_models = bind_rows(cor_abs_diff_all_models, calc_corr(conditional_unseen_data_samples_cov, cov_singular_value_df_merged_wide_kappa_subset, cov_cols, dataset_info_str, 'conditional_unseen_data', use_conditional_str))
    
  }
}


calc_kl_divergence = function(prediction, ground_truth, columns, dataset_info_str, samples_conditional_str, model_conditional_str)
{
  triplet_idx = list(c(1,2,3), c(1,3,4), c(2,3,4))
  num_neighbors = 10
  
  out = data.frame()
  
  for (i in 1:(length(columns)-1))
  {
    for (j in ((i+1):length(columns)))
    {
      col_pair_str = paste(columns[i], columns[j], sep='-')
      col_pair_str = gsub("mean", "m", col_pair_str)
      col_pair_str = gsub("cov", "c", col_pair_str)
      
      
      kl_div = KL.divergence(as.matrix(prediction[,columns[c(i,j)]]), as.matrix(ground_truth[,columns[c(i,j)]]), k=num_neighbors)
      out = bind_rows(out, data.frame(KL_div = kl_div[num_neighbors], col_pair = col_pair_str))
    }
  }
  
  
  for (i in 1:length(triplet_idx))
  {
    col1 = columns[triplet_idx[[i]][1]]
    col2 = columns[triplet_idx[[i]][2]]
    col3 = columns[triplet_idx[[i]][3]]

    triplet_pair_str = paste(col1, col2, col3, sep='-')
    triplet_pair_str = gsub("mean", "m", triplet_pair_str)
    triplet_pair_str = gsub("cov", "c", triplet_pair_str)  
    
    kl_div = KL.divergence(as.matrix(prediction[,columns[triplet_idx[[i]]]]), as.matrix(ground_truth[,columns[triplet_idx[[i]]]]), k=num_neighbors)
    out = bind_rows(out, data.frame(KL_div = kl_div[num_neighbors], col_pair = triplet_pair_str))
  }
  
  out$dataset_info = dataset_info_str
  out$samples_conditional_str = samples_conditional_str
  out$model_conditional_str = model_conditional_str
  
  return(out)
}



calc_corr = function(prediction, ground_truth, columns, dataset_info_str, samples_conditional_str, model_conditional_str)
{
  cor_abs_diff = abs(cor(prediction[,columns], method='spearman') - cor(ground_truth[,columns], method='spearman'))
  cor_abs_diff = melt(replace(cor_abs_diff, lower.tri(cor_abs_diff, TRUE), NA), na.rm = TRUE)
  cor_abs_diff$dataset_info = dataset_info_str
  cor_abs_diff$samples_conditional_str = samples_conditional_str
  cor_abs_diff$model_conditional_str = model_conditional_str
  
  colnames(cor_abs_diff)[3] = 'cor_abs_diff'
  
  return(cor_abs_diff)
}



plot_save_dir = '/home/ishan/Lasker/idr_design/plots/test_losses'
dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)

print(ggplot(test_losses_final_epoch_all_models, aes(x = as.factor(dataset_info), y = test_losses_final_epoch, fill = model_conditional_str)) +
        geom_violin(alpha=.5) + 
        labs(x ="Dataset Partition", y = "NLL Loss") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('test_losses', '.tiff'), sep = '/'), width = 6, height = 4)


print(ggplot(test_losses_final_epoch_all_models, aes(x = as.factor(dataset_info), y = test_losses_final_epoch, fill = model_conditional_str)) +
        geom_violin(alpha=.5) + 
        labs(x ="Dataset Partition", y = "NLL Loss") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('test_losses', '.png'), sep = '/'), width = 6, height = 4)




plot_save_dir = '/home/ishan/Lasker/idr_design/plots/kll_div'
dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)

kl_divergence_all_models_subset = kl_divergence_all_models %>% filter(model_conditional_str == 'use_conditional=True') %>% as.data.frame()
kl_divergence_all_models_subset$col_pair = factor(kl_divergence_all_models_subset$col_pair, levels = unique(kl_divergence_all_models_subset$col_pair))
kl_divergence_all_models_subset = kl_divergence_all_models_subset %>% filter(dataset_info != 'rand_split') %>% as.data.frame()
kl_divergence_all_models_subset$dataset_info = factor(kl_divergence_all_models_subset$dataset_info, levels = c('kappa_0.0_0.1', 'kappa_0.1_0.2', 'kappa_0.2_0.3', 'kappa_0.3_0.4', 'kappa_0.4_0.5', 'kappa_0.5_0.6', 'kappa_0.6_0.7', 'kappa_0.7_0.8', 'kappa_0.8_0.9', 'kappa_0.0_0.2', 'kappa_0.2_0.4', 'kappa_0.4_0.6', 'kappa_0.6_0.8')) 

kl_divergence_all_models_subset_mean = kl_divergence_all_models_subset %>% filter(grepl('_m', col_pair)) %>% as.data.frame()
kl_divergence_all_models_subset_cov = kl_divergence_all_models_subset %>% filter(grepl('_c', col_pair)) %>% as.data.frame()


print(ggplot(kl_divergence_all_models_subset_mean, aes(x = col_pair, y = KL_div)) +
        geom_point(aes(color = samples_conditional_str), position = position_jitterdodge()) +
        labs(x ="Singular Value Combination", y = "KL Divergence") +
        facet_wrap(vars(dataset_info)) + 
        theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('kl_div_mean', '.tiff'), sep = '/'), width = 12, height = 8)


print(ggplot(kl_divergence_all_models_subset_mean, aes(x = col_pair, y = KL_div)) +
        geom_point(aes(color = samples_conditional_str), position = position_jitterdodge()) +
        labs(x ="Singular Value Combination", y = "KL Divergence") +
        facet_wrap(vars(dataset_info)) + 
        theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('kl_div_mean', '.png'), sep = '/'), width = 12, height = 8)


print(ggplot(kl_divergence_all_models_subset_cov, aes(x = col_pair, y = KL_div)) +
        geom_point(aes(color = samples_conditional_str), position = position_jitterdodge()) +
        labs(x ="Singular Value Combination", y = "KL Divergence") +
        facet_wrap(vars(dataset_info)) + 
        theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('kl_div_cov', '.tiff'), sep = '/'), width = 12, height = 8)


print(ggplot(kl_divergence_all_models_subset_cov, aes(x = col_pair, y = KL_div)) +
        geom_point(aes(color = samples_conditional_str), position = position_jitterdodge()) +
        labs(x ="Singular Value Combination", y = "KL Divergence") +
        facet_wrap(vars(dataset_info)) + 
        theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('kl_div_cov', '.png'), sep = '/'), width = 12, height = 8)




plot_save_dir = '/home/ishan/Lasker/idr_design/plots/kll_div'
dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)

kl_divergence_all_models_subset = cor_abs_diff_all_models %>% filter(model_conditional_str == 'use_conditional=True') %>% as.data.frame()
kl_divergence_all_models_subset$col_pair = factor(kl_divergence_all_models_subset$col_pair, levels = unique(kl_divergence_all_models_subset$col_pair))
kl_divergence_all_models_subset = kl_divergence_all_models_subset %>% filter(dataset_info != 'rand_split') %>% as.data.frame()
kl_divergence_all_models_subset$dataset_info = factor(kl_divergence_all_models_subset$dataset_info, levels = c('kappa_0.0_0.1', 'kappa_0.1_0.2', 'kappa_0.2_0.3', 'kappa_0.3_0.4', 'kappa_0.4_0.5', 'kappa_0.5_0.6', 'kappa_0.6_0.7', 'kappa_0.7_0.8', 'kappa_0.8_0.9', 'kappa_0.0_0.2', 'kappa_0.2_0.4', 'kappa_0.4_0.6', 'kappa_0.6_0.8')) 

kl_divergence_all_models_subset_mean = kl_divergence_all_models_subset %>% filter(grepl('_m', col_pair)) %>% as.data.frame()
kl_divergence_all_models_subset_cov = kl_divergence_all_models_subset %>% filter(grepl('_c', col_pair)) %>% as.data.frame()


print(ggplot(kl_divergence_all_models_subset_mean, aes(x = col_pair, y = KL_div)) +
        geom_point(aes(color = samples_conditional_str), position = position_jitterdodge()) +
        labs(x ="Singular Value Combination", y = "KL Divergence") +
        facet_wrap(vars(dataset_info)) + 
        theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('kl_div_mean', '.tiff'), sep = '/'), width = 12, height = 8)


print(ggplot(kl_divergence_all_models_subset_mean, aes(x = col_pair, y = KL_div)) +
        geom_point(aes(color = samples_conditional_str), position = position_jitterdodge()) +
        labs(x ="Singular Value Combination", y = "KL Divergence") +
        facet_wrap(vars(dataset_info)) + 
        theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('kl_div_mean', '.png'), sep = '/'), width = 12, height = 8)


print(ggplot(kl_divergence_all_models_subset_cov, aes(x = col_pair, y = KL_div)) +
        geom_point(aes(color = samples_conditional_str), position = position_jitterdodge()) +
        labs(x ="Singular Value Combination", y = "KL Divergence") +
        facet_wrap(vars(dataset_info)) + 
        theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('kl_div_cov', '.tiff'), sep = '/'), width = 12, height = 8)


print(ggplot(kl_divergence_all_models_subset_cov, aes(x = col_pair, y = KL_div)) +
        geom_point(aes(color = samples_conditional_str), position = position_jitterdodge()) +
        labs(x ="Singular Value Combination", y = "KL Divergence") +
        facet_wrap(vars(dataset_info)) + 
        theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('kl_div_cov', '.png'), sep = '/'), width = 12, height = 8)



###kl boxplot###


plot_save_dir = '/home/ishan/Lasker/idr_design/plots/kll_div'
dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)


print(ggplot(kl_divergence_all_models_subset_mean, aes(x = col_pair, y = KL_div, fill = samples_conditional_str)) +
        geom_boxplot(alpha=.5) +
        labs(x ="Singular Value Combination", y = "KL Divergence") +
        theme(axis.text.x = element_text(angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('kl_div_mean_boxplot', '.png'), sep = '/'), width = 8, height = 4)

print(ggplot(kl_divergence_all_models_subset_cov, aes(x = col_pair, y = KL_div, fill = samples_conditional_str)) +
        geom_boxplot(alpha=.5) +
        labs(x ="Singular Value Combination", y = "KL Divergence") +
        theme(axis.text.x = element_text(angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('kl_div_cov_boxplot', '.png'), sep = '/'), width = 8, height = 4)








######


plot_save_dir = '/home/ishan/Lasker/idr_design/plots/corr_error'
dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)

colnames(cor_abs_diff_all_models)[1:2] = c('SV_1', 'SV_2')

cor_abs_diff_all_models_subset = cor_abs_diff_all_models %>% filter(model_conditional_str == 'use_conditional=True') %>% as.data.frame()

cor_abs_diff_all_models_subset$SV_1 = gsub("mean", "m", cor_abs_diff_all_models_subset$SV_1)
cor_abs_diff_all_models_subset$SV_1 = gsub("cov", "c", cor_abs_diff_all_models_subset$SV_1)
cor_abs_diff_all_models_subset$SV_2 = gsub("mean", "m", cor_abs_diff_all_models_subset$SV_2)
cor_abs_diff_all_models_subset$SV_2 = gsub("cov", "c", cor_abs_diff_all_models_subset$SV_2)

cor_abs_diff_all_models_subset$SV_1 = factor(cor_abs_diff_all_models_subset$SV_1, levels = c('S1_m', 'S2_m', 'S3_m', 'S1_c', 'S2_c', 'S3_c'))
cor_abs_diff_all_models_subset$SV_2 = factor(cor_abs_diff_all_models_subset$SV_2, levels = c('S2_m', 'S3_m', 'S4_m', 'S2_c', 'S3_c', 'S4_c'))


cor_abs_diff_all_models_subset = cor_abs_diff_all_models_subset %>% filter(dataset_info != 'rand_split') %>% as.data.frame()
cor_abs_diff_all_models_subset$dataset_info = factor(cor_abs_diff_all_models_subset$dataset_info, levels = c('kappa_0.0_0.1', 'kappa_0.1_0.2', 'kappa_0.2_0.3', 'kappa_0.3_0.4', 'kappa_0.4_0.5', 'kappa_0.5_0.6', 'kappa_0.6_0.7', 'kappa_0.7_0.8', 'kappa_0.8_0.9', 'kappa_0.0_0.2', 'kappa_0.2_0.4', 'kappa_0.4_0.6', 'kappa_0.6_0.8')) 

#cor_abs_diff_all_models_subset_mean = cor_abs_diff_all_models_subset %>% filter(grepl('_m', col_pair)) %>% as.data.frame()
#cor_abs_diff_all_models_subset_cov = cor_abs_diff_all_models_subset %>% filter(grepl('_c', col_pair)) %>% as.data.frame()



cor_abs_diff_all_models_subset_conditional_all = cor_abs_diff_all_models_subset %>% filter(samples_conditional_str == 'conditional_all_data') %>% as.data.frame()
print(ggplot(cor_abs_diff_all_models_subset_conditional_all, aes(x = SV_1, y = SV_2, fill = cor_abs_diff)) + 
  geom_tile() + scale_fill_gradientn(colours = rainbow(5), limits=c(0, 1), breaks=seq(0,1,by=0.2))  + 
  labs(title="Conditional, All Data") + 
  facet_wrap(vars(dataset_info)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('corr_error_conditional_all', '.png'), sep = '/'), width = 8, height = 8)


print(ggplot(cor_abs_diff_all_models_subset_conditional_all, aes(x = SV_1, y = SV_2, fill = cor_abs_diff)) + 
        geom_tile() + scale_fill_gradientn(colours = rainbow(5), limits=c(0, 1), breaks=seq(0,1,by=0.2))  + 
        labs(title="Conditional, All Data") + 
        facet_wrap(vars(dataset_info)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('corr_error_conditional_all', '.tiff'), sep = '/'), width = 8, height = 8)



cor_abs_diff_all_models_subset_conditional_unseen_data = cor_abs_diff_all_models_subset %>% filter(samples_conditional_str == 'conditional_unseen_data') %>% as.data.frame()
print(ggplot(cor_abs_diff_all_models_subset_conditional_unseen_data, aes(x = SV_1, y = SV_2, fill = cor_abs_diff)) + 
        geom_tile() + scale_fill_gradientn(colours = rainbow(5), limits=c(0, 1), breaks=seq(0,1,by=0.2))  + 
        labs(title="Conditional, Unseen Data") + 
  facet_wrap(vars(dataset_info)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('corr_error_conditional_unseen', '.png'), sep = '/'), width = 8, height = 8)

print(ggplot(cor_abs_diff_all_models_subset_conditional_unseen_data, aes(x = SV_1, y = SV_2, fill = cor_abs_diff)) + 
        geom_tile() + scale_fill_gradientn(colours = rainbow(5), limits=c(0, 1), breaks=seq(0,1,by=0.2))  + 
        labs(title="Conditional, Unseen Data") + 
        facet_wrap(vars(dataset_info)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('corr_error_conditional_unseen', '.tiff'), sep = '/'), width = 8, height = 8)



cor_abs_diff_all_models_subset_unconditional = cor_abs_diff_all_models_subset %>% filter(samples_conditional_str == 'unconditional') %>% as.data.frame()
print(ggplot(cor_abs_diff_all_models_subset_unconditional, aes(x = SV_1, y = SV_2, fill = cor_abs_diff)) + 
        geom_tile() + scale_fill_gradientn(colours = rainbow(5), limits=c(0, 1.5), breaks=seq(0,1.5,by=0.2))  + 
        labs(title="Unconditional, All Data") + 
        facet_wrap(vars(dataset_info)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('corr_error_unconditional', '.png'), sep = '/'), width = 8, height = 8)

print(ggplot(cor_abs_diff_all_models_subset_unconditional, aes(x = SV_1, y = SV_2, fill = cor_abs_diff)) + 
        geom_tile() + scale_fill_gradientn(colours = rainbow(5), limits=c(0, 1.5), breaks=seq(0,1.5,by=0.2))  + 
        labs(title="Unconditional, All Data") + 
        facet_wrap(vars(dataset_info)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('corr_error_unconditional', '.tiff'), sep = '/'), width = 8, height = 8)


#########boxplot##########

cor_abs_diff_all_models_subset$col_pair = paste(cor_abs_diff_all_models_subset$SV_1, cor_abs_diff_all_models_subset$SV_2, sep='-')
cor_abs_diff_all_models_subset$col_pair = factor(cor_abs_diff_all_models_subset$col_pair, levels = unique(cor_abs_diff_all_models_subset$col_pair))

plot_save_dir = '/home/ishan/Lasker/idr_design/plots/corr_error'
dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)

print(ggplot(cor_abs_diff_all_models_subset, aes(x = col_pair, y = cor_abs_diff, fill = samples_conditional_str)) +
        geom_boxplot(alpha=.5) +
        labs(x ="Singular Value Combination", y = "|Spearman Correlation Diff|") +
        theme(axis.text.x = element_text(angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('corr_boxplot', '.png'), sep = '/'), width = 8, height = 4)


print(ggplot(cor_abs_diff_all_models_subset, aes(x = col_pair, y = cor_abs_diff, fill = samples_conditional_str)) +
        geom_boxplot(alpha=.5) +
        labs(x ="Singular Value Combination", y = "|Spearman Correlation Diff|") +
        theme(axis.text.x = element_text(angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('corr_boxplot', '.tiff'), sep = '/'), width = 8, height = 4)








##for reconstruction: show mean_ij and covariance 




