library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS)
library(viridis)
library(reshape2)
library(GGally)
library(FNN)


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


mean_singular_value_df_merged_wide$dataset_info = 'original'
cov_singular_value_df_merged_wide$dataset_info = 'original'
mean_cov_singular_value_df_merged_wide$dataset_info = 'original'

mean_singular_value_df_merged_wide = subset(mean_singular_value_df_merged_wide, select = -seq_id)
cov_singular_value_df_merged_wide = subset(cov_singular_value_df_merged_wide, select = -seq_id)
mean_cov_singular_value_df_merged_wide = subset(mean_cov_singular_value_df_merged_wide, select = -seq_id)






mean_singular_value_df_hydrophobic_variants = read.csv('/home/ishan/Lasker/mpipi/ood_constructs/hydrophobic_variants/output_data/all_samples/mean_ij_singular_value_data.csv')
cov_singular_value_df_hydrophobic_variants = read.csv('/home/ishan/Lasker/mpipi/ood_constructs/hydrophobic_variants/output_data/all_samples/cov_ij_singular_value_data.csv')

mean_singular_value_df_kappa_variants_ood = read.csv('/home/ishan/Lasker/mpipi/ood_constructs/kappa_variants/output_data/all_samples/mean_ij_singular_value_data.csv')
cov_singular_value_df_kappa_variants_ood = read.csv('/home/ishan/Lasker/mpipi/ood_constructs/kappa_variants/output_data/all_samples/cov_ij_singular_value_data.csv')

mean_singular_value_df_hydrophobic_variants$seq_id = mean_singular_value_df_hydrophobic_variants$full_seq
cov_singular_value_df_hydrophobic_variants$seq_id = cov_singular_value_df_hydrophobic_variants$full_seq

mean_singular_value_df_kappa_variants_ood$seq_id = mean_singular_value_df_kappa_variants_ood$full_seq
cov_singular_value_df_kappa_variants_ood$seq_id = cov_singular_value_df_kappa_variants_ood$full_seq

mean_singular_value_df_hydrophobic_variants = mean_singular_value_df_hydrophobic_variants[,c(1,2,3,4,7,8,10)]
cov_singular_value_df_hydrophobic_variants = cov_singular_value_df_hydrophobic_variants[,c(1,2,3,4,7,8,10)]

mean_singular_value_df_kappa_variants_ood = mean_singular_value_df_kappa_variants_ood[,c(1,2,3,4,7,8,9)]
cov_singular_value_df_kappa_variants_ood = cov_singular_value_df_kappa_variants_ood[,c(1,2,3,4,7,8,9)]




#convert to wide
mean_singular_value_df_hydrophobic_variants$num_singular_values = paste0('S', mean_singular_value_df_hydrophobic_variants$num_singular_values, '_mean')
mean_singular_value_df_hydrophobic_variants = mean_singular_value_df_hydrophobic_variants[,c('seq_id', 'ncpr', 'kappa', 'singular_value', 'num_singular_values')]
mean_singular_value_df_hydrophobic_variants_wide = spread(mean_singular_value_df_hydrophobic_variants, key = num_singular_values, value = singular_value) %>% as.data.frame()
mean_singular_value_df_hydrophobic_variants_wide[mean_singular_value_df_hydrophobic_variants_wide$kappa == -1, 'kappa'] = 0
mean_singular_value_df_hydrophobic_variants_wide = mean_singular_value_df_hydrophobic_variants_wide[,c('seq_id', 'ncpr', 'kappa', mean_singular_value_df_hydrophobic_variants$num_singular_values[1:24])]

cov_singular_value_df_hydrophobic_variants$num_singular_values = paste0('S', cov_singular_value_df_hydrophobic_variants$num_singular_values, '_cov')
cov_singular_value_df_hydrophobic_variants = cov_singular_value_df_hydrophobic_variants[,c('seq_id', 'ncpr', 'kappa', 'singular_value', 'num_singular_values')]
cov_singular_value_df_hydrophobic_variants_wide = spread(cov_singular_value_df_hydrophobic_variants[,c('seq_id', 'ncpr', 'kappa', 'singular_value', 'num_singular_values')], key = num_singular_values, value = singular_value) %>% as.data.frame()
cov_singular_value_df_hydrophobic_variants_wide[cov_singular_value_df_hydrophobic_variants_wide$kappa == -1, 'kappa'] = 0
cov_singular_value_df_hydrophobic_variants_wide = cov_singular_value_df_hydrophobic_variants_wide[,c('seq_id', 'ncpr', 'kappa', cov_singular_value_df_hydrophobic_variants$num_singular_values[1:276])]


tmp1 = mean_singular_value_df_hydrophobic_variants_wide[,1:7]
tmp2 = cov_singular_value_df_hydrophobic_variants_wide[,c(1,4:7)]
mean_cov_singular_value_df_hydrophobic_variants_wide = left_join(tmp1, tmp2, by = 'seq_id') %>% as.data.frame()


mean_singular_value_df_hydrophobic_variants_wide$dataset_info = 'hydrophobic'
cov_singular_value_df_hydrophobic_variants_wide$dataset_info = 'hydrophobic'
mean_cov_singular_value_df_hydrophobic_variants_wide$dataset_info = 'hydrophobic'

mean_singular_value_df_hydrophobic_variants_wide = subset(mean_singular_value_df_hydrophobic_variants_wide, select = -seq_id)
cov_singular_value_df_hydrophobic_variants_wide = subset(cov_singular_value_df_hydrophobic_variants_wide, select = -seq_id)
mean_cov_singular_value_df_hydrophobic_variants_wide = subset(mean_cov_singular_value_df_hydrophobic_variants_wide, select = -seq_id)



#convert to wide
mean_singular_value_df_kappa_variants_ood$num_singular_values = paste0('S', mean_singular_value_df_kappa_variants_ood$num_singular_values, '_mean')
mean_singular_value_df_kappa_variants_ood = mean_singular_value_df_kappa_variants_ood[,c('seq_id', 'ncpr', 'kappa', 'singular_value', 'num_singular_values')]
mean_singular_value_df_kappa_variants_ood_wide = spread(mean_singular_value_df_kappa_variants_ood, key = num_singular_values, value = singular_value) %>% as.data.frame()
mean_singular_value_df_kappa_variants_ood_wide[mean_singular_value_df_kappa_variants_ood_wide$kappa == -1, 'kappa'] = 0
mean_singular_value_df_kappa_variants_ood_wide = mean_singular_value_df_kappa_variants_ood_wide[,c('seq_id', 'ncpr', 'kappa', mean_singular_value_df_kappa_variants_ood$num_singular_values[1:24])]

cov_singular_value_df_kappa_variants_ood$num_singular_values = paste0('S', cov_singular_value_df_kappa_variants_ood$num_singular_values, '_cov')
cov_singular_value_df_kappa_variants_ood = cov_singular_value_df_kappa_variants_ood[,c('seq_id', 'ncpr', 'kappa', 'singular_value', 'num_singular_values')]
cov_singular_value_df_kappa_variants_ood_wide = spread(cov_singular_value_df_kappa_variants_ood[,c('seq_id', 'ncpr', 'kappa', 'singular_value', 'num_singular_values')], key = num_singular_values, value = singular_value) %>% as.data.frame()
cov_singular_value_df_kappa_variants_ood_wide[cov_singular_value_df_kappa_variants_ood_wide$kappa == -1, 'kappa'] = 0
cov_singular_value_df_kappa_variants_ood_wide = cov_singular_value_df_kappa_variants_ood_wide[,c('seq_id', 'ncpr', 'kappa', cov_singular_value_df_kappa_variants_ood$num_singular_values[1:276])]


tmp1 = mean_singular_value_df_kappa_variants_ood_wide[,1:7]
tmp2 = cov_singular_value_df_kappa_variants_ood_wide[,c(1,4:7)]
mean_cov_singular_value_df_kappa_variants_ood_wide = left_join(tmp1, tmp2, by = 'seq_id') %>% as.data.frame()


mean_singular_value_df_kappa_variants_ood_wide$dataset_info = 'kappa'
cov_singular_value_df_kappa_variants_ood_wide$dataset_info = 'kappa'
mean_cov_singular_value_df_kappa_variants_ood_wide$dataset_info = 'kappa'

mean_singular_value_df_kappa_variants_ood_wide = subset(mean_singular_value_df_kappa_variants_ood_wide, select = -seq_id)
cov_singular_value_df_kappa_variants_ood_wide = subset(cov_singular_value_df_kappa_variants_ood_wide, select = -seq_id)
mean_cov_singular_value_df_kappa_variants_ood_wide = subset(mean_cov_singular_value_df_kappa_variants_ood_wide, select = -seq_id)


######

mean_singular_value_df_all = bind_rows(mean_singular_value_df_merged_wide, mean_singular_value_df_hydrophobic_variants_wide, mean_singular_value_df_kappa_variants_ood_wide)
cov_singular_value_df_all = bind_rows(cov_singular_value_df_merged_wide, cov_singular_value_df_hydrophobic_variants_wide, cov_singular_value_df_kappa_variants_ood_wide)
mean_cov_singular_value_df_all = bind_rows(mean_cov_singular_value_df_merged_wide, mean_cov_singular_value_df_hydrophobic_variants_wide, mean_cov_singular_value_df_kappa_variants_ood_wide)

mean_singular_value_df_all = mean_singular_value_df_all[,c('ncpr', 'kappa', 'S1_mean', 'S2_mean', 'S3_mean', 'S4_mean', 'S5_mean', 'S6_mean', 'S7_mean', 'S8_mean', 'dataset_info')]
cov_singular_value_df_all = cov_singular_value_df_all[,c('ncpr', 'kappa', 'S1_cov', 'S2_cov', 'S3_cov', 'S4_cov', 'S5_cov', 'S6_cov', 'S7_cov', 'S8_cov', 'dataset_info')]
mean_cov_singular_value_df_all = mean_cov_singular_value_df_all[,c('ncpr', 'kappa', 'S1_mean', 'S2_mean', 'S3_mean', 'S4_mean', 'S1_cov', 'S2_cov', 'S3_cov', 'S4_cov', 'dataset_info')]

mean_singular_value_df_all$dataset_info = factor(mean_singular_value_df_all$dataset_info, levels = c('original', 'kappa', 'hydrophobic'))
mean_cov_singular_value_df_all$dataset_info = factor(mean_cov_singular_value_df_all$dataset_info, levels = c('original', 'kappa', 'hydrophobic'))

######

plot_save_dir = paste('/home/ishan/Lasker/idr_design/plots/singular_value_pairwise_distribution/ood_constructs/hydrophobic_variants', sep='/')
dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)

plot_name = 'mean_ij_singular_values_distribution_pairwise'
pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
g = (ggpairs(mean_singular_value_df_all %>%  filter(dataset_info != 'kappa') %>% as.data.frame(), aes(colour = dataset_info, alpha = .4), upper = list(continuous = "density")) +
       theme(axis.text.x = element_text(angle = 45, hjust = 1)))
print(g, progress=F)
while (!is.null(dev.list()))  { dev.off() }

plot_name = 'cov_ij_singular_values_distribution_pairwise'
pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
g = (ggpairs(cov_singular_value_df_all %>%  filter(dataset_info != 'kappa') %>% as.data.frame(), aes(colour = dataset_info, alpha = .4), upper = list(continuous = "density")) +
       theme(axis.text.x = element_text(angle = 45, hjust = 1)))
print(g, progress=F)
while (!is.null(dev.list()))  { dev.off() }

plot_name = 'mean_cov_ij_singular_values_distribution_pairwise'
pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
g = (ggpairs(mean_cov_singular_value_df_all %>%  filter(dataset_info != 'kappa') %>% as.data.frame(), aes(colour = dataset_info, alpha = .4), upper = list(continuous = "density")) +
       theme(axis.text.x = element_text(angle = 45, hjust = 1)))
print(g, progress=F)
while (!is.null(dev.list()))  { dev.off() }


######

plot_save_dir = paste('/home/ishan/Lasker/idr_design/plots/singular_value_pairwise_distribution/ood_constructs/kappa_variants', sep='/')
dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)

plot_name = 'mean_ij_singular_values_distribution_pairwise'
pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
g = (ggpairs(mean_singular_value_df_all %>%  filter(dataset_info != 'hydrophobic') %>% as.data.frame(), aes(colour = dataset_info, alpha = .4), upper = list(continuous = "density")) +
       theme(axis.text.x = element_text(angle = 45, hjust = 1)))
print(g, progress=F)
while (!is.null(dev.list()))  { dev.off() }

plot_name = 'cov_ij_singular_values_distribution_pairwise'
pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
g = (ggpairs(cov_singular_value_df_all %>%  filter(dataset_info != 'hydrophobic') %>% as.data.frame(), aes(colour = dataset_info, alpha = .4), upper = list(continuous = "density")) +
       theme(axis.text.x = element_text(angle = 45, hjust = 1)))
print(g, progress=F)
while (!is.null(dev.list()))  { dev.off() }

plot_name = 'mean_cov_ij_singular_values_distribution_pairwise'
pdf(paste(plot_save_dir, paste0(plot_name, '.pdf'), sep = '/'), width=12, height=12)
g = (ggpairs(mean_cov_singular_value_df_all %>%  filter(dataset_info != 'hydrophobic') %>% as.data.frame(), aes(colour = dataset_info, alpha = .4), upper = list(continuous = "density")) +
       theme(axis.text.x = element_text(angle = 45, hjust = 1)))
print(g, progress=F)
while (!is.null(dev.list()))  { dev.off() }


######





calc_kl_divergence = function(prediction, ground_truth, columns, dataset_info_str)
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

  
  return(out)
}


calc_corr = function(prediction, ground_truth, columns, dataset_info_str)
{
  cor_abs_diff = abs(cor(prediction[,columns], method='spearman') - cor(ground_truth[,columns], method='spearman'))
  cor_abs_diff = melt(replace(cor_abs_diff, lower.tri(cor_abs_diff, TRUE), NA), na.rm = TRUE)
  cor_abs_diff$dataset_info = dataset_info_str
  
  colnames(cor_abs_diff)[3] = 'cor_abs_diff'
  
  return(cor_abs_diff)
}


######

kl_divergence_all = data.frame()
cor_abs_diff_all = data.frame()

mean_cols = c('S1_mean', 'S2_mean', 'S3_mean', 'S4_mean')
cov_cols = c('S1_cov', 'S2_cov', 'S3_cov', 'S4_cov')

kl_divergence_all = bind_rows(kl_divergence_all, calc_kl_divergence(mean_singular_value_df_all %>%  filter(dataset_info == 'hydrophobic') %>% as.data.frame(), mean_singular_value_df_merged_wide, mean_cols, 'hydrophobic'))
kl_divergence_all = bind_rows(kl_divergence_all, calc_kl_divergence(cov_singular_value_df_all %>%  filter(dataset_info == 'hydrophobic') %>% as.data.frame(), cov_singular_value_df_merged_wide, cov_cols, 'hydrophobic'))

cor_abs_diff_all = bind_rows(cor_abs_diff_all, calc_corr(mean_singular_value_df_all %>%  filter(dataset_info == 'hydrophobic') %>% as.data.frame(), mean_singular_value_df_merged_wide, mean_cols, 'hydrophobic'))
cor_abs_diff_all = bind_rows(cor_abs_diff_all, calc_corr(cov_singular_value_df_all %>%  filter(dataset_info == 'hydrophobic') %>% as.data.frame(), cov_singular_value_df_merged_wide, cov_cols, 'hydrophobic'))

kl_divergence_all = bind_rows(kl_divergence_all, calc_kl_divergence(mean_singular_value_df_all %>%  filter(dataset_info == 'kappa') %>% as.data.frame(), mean_singular_value_df_merged_wide, mean_cols, 'kappa'))
kl_divergence_all = bind_rows(kl_divergence_all, calc_kl_divergence(cov_singular_value_df_all %>%  filter(dataset_info == 'kappa') %>% as.data.frame(), cov_singular_value_df_merged_wide, cov_cols, 'kappa'))

cor_abs_diff_all = bind_rows(cor_abs_diff_all, calc_corr(mean_singular_value_df_all %>%  filter(dataset_info == 'kappa') %>% as.data.frame(), mean_singular_value_df_merged_wide, mean_cols, 'kappa'))
cor_abs_diff_all = bind_rows(cor_abs_diff_all, calc_corr(cov_singular_value_df_all %>%  filter(dataset_info == 'kappa') %>% as.data.frame(), cov_singular_value_df_merged_wide, cov_cols, 'kappa'))

######



plot_save_dir = '/home/ishan/Lasker/idr_design/plots/kll_div/ood'
dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)

kl_divergence_all$col_pair = factor(kl_divergence_all$col_pair, levels = unique(kl_divergence_all$col_pair))
kl_divergence_all$dataset_info = factor(kl_divergence_all$dataset_info, levels = c('original', 'hydrophobic', 'kappa'))

kl_divergence_all_mean = kl_divergence_all %>% filter(grepl('_m', col_pair)) %>% as.data.frame()
kl_divergence_all_cov = kl_divergence_all %>% filter(grepl('_c', col_pair)) %>% as.data.frame()


print(ggplot(kl_divergence_all_mean, aes(x = col_pair, y = KL_div, fill = dataset_info)) +
        geom_bar(stat='identity', position=position_dodge()) +
        labs(x ="Singular Value Combination", y = "KL Divergence") +
        theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('kl_div_mean', '.tiff'), sep = '/'), width = 6, height = 4)


print(ggplot(kl_divergence_all_mean, aes(x = col_pair, y = KL_div, fill = dataset_info)) +
        geom_bar(stat='identity', position=position_dodge()) +
        labs(x ="Singular Value Combination", y = "KL Divergence") +
        theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('kl_div_mean', '.png'), sep = '/'), width = 6, height = 4)


print(ggplot(kl_divergence_all_cov, aes(x = col_pair, y = KL_div, fill = dataset_info)) +
        geom_bar(stat='identity', position=position_dodge()) +
        labs(x ="Singular Value Combination", y = "KL Divergence") +
        theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('kl_div_cov', '.tiff'), sep = '/'), width = 6, height = 4)


print(ggplot(kl_divergence_all_cov, aes(x = col_pair, y = KL_div, fill = dataset_info)) +
        geom_bar(stat='identity', position=position_dodge()) +
        labs(x ="Singular Value Combination", y = "KL Divergence") +
        theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('kl_div_cov', '.png'), sep = '/'), width = 6, height = 4)




######


plot_save_dir = '/home/ishan/Lasker/idr_design/plots/corr_error/ood'
dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)

colnames(cor_abs_diff_all)[1:2] = c('SV_1', 'SV_2')

cor_abs_diff_all$SV_1 = gsub("mean", "m", cor_abs_diff_all$SV_1)
cor_abs_diff_all$SV_1 = gsub("cov", "c", cor_abs_diff_all$SV_1)
cor_abs_diff_all$SV_2 = gsub("mean", "m", cor_abs_diff_all$SV_2)
cor_abs_diff_all$SV_2 = gsub("cov", "c", cor_abs_diff_all$SV_2)

cor_abs_diff_all$SV_1 = factor(cor_abs_diff_all$SV_1, levels = c('S1_m', 'S2_m', 'S3_m', 'S1_c', 'S2_c', 'S3_c'))
cor_abs_diff_all$SV_2 = factor(cor_abs_diff_all$SV_2, levels = c('S2_m', 'S3_m', 'S4_m', 'S2_c', 'S3_c', 'S4_c'))


cor_abs_diff_all$dataset_info = factor(cor_abs_diff_all$dataset_info, levels = c('original', 'hydrophobic', 'kappa')) 

print(ggplot(cor_abs_diff_all, aes(x = SV_1, y = SV_2, fill = cor_abs_diff)) + 
        geom_tile() + scale_fill_gradientn(colours = rainbow(5), limits=c(0, 1), breaks=seq(0,1,by=0.2))  + 
        facet_wrap(vars(dataset_info)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('corr_error', '.png'), sep = '/'), width = 6, height = 4)


print(ggplot(cor_abs_diff_all, aes(x = SV_1, y = SV_2, fill = cor_abs_diff)) + 
        geom_tile() + scale_fill_gradientn(colours = rainbow(5), limits=c(0, 1), breaks=seq(0,1,by=0.2))  + 
        facet_wrap(vars(dataset_info)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(paste(plot_save_dir, paste0('corr_error', '.tiff'), sep = '/'), width = 6, height = 4)

