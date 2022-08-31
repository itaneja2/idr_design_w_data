library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)


samples = read.csv('/home/ishan/Lasker/idr_design/gen_samples/umap.csv')

#samples = read.csv('/home/ishan/Lasker/idr_design/gen_samples/made_baseline.csv')


home_dir = '/home/ishan'


plot_name = 'mean_ratio_ij'
print(ggplot(samples, aes(x = atom_id1, y = atom_id2, fill = mean)) + 
        geom_tile() + scale_fill_viridis(discrete = FALSE) + 
        facet_wrap(vars(seq)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(paste(plot_save_dir, paste0(plot_name, '.png'), sep = '/'), dpi=300)








colnames(samples) = sub("^X", "", colnames(samples))
samples_long = samples %>% gather(atom_pair_id, mean, colnames(samples)[1:(ncol(samples)-2)]) %>% as.data.frame()
samples_long = samples_long %>% separate(col = atom_pair_id, into = c("atom_id1", "atom_id2"), sep = "_", remove=FALSE) %>% as.data.frame()

for (curr_seq in sort(unique(samples_long$nearest_seq)))
{
  print(curr_seq)
  tmp = ratio_mean_sd_long_all %>% filter(seq == curr_seq) %>% as.data.frame()
  tmp = tmp[,c(1,2,3,4)]
  
  tmp$nearest_seq = curr_seq
  tmp$nearest_seq_dist = 0 
  tmp = tmp[,c(5,6,1,2,3,4)]
  
  curr_plot_data = bind_rows(tmp, samples_long %>% filter(nearest_seq == curr_seq) %>% as.data.frame())
  
  curr_plot_data$seq_w_dist = paste(curr_plot_data$nearest_seq, curr_plot_data$nearest_seq_dist, sep = '_')
  curr_plot_data$seq_w_dist = factor(curr_plot_data$seq_w_dist, levels = sort(unique(curr_plot_data$seq_w_dist)))
  
  curr_plot_data$atom_id1 = factor(curr_plot_data$atom_id1, levels = sort(unique(as.numeric(curr_plot_data$atom_id1))))
  curr_plot_data$atom_id2 = factor(curr_plot_data$atom_id2, levels = sort(unique(as.numeric(curr_plot_data$atom_id2))))
  

  plot_save_dir = paste('/home/ishan/Lasker/idr_design/plots/umap', curr_seq, sep ='/')
  dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)
  
  plot_name = 'mean_ratio_ij'
  print(ggplot(curr_plot_data, aes(x = atom_id1, y = atom_id2, fill = mean)) + 
          geom_tile() + scale_fill_viridis(discrete = FALSE) + 
          facet_wrap(vars(seq_w_dist)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  ggsave(paste(plot_save_dir, paste0(plot_name, '.png'), sep = '/'), dpi=300)
  
}


