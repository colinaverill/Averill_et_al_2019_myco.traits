#paths for mycorrhizal trait analysis project.
#main data directory.----
#user must change this to path where raw data is hosted.
data.dir <- '/fs/data3/caverill/myc_traits/'

#Input data files to replicate analysis.----
#phylogeny file.
phylogeny_raw.path <- paste0(data.dir,'colin_2018-12--2.tre')
#intra- and interspecific analysis files.
intra_specific_analysis_data.path <- paste0(data.dir,'intra_specific_for_analysis.rds')
inter_specific_analysis_data.path <- paste0(data.dir,'inter_specific_for_analysis.rds')
                   
#Analysis output.----
dir <- paste0(data.dir,'analysis_output/')
system(paste0('mkdir -p ',dir))
                                 variance_decomp_output.path <- paste0(dir,'variance_decomposition.rds')
pgls.glmm_myc.biome3_interaction_no.selection_DECIDUOUS.path <- paste0(dir, 'pgls.glmm_myc.biome3_interaction_no.selection_DECIDUOUS.rds')
  lm.glmm_myc.biome3_interaction_no.selection_DECIDUOUS.path <- paste0(dir,   'lm.glmm_myc.biome3_interaction_no.selection.rds')
#phylo estimated traits analysis and output.
  phylo_estimated_traits.path <- paste0(dir,'phylo_estimated_traits.rds')
     phy_est_models_data.path <- paste0(dir,'phy_est_models_data.path')

#Figure and table paths.----
dir <- 'figures/'
            var_decomp_figure.path <- paste0(dir,'variance_decomp.png')
            sample_map_figure.path <- paste0(dir,'global_sampling_map.png')
             phylogeny_figure.path <- paste0(dir,'phylogeny_figure.png')
       lm_pgls_effects_figure.path <- paste0(dir,'lm_vs_pgls_effects_figure.png')
   pgls_model_parameter_table.path <- paste0(dir,'pgls_model_parameter_table.csv')
phylo_estimated_traits_figure.path <- paste0(dir,'phylo_est_traits.png')
         lat_myco_trait_means.path <- paste0(dir,'lat_myco_trait_means.png')
                trait_N_table.path <- paste0(dir,'trait_N_table.csv')
    Supplementary_Data_File_1.path <- paste0(dir,'Supplementary_Data_File_1.csv')
