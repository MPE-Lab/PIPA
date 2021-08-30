# PIPA
Prognosis-informed Phenotype Assignment (PIPA): A Novel Method for Unsupervised Discovery of Cell Phenotypes with Prognostic Significance in the Tumor Microenvironment

## Initialization
````
# Package installation
# devtools::install_github('MPE-Lab/PIPA')
package_root_dir <- 'C:/Desktop'
install.packages(file.path(package_root_dir,'PIPA.tar.gz'), repos = NULL, type = 'source')
library(PIPA)

# Set paths
root_dir <- "C:/Desktop/";
results_folder <- 'PIPA_analysis'
output_dir <- file.path(root_dir, results_folder)
dir.create(output_dir)

# Set parameters
Lasso_run_no<- 10 #min = 2, default = 10
min_cluster_size <- 5 #please choose according the problem, this is ONLY for demonstration purpose
````

## Prepare data for feature selection
````
# prepare cell feature data
cell_data_sub <- PIPA::cell_data
cell_data_sub[,c("celltype","core_ids","cell_ids")] <- NULL
# prepare surv data
surv_data <- patient_data[, c("cens","time")]
# feature selection
step1_dir <- output_dir
featureSelection_wrapper(data=cell_data_sub, surv_data=surv_data,
                         log10_featureNms=c('feature6','feature7'),  
                         no_of_runs=Lasso_run_no,
                         output_dir=step1_dir,
                         min_cluster_size= min_cluster_size)
````

## Build phenotyping model using each (or the user-selected) cell feature subsets determined above
````
# get folder names for each candidate feature subsets
feature_subset_fnms <- list.files(path = step1_dir, pattern = "^detection_rate")
# MANUAL selection (exclusion) for feature subsets
feature_subset_fnms <- setdiff(feature_subset_fnms, 'detection_rate0.9')
 
# loop over individual candidate feature subsets
for (ff in feature_subset_fnms){
  subtype_subdir <- file.path(step1_dir, ff)
  if(!dir.exists(subtype_subdir))next
  
  # build cell phenotyping model (using previously selected training cells)
  model_dir <- file.path(subtype_subdir, 'model')
  dir.create(model_dir)
  
  cell_data <- PIPA::cell_data
  cell_phenotyping_model(cell_data=cell_data,log10_featureNms=log10_featureNms,
                         input_dir=subtype_subdir,output_dir=model_dir,
                         batch_size=500, 
                         consensus_rep=2,
                         no_of_random_merging=5,
                         max_batch=5,
                         cell_batch_clustering = TRUE,
                         cell_batch_merging=TRUE)
  
  # cell phenotype assignment of all cell of interest
  pheno_assignment_dir <- file.path(subtype_subdir, 'pheno_assignment')
  dir.create(pheno_assignment_dir)
  
  cell_phenotype_assignment(cell_data=cell_data,log10_featureNms=log10_featureNms, 
                            selected_feature_dir=subtype_subdir,
                            cell_model_input_dir=model_dir,
                            output_dir=pheno_assignment_dir)
  
  # refine phenotypes by merging based on prognostic significance
  
  # refine phenotypes: calculate phenotype density
  calculate_pheno_density(area_data=area_data,  
                          input_dir=pheno_assignment_dir, output_dir=pheno_assignment_dir)
  # refine phenotypes: cell density-based survival analysis
  final_pheno_dir <- file.path(subtype_subdir, 'final_pheno')
  dir.create(final_pheno_dir)
  
  colnames(patient_data)
  density_surv_analysis(surv_data= patient_data,  
                        input_dir=pheno_assignment_dir, output_dir=final_pheno_dir, 
                        density_fnm= 'cell_phenotype_density.RData',
                        KM_xlab= 'survivals')
  # refine phenotypes: merging
  cutoff_pval <- 0.05 # this is only for demonstration, please use cutoff 0.005
  merged_pheno_dir <- file.path(final_pheno_dir, paste0('merged_pheno_P',cutoff_pval))
  dir.create(merged_pheno_dir)
  out <- refine_pheno_byMerging(surv_data= patient_data, 
                                KM_xlab='survivals', output_dir=merged_pheno_dir,
                                preMerge_density_dir =pheno_assignment_dir, 
                                preMerge_pheno_dir=pheno_assignment_dir,
                                uni_Cox_dir=file.path(final_pheno_dir,'density_surv/ptrend'), 
                                uni_Cox_fnm='uniCoxPH_summary.txt',uni_P_cutoff=cutoff_pval)
  if(out==(-1))next
  
  
  # Post-PIPA analysis:: calculate cell density of PP/IP/GP phenotypes
  # compute phenotype density
  calculate_pheno_density(area_data=area_data,  
                          input_dir=merged_pheno_dir, output_dir=merged_pheno_dir)
  
  # Post-PIPA analysis:: phenotype density-based survival analysis
  density_surv_analysis(surv_data= patient_data,  
                        input_dir=merged_pheno_dir, output_dir=merged_pheno_dir, 
                        density_fnm= 'cell_phenotype_density.RData',
                        KM_xlab= 'Survivals')
  
  # Post-PIPA analysis:: cell composition survival analysis using PhenoGraph
  composition_surv_analysis(surv_data= patient_data,  
                            density_data = raw_dens_data,
                            min_cluster_size = min_cluster_size, 
                            phenotype_dir=merged_pheno_dir, output_dir=merged_pheno_dir,
                            neighbor_size_byFrac= seq(0.05,0.25,by=0.15),
                            KM_xlab='Survivals',area_split=TRUE)
  
  # =========================
  # Post-PIPA visualization 
  # =========================
  selected_feature <- read.csv(file = file.path(subtype_subdir,'selectedCellFeatures.csv'), as.is = TRUE)
  
  # heat map - distribution of PIPA-selected cell features by phenotypes
  cell_data_sub <- PIPA::cell_data[, c("cell_ids",selected_feature$features)]
  Viz_feature_boxplot_byPheno(cell_data= cell_data_sub,plot_bnm='PIPA_features',outlier_removal=TRUE,
                              pheno_levels = c('PP','IP','GP'),
                              pheno_dir=merged_pheno_dir, output_dir=merged_pheno_dir,
                              ncol = 4, plot_height = 16, plot_width = 15)
  # box plot - distribution of PIPA-selected cell features by phenotypes
  cell_data_sub <- PIPA::cell_data[, c("cell_ids","tumor_ids",selected_feature$features)]
  Viz_feature_heatmap_byPheno(cell_data= cell_data_sub,
                              pheno_dir=merged_pheno_dir, output_dir=merged_pheno_dir,
                              pheno_levels = c('PP','IP','GP'),
                              column_scale=TRUE)
  # UMAP - relatedness of PIPA-derived phenotypes in the PIPA-selected feature space, labelled by phenotypes & features
  cell_data_sub <- PIPA::cell_data[, c("cell_ids",selected_feature$features)]
  Viz_feature_UMAP(cell_data= cell_data_sub, downspl_size=100,
                   pheno_dir=merged_pheno_dir, output_dir=merged_pheno_dir,
                   log10_featureNms=log10_featureNms,
                   plot_height = 8, plot_width = 8)
} 

````
## Visualization for comparing parameter settings e.g. (i) detection rates for optimal feature selection (ii) batch_size/max_batch/no_of_random_merging 
````
# Viz: (a) cell count and (b) cell detection of patients/tumors, for each phenotypes
param_bnm='detection_rate'
output_dir2 <- file.path(step1_dir, paste0("Viz_all_",param_bnm)); dir.create(output_dir2)
cutoff_pval <- 0.05 # this is only for demonstration, please use cutoff 0.005
pheno_dir_path <- paste0('final_pheno/merged_pheno_P',cutoff_pval)
Viz_cellcountNdetection_byParam(root_dir=step1_dir,
                                output_dir= output_dir2,
                                pheno_dir_path=pheno_dir_path,
                                param_bnm=param_bnm,
                                tissue_area_list=c('T','S'))

# Forest plot summarizing Cox PH regression analysis results based on density quartiles of PP, IP, GP
Viz_foresPlot_denQ_CoxPH(root_dir=step1_dir, output_dir= output_dir2,
                         pheno_dir_path=pheno_dir_path,
                         param_bnm='detection_rate',
                         Cox_result_folderNm='density_surv/cat')

````
