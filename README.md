# PIPA
Prognosis-informed Phenotype Assignment (PIPA): A Novel Method for Unsupervised Discovery of Cell Phenotypes with Prognostic Significance in the Tumor Microenvironment

### package installation
devtools::install_github('MPE-Lab/PIPA')
library(PIPA)

### Set paths
root_dir <- "C:/Desktop/"
results_folder <- 'PIPA_analysis'
output_dir <- file.path(root_dir, results_folder)
dir.create(output_dir)

### set parameters
Lasso_run_no<- 3 #min 2 repeats
min_cluster_size <- 5 #please choose according the problem, this is ONLY for demonstration purpose

## prepare data for feature selection
### prepare cell feature data
cell_data_sub <- PIPA::cell_data
cell_data_sub[,c("celltype","core_ids","cell_ids")] <- NULL
### prepare surv data
surv_data <- patient_data[, c("cens","time")]
### feature selection
step1_dir <- output_dir
featureSelection_wrapper(data=cell_data_sub, surv_data=surv_data,
                         log10_featureNms=c('feature6','feature7'),  
                         no_of_runs=Lasso_run_no,
                         output_dir=step1_dir,
                         min_cluster_size= min_cluster_size)
                         
## 
