library(PEcAn.all)
library(PEcAn.SIPNET)
library(PEcAn.uncertainty)
library(PEcAn.settings)
library(future)
library(furrr)
library(data.table)

setwd("/projectnb/dietzelab/guYANG/pecan/")
## Load new PFT classification
newpft <- fread("/projectnb/dietzelab/guYANG/SIPNET_Model_Calibration/Final_PFT_assignment_v4/final_8000_sites_with_final_pft_v4.csv")
newpft <- newpft[, .(
  site = index,
  pft  = final_pft
)]
# fwrite(newpft,"/projectnb/dietzelab/guYANG/SIPNET_Model_Calibration/newpft.csv")
new_path <- "/projectnb/dietzelab/guYANG/SIPNET_Model_Calibration/newpft.csv"
# site_names <- grep("^site\\.", names(settings$run), value = TRUE)
# settings$run[site_names] <- lapply(
#   settings$run[site_names],
#   function(site_setting) {
#     site_setting$inputs$pft.site$path <- new_path
#     site_setting
#   }
# )

# # Old PFT
# load("/projectnb/dietzelab/guYANG/pecan/pecan.Rdata")
# # New PFT
# load("/projectnb/dietzelab/guYANG/pecan/pecan_flux_newpft.Rdata")
# New PFT with updated ERA clim with PFT specific tau
load("/projectnb/dietzelab/guYANG/pecan/pecan_flux_new_with_changed_clim.Rdata")
# # Old PFT samples
# load("/projectnb/dietzelab/guYANG/pecan/samples.Rdata")
# New PFT samples
load("/projectnb/dietzelab/guYANG/pecan/ensemble.samples_16PFT_100ens.RData")
load("/projectnb/dietzelab/guYANG/pecan/output/6obs_monthly/sda.output144.Rdata")
source("/projectnb/dietzelab/guYANG/pecan/pecan/local/diagnose_functions_monthly.R")
settings$outdir      <- "/projectnb/dietzelab/guYANG/pecan/pft_clim"
settings$rundir      <- file.path(settings$outdir, "run")
settings$modeloutdir <- file.path(settings$outdir, "out")
settings$host$rundir <- file.path(settings$outdir, "run")
settings$host$outdir <- file.path(settings$outdir, "out")
settings$host$folder <- file.path(settings$outdir, "out")
settings$ensemble$size <- 25
settings$state.data.assimilation$adjustment <- "FALSE"
settings$host$prerun <- "module load R/4.4.0"
settings$state.data.assimilation$q.type <- "vector"
settings$state.data.assimilation$aqq.Init <- "1"
settings$state.data.assimilation$bqq.Init <- "1"
## Fix the multi output in one timestep bug
settings$model$jobtemplate <- "/projectnb/dietzelab/guYANG/pecan/runners/test7/sipnet_template.job"
# Load the selected sites
load("/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/sda_idx.Rdata")
all_ids  <- vapply(settings, \(s) as.character(s$run$site$id), "")
# Sub settings for this run
# keep_ids <- c("341")
settings <- settings[all_ids %in% keep_ids]
settings <- PEcAn.settings::as.MultiSettings(settings)
# setup the batch job settings.
general.job <- list(cores = 28, folder.num = 80)
batch.settings = structure(list(
  general.job = general.job,
  qsub.cmd = "qsub -l h_rt=24:00:00 -l mem_per_core=4G -l buyin -pe omp @CORES@ -V -N @NAME@ -o @STDOUT@ -e @STDERR@ -S /bin/bash"
))
settings$state.data.assimilation$batch.settings <- batch.settings
# update settings with the actual PFTs.
settings <- PEcAn.settings::prepare.settings(settings)
settings$state.data.assimilation$start.date <- "2012-07-15 00:00:00"
settings$state.data.assimilation$end.date   <- "2024-07-15 23:59:59"

### Generated updated clim file for SoilT and VPDsoil
# source("/projectnb/dietzelab/guYANG/pecan/pecan/local/met2model.SIPNET_tau.R")
# 
# tau15_manifest <- generate_tau_clims_for_lookup(
#   lookup = lookup,
# 
#   input_root =
#     "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/ERA5_2012_2024",
# 
#   output_root =
#     "/projectnb/dietzelab/guYANG/pecan/modified_met/tau15_lookup",
# 
#   members = 1:10,
# 
#   tau_days = 15,
# 
#   n_cores = 14,
# 
#   overwrite = TRUE,
# 
#   verbose = TRUE
# )
# 
##### Parellel Run with SIPNET model
res <- run_sipnet_only_parallel(
  settings = settings,
  ensemble.samples = ensemble.samples,
  cores = 28,
  overwrite_outdir = TRUE,
  run_model = TRUE
)
