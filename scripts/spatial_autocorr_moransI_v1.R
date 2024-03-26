# Run Spatial Autocorrelation using Moran's I 
### make sure all libs are installed a priori!

## load libs ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(scater)
  library(Seurat)
  library(dplyr)
  library(magrittr)
  library(png)
  library(cowplot)
  library(jpeg)
  library(scales)    
  library(patchwork)
  library(Matrix)
  library(progressr)
  library(BiocParallel)    
})

# set R options, paths, object params, names, functions, etc..
message("Setting options, paths and object params..")
Sys.setenv(LANG = "en")

#' To get metadata df from sfe object
#' @param object SFE or Seurat object
#' @noRd
getMeta <- function(object = NULL) {
  if (class(object) == "Seurat") {
    #is(object, "Seurat")
    # return Seurat metadata
    return(methods::slot(object, name = "meta.data"))
  } else {
    #return(colData(object) |> slot(name = "listData") |> as.data.frame.list())
    return(colData(object) |> 
             methods::slot(name = "listData") |> 
             as.data.frame.list())
  }
}

# set root dir
dir_working <- "./spatial_techs_compare/segmentation_based/"
#dir_path <- file.path(dir_working, "objects")
setwd(dir_working)

obj.bg <- 
  readRDS(file = "./spatial_techs_compare/segmentation_based/objects/vz.re.xen_matched_bg.rds")
obj.bg
obj.bg %>% getMeta %>% str

# Subset to keep only 2 samples - mb295 and mb266 ----
if (TRUE) {
  cells.use <- 
    obj.bg %>%
    getMeta %>%  
    slice(pull(., samples) %>% grep("295|266", .)) %>% rownames
  
  # subset obj and FOVs ----
  # load modified subset function
  source("./scripts/subset_obj_seurat_v2.R")
  obj.bg %<>% subset_opt(cells = cells.use)
  
  # Update metadata ----
  obj.bg@meta.data %<>%
    mutate(samples = stringi::stri_replace_all_regex(samples, 
                                                     pattern = "region[0-1].|w(.)[^:][0-9].", 
                                                     replacement = ""))
  }


# moran's I (global) -> compute per sample section
# set assay to use
assay_use <- "MultiTech"
message("Assay: ", assay_use, " will be used for count matrix")
count_matrix <- GetAssayData(obj.bg, assay = assay_use, slot = "counts")

# moran's I (global) -> compute per sample section/FOV
start.time <- Sys.time()
pbapply::pboptions(type = "timer", style = 1, char = "=")
obj.moransI.fast <-
  pbapply::pblapply(obj.bg %>% 
                      Images() %>% seq, 
                    function(im) { 
                      BiocParallel::bplapply(obj.bg %>% rownames %>% seq, function(i) {
                        msi <-  
                          moranfast::moranfast(x = 
                                                 # count matrix per sample
                                                 # using raw or normalized assay
                                                 count_matrix %>%  
                                                 .[i, match(obj.bg[[Images(obj.bg)[im]]] %>% Cells(),
                                                            obj.bg %>% Cells()) %>% na.omit()],
                                               c1 = obj.bg[[Images(obj.bg)[im]]][["centroids"]] %>% 
                                                 GetTissueCoordinates() %>% pull(x),
                                               c2 = obj.bg[[Images(obj.bg)[im]]][["centroids"]] %>% 
                                                 GetTissueCoordinates() %>% pull(y),
                                               alternative = "greater")   
                        # add total counts of a gene in all cells - per sample section/FOV   
                        sum.counts <- 
                          count_matrix %>% .[i, match(obj.bg[[Images(obj.bg)[im]]] %>% Cells(),
                                                      obj.bg %>% Cells()) %>% na.omit()] %>% sum
                        msi$sum.counts <- sum.counts
                        return(msi)
                      }, BPPARAM = BiocParallel::MulticoreParam(20, tasks = 20L,
                                                                force.GC = FALSE,
                                                                progressbar = TRUE)
                      ) %>%
                        # convert/make df data
                        data.table::rbindlist(.) %>%
                        mutate(samples = Images(obj.bg)[im],
                               genes = rownames(obj.bg),
                               genes_group = 
                                 case_when(
                                   grepl(pattern = "Blank-|FP |NegControlProbe|UnassignedCodeword", x = rownames(obj.bg)) ~ "Background",
                                   !grepl(pattern = "Blank-|FP |NegControlProbe|UnassignedCodeword", x = rownames(obj.bg)) ~ "Target"))
                    })

end.time <- Sys.time()
end.time - start.time

# update and make a dataframe
obj.moransI.fast %<>%
  data.table::rbindlist() %>%
  mutate(samples = stringi::stri_replace_all_regex(samples, 
                                                   pattern = "region[0-1].|w[4,5,7]a..|", 
                                                   replacement = ""),
         spatial_tech = stringi::stri_replace_all_regex(samples, 
                                                        pattern = "region[0-1].|w(.)[^:][0-9].|mb2[6-9][3-9].", 
                                                        replacement = ""),
         spatial_statistic = "Morans_I") %>% 
  mutate_if(is.character,
            ~ stringr::str_replace_all(., c("Resolve" = "MC", "Vizgen" = "MERFISH"))) %>% 
  arrange(observed %>% desc)

# export as table .csv 
data.table::fwrite(obj.moransI.fast, 
                   file = "./spatial_techs_compare/fig4/panel_C/spatial_autocorr_moransI.csv")
