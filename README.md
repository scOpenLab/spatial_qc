## This repo contains scripts used for spatial QC steps and can also be used for specific analysis parts
`./scripts` contains the R scripts, details can be found within each script documentaitons
  - `crop_seurat_v1.R` allows user to crop (or keep) Seurat object FOV part or a dataframe
    - discussion can be found on [Seurat issue here](https://github.com/satijalab/seurat/issues/8457)
  - `spatial_autocorr_moransI_v1.R` this runs spatial autocorrelation (Moran's I) using [`moranfast` R package](https://github.com/mcooper/moranfast)
  - `subset_obj_seurat_v2.R` optimized subset function that can subset Seurat object FOVs for when > 1 FOVs or samples are present
    - see [this discussion](https://github.com/satijalab/seurat/issues/7462) and [this issue](https://github.com/satijalab/seurat/issues/6409)   
