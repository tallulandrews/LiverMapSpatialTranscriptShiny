
require(Seurat)
slice <- "C73_D1"
set.seed(101)

Normal_obj <- Load10X_Spatial(paste("./", slice, sep=""), 
		filename="filtered_feature_bc_matrix.h5", 
		slice=slice)

Normal_obj@meta.data$orig.ident <- rep(slice, ncol(Normal_obj))


varimax_obj <- readRDS(paste("C:/Users/tandrews/Documents/UHNSonya/Spatial/", slice, "/", slice, "_varimax_Spatial_Seurat.rds", sep=""))

varimax_component <- 1
zonation_score <- varimax_obj@meta.data[,paste("RotPC", varimax_component, sep="_")]

Normal_obj@meta.data$zonation_score <- zonation_score

SpatialFeaturePlot(Normal_obj, features="zonation_score")

saveRDS(Normal_obj, paste("C:/Users/tandrews/Documents/UHNSonya/Spatial/HealthyLiverShinyApp/data/", slice, "_raw.rds", sep=""))