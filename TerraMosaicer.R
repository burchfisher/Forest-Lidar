# ------------------------------------------------------------------------
# 🧩  DTM / MULTIBAND METRIC MOSAIC BUILDER
# Created by Burch Fisher 10/15/25
# ------------------------------------------------------------------------
# What this script does
# 1) Reads all GeoTIFF tiles from a given folder (e.g., per-tile DTMs or
#    per-tile metric stacks from the LiDAR pipeline).
# 2) Automatically handles both single-band (e.g., DTMs, CHMs) and
#    multi-band rasters (e.g., pixel-metrics outputs with many variables).
# 3) Checks that all rasters share the same number of bands and layer names.
# 4) Merges all tiles into a single seamless raster mosaic using the chosen
#    function for overlapping areas (e.g., `mean`, `min`, `max`).
# 5) Writes the resulting mosaic to disk as a GeoTIFF.
#
# Notes:
# • All input rasters must share CRS, resolution, and alignment.
# • The merge function defines how overlapping pixels are resolved:
#      - `mean`  → smooths overlaps (best for continuous data)
#      - `min`   → keeps lowest values (e.g., DEMs with edge spikes)
#      - `max`   → keeps highest values (e.g., canopy surfaces)
# • For very large mosaics, you can stream directly to disk by adding:
#      `filename = out_tif, overwrite = TRUE` inside `mosaic()` at line 59.
# ------------------------------------------------------------------------

library(terra)

# ----------------------------
# CONFIG
# ----------------------------
# Folder containing GeoTIFF tiles
in_dir  <- "/path/to/input/files"

# Output mosaic path
out_tif <- "/path/to/output/file/delivery_area_epsg6346_DTM_3m.tif"

# ----------------------------
# PROCESS
# ----------------------------
# Find all GeoTIFFs
tif_files <- list.files(in_dir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(tif_files) > 0)

# Read as SpatRaster objects
rasters <- lapply(tif_files, rast)

# Verify all rasters have the same number of layers
nlayers_each <- sapply(rasters, nlyr)
if (length(unique(nlayers_each)) > 1) {
  stop("❌ Rasters have different numbers of layers — cannot mosaic safely.")
}

# Warn if layer names differ
layer_names <- sapply(rasters, function(r) paste(names(r), collapse = ","))
if (length(unique(layer_names)) > 1) {
  warning("⚠️ Layer names differ among rasters — bands may not align perfectly.")
}

# Mosaic all bands together using mean (change to min/max if preferred)
mosaic_rast <- do.call(mosaic, c(rasters, fun = mean))
# mosaic_rast <- do.call(mosaic, c(rasters, fun = mean, filename = out_tif, overwrite = TRUE))

# Write to disk
writeRaster(mosaic_rast, out_tif, overwrite = TRUE)

cat("✅ Mosaic written to:", out_tif, "\n")
cat("📊 Bands in output:", paste(names(mosaic_rast), collapse = ", "), "\n")
