# install.packages(c("lidR","terra"))  # if needed
library(lidR)
library(terra)

# ------------------------------------------------------------------------
# 🏔️  LIDAR DTM BUILDER (PER TILE) — unit-aware + skips all-water tiles
# Created by Burch Fisher 10/15/25
# ------------------------------------------------------------------------
# What this does
# 1) Loads a folder of LAS/LAZ tiles into a LAScatalog.
# 2) Filters to **class 2 (ground)** only.
# 3) Processes **per source file** with a **buffer** for smooth edges.
# 4) Builds per-tile DTMs at your requested resolution (given in **meters**).
#    • Converts the meter inputs to the dataset’s **native XY units** (m or US-ft).
#    • Tries TIN first; if TIN fails for a tile/chunk, automatically retries with IDW.
# 5) **Skips tiles/chunks with no ground points** (e.g., all water) so the run never crashes.
# 6) Writes one GeoTIFF per input tile using `{*}_DTM_<res_m>m.tif`.
#
# Notes:
# • Set `res_m` and `buffer_m` in **meters**; the script converts them to native units.
# • With chunk_size = 0, each source file is processed as a single chunk (+ buffer).
# ------------------------------------------------------------------------

# ----------------------------
# CONFIG (set in METERS)
# ----------------------------
laz_dir  <- "/path/to/input/LAZ/"
out_dir  <- "/path/to/output/DTM"
res_m    <- 3
buffer_m <- 50      # try 50–100 m for TIN in steep terrain
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helpers: unit detection & conversion
# ----------------------------
USFT_TO_M <- 1200/3937         # 0.3048006096 m
M_TO_USFT <- 1 / USFT_TO_M

.get_proj_text <- function(x) {
  # Safely extract CRS text from LAS/LAScatalog
  p <- tryCatch(projection(x), error = function(e) NULL)
  if (is.null(p)) "" else as.character(p)
}

detect_us_survey_foot <- function(x) {
  w <- .get_proj_text(x)
  if (!nzchar(w)) return(NA)
  if (grepl("US.?survey.?foot|Foot_US|US[-_ ]?ft|US[-_ ]?Survey[-_ ]?Foot", w, ignore.case = TRUE)) return(TRUE)
  if (grepl("metre|meter|\\bunit\\s*=\\s*metre", w, ignore.case = TRUE)) return(FALSE)
  NA
}

meters_to_native <- function(x_meters, las_or_ctg) {
  det <- detect_us_survey_foot(las_or_ctg)
  is_usft <- if (is.na(det)) FALSE else det
  if (is_usft) x_meters * M_TO_USFT else x_meters
}

# ----------------------------
# Catalog setup (per-tile, buffered)
# ----------------------------
ctg <- readLAScatalog(laz_dir)

opt_filter(ctg)           <- "-keep_class 2"   # ground only
opt_chunk_size(ctg)       <- 0                 # process per SOURCE FILE
opt_chunk_buffer(ctg)     <- meters_to_native(buffer_m, ctg)  # convert meters -> native XY
opt_laz_compression(ctg)  <- TRUE
opt_progress(ctg)         <- TRUE

# Output file template uses the requested metric resolution (meters) in the name
opt_output_files(ctg)     <- file.path(out_dir, sprintf("{*}_DTM_%dm", res_m))

# Convert resolution (meters -> native XY units) for the actual computation
res_native <- meters_to_native(res_m, ctg)

# Warn if CRS missing (we’ll assume meters for XY if so)
proj_txt <- .get_proj_text(ctg)
if (!nzchar(proj_txt)) {
  message("⚠️  Catalog has no CRS. Assuming meters for XY units. ",
          "If data are in US survey feet, results will use meter values as given.")
}

# ----------------------------
# Per-chunk DTM builder with a guard
# ----------------------------
build_dtm_chunk <- function(las, res_native) {
  # Because we filtered to class 2 at read, total points == ground points.
  if (npoints(las) == 0) {
    e <- tryCatch(lidR::lasextent(las), error = function(...) NULL)
    if (!is.null(e)) {
      message(sprintf("⛵ No ground points in chunk (%.1f,%.1f)-(%.1f,%.1f). Skipping.",
                      e@xmin, e@ymin, e@ymax, e@ymax))
    } else {
      message("⛵ No ground points in chunk. Skipping.")
    }
    return(NULL)
  }
  
  # Try TIN, fall back to IDW if needed (per chunk/tile)
  r <- tryCatch(
    rasterize_terrain(las, res = res_native, algorithm = tin(), na.rm = TRUE),
    error = function(e) {
      message("TIN failed here: ", conditionMessage(e), " → retrying with IDW")
      rasterize_terrain(las, res = res_native, algorithm = knnidw(k = 10, p = 2), na.rm = TRUE)
    }
  )
  
  if (!inherits(r, "SpatRaster")) {
    message("⚠️  Could not build DTM for this chunk even with IDW. Skipping.")
    return(NULL)
  }
  
  r
}

# ----------------------------
# Run (writes one DTM per source tile)
# ----------------------------
invisible(
  catalog_map(
    ctg,
    build_dtm_chunk,
    res_native = res_native
  )
)

message(sprintf("✅ Done. Per-tile DTMs written to: %s (cell size: %g %s)",
                normalizePath(out_dir),
                res_m,
                "meters (file names reflect your requested metric resolution)"))
