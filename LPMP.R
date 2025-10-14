# ------------------------------------------------------------------------
# 🦙 LIDAR PIXEL METRICS PIPELINE
#    (All inputs in METERS; optional global DTM; lidR ≥ 4.1.x)
# Created by Burch Fisher 10/14/25
# ------------------------------------------------------------------------

# What this script does
# 1) Reads a FOLDER of LAS/LAZ tiles into a LAScatalog (no '*.laz' glob).
# 2) Filters to classes 1–5 and drops withheld/overlap at read.
# 3) Normalizes heights either:
#      • useGlobalDTM = TRUE  → build ONE DTM for the whole area (SAVED), then normalize
#      • useGlobalDTM = FALSE → per-tile TIN normalization in memory (no DTM on disk)
# 4) Removes exact duplicates (and optional near-duplicates).
# 5) Computes per-pixel metrics; writes one multilayer GeoTIFF per tile
#    with band names matching metric names (CamelCase where requested).
# 6) **Automatically skips tiles with no ground points (e.g., all water) to avoid errors.**
#
# Units policy
# • You specify ALL sizes in METERS (px_res_m, dtm_res_m, chunk_size_m, etc.).
# • The script converts those to native XY units (feet/meters) where lidR needs it.
# • If input data are US survey feet, Z is converted to METERS inside myMetrics(),
#   so all height metrics are in meters in the outputs.
# • Output rasters keep the input CRS; reproject later if you need XY in meters.
#
# macOS/RStudio note
# • CRAN binary of lidR on macOS is often built w/o OpenMP → internal threads disabled.
#   We set threads to 1 to avoid warnings. You can still parallelize across areas in R.
# ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(lidR)      # ≥ 4.1.x
  library(lubridate)
  library(terra)
})

# ============================================================================
# 1) USER SETTINGS (ALL IN METERS)
# ============================================================================

# --- Input/Output Paths ---
# IMPORTANT: in_path must be a DIRECTORY containing .las/.laz files, not a single file
in_path        <- "/path/to/delivery_area/LAZ"  # Source directory with LiDAR tiles
out_dir        <- "/path/to/outdir"              # Destination for processed outputs
out_prefix     <- "metrics"                     # Prefix added to output filenames

# --- Normalization Method ---
# Controls how ground height normalization is performed:
#   TRUE  = Create a single DTM for the entire dataset, save to disk, and reuse for all tiles
#   FALSE = Generate a temporary TIN for each tile independently (no DTM saved)
useGlobalDTM   <- FALSE

# --- Spatial Resolution ---
px_res_m       <- 5    # Pixel size for output raster metrics (meters)
dtm_res_m      <- 3     # DTM resolution (meters) - only used when useGlobalDTM = TRUE

# --- Tile Processing Strategy ---
# Controls whether large tiles are subdivided during processing:
#   0  = Process each input tile as a single chunk (no subdivision)
#   >0 = Split tiles into sub-chunks of this size (meters) for memory efficiency
chunk_size_m   <- 0
chunk_buffer_m <- 30    # Buffer around each chunk (meters) - ignored when chunk_size_m = 0

# --- Performance ---
n_threads      <- 4     # Number of parallel threads (set to 1 on macOS CRAN builds without OpenMP)

# --- Point Cloud Filtering (Optional) ---
# Remove near-duplicate points before processing:
use_eps_dedup  <- FALSE   # Enable/disable epsilon-based deduplication
eps_xy_m       <- 0.01    # Horizontal tolerance (meters) - points closer than this are duplicates
eps_z_m        <- 0.01    # Vertical tolerance (meters)

# --- Coordinate System ---
# Override automatic unit detection from CRS:
#   "auto"  = Detect from coordinate system metadata
#   "us-ft" = Force US survey feet
#   "m"     = Force meters
force_units    <- "auto"

# --- Initialize Output Directory ---
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------- 2) UNIT HELPERS --------------------------------------------------
USFT_TO_M <- 1200/3937         # 0.3048006096 m
M_TO_USFT <- 1 / USFT_TO_M

.get_proj_text <- function(x) {
  # Safely extract CRS WKT/PROJ text from LAS/LAScatalog
  tryCatch({
    p <- projection(x)
    if (is.null(p)) "" else as.character(p)
  }, error = function(e) "")
}

detect_us_survey_foot <- function(x) {
  # Read CRS text and detect US survey foot
  w <- .get_proj_text(x)
  if (!nzchar(w)) return(NA)
  if (grepl("US.?survey.?foot|Foot_US|US[-_ ]?ft|US[-_ ]?Survey[-_ ]?Foot", w, ignore.case = TRUE)) return(TRUE)
  if (grepl("metre|meter|\\bunit\\s*=\\s*metre", w, ignore.case = TRUE)) return(FALSE)
  NA
}

meters_to_native <- function(x_meters, las_or_ctg) {
  # Convert a value in meters to the LAS/CTG’s native XY units
  is_usft <- switch(
    force_units,
    "us-ft" = TRUE,
    "m"     = FALSE,
    "auto"  = {
      det <- detect_us_survey_foot(las_or_ctg)
      if (is.na(det)) FALSE else det
    }
  )
  if (is_usft) x_meters * M_TO_USFT else x_meters
}

# -------- 3) METRICS FUNCTION (Z→meters if needed) -------------------------
myMetrics <- function(z, rn, i, a, g, c) {
  z  <- z[!is.na(z)]
  nz <- length(z)
  
  angle_missing <- is.null(a) || length(a) == 0 || all(is.na(a))
  safe_mean <- function(x) if (length(x) == 0) NA_real_ else mean(x, na.rm = TRUE)
  safe_max  <- function(x) if (length(x) == 0) NA_real_ else suppressWarnings(max(x, na.rm = TRUE))
  safe_min  <- function(x) if (length(x) == 0) NA_real_ else suppressWarnings(min(x, na.rm = TRUE))
  q         <- function(p) if (nz == 0) NA_real_ else as.numeric(quantile(z, probs = p, na.rm = TRUE, names = FALSE))
  pct_bin   <- function(lo, hi) if (nz == 0) NA_real_ else 100 * sum(z >= lo & z < hi, na.rm = TRUE) / nz
  
  has_g <- !(is.null(g) || length(g) == 0 || all(is.na(g)))
  maxg  <- if (has_g) safe_max(g) else NA_real_
  ming  <- if (has_g) safe_min(g) else NA_real_
  
  list(
    NoPoints      = nz,
    meanIntensity = safe_mean(i),
    meanAngle     = if (angle_missing) NA_real_ else safe_mean(a),
    maxAngle      = if (angle_missing) NA_real_ else safe_max(a),
    minAngle      = if (angle_missing) NA_real_ else safe_min(a),
    
    iqr     = if (nz == 0) NA_real_ else IQR(z, na.rm = TRUE),
    meanHt  = safe_mean(z),
    maxHt   = safe_max(z),
    
    p10 = q(0.10), p20 = q(0.20), p25 = q(0.25), p30 = q(0.30),
    p40 = q(0.40), p50 = q(0.50), p60 = q(0.60), p70 = q(0.70),
    p75 = q(0.75), p80 = q(0.80), p90 = q(0.90), p95 = q(0.95), p99 = q(0.99),
    
    v1  = pct_bin(0,1), v2=pct_bin(1,2), v3=pct_bin(2,3), v4=pct_bin(3,4),
    v5  = pct_bin(4,5), v6=pct_bin(5,6), v9=pct_bin(6,9), v12=pct_bin(9,12),
    v15 = pct_bin(12,15), v18=pct_bin(15,18), v21=pct_bin(18,21),
    v24 = pct_bin(21,24), v27=pct_bin(24,27), v30=pct_bin(27,30),
    v33 = pct_bin(30,33), v36=pct_bin(33,36), v39=pct_bin(36,39),
    v42 = pct_bin(39,42), v45=pct_bin(42,45), v48=pct_bin(45,48), v51=pct_bin(48,51),
    
    maxTimeInt = if (!has_g) NA_integer_ else {
      as.integer(gsub("-", "", substr(lubridate::as_datetime(1315576000 + maxg), 3, 10)))
    },
    minTimeInt = if (!has_g) NA_integer_ else {
      as.integer(gsub("-", "", substr(lubridate::as_datetime(1315576000 + ming), 3, 10)))
    },
    diffTime = if (!has_g) NA_real_ else (maxg - ming) / 86400
  )
}


# -------- 4) DEDUP HELPER (tolerances passed in METERS) --------------------
dedup_eps <- function(las, eps_xy_m, eps_z_m) {
  if (npoints(las) == 0) return(las)
  eps_xy_native <- meters_to_native(eps_xy_m, las)
  eps_z_native  <- meters_to_native(eps_z_m,  las)
  d <- las@data
  key <- paste(
    round(d$X / eps_xy_native) * eps_xy_native,
    round(d$Y / eps_xy_native) * eps_xy_native,
    round(d$Z / eps_z_native)  * eps_z_native
  )
  keep <- !duplicated(key)
  lasfilter(las, keep)
}

# -------- 5) GROUND-CHECK HELPER ------------------------------------------
has_ground_points <- function(las) {
  # class 2 = ground
  !is.null(las$Classification) && any(las$Classification == 2L, na.rm = TRUE)
}

# -------- 6) BUILD CATALOG -------------------------------------------------
stopifnot(dir.exists(in_path))
ctg <- readLAScatalog(in_path)

# Keep classes 1–5; drop withheld/overlap; select attributes
opt_filter(ctg)           <- "-keep_class 1 2 3 4 5 -drop_withheld -drop_overlap"
opt_select(ctg)           <- "*"
opt_independent_files(ctg) <- TRUE

# Chunking controls (convert meter inputs to native XY units); 0 = per tile
opt_chunk_size(ctg)   <- if (chunk_size_m == 0) 0 else meters_to_native(chunk_size_m, ctg)
opt_chunk_buffer(ctg) <- meters_to_native(chunk_buffer_m, ctg)

# Output template for per-tile raster results
opt_output_files(ctg) <- file.path(out_dir, paste0(out_prefix, "_{*}"))

# Threads (mac CRAN usually no OpenMP → keep 1)
try({ set_lidr_threads(n_threads) }, silent = TRUE)

# Warn if CRS is missing
proj_txt <- .get_proj_text(ctg)
if (!nzchar(proj_txt) && force_units == "auto") {
  warning("⚠️  Catalog has no CRS. Assuming meters for XY+Z. If data are in US ft, set force_units='us-ft'.")
}

# -------- 7) (OPTIONAL) BUILD & SAVE GLOBAL DTM ----------------------------
# If useGlobalDTM = TRUE, build one TIN-based DTM and save to out_dir/<prefix>_DTM.tif
if (isTRUE(useGlobalDTM)) {
  dtm_res_native <- meters_to_native(dtm_res_m, ctg)
  message(sprintf("🛠️  Building global DTM at %.2f m resolution...", dtm_res_m))
  dtm <- rasterize_terrain(ctg, res = dtm_res_native, algorithm = tin())
  
  dtm_path <- file.path(out_dir, paste0(out_prefix, "_DTM.tif"))
  writeRaster(dtm, dtm_path, overwrite = TRUE)
  message(sprintf("💾 Saved DTM → %s", dtm_path))
} else {
  dtm <- NULL
}

# -------- 8) PER-TILE PROCESSING FUNCTION ----------------------------------
process_tile <- function(las, res_native, dtm, useGlobalDTM) {
  if (npoints(las) == 0) return(NULL)
  
  # --- If using per-tile TIN, ensure the tile actually has ground points ---
  if (!isTRUE(useGlobalDTM) && !has_ground_points(las)) {
    message("⛵ Tile has no class-2 ground points (likely all water). Skipping.")
    return(NULL)
  }
  
  # --- Normalize heights (robust to failures) ---
  las <- tryCatch({
    if (isTRUE(useGlobalDTM) && !is.null(dtm)) {
      normalize_height(las, dtm, na.rm = TRUE, method = "bilinear")
    } else {
      normalize_height(las, tin(), na.rm = TRUE)
    }
  }, error = function(e) {
    message("⚠️ normalize_height failed on this tile: ", conditionMessage(e), " — skipping tile.")
    return(NULL)
  })
  
  if (is.null(las) || npoints(las) == 0) {
    message("ℹ️ Tile has zero points after normalization (all NA ground?); skipping.")
    return(NULL)
  }
  
  # --- De-duplication ---
  las <- filter_duplicates(las)
  if (use_eps_dedup) las <- dedup_eps(las, eps_xy_m, eps_z_m)
  
  # --- If Z is in US survey feet, convert Z -> meters (post-normalization) ---
  detected <- detect_us_survey_foot(las)
  is_usft_flag <- switch(
    force_units,
    "us-ft" = TRUE, "m" = FALSE,
    "auto"  = if (is.na(detected)) FALSE else detected
  )
  if (isTRUE(is_usft_flag)) {
    las@data$Z <- las@data$Z * USFT_TO_M
  }
  
  # --- Ensure ScanAngle exists if only ScanAngleRank is present ---
  if (is.null(las$ScanAngle) && !is.null(las$ScanAngleRank)) {
    las@data$ScanAngle <- las@data$ScanAngleRank
  }
  
  # --- Pixel metrics (heights guaranteed meters) ---
  pixel_metrics(
    las,
    ~ myMetrics(
      z  = Z,
      rn = ReturnNumber,
      i  = Intensity,
      a  = ScanAngle,
      g  = gpstime,
      c  = Classification
    ),
    res = res_native
  )
}

# -------- 9) RUN PIPELINE --------------------------------------------------
px_res_native <- meters_to_native(px_res_m, ctg)

message(sprintf(
  "🚀 Computing pixel metrics at %.2f m res (%s)...",
  px_res_m, if (isTRUE(useGlobalDTM)) "global DTM" else "per-tile TIN"
))

rlist <- catalog_map(
  ctg,
  process_tile,
  res_native   = px_res_native,
  dtm          = dtm,
  useGlobalDTM = useGlobalDTM
)

message(sprintf(
  "✅ Done! Per-tile GeoTIFFs → %s\n   Pixel size: %.2f m%s",
  normalizePath(out_dir),
  px_res_m,
  if (isTRUE(useGlobalDTM)) sprintf(" | DTM saved: %s", file.path(out_dir, paste0(out_prefix, "_DTM.tif"))) else ""
))

# ------------------------------------------------------------------------
# End of Script
# ------------------------------------------------------------------------
