library(macrosheds)
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(httr)

sites <- c('ALBION', 'GREEN4', 'NAVAJO', 'w8', 'w9', 'W-9')

sheds <- ms_load_spatial_product(
  macrosheds_root = '~/ssd2/macrosheds_stuff/ms_test/',
  spatial_product = 'ws_boundary',
  site_codes = sites
)

# ---- configuration ----

nadp_dir <- 'nadp'
dir.create(nadp_dir, showWarnings = FALSE)

# NADP NTN annual raster base URL
# Variables: pH (lab pH), Cond (conductivity, uS/cm)
# Files are named like: pH_dep_2022.tif, Cond_dep_2022.tif
# Available from NADP's gridded data: https://nadp.slh.wisc.edu/maps-data/ntn-gradient-maps/

nadp_vars <- c('pH', 'Cond')
years <- 1985:2023  # adjust range as needed

# Base URL pattern for NADP NTN concentration grids (annual)
# NADP provides these as GeoTIFFs via their data portal
base_url <- 'https://nadp.slh.wisc.edu/datalib/ntn/grids/annual'

# ---- prepare watershed geometries ----

# Reproject watersheds to WGS84 for consistency with NADP rasters
sheds_wgs84 <- st_transform(sheds, 4326)

# Get a bounding box with some buffer for cropping
sheds_bbox <- st_bbox(sheds_wgs84)

# ---- download NADP rasters ----

download_nadp_raster <- function(variable, year, dest_dir) {

    # NADP NTN annual concentration grids URL pattern
    # pH: precipitation-weighted mean pH
    # Cond: precipitation-weighted mean conductivity (uS/cm)
    fname <- paste0(variable, '_', year, '.tif')
    dest_file <- file.path(dest_dir, fname)

    if(file.exists(dest_file)) {
        message('Already downloaded: ', fname)
        return(dest_file)
    }

    # Try common NADP URL patterns
    urls_to_try <- c(
        paste0(base_url, '/', variable, '/', year, '/', fname),
        paste0(base_url, '/', variable, '/', fname),
        paste0(base_url, '/', year, '/', fname)
    )

    for(url in urls_to_try) {
        resp <- try(GET(url, write_disk(dest_file, overwrite = TRUE),
                        timeout(60)), silent = TRUE)

        if(! inherits(resp, 'try-error') && status_code(resp) == 200) {
            # Verify it's a valid raster
            test <- try(rast(dest_file), silent = TRUE)
            if(! inherits(test, 'try-error')) {
                message('Downloaded: ', fname)
                return(dest_file)
            }
        }

        # Clean up failed download
        if(file.exists(dest_file)) file.remove(dest_file)
    }

    message('Could not download: ', variable, ' ', year)
    return(NA_character_)
}

message('Downloading NADP NTN annual rasters...')
message('If automatic download fails, manually download rasters from:')
message('  https://nadp.slh.wisc.edu/maps-data/ntn-gradient-maps/')
message('Place .tif files in the nadp/ directory as {Variable}_{Year}.tif')

download_log <- expand.grid(variable = nadp_vars, year = years,
                            stringsAsFactors = FALSE) %>%
    as_tibble()

download_log$filepath <- mapply(
    download_nadp_raster,
    variable = download_log$variable,
    year = download_log$year,
    dest_dir = nadp_dir
)

download_log <- download_log %>%
    filter(! is.na(filepath))

if(nrow(download_log) == 0) {
    stop('No NADP rasters were downloaded. Please download manually from:\n',
         '  https://nadp.slh.wisc.edu/maps-data/ntn-gradient-maps/\n',
         'and place them in nadp/ as {Variable}_{Year}.tif\n',
         'Then re-run this script.')
}

message(nrow(download_log), ' rasters available for processing.')

# ---- also check for any manually placed rasters ----

manual_tifs <- list.files(nadp_dir, pattern = '\\.tif$', full.names = TRUE)
if(length(manual_tifs) > 0) {
    # Parse variable and year from filenames
    manual_info <- tibble(filepath = manual_tifs) %>%
        mutate(
            basename = tools::file_path_sans_ext(basename(filepath)),
            variable = sub('_[0-9]{4}$', '', basename),
            year = as.integer(sub('.*_', '', basename))
        ) %>%
        filter(variable %in% nadp_vars, ! is.na(year)) %>%
        select(variable, year, filepath)

    download_log <- bind_rows(download_log, manual_info) %>%
        distinct(variable, year, .keep_all = TRUE)
}

# ---- extract raster values for each watershed ----

extract_nadp_for_sheds <- function(raster_path, sheds_sf) {

    r <- try(rast(raster_path), silent = TRUE)
    if(inherits(r, 'try-error')) return(NULL)

    # Reproject watersheds to match raster CRS
    sheds_reproj <- st_transform(sheds_sf, crs(r))

    # Crop raster to watershed extent (with buffer) to speed extraction
    sheds_extent <- ext(vect(sheds_reproj)) + 0.5  # ~0.5 degree buffer
    r_crop <- try(crop(r, sheds_extent), silent = TRUE)
    if(inherits(r_crop, 'try-error')) r_crop <- r

    # Extract area-weighted mean for each polygon
    vals <- exact_extract(r_crop, sheds_reproj, fun = 'mean')

    if(is.null(vals)) {
        # Fallback to terra::extract
        vals <- terra::extract(r_crop, vect(sheds_reproj), fun = mean,
                               na.rm = TRUE, weights = TRUE)
        vals <- vals[, 2]
    }

    return(vals)
}

# Use exactextractr if available, otherwise terra::extract
use_exactextractr <- requireNamespace('exactextractr', quietly = TRUE)

if(use_exactextractr) {
    library(exactextractr)
    message('Using exactextractr for weighted extraction.')
} else {
    message('Using terra::extract. Install exactextractr for better accuracy.')
}

results <- list()

for(i in seq_len(nrow(download_log))) {

    row <- download_log[i, ]
    message('Processing: ', row$variable, ' ', row$year)

    r <- try(rast(row$filepath), silent = TRUE)
    if(inherits(r, 'try-error')) {
        message('  Skipping (invalid raster)')
        next
    }

    # Reproject watersheds to match raster CRS
    sheds_reproj <- st_transform(sheds_wgs84, crs(r))

    # Crop to reduce memory
    sheds_ext <- ext(vect(sheds_reproj))
    sheds_ext_buf <- sheds_ext + 0.5
    r_crop <- try(crop(r, sheds_ext_buf), silent = TRUE)
    if(inherits(r_crop, 'try-error')) r_crop <- r

    if(use_exactextractr) {
        vals <- exact_extract(r_crop, sheds_reproj, fun = 'mean')
    } else {
        ex <- terra::extract(r_crop, vect(sheds_reproj), fun = mean,
                             na.rm = TRUE, weights = TRUE)
        vals <- ex[, 2]
    }

    res <- tibble(
        site_code = sheds_wgs84$site_code,
        variable = row$variable,
        year = row$year,
        value = vals
    )

    results[[length(results) + 1]] <- res
}

# ---- combine and format results ----

nadp_summary <- bind_rows(results)

# Add descriptive labels
nadp_summary <- nadp_summary %>%
    mutate(
        variable_description = case_when(
            variable == 'pH'   ~ 'Precipitation-weighted mean pH',
            variable == 'Cond' ~ 'Precipitation-weighted mean conductivity (uS/cm)',
            TRUE ~ variable
        )
    )

# Pivot to wide format (one column per variable)
nadp_wide <- nadp_summary %>%
    select(site_code, year, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    arrange(site_code, year)

# ---- save results ----

write.csv(nadp_summary, 'nadp_watershed_summary_long.csv', row.names = FALSE)
write.csv(nadp_wide, 'nadp_watershed_summary_wide.csv', row.names = FALSE)

message('\n---- Summary ----')
message('Sites: ', paste(unique(nadp_summary$site_code), collapse = ', '))
message('Variables: ', paste(unique(nadp_summary$variable), collapse = ', '))
message('Years: ', min(nadp_summary$year), '-', max(nadp_summary$year))
message('Results saved to:')
message('  nadp_watershed_summary_long.csv')
message('  nadp_watershed_summary_wide.csv')

# Print a quick summary
nadp_summary %>%
    group_by(site_code, variable) %>%
    summarize(
        n_years = n(),
        mean_value = mean(value, na.rm = TRUE),
        min_value = min(value, na.rm = TRUE),
        max_value = max(value, na.rm = TRUE),
        .groups = 'drop'
    ) %>%
    print(n = Inf)

