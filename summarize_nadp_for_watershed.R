library(macrosheds)

sites <- c('ALBION', 'GREEN4', 'NAVAJO', 'w8', 'w9', 'W-9')

sheds <- ms_load_spatial_product(
  macrosheds_root = '~/ssd2/macrosheds_stuff/ms_test/',
  spatial_product = 'ws_boundary',
  site_codes = sites
)

