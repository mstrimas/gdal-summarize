library(raster)
library(prioritizr)
library(rasterVis)
library(ebirdst)
set.seed(1)

# small example ----

# 2 files with 2 bands each with dimensions 3x3
# filled with count data, including zeros and missing values
nr <- 3
nc <- 3
for (i in 1:2) {
  r1 <- raster(matrix(rpois(nr * nc, lambda = 1), nrow = nr))
  r2 <- raster(matrix(rpois(nr * nc, lambda = 10), nrow = nr))
  
  # add some zeros
  s <- sample(ncell(r1), 5)
  r1[s[1:2]] <- 0
  r2[s[2:3]] <- 0
  r1[s[3:4]] <- NA
  r2[s[4:5]] <- NA
  s <- stack(r1, r2)
  
  # save
  writeRaster(s, paste0("data/small", i, ".tif"), overwrite = TRUE)
}

# large example ----

# 2 files with 9 bands each with dimensions 250x250
# filled with simulate species occupancy data
r <- raster(ncol = 1000, nrow = 1000, xmn = 0, xmx = 1, ymn = 0, ymx = 1)
values(r) <- 1
sim1 <- simulate_species(r, 9)
# add some zeros and missing data
sim1[sim1 < quantile(sim1, 0.1)[, 1]] <- 0
sim1[sim1 > quantile(sim1, 0.9)[, 1]] <- NA
sim2 <- simulate_species(r, 9)
# add some zeros and missing data
sim2[sim2 < quantile(sim2, 0.1)[, 1]] <- 0
sim2[sim2 > quantile(sim2, 0.9)[, 1]] <- NA
writeRaster(sim1, "data/large1.tif", overwrite = TRUE)
writeRaster(sim2, "data/large2.tif", overwrite = TRUE)

# ebirdst ----

# download status and trends data
abd <- get_species_path("woothr") %>% 
  load_raster("abundance", .)
# subset to the breeding season
dts <- parse_raster_dates(abd)
r <- abd[[which(dts >= as.Date("2018-05-24") & dts <= as.Date("2018-09-07"))]]
r <- writeRaster(r, "data/woothr.tif", options = "COMPRESS=LZW")