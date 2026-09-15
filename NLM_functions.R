### NLM Functions ###
# A simple random NLM: every cell is drawn independently.
nlm_simple <- function(ncol, nrow, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- rast(ncols = ncol, nrows = nrow,
            xmin = 0, xmax = ncol, ymin = 0, ymax = nrow)
  values(x) <- runif(ncell(x))
  names(x) <- "value"
  x
}

# A multiscale, fractal-like NLM. Larger H values place more weight on
# broad-scale variation and therefore produce stronger spatial autocorrelation.
nlm_fractal <- function(ncol, nrow, H = 0.5, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  noise <- ambient::noise_perlin(
    dim = c(nrow, ncol),
    frequency = (1 - H) / 10,
    octaves = 8,
    lacunarity = 3
  )
  
  raster_fractal <- terra::rast(noise)
  terra::ext(raster_fractal) <- c(0, ncol, 0, nrow)
  
  minimum <- terra::global(
    raster_fractal, "min", na.rm = TRUE
  )[1, 1]
  
  maximum <- terra::global(
    raster_fractal, "max", na.rm = TRUE
  )[1, 1]
  
  raster_fractal <-
    (raster_fractal - minimum) / (maximum - minimum)
  
  names(raster_fractal) <- paste0("H_", H)
  
  return(raster_fractal)
}

# Estimate an empirical semivariogram from a random sample of raster cells.
# For the NLMs below, distance is measured in cell-width units.
semivariogram <- function(x, sample_size = 800,
                          n_bins = 18, seed = 1) {
  
  set.seed(seed)
  
  pts <- terra::spatSample(
    x,
    size = sample_size,
    method = "random",
    na.rm = TRUE,
    xy = TRUE,
    as.df = TRUE
  )
  
  z <- pts[[ncol(pts)]]
  
  pair_distance <- as.vector(
    dist(pts[, c("x", "y")])
  )
  
  # Divide by the variance so a sill is approximately 1.
  pair_semivariance <- 
    0.5 * as.vector(dist(z))^2 / var(z)
  
  maximum_distance <- quantile(
    pair_distance,
    probs = 0.75,
    names = FALSE
  )
  
  breaks <- seq(
    0,
    maximum_distance,
    length.out = n_bins + 1
  )
  
  bin <- cut(
    pair_distance,
    breaks = breaks,
    include.lowest = TRUE,
    labels = FALSE
  )
  
  data.frame(
    distance = as.numeric(
      tapply(pair_distance, bin, mean)
    ),
    semivariance = as.numeric(
      tapply(pair_semivariance, bin, mean)
    )
  )
}

plot_variogram <- function(v, main = "Semivariogram") {
  plot(v$distance, v$semivariance, type = "b", pch = 16,
       xlab = "Distance", ylab = "Semivariance", main = main)
}