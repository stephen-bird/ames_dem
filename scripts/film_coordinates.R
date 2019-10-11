
# This script calculates the film coordinates of the fiducial marks based on the calibration report.


# This is a Packrat project.

library(here)
library(dplyr)
source("scripts/pp_functions.R")


film_coord <-
        data.frame(
                fids = c("1", "2", "3", "4"),
                x = c(-105.993, 106.001, -106.015, 106.015),
                y = c(-106.006, 106.002, 105.982, -106.006)
        )


fid_dist <-
        data.frame(
                fids = c("1_2", "1_3", "1_4", "2_3", "2_4", "3_4"),
                # For example, 1_2 is distance between fiducial 1 and 2
                dist_mm = c(299.814, 211.988, 212.015, 212.008, 212.007, 299.826)
        )

# Enter the error in the film measurements as stated in the report. If none is
# given, estimate a reasonable value:
dist_err <- 0.003 # Stated error between two points
point_err <- sqrt((dist_err ^ 2) / 2) # Error of placing a single point
# Enter the angle between lines intersecting the fiducial marks
fid_angle <- 90 + 13 / 3600

# ------------------------------------------------------------------------------

# How many samples to generate?
n <- 10000

# Generate normal random samples of each x,y coordinate as given in the
# calibration report using the stated measurement error (e.g. sfilm_1x = sample
# film, fiducial 1, x):
set.seed(42)
sfilm_1x <- rnorm(n, mean = film_coord$x[1], sd = point_err)
sfilm_2x <- rnorm(n, mean = film_coord$x[2], sd = point_err)
sfilm_3x <- rnorm(n, mean = film_coord$x[3], sd = point_err)
sfilm_4x <- rnorm(n, mean = film_coord$x[4], sd = point_err)
sfilm_1y <- rnorm(n, mean = film_coord$y[1], sd = point_err)
sfilm_2y <- rnorm(n, mean = film_coord$y[2], sd = point_err)
sfilm_3y <- rnorm(n, mean = film_coord$y[3], sd = point_err)
sfilm_4y <- rnorm(n, mean = film_coord$y[4], sd = point_err)

# Gather into a df:
sfilm_coord <-
        data.frame(sfilm_1x,
                   sfilm_2x,
                   sfilm_3x,
                   sfilm_4x,
                   sfilm_1y,
                   sfilm_2y,
                   sfilm_3y,
                   sfilm_4y)


# Make distance measurements among all fiducial points using the sample coordinates generated above:
sfilm_coord$sdist_12 <- hypo(sfilm_coord$sfilm_1x,sfilm_coord$sfilm_1y,sfilm_coord$sfilm_2x,sfilm_coord$sfilm_2y)
sfilm_coord$sdist_13 <- hypo(sfilm_coord$sfilm_1x,sfilm_coord$sfilm_1y,sfilm_coord$sfilm_3x,sfilm_coord$sfilm_3y)
sfilm_coord$sdist_14 <- hypo(sfilm_coord$sfilm_1x,sfilm_coord$sfilm_1y,sfilm_coord$sfilm_4x,sfilm_coord$sfilm_4y)
sfilm_coord$sdist_23 <- hypo(sfilm_coord$sfilm_2x,sfilm_coord$sfilm_2y,sfilm_coord$sfilm_3x,sfilm_coord$sfilm_3y)
sfilm_coord$sdist_24 <- hypo(sfilm_coord$sfilm_2x,sfilm_coord$sfilm_2y,sfilm_coord$sfilm_4x,sfilm_coord$sfilm_4y)
sfilm_coord$sdist_34 <- hypo(sfilm_coord$sfilm_3x,sfilm_coord$sfilm_3y,sfilm_coord$sfilm_4x,sfilm_coord$sfilm_4y)

# Generate normal random samples based on film measured distances and the calibrated measurement error:
set.seed(24)
sfilm_coord$fdist_12 <- rnorm(n, mean = fid_dist$dist_mm[1], sd = 0.003)
sfilm_coord$fdist_13 <- rnorm(n, mean = fid_dist$dist_mm[2], sd = 0.003)
sfilm_coord$fdist_14 <- rnorm(n, mean = fid_dist$dist_mm[3], sd = 0.003)
sfilm_coord$fdist_23 <- rnorm(n, mean = fid_dist$dist_mm[4], sd = 0.003)
sfilm_coord$fdist_24 <- rnorm(n, mean = fid_dist$dist_mm[5], sd = 0.003)
sfilm_coord$fdist_34 <- rnorm(n, mean = fid_dist$dist_mm[6], sd = 0.003)

# Calculate residuals^2 between distances measured on the film and distances sampled above:
sfilm_coord$res_12 <- (sfilm_coord$sdist_12 - sfilm_coord$fdist_12) ^ 2
sfilm_coord$res_13 <- (sfilm_coord$sdist_13 - sfilm_coord$fdist_13) ^ 2
sfilm_coord$res_14 <- (sfilm_coord$sdist_14 - sfilm_coord$fdist_14) ^ 2
sfilm_coord$res_23 <- (sfilm_coord$sdist_23 - sfilm_coord$fdist_23) ^ 2
sfilm_coord$res_24 <- (sfilm_coord$sdist_24 - sfilm_coord$fdist_24) ^ 2
sfilm_coord$res_34 <- (sfilm_coord$sdist_34 - sfilm_coord$fdist_34) ^ 2

# Remove any rows that exceedes the angle threhold between intersecting fiducial lines (remove if > 1 SD from calibration reported value):
sfilm_coord$angle <- 180 - intersect.angle(sfilm_coord)
angle_error <- sd(sfilm_coord$angle)
sfilm_coord <- filter(sfilm_coord,angle > fid_angle - angle_error & angle < fid_angle + angle_error)

# Minimize RSS and select the best row of sample data:
sfilm_coord$rss <- apply(select(sfilm_coord,res_12,res_13,res_14,res_23,res_24,res_34), 1, sum)
best_row <- which(sfilm_coord$rss == min(sfilm_coord$rss), arr.ind = TRUE)
# Subset the image_fids df by the best row:
best_data <- sfilm_coord[best_row,]

# Collect results:
derived_coord <- data.frame(x_revised = c(best_data$sfilm_1x,best_data$sfilm_2x,best_data$sfilm_3x,best_data$sfilm_4x),
           y_revised = c(best_data$sfilm_1y,best_data$sfilm_2y,best_data$sfilm_3y,best_data$sfilm_4y))
film_coord <- cbind(film_coord,derived_coord)

# Derive residuals:
film_coord$xres <- film_coord$x - film_coord$x_revised
film_coord$yres <- film_coord$y - film_coord$y_revised

# Caluclate RMS of final solution:
film_rms <- sqrt(mean(c(film_coord$xres^2,film_coord$yres^2)))
film_rms






