
# Projet Name: Frame Camera Calibration
# Script Name: 1_film_coordinates

# Stephen Bird, 2019-11-04

# This project takes a camera calibration report and a scanned airphoto(s) with
# fiducial marks and determines the image resolution, focal length, and principal
# point all in units of pixels.

# This script calculates the film coordinates of the fiducial marks based on the
# calibration report.


# Load libraries, enter data & set parameters ---------------------------------

# This is a Packrat project.

library(here)
library(dplyr)
source("scripts/pp_functions.R")

# Film coordinates of the fiducial marks (mm):
film_coord <-
        data.frame(
                fids = c("1", "2", "3", "4"),
                x = c(-105.993, 106.001,-106.015, 106.015),
                y = c(-106.006, 106.002, 105.982,-106.006)
        )

# Film distance between fiducial marks (mm):
fid_dist <-
        data.frame(
                fids = c("1_2", "1_3", "1_4", "2_3", "2_4", "3_4"),
                # For example, 1_2 is distance between fiducial 1 and 2
                dist_mm = c(299.814, 211.988, 212.015, 212.008, 212.007, 299.826)
        )

# Angle between lines intersecting the fiducial marks:
fid_angle <- 90 + 13 / 3600

# Enter the error in the film measurements as stated in the report. If none is
# given, estimate a reasonable value:
dist_err <- 0.003 # Stated error between two points
point_err <-
        sqrt((dist_err ^ 2) / 2) # Error of placing a single point

ipp_cal <- c(-0.003, -0.009) # Indicated principal point (i.e. intersetion of fiducials) (mm)
ppa <- c(0, 0) # Principal point of autocollimation (mm)
cpp_cal <- c(0.0, 0.007) # Calibrated principal point (mm)
focal_cal <- 153.228 # Focal length (mm)

# How many samples to generate?
n <- 10000


# Fit the film-based fiducial coordinates --------------------------------------

# Refine the measured fiducial coordinates and distances among fiducials as
# stated in the calibration report.

# Generate normal random samples of each x,y coordinate as given in the
# calibration report using the stated measurement error (e.g. sfilm_1x = sample
# film, fiducial 1, x) (the "s" prefix denotes "sample"):
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

# Make distance measurements among all fiducial points using the sample coordinates generated above. Use the "hypo" function (1):
sfilm_coord$sdist_12 <- hypo(sfilm_coord$sfilm_1x,sfilm_coord$sfilm_1y,sfilm_coord$sfilm_2x,sfilm_coord$sfilm_2y)
sfilm_coord$sdist_13 <- hypo(sfilm_coord$sfilm_1x,sfilm_coord$sfilm_1y,sfilm_coord$sfilm_3x,sfilm_coord$sfilm_3y)
sfilm_coord$sdist_14 <- hypo(sfilm_coord$sfilm_1x,sfilm_coord$sfilm_1y,sfilm_coord$sfilm_4x,sfilm_coord$sfilm_4y)
sfilm_coord$sdist_23 <- hypo(sfilm_coord$sfilm_2x,sfilm_coord$sfilm_2y,sfilm_coord$sfilm_3x,sfilm_coord$sfilm_3y)
sfilm_coord$sdist_24 <- hypo(sfilm_coord$sfilm_2x,sfilm_coord$sfilm_2y,sfilm_coord$sfilm_4x,sfilm_coord$sfilm_4y)
sfilm_coord$sdist_34 <- hypo(sfilm_coord$sfilm_3x,sfilm_coord$sfilm_3y,sfilm_coord$sfilm_4x,sfilm_coord$sfilm_4y)

# Generate normal random samples based on film measured distances and the calibrated measurement error:
set.seed(24)
sfilm_coord$fdist_12 <- rnorm(n, mean = fid_dist$dist_mm[1], sd = dist_err)
sfilm_coord$fdist_13 <- rnorm(n, mean = fid_dist$dist_mm[2], sd = dist_err)
sfilm_coord$fdist_14 <- rnorm(n, mean = fid_dist$dist_mm[3], sd = dist_err)
sfilm_coord$fdist_23 <- rnorm(n, mean = fid_dist$dist_mm[4], sd = dist_err)
sfilm_coord$fdist_24 <- rnorm(n, mean = fid_dist$dist_mm[5], sd = dist_err)
sfilm_coord$fdist_34 <- rnorm(n, mean = fid_dist$dist_mm[6], sd = dist_err)

# Calculate residuals^2 between distances measured on the film and distances sampled above:
sfilm_coord$res_12 <- (sfilm_coord$sdist_12 - sfilm_coord$fdist_12) ^ 2
sfilm_coord$res_13 <- (sfilm_coord$sdist_13 - sfilm_coord$fdist_13) ^ 2
sfilm_coord$res_14 <- (sfilm_coord$sdist_14 - sfilm_coord$fdist_14) ^ 2
sfilm_coord$res_23 <- (sfilm_coord$sdist_23 - sfilm_coord$fdist_23) ^ 2
sfilm_coord$res_24 <- (sfilm_coord$sdist_24 - sfilm_coord$fdist_24) ^ 2
sfilm_coord$res_34 <- (sfilm_coord$sdist_34 - sfilm_coord$fdist_34) ^ 2

# Remove any rows that exceeds the angle threshold between intersecting fiducial lines (remove if > 1 SD from calibration reported value):
sfilm_coord$angle <- 180 - intersect.angle(sfilm_coord)
angle_error <- sd(sfilm_coord$angle)
sfilm_coord <- filter(sfilm_coord,angle > fid_angle - angle_error & angle < fid_angle + angle_error)

# Minimize RSS and select the best row of sample data:
sfilm_coord$rss <-
        apply(select(sfilm_coord, res_12, res_13, res_14, res_23, res_24, res_34),
              1,
              sum)
best_row <-
        which(sfilm_coord$rss == min(sfilm_coord$rss), arr.ind = TRUE)
# Subset the image_fids df by the best row:
best_data <- sfilm_coord[best_row,]

# Collect results:
derived_coord <-
        data.frame(
                x_revised = c(
                        best_data$sfilm_1x,
                        best_data$sfilm_2x,
                        best_data$sfilm_3x,
                        best_data$sfilm_4x
                ),
                y_revised = c(
                        best_data$sfilm_1y,
                        best_data$sfilm_2y,
                        best_data$sfilm_3y,
                        best_data$sfilm_4y
                )
        )
film_coord <- cbind(film_coord, derived_coord)

# Derive residuals:
film_coord$xres <- film_coord$x - film_coord$x_revised
film_coord$yres <- film_coord$y - film_coord$y_revised

# Calculate RMS of final solution:
film_rms <- sqrt(mean(c(film_coord$xres^2,film_coord$yres^2)))
# RMS in mm derived among the x,y fiducial coordinates stated in the calibration
# report and those adjusted here based on minimizing the RSS (units of mm):
print(film_rms) 


# Generate n samples -----------------------------------------------------------

# Generate a random sample of film fiducial coordinates based on the measurement
# error and the refined coordinates taken from the calibration report:
fid_film <- select(film_coord, x_revised, y_revised)
names(fid_film) <- c("x", "y") # Column names much be in this form
fid_film$x_err <- rep(point_err, 4)
fid_film$y_err <- rep(point_err, 4)
# Generate a sample of coordinates for each fiducial mark based on the
# measurement error. Use the "fiducial.samples" function (2):
film_fids <- fiducial.samples(fid_film, n)
# Calculate the angle between crossing lines that connect fiducials. Use the "intersect.angle" function (3):
film_fids$m <- 180 - intersect.angle(film_fids[, 1:8])
intersect_mean <- mean(film_fids$m)
# Error in the two intersecting fiducial lines derived from the random sample:
intersect_err <- sd(film_fids$m)
ifelse(
        abs(intersect_mean - fid_angle) <= intersect_err,
        "Good result -- the measured angle is within one standard deviation of the mean angle.",
        "Poor result -- the measured angle is greater than one standard deviation of the mean angle."
)

# A data frame of plausible coordinates as stated in the calibration report (in units of mm):
head(film_fids)

