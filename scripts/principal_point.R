
# This is a Packrat project.

library(here)
library(dplyr)
source("scripts/pp_functions.R")


# Enter data from the camera calibration report and a given image --------------

# Fid = fiducial
# Numbers refer to the numbering system in the report (fiducial 1, 2, 3, and 4...usually)
# Sub_x = x coordinate, sub_y = y coordinate
# Px = units of pixels, mm = units of millimeters


image_no <- c("bcc04050_041")
cal_focal <- 153.228 # Calibrated focal length in mm

# Enter the number of pixels to the fiducial mark and the film coordinates for
# the same mark in mm. Use a pixel/image coordinate system for the pixels (0,0
# at top left). Note that the pixel coordinate will later be refined...at this
# point, just make sure the coordinate is +/- 1 pixel from the true position.
# Note also that the origin of the image and film coordinate systems do not need
# to match.
meas_coord <-
        data.frame(
                # Fiducial number
                fid = c(1, 2, 3, 4),
                # x-coordinate in pixels
                px_x = c(15943, 801, 15947, 798),
                # y-coordinate in pixels
                px_y = c(15985, 836, 839, 15978),
                # x-coordinate in mm
                mm_x = c(-105.993, 106.001,-106.015, 106.015),
                # y-coordinate in mm
                mm_y = c(-106.006, 106.002, 105.982,-106.006) 
        )

image_dimensions <- c(16730,16816) # Width and hieght in units of pixels
ipp_cal <- c(-.003,-.009) # Indicated principal point in mm
ppa_cal <- c(0,0) # Principal point of autocollimination
cpp_cal <- c(0.000,0.007) # Calibrated principal point

# Enter the distance between fiducails in mm:
dist2fid_mm <-
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

# Enter the DN's for each pixel in a 9x9 window centered at the pixel at the
# percieved center of the fiducial. Start in the top right position in the
# matrix (assuming we're viewing the image from the focal plane). If the image
# is in RGB, do all three channels before moving to the next cell. Filter the
# fiducial by Outliers, Despekle, and Gaussian Blur using the Fiji tool.
fid1_dn <- matrix(c(106,74,65,111,77,67,101,68,59, # Row 1, rgb,rgb,...
                    108,77,66,113,80,69,102,71,60, # Row 2
                    95,67,57,99,70,60,90,62,52) # Row 3
                   ,nrow = 3)

fid2_dn <- matrix(c(127,84,74,141,97,85,139,95,81, # Row 1
                    138,94,83,154,107,95,151,104,91, # Row 2
                    133,88,78,148,101,89,146,97,86) # Row 3
                  ,nrow = 3)

fid3_dn <- matrix(c(188,123,102,188,125,102,174,110,88, # Row 1
                    194,134,112,195,137,113,180,120,98, # Row 2
                    187,128,107,187,130,108,174,114,93) # Row 3
                   ,nrow = 3)

fid4_dn <- matrix(c(181,123,104,185,129,108,174,118,96, # Row 1
                    193,135,114,197,142,118,185,129,105, # Row 2
                    190,129,108,193,135,112,181,123,99) # Row 33
                  ,nrow = 3)


# Find the coordinates of the fiducial points on the film ----------------------

fid_film <- select(meas_coord, mm_x, mm_y)
names(fid_film) <- c("x", "y") # Column names much be in this form
fid_film$x_err <- rep(point_err, 4)
fid_film$y_err <- rep(point_err, 4)
film_fids <- fiducial.samples(fid_film, 1000000)
film_fids$m <- 180 - intersect.angle(film_fids[, 1:8])
# Intersecting angle between fiducial lines derived from the random sample"
intersect_mean <- mean(film_fids$m)
# Error in the two intersecting fiducial lines derived from the random sample:
intersect_err <- sd(film_fids$m)
ifelse(
        abs(intersect_mean - fid_angle) <= intersect_err,
        "Good result -- the measured angle is within one standard deviation of the mean angle.",
        "Poor result -- the measured angle is greater than one standard deviation of the mean angle."
)


# Find the coordinates of the fiducial points on the image ---------------------

# Important: the center row and center column must both be 2. If not, re-center
# the meas_coord matrix.

# Fiducial 1:
fid1_sub <-fid.coord(fid1_dn,meas_coord,1)

# Fiducial 2:
fid2_sub <- fid.coord(fid2_dn,meas_coord,2)

# Fiducial 3:
fid3_sub <- fid.coord(fid3_dn,meas_coord,3)

# Fiducial 4:
fid4_sub <- fid.coord(fid4_dn,meas_coord,4)

# These are the image coordinates of the fiducial marks. They are registered at
# the sub-pixel level and are within +/- 0.5 pixels to their true location.
fid_sub <- rbind(fid1_sub,fid2_sub,fid3_sub,fid4_sub)
fid_sub <- as.data.frame(fid_sub)
print(fid_sub)

# Generate normal random samples of the sub-pixel fiducial location based on the
# sub-pixel estimate and the uncertainty:
image_fids <- fiducial.samples(fid_sub,1000000)
image_fids$m <- 180 - intersect.angle(image_fids[,1:8])


# Calculate the scan resolution of the image -----------------------------------

# Compute the pixel resolution based on the sample of 10,000 fiducial pixel
# distances and lab distances:
image_fids$res_12 <- film_fids$dist_12 / image_fids$dist_12
image_fids$res_13 <- film_fids$dist_13 / image_fids$dist_13
image_fids$res_14 <- film_fids$dist_14 / image_fids$dist_14
image_fids$res_23 <- film_fids$dist_23 / image_fids$dist_23
image_fids$res_24 <- film_fids$dist_24 / image_fids$dist_24
image_fids$res_34 <- film_fids$dist_34 / image_fids$dist_34

# Calculate the mean resolution for each sample row and correspoinding standard
# deviation. Select the smallest standard deviation as the best fit (i.e. the
# distnaces bewtween all fiducials is the ~same resolution):
image_fids$mean_res <- apply(select(image_fids,res_12,res_13,res_14,res_23,res_24,res_34), 1, mean)
image_fids$sd_res <- apply(select(image_fids,res_12,res_13,res_14,res_23,res_24,res_34), 1, sd)
# Calculate the weighted average based on nucertainty:
image_fids$w <- 1 / image_fids$sd_res
image_fids$mean_w <- image_fids$mean_res * image_fids$w
mean_w <- apply(select(image_fids,mean_w), 2, sum)
w <- apply(select(image_fids,w), 2, sum)
mean_res <- mean_w / w # This is the best estimate of the pixel resolution for this image
sd_res <- 1 / sqrt(w) # This is the uncertainty in the best estimate of the pixel resolution for this image


# Fit the coordinates by minimizing the sum of squares -------------------------

# Translate the film distance between fiducials from mm to pixels based on the
# mean image resolution:
film_fids$pxdist_12 <- film_fids$dist_12 / mean_res
film_fids$pxdist_13 <- film_fids$dist_13 / mean_res
film_fids$pxdist_14 <- film_fids$dist_14 / mean_res
film_fids$pxdist_23 <- film_fids$dist_23 / mean_res
film_fids$pxdist_24 <- film_fids$dist_24 / mean_res
film_fids$pxdist_34 <- film_fids$dist_34 / mean_res

# Calculate the residuals between the image-pixel distances and the film-pixel
# distances between fiducials. Residuals are Internally Studentized (divide by
# sd) so they can be scaled and compared to one another (allowing the intersect
# angle to weigh equally in the RSS calc):
sqres_12 <- resid.calc(image_fids$dist_12,film_fids$pxdist_12)
sqres_13 <- resid.calc(image_fids$dist_13,film_fids$pxdist_13)
sqres_14 <- resid.calc(image_fids$dist_14,film_fids$pxdist_14)
sqres_23 <- resid.calc(image_fids$dist_23,film_fids$pxdist_23)
sqres_24 <- resid.calc(image_fids$dist_24,film_fids$pxdist_24)
sqres_34 <- resid.calc(image_fids$dist_34,film_fids$pxdist_34)
sqres_m <- resid.calc(image_fids$m,film_fids$m)

# Gather results into a df
resids <- data.frame(sqres_12,sqres_13,sqres_14,sqres_23,sqres_24,sqres_34,sqres_m)
# Calculate the RSS:
resids$rss <- apply(resids, 1, sum)
# Select the row with the smallest RSS:
best_row <- which(resids$rss == min(resids$rss), arr.ind = TRUE)
# Subset the image_fids df by the best row:
best_data <- image_fids[best_row,]

# Calculate the RMS of the fit among measurements between fiducials. The results should be at the sub-pixel level:
best_image_fids <- select(image_fids[best_row,],dist_12,dist_13,dist_14,dist_23,dist_24,dist_34)
best_film_fids <- select(film_fids[best_row,],pxdist_12,pxdist_13,pxdist_14,pxdist_23,pxdist_24,pxdist_34)
rms_fids <- rms.calc(best_image_fids,best_film_fids)
print(rms_fids) # RMS among fiducials in units of pixels

# Translate the pixel coordinates on the image to minimize the residuals between
# the original sub-pixel estimates and the Monte Carlo sample estimate found in
# the "best_data" df:
fid_sub$best_x <- t(best_data[1:4])
fid_fitx <- as.matrix(fid_sub[,c(1,7)],ncol=2)
colnames(fid_fitx) <- c("measured","derived") # Rename columns
fid_fitx <- as.data.frame(fid_fitx)
x_trans <- optimize(min.x,c(-1,1),data = fid_fitx)
fid_fitx$trans <- fid_fitx$derived - x_trans$minimum

fid_sub$best_y <- t(best_data[5:8])
fid_fity <- as.matrix(fid_sub[,c(2,8)],ncol=2)
colnames(fid_fity) <- c("measured","derived") # Rename columns
fid_fity <- as.data.frame(fid_fity)
y_trans <- optimize(min.x,c(-1,1),data = fid_fity)
fid_fity$trans <- fid_fity$derived - y_trans$minimum
# Final calibrated fiducial mark coordinates (units of pixels):
fid_cal <- data.frame(x = fid_fitx$trans,y = fid_fity$trans) 
rms_imagefids <- rms.calc(fid_cal,fid_sub[,1:2])
print(rms_imagefids)


# Find the indicated principal point (i.e. where the fiducials intersect) ------

# Slope and intercept of line 1 - 2
rise1 <- fid_cal$y[1] - fid_cal$y[2]
run1 <- fid_cal$x[1] - fid_cal$x[2]
slope1 <- rise1/run1
inter1 <- fid_cal$y[1] - slope1 * fid_cal$x[1]

# Slope and intercept of line 3 - 4
rise2 <- fid_cal$y[3] - fid_cal$y[4]
run2 <- fid_cal$x[3] - fid_cal$x[4]
slope2 <- rise2/run2
inter2 <- fid_cal$y[3] - slope2 * fid_cal$x[3]

# Intersection point or "indicated principal point":

icent <- image_dimensions/2 # Find the center pixel of the image
ipp <- line.sect(slope1,slope2,inter1,inter2) # In an image coordinate system with (0,0) at top left
icent2ipp <- icent - ipp # This is the location of the IPP relative to the image center
icent2ipp[1] <- -1 * icent2ipp[1] # Change the sign of the x-coord since the origin is now in the image center with the y-axis point up
ipp2ppa <- icent2ipp - ipp_cal / mean_res # PPA relative to the image center
ppa2cpp <- ipp2ppa + cpp_cal / mean_res # Calibrated principal point
image_1 <- list(image_no,ppa2cpp,mean_res,sd_res,cal_focal) # Gather the results
print(image_1)











