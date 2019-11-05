
# Projet Name: Frame Camera Calibration
# Script Name: 3_image_fid_fit

# Stephen Bird, 2019-11-04

# This script determines the pixel coordinates of the principal point for each
# image. Replicate for each image in the block.


# Fit the coordinates by minimizing the sum of squares -------------------------

# Translate the film distance between fiducials from mm to pixels based on the
# mean image resolution:
film_fids$pxdist_12 <- film_fids$dist_12 / mission_res
film_fids$pxdist_13 <- film_fids$dist_13 / mission_res
film_fids$pxdist_14 <- film_fids$dist_14 / mission_res
film_fids$pxdist_23 <- film_fids$dist_23 / mission_res
film_fids$pxdist_24 <- film_fids$dist_24 / mission_res
film_fids$pxdist_34 <- film_fids$dist_34 / mission_res


# Image bcc04050_043 -----------------------------------------------------------

# Calculate the residuals between the image-pixel distances and the film-pixel
# distances between from a sample of n distance measurements. Select the
# combination of distances with the smallest RMSE. Then, georeference the film
# coordinates to the image coordinates, and translate the first fiducial to
# match coordinates (film vs. image). Use the "fit.fids" function (16):
fitted_fids <- fit.fids(bcc04050_043,n)
print(fitted_fids)
best_image_fids <- fitted_fids[[2]] # Best set of image fids that match the film fids
film_coord <- fitted_fids[[1]] # Translation based on fid 1
rms_fids <- fitted_fids[[3]] # RMS for distance measurements among fiducials (pixels)

# Once the fids are translated, find the rotation angle that minimizes the RMSE
# among the corresponding film and image measurements (rotate the coordinates
# about Fid 1):
optimized_rotation <- optimize(rot, c(-2 * pi, 2 * pi), best_image_fids)
rms.xy <- optimized_rotation$objective  #report this value for planimetric error (units of pixels...should be < 1)
print(rms.xy) # Fit among pixels
x <- optimized_rotation$minimum # Extract the target angle
# Rotate the data:
film_coord$rotx <- film_coord$x_trans[1] + (film_coord$x_trans-film_coord$x_trans[1]) * cos(x) - (film_coord$y_trans-film_coord$y_trans[1]) * sin(x) #rotate easting by x radians
film_coord$roty <- film_coord$y_trans[1] + (film_coord$x_trans-film_coord$x_trans[1]) * sin(x) + (film_coord$y_trans-film_coord$y_trans[1]) * cos(x) #rotate northing by x radians

# Optimize the fit of the rotated x fid coordinate by minimizing the difference between the best image coordinates and the rotated film coordinates:
image_axis <- best_image_fids$x
film_axis <- film_coord$rotx
adjustment <- optimize(fid.rms,c(-1,1))
film_coord$optim_x <- adjustment$minimum + film_coord$rotx

# Optimize the fit of the rotated y fid coordinate by minimizing the difference between the best image coordinates and the rotated film coordinates:
image_axis <- best_image_fids$y
film_axis <- film_coord$roty
adjustment <- optimize(fid.rms,c(-1,1))
film_coord$optim_y <- adjustment$minimum + film_coord$roty

# Calculate error in the fit among image pixel centers and those predicted by the calibration measurements (should be sub-pixel:
x_rms <- rms.calc(best_image_fids$x,film_coord$optim_x)
y_rms <- rms.calc(best_image_fids$y,film_coord$optim_y)
rms_imagefids <- data.frame(x_rms,y_rms)


# Find the indicated principal point (i.e. where the fiducials intersect) ------

fid_cal <- data.frame(x = film_coord$optim_x,y = film_coord$optim_y)
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
image_dimensions_sub <- image_dimensions[1,2:3]

icent <- image_dimensions_sub/2 # Find the center pixel of the image
ipp <- line.sect(slope1,slope2,inter1,inter2) # In an image coordinate system with (0,0) at top left
icent2ipp <- icent - ipp # This is the location of the IPP relative to the image center
icent2ipp[1] <- -1 * icent2ipp[1] # Change the sign of the x-coord since the origin is now in the image center with the y-axis point up
ipp2ppa <- icent2ipp - ipp_cal / mission_res # PPA relative to the image center
ppa2cpp <- ipp2ppa + cpp_cal / mission_res # Calibrated principal point
focal_cal_px <- focal_cal / mission_res
final_results <- list(image_no[1],fid_cal,ppa2cpp,mission_res,focal_cal_px,rms_fids,rms_imagefids) # Gather the results
names(final_results) <- c("Images","Image coordinates of the fiducial marks (pixels)","Calibrated principal point","Image resolution (mm)","focal length (pixels)","RMS of the fids","Fiducial placement RMS")
print(final_results)

bcc04050_043 <- append(bcc04050_043,final_results)


# Image bcc04050_042 -----------------------------------------------------------

# Calculate the residuals between the image-pixel distances and the film-pixel
# distances between from a sample of n distance measurements. Select the
# combination of distances with the smallest RMSE. Then, georeference the film
# coordinates to the image coordinates, and translate the first fiducial to
# match coordinates (film vs. image). Use the "fit.fids" function (16):
fitted_fids <- fit.fids(bcc04050_042,n)
print(fitted_fids)
best_image_fids <- fitted_fids[[2]] # Best set of image fids that match the film fids
film_coord <- fitted_fids[[1]] # Translation based on fid 1
rms_fids <- fitted_fids[[3]] # RMS for distance measurements among fiducials (pixels)

# Once the fids are translated, find the rotation angle that minimizes the RMSE
# among the corresponding film and image measurements (rotate the coordinates
# about Fid 1):
optimized_rotation <- optimize(rot, c(-2 * pi, 2 * pi), best_image_fids)
rms.xy <- optimized_rotation$objective  #report this value for planimetric error (units of pixels...should be < 1)
print(rms.xy) # Fit among pixels
x <- optimized_rotation$minimum # Extract the target angle
# Rotate the data:
film_coord$rotx <- film_coord$x_trans[1] + (film_coord$x_trans-film_coord$x_trans[1]) * cos(x) - (film_coord$y_trans-film_coord$y_trans[1]) * sin(x) #rotate easting by x radians
film_coord$roty <- film_coord$y_trans[1] + (film_coord$x_trans-film_coord$x_trans[1]) * sin(x) + (film_coord$y_trans-film_coord$y_trans[1]) * cos(x) #rotate northing by x radians

# Optimize the fit of the rotated x fid coordinate by minimizing the difference between the best image coordinates and the rotated film coordinates:
image_axis <- best_image_fids$x
film_axis <- film_coord$rotx
adjustment <- optimize(fid.rms,c(-1,1))
film_coord$optim_x <- adjustment$minimum + film_coord$rotx

# Optimize the fit of the rotated y fid coordinate by minimizing the difference between the best image coordinates and the rotated film coordinates:
image_axis <- best_image_fids$y
film_axis <- film_coord$roty
adjustment <- optimize(fid.rms,c(-1,1))
film_coord$optim_y <- adjustment$minimum + film_coord$roty

# Calculate error in the fit among image pixel centers and those predicted by the calibration measurements (should be sub-pixel:
x_rms <- rms.calc(best_image_fids$x,film_coord$optim_x)
y_rms <- rms.calc(best_image_fids$y,film_coord$optim_y)
rms_imagefids <- data.frame(x_rms,y_rms)


# Find the indicated principal point (i.e. where the fiducials intersect) ------

fid_cal <- data.frame(x = film_coord$optim_x,y = film_coord$optim_y)
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
image_dimensions_sub <- image_dimensions[1,2:3]

icent <- image_dimensions_sub/2 # Find the center pixel of the image
ipp <- line.sect(slope1,slope2,inter1,inter2) # In an image coordinate system with (0,0) at top left
icent2ipp <- icent - ipp # This is the location of the IPP relative to the image center
icent2ipp[1] <- -1 * icent2ipp[1] # Change the sign of the x-coord since the origin is now in the image center with the y-axis point up
ipp2ppa <- icent2ipp - ipp_cal / mission_res # PPA relative to the image center
ppa2cpp <- ipp2ppa + cpp_cal / mission_res # Calibrated principal point
focal_cal_px <- focal_cal / mission_res
final_results <- list(image_no[2],fid_cal,ppa2cpp,mission_res,focal_cal_px,rms_fids,rms_imagefids) # Gather the results
names(final_results) <- c("Images","Image coordinates of the fiducial marks (pixels)","Calibrated principal point","Image resolution (mm)","focal length (pixels)","RMS of the fids","Fiducial placement RMS")
print(final_results)

bcc04050_042 <- append(bcc04050_042,final_results)


# Image bcc04050_041 -----------------------------------------------------------

# Calculate the residuals between the image-pixel distances and the film-pixel
# distances between from a sample of n distance measurements. Select the
# combination of distances with the smallest RMSE. Then, georeference the film
# coordinates to the image coordinates, and translate the first fiducial to
# match coordinates (film vs. image). Use the "fit.fids" function (16):
fitted_fids <- fit.fids(bcc04050_041,n)
print(fitted_fids)
best_image_fids <- fitted_fids[[2]] # Best set of image fids that match the film fids
film_coord <- fitted_fids[[1]] # Translation based on fid 1
rms_fids <- fitted_fids[[3]] # RMS for distance measurements among fiducials (pixels)

# Once the fids are translated, find the rotation angle that minimizes the RMSE
# among the corresponding film and image measurements (rotate the coordinates
# about Fid 1):
optimized_rotation <- optimize(rot, c(-2 * pi, 2 * pi), best_image_fids)
rms.xy <- optimized_rotation$objective  #report this value for planimetric error (units of pixels...should be < 1)
print(rms.xy) # Fit among pixels
x <- optimized_rotation$minimum # Extract the target angle
# Rotate the data:
film_coord$rotx <- film_coord$x_trans[1] + (film_coord$x_trans-film_coord$x_trans[1]) * cos(x) - (film_coord$y_trans-film_coord$y_trans[1]) * sin(x) #rotate easting by x radians
film_coord$roty <- film_coord$y_trans[1] + (film_coord$x_trans-film_coord$x_trans[1]) * sin(x) + (film_coord$y_trans-film_coord$y_trans[1]) * cos(x) #rotate northing by x radians

# Optimize the fit of the rotated x fid coordinate by minimizing the difference between the best image coordinates and the rotated film coordinates:
image_axis <- best_image_fids$x
film_axis <- film_coord$rotx
adjustment <- optimize(fid.rms,c(-1,1))
film_coord$optim_x <- adjustment$minimum + film_coord$rotx

# Optimize the fit of the rotated y fid coordinate by minimizing the difference between the best image coordinates and the rotated film coordinates:
image_axis <- best_image_fids$y
film_axis <- film_coord$roty
adjustment <- optimize(fid.rms,c(-1,1))
film_coord$optim_y <- adjustment$minimum + film_coord$roty

# Calculate error in the fit among image pixel centers and those predicted by the calibration measurements (should be sub-pixel:
x_rms <- rms.calc(best_image_fids$x,film_coord$optim_x)
y_rms <- rms.calc(best_image_fids$y,film_coord$optim_y)
rms_imagefids <- data.frame(x_rms,y_rms)


# Find the indicated principal point (i.e. where the fiducials intersect) ------

fid_cal <- data.frame(x = film_coord$optim_x,y = film_coord$optim_y)
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
image_dimensions_sub <- image_dimensions[1,2:3]

icent <- image_dimensions_sub/2 # Find the center pixel of the image
ipp <- line.sect(slope1,slope2,inter1,inter2) # In an image coordinate system with (0,0) at top left
icent2ipp <- icent - ipp # This is the location of the IPP relative to the image center
icent2ipp[1] <- -1 * icent2ipp[1] # Change the sign of the x-coord since the origin is now in the image center with the y-axis point up
ipp2ppa <- icent2ipp - ipp_cal / mission_res # PPA relative to the image center
ppa2cpp <- ipp2ppa + cpp_cal / mission_res # Calibrated principal point
focal_cal_px <- focal_cal / mission_res
final_results <- list(image_no[3],fid_cal,ppa2cpp,mission_res,focal_cal_px,rms_fids,rms_imagefids) # Gather the results
names(final_results) <- c("Images","Image coordinates of the fiducial marks (pixels)","Calibrated principal point","Image resolution (mm)","focal length (pixels)","RMS of the fids","Fiducial placement RMS")
print(final_results)

bcc04050_041 <- append(bcc04050_041,final_results)



