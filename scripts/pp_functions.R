
# Load functions ---------------------------------------------------------------


# 1. Function to calculate the hypotonuse based on xy coordinates of two points:
hypo <- function(x1, y1, x2, y2) {
        sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2)
}


# 2 Function to generate a sample of coordinates for each fiducial mark based on the measurement error:
fiducial.samples <- function(fid_coord,n) {
        # Generate samples based on the fiducial mark uncertainty. Generate normal
        # random samples (n = 10,000) for each sub-pixel estimatation of the image
        # fiducial center by x and y coordinate.
        set.seed(42)
        x1 <- rnorm(n, mean = fid_coord$x[1], sd = fid_coord$x_err[1])
        x2 <- rnorm(n, mean = fid_coord$x[2], sd = fid_coord$x_err[2])
        x3 <- rnorm(n, mean = fid_coord$x[3], sd = fid_coord$x_err[3])
        x4 <- rnorm(n, mean = fid_coord$x[4], sd = fid_coord$x_err[4])
        y1 <- rnorm(n, mean = fid_coord$y[1], sd = fid_coord$y_err[1])
        y2 <- rnorm(n, mean = fid_coord$y[2], sd = fid_coord$y_err[2])
        y3 <- rnorm(n, mean = fid_coord$y[3], sd = fid_coord$y_err[3])
        y4 <- rnorm(n, mean = fid_coord$y[4], sd = fid_coord$y_err[4])
        
        # Calcualte the distance between fiducials based on the 10,000 random sample of
        # x,y coordinates:
        dist_12 <- hypo(x1,y1,x2,y2)
        dist_13 <- hypo(x1,y1,x3,y3)
        dist_14 <- hypo(x1,y1,x4,y4)
        dist_23 <- hypo(x2,y2,x3,y3)
        dist_24 <- hypo(x2,y2,x4,y4)
        dist_34 <- hypo(x3,y3,x4,y4)
        
        # Gather results into df:
        fids <- data.frame(x1, x2, x3, x4, y1, y2, y3, y4, dist_12,dist_13,dist_14,dist_23,dist_24,dist_34)
        fids
        
}


# 3. Function to calcualte the angle between crossing lines that connect fiducials based on the n random sample of
# x,y coordinates:
intersect.angle <- function(fid_pts){
        
        # Input includes 8 columns of data: columns 1:4 are the x coordinates of the fiducials, columns 5:8 are the corresponding y coordinates.
        slope_line12 <- ((fid_pts[,5] - fid_pts[,6]) / (fid_pts[,1] - fid_pts[,2])) # Slope of line 1_2
        slope_line34 <- ((fid_pts[,7] - fid_pts[,8]) / (fid_pts[,3] - fid_pts[,4])) # Slope of line 3_4
        intersect_ang <- ((180 / pi) * atan((slope_line12 - slope_line34) / (1 + slope_line12 * slope_line34))) # Intersection angle line 12 vs line 34
        intersect_ang
        
}



# 4. Function calculates image resolution:
mission.res <- function(fid1_pxy,fid2_pxy,fid3_pxy,fid4_pxy,film_coord,film_fids,n){
        
        # fid1_pxy = c(x_coord, y_coord)...
        # film_coord = measured & revised pixel coordinates
        # film_fids = sampled distances
        # n = sample size
        
        # Create a df with the fiducial ID, fiducial coordinates in pixels, and fiducial coordinates in mm:
        meas_coord <- data.frame(fid = c(1,2,3,4), 
                                 px_x = c(fid1_pxy[1], fid2_pxy[1], fid3_pxy[1], fid4_pxy[1]),
                                 px_y = c(fid1_pxy[2], fid2_pxy[2], fid3_pxy[2], fid4_pxy[2]),
                                 mm_x = c(film_coord$x_revised[1], film_coord$x_revised[2], film_coord$x_revised[3], film_coord$x_revised[4]),
                                 mm_y = c(film_coord$y_revised[1], film_coord$y_revised[2], film_coord$y_revised[3], film_coord$y_revised[4])
        )
        
        # Make an initial estimate of the sub-pixel coordinates of each pixel:
        fid1_sub <-fid.coord(fid1_dn,meas_coord,1)
        fid2_sub <-fid.coord(fid2_dn,meas_coord,2)
        fid3_sub <-fid.coord(fid3_dn,meas_coord,3)
        fid4_sub <-fid.coord(fid4_dn,meas_coord,4)
        
        # These are the image coordinates of the fiducial marks. They are registered at
        # the sub-pixel level and are within +/- 0.5 pixels to their true location.
        fid_sub <- rbind(fid1_sub,fid2_sub,fid3_sub,fid4_sub)
        fid_sub <- as.data.frame(fid_sub)
        writeLines("\n") # Add a carrage return for readability of the ouput
        print(fid_sub)
        
        # Generate normal random samples of the sub-pixel fiducial location based on the
        # sub-pixel estimate and the uncertainty:
        image_fids <- fiducial.samples(fid_sub,n)
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
        writeLines("\n") # Add a carrage return for readability of the ouput
        print(mean_res)
        
        list(fid_sub,mean_res,sd_res)
        
}


# 5. Function to calculate the weigthed mean:
weighted.mean <- function(mean_res,sd_res){
        
        w <- 1/sd_res^2
        R <- sum(mean_res * w) / sum(w)
        R
        
}


# 6. Function to calculate the residuals between the image-pixel distances and
# the film-pixel distances between fiducials. Residuals are Internally
# Studentized (divide by sd) so they can be scaled and compared to one another:
resid.calc <- function(image_dist,film_pix_dist){
        
        resid <- image_dist - film_pix_dist
        resid_sd <- sd(resid)
        resid_st <- resid/resid_sd
        resid_st2 <- resid_st^2
        resid_st2
}


# 7. Fucntion to determine if a 9x9 matrix contains the maximum pixel value in the center of the matrix:
subfid.mx <- function(channel,fid_dn){
        # This function takes a matrix of pixel DN's for each channel and returns a
        # matrix for a given channel (1 = red, 2 = blue, 3 = green).
        chan_dn <- fid_dn[channel,] # Subset the matrix by channel
        chan_dn <- matrix(chan_dn,ncol=3) # Rebuild the original (but transposed) matrix by row and column
        t(chan_dn) # Transpose to the original matrix by channel
        
}


# 8. Function to calculate the sub-pixel center of a line made-up of three
# points. Fucntion assumes the center coordinate is the middle: pa = dn at
# pixel 1; pb = dn at pixel 2; c = dn at pixel 3:
sub.px <- function(pa,pb,pc) {
        
        0.5*(pa - pc) / (pa-2*pb+pc)
        
}


# 9. Function to find the coordinates of two intersecting lines:
line.sect <- function(slope1,slope2,inter1,inter2){
        slope_diff <- slope1 - slope2
        intercept_sum <- inter1 + inter2
        solved_x <- intercept_sum/slope_diff
        solved_y <- slope1 * solved_x + inter1
        point <- c(solved_x,solved_y)
        names(point) <- c("x","y")
        point
}



# 10. Function to find the coordinates of a fiducial point. Inputs: fid_dn = 9x3 matrix of RGB data:
fid.coord <- function(fid_dn,meas_coord,fid_no){
        
        # Sum the DN's for each channel (e.g. red + green + blue) and save as a matrix:
        fid_mx <- subfid.mx(1,fid_dn) + subfid.mx(2,fid_dn) + subfid.mx(3,fid_dn)
        # Find the row and column with the maximum DN:
        max_coord <- which(fid_mx == max(fid_mx), arr.ind = TRUE)
        names(max_coord) <- c("center row","center col")
        # Subset the max row and column:
        max_row <- fid_mx[max_coord[1],]
        max_col <- fid_mx[,max_coord[2]]
        center <- c(max_row,max_col)
        
        # Find the sub-pixel center for each row and column, relative to the "max_coord" matrix cell:
        row_fid <- sub.px(max_row[1],max_row[2],max_row[3])
        col_fid <- sub.px(max_col[1],max_col[2],max_col[3])
        
        # Need to give the relative position a real pixel coordinate: 
        subpix_x <- ifelse(max_coord[1] == 1, meas_coord[fid_no,2] - 1 + row_fid,
                           ifelse(max_coord[1] == 2, meas_coord[fid_no,2] + row_fid,
                                  meas_coord[fid_no,2] + 1 + row_fid))
        
        subpix_y <- ifelse(max_coord[2] == 1, meas_coord[fid_no,3] - 1 + col_fid,
                           ifelse(max_coord[2] == 2, meas_coord[fid_no,3] + col_fid,
                                  meas_coord[fid_no,3] + 1 + col_fid))
        coord <- c(subpix_x,subpix_y,0.5,0.5)
        names(coord) <- c("x","y","x_err","y_err")
        coord_summary <- c(coord,max_coord)
        coord_summary
        
}


# 11. Function to calucalte the RMSE between the image and film fiducuals:
rms.calc <- function(best_image_fids,best_film_fids){
        
        resid <- best_image_fids - best_film_fids
        resid2 <- resid^2
        mean_resid2 <- mean(resid2)
        sqrt(mean_resid2)
        
}


# 12. Function to minimize the shift between measured sub-pixel centers and those dervied from the distance between fiducials using a simple translation:
min.x <- function(data, par) {
        
        with(data, abs(sum(
                
                (derived - par) - measured)))
}


# 13. Function to calcualte the angle between crossing lines that connect fiducials:
single.intersect.angle <- function(x1,y1,x2,y2,x3,y3,x4,y4){
        
        # Input includes 8 columns of data: columns 1:4 are the x coordinates of the fiducials, columns 5:8 are the corresponding y coordinates.
        slope_line12 <- ((y1 - y2) / (x1 - x2)) # Slope of line 1_2
        slope_line34 <- ((y3 - y4) / (x3 - x4))  # Slope of line 3_4
        intersect_ang <- ((180 / pi) * atan((slope_line12 - slope_line34) / (1 + slope_line12 * slope_line34))) # Intersection angle line 12 vs line 34
        intersect_ang
        
}


# 14. Fucntion to rotate coordinates and  minimize the RMS error (field vs control data)

rot <- function(x,best_image_fids) {
        rotated_east <-
                film_coord$x_trans[1] + (film_coord$x_trans - film_coord$x_trans[1]) * cos(x) - (film_coord$y_trans -
                                                                                                         film_coord$x_trans[1]) * sin(x) #rotate easting by x radians
        rotated_north <-
                film_coord$y_trans[1] + (film_coord$x_trans - film_coord$x_trans[1]) * sin(x) + (film_coord$y_trans -
                                                                                                         film_coord$y_trans[1]) * cos(x) #rotate northing by x radians
        
        rms_xy <-
                (mean((rotated_east - best_image_fids$x) ^ 2) +  mean((rotated_north - best_image_fids$y) ^
                                                                              2)) ^ 0.5 #RMS error
}


# 15. Function to calcualte the RMS error between the 

fid.rms <- function(adj){
        # Calculate the residuals:
        residuals <- image_axis - film_axis
        # Add the optimization adjustment:
        adj_residuals <- adj + residuals
        # Square the residuals
        sq_residuals <- adj_residuals^2
        # Average the residuals and take the sqrt:
        sqrt(mean(sq_residuals))
}


# 16. Function to calculate the image coordinates of the fiducials

fit.fids <- function(image,n){
        
        image_fids <- fiducial.samples(image[[1]],n)
        
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
        #sqres_m <- resid.calc(image_fids$m,film_fids$m)
        
        # Gather results into a df
        resids <- data.frame(sqres_12,sqres_13,sqres_14,sqres_23,sqres_24,sqres_34)
        # Calculate the RSS:
        resids$rss <- apply(resids, 1, sum)
        # Select the row with the smallest RSS:
        best_row <- which(resids$rss == min(resids$rss), arr.ind = TRUE)
        # Subset the image_fids df by the best row:
        best_data <- image_fids[best_row,]
        
        # Calculate the RMS of the fit among measurements between fiducials. The results should be at the sub-pixel level:
        best_image_dist <- select(image_fids[best_row,],dist_12,dist_13,dist_14,dist_23,dist_24,dist_34)
        best_film_fids <- select(film_fids[best_row,],pxdist_12,pxdist_13,pxdist_14,pxdist_23,pxdist_24,pxdist_34)
        rms_fids <- rms.calc(best_image_dist,best_film_fids)
        print(rms_fids) # RMS among fiducials in units of pixels
        
        
        # Georeference the film coordinates to the image coordinatess
        
        # First extract the best image fiducial coordinates:
        x_image_fids <- t(select(image_fids[best_row,],x1,x2,x3,x4))
        y_image_fids <- t(select(image_fids[best_row,],y1,y2,y3,y4))
        best_image_fids <- data.frame(x_image_fids,y_image_fids)
        colnames(best_image_fids) <- c("x","y")
        rownames(best_image_fids) <- NULL
        
        # Then extract the best film fiducial coordinates:
        film_coordx <- t(select(film_fids[best_row,],x1,x2,x3,x4) / mission_res) # Convert the film "x" best row into units of pixels
        film_coordy <- t(select(film_fids[best_row,],y1,y2,y3,y4) / mission_res) # Convert the film "y" best row into units of pixels
        film_coord <- data.frame(x = film_coordx, y = film_coordy)
        colnames(film_coord) <- c("x","y")
        rownames(film_coord) <- NULL
        # Reflect the film coordinates across the y-axis since we're viewing the film from the back:
        film_coord$x <- film_coord$x * -1
        film_coord$y <- film_coord$y * -1
        # Calcualte translation factor for the film coordinates so that Fid 1 on film = Fid 1 on the image:
        x_trans <- film_coord[1,1] - best_image_fids[1,1]
        y_trans <- film_coord[1,2] - best_image_fids[1,2]
        # Translate the fim coordinates
        film_coord$x_trans <- film_coord$x - x_trans
        film_coord$y_trans <- film_coord$y - y_trans
        list(film_coord,best_image_fids,rms_fids)
        
}


