
# Load functions ---------------------------------------------------------------

# Determine if the 9x9 matrix contains the maximum pixel value in the center of the matrix:


# This function takes a matrix of pixel DN's for each channel and returns a
# matrix for a given channel (1 = red, 2 = blue, 3 = green).

subfid.mx <- function(channel,fid_dn){
        
        chan_dn <- fid_dn[channel,] # Subset the matrix by channel
        chan_dn <- matrix(chan_dn,ncol=3) # Rebuild the original (but transposed) matrix by row and column
        t(chan_dn) # Transpose to the original matrix by channel
        
}


# Function to calculate the sub-pixel center of a line made-up of three points.
# Fucntion assumes the center coordinate is the middle.
# pa = dn at pixel 1; pb = dn at pixel 2; c = dn at pixel 3

sub.px <- function(pa,pb,pc) {
        
        0.5*(pa - pc) / (pa-2*pb+pc)
        
}


# Function to find the coordinates of two intersecting lines.
line.sect <- function(slope1,slope2,inter1,inter2){
        slope_diff <- slope1 - slope2
        intercept_sum <- inter1 + inter2
        solved_x <- intercept_sum/slope_diff
        solved_y <- slope1 * solved_x + inter1
        point <- c(solved_x,solved_y)
        names(point) <- c("x","y")
        point
}





# Function to find the coordinates of a fiducial point. Inputs: fid_dn = 9x3 matrix of RGB data
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


# Function to calcualte the angle between crossing lines that connect fiducials based on the 10,000 random sample of
# x,y coordinates:
intersect.angle <- function(fid_pts){
        
        # Input includes 8 columns of data: columns 1:4 are the x coordinates of the fiducials, columns 5:8 are the corresponding y coordinates.
        slope_line12 <- ((fid_pts[,5] - fid_pts[,6]) / (fid_pts[,1] - fid_pts[,2])) # Slope of line 1_2
        slope_line34 <- ((fid_pts[,7] - fid_pts[,8]) / (fid_pts[,3] - fid_pts[,4])) # Slope of line 3_4
        intersect_ang <- ((180 / pi) * atan((slope_line12 - slope_line34) / (1 + slope_line12 * slope_line34))) # Intersection angle line 12 vs line 34
        intersect_ang
        
}


hypo <- function(x1, y1, x2, y2) {
        sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2)
}


fiducial.samples <- function(fid_coord,n) {
        # Generate samples based on the fiducial mark uncertainty. Generate normal
        # random samples (n = 10,000) for each sub-pixel estimatation of the image
        # fiducial center (+/- 0.5 pixels) by x and y coordinate.
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


# Calculate the residuals between the image-pixel distances and the film-pixel distances between fiducials. Residuals are Internally Studentized (divide by sd) so they can be scaled and compared to one another:


resid.calc <- function(image_dist,film_pix_dist){
        
        resid <- image_dist - film_pix_dist
        resid_sd <- sd(resid)
        resid_st <- resid/resid_sd
        resid_st2 <- resid_st^2
        resid_st2
}

rms.calc <- function(best_image_fids,best_film_fids){
        
        resid <- best_image_fids - best_film_fids
        resid2 <- resid^2
        sum_resid2 <- sum(resid2)
        sqrt(sum_resid2)
        
}


# Function to minimize the shift between measured sub-pixel centers and those dervied from the distance between fiducials using a simple translation.
min.x <- function(data, par) {
        
        with(data, abs(sum(
                
                (derived - par) - measured)))
}

