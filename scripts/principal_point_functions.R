
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


# Function to find the coordinates of two intersecting lines. The inputs include
# linear models ("lm1" and "lm2") fitted through the horizontal and vertical
# sub-pixel centers (the order doens't matter).
# line.x <- function(lm1,lm2){
#         slope_diff <- lm1$coefficients[2] - lm2$coefficients[2]
#         intercept_diff <- lm2$coefficients[1] - lm1$coefficients[1]
#         solved_x <- intercept_diff/slope_diff
#         solved_y <- lm1$coefficients[2] * solved_x + lm1$coefficients[1]
#         point <- c(solved_x,solved_y)
#         names(point) <- c("x","y")
#         point
# }


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
        intersect_ang <- abs((180 / pi) * atan((slope_line12 - slope_line34) / (1 + slope_line12 * slope_line34))) # Intersection angle line 12 vs line 34
        intersect_ang
        