
# This is a Packrat project.

library(here)






# Enter data from the camera calibration report and a given image --------------

# Fid = fiducial
# Numbers refer to the numbering system in the report (fiducial 1, 2, 3, and 4...usually)
# Sub_x = x coordinate, sub_y = y coordinate
# Px = units of pixels, mm = units of millimeters

# Enter the number of pixels to the fiducial mark and the film coordinates for
# the same mark in mm. Use a pixel/image coordinate system for the pixels (0,0
# at top left). Note that the pixel coordinate will later be refined...at this
# point, just make sure the coordinate is +/- 1 pixel from the true position.
# Note also that the origin of the image and film coordinate systems do not need
# to match.
pix_coord <-
        data.frame(
                fid = c(1, 2, 3, 4), # Fiducial number
                px_x = c(798, 15947, 801, 15942), # x-coordinate in pixels
                px_y = c(15978, 839, 836, 15985), # y-coordinate in pixels
                mm_x = c(-105.993, 106.001, -106.015, 106.015), # x-coordinate in mm
                mm_y = c(-106.006, 106.002, 105.982, -106.006) # y-coordinate in mm
        )

# Enter the distance between fiducails in mm:
dist2fid_mm <-
        data.frame(
                fids = c("1_2", "1_3", "1_4", "2_3", "2_4", "3_4"), # For example, 1_2 is distance between fiducial 1 and 2
                dist_mm = c(299.814, 211.988, 212.015, 212.008, 212.007, 299.826)
        )

# Enter the error in the film measurements as stated in the report. If none is given, estimate a reasonable value:
dist_err <- 0.003

# Enter the DN's for each pixel in a 9x9 window centered at the pixel at the
# percieved center of the fiducial. Start in the top left position in the
# matrix. If the image is in RGB, do all three channels before moving to the
# next cell.

fid1_dn <- matrix(c(255,249,222,253,252,231,255,250,224, # Row 1, rgb,rgb,...
                     255,255,234,255,255,243,255,251,228, # Row 2
                     255,250,223,255,255,234,255,250,220) # Row 3
                   ,nrow = 3)

fid2_dn <- matrix(c(255,248,212,255,252,230,255,255,233, # Row 1m
                    255,241,204,255,252,238,248,246,231, # Row 2
                    255,232,201,255,246,230,255,246,228) # Row 3
                  ,nrow = 3)

fid3_dn <- matrix(c(255,242,219,255,248,219,255,249,216, # Row 1
                     255,245,236,255,253,246,255,245,236, # Row 2
                     255,241,235,255,244,241,255,243,239) # Row 3
                   ,nrow = 3)

# fid4_dn <- matrix(c(255,242,224,255,251,225,255,238,200, # Row 1
#                     255,250,228,255,241,221,255,243,201, # Row 2
#                     255,230,207,255,232,203,255,206,169) # Row 3
#                   ,nrow = 3)

fid4_dn <- matrix(c(255,216,198,255,242,224,255,251,225, # Row 1
                    255,219,200,255,250,228,255,241,221, # Row 2
                    235,184,165,255,230,207,255,232,203) # Row 3
                  ,nrow = 3)



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


fid_dn <- fid3_dn
fid_no <- 3

# Function to find the coordinates of a fiducial point. Inputs: fid_dn = 9x3 matrix of RGB data
fid.coord <- function(fid_dn,pix_coord,fid_no){
        
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
        subpix_x <- ifelse(max_coord[1] == 1, pix_coord[fid_no,2] - 1 + row_fid,
                           ifelse(max_coord[1] == 2, pix_coord[fid_no,2] + row_fid,
                                  pix_coord[fid_no,2] + 1 + row_fid))
        
        subpix_y <- ifelse(max_coord[2] == 1, pix_coord[fid_no,3] - 1 + col_fid,
                           ifelse(max_coord[2] == 2, pix_coord[fid_no,3] + col_fid,
                                  pix_coord[fid_no,3] + 1 + col_fid))
        coord <- c(subpix_x,subpix_y,0.5,0.5)
        names(coord) <- c("x","y","x_err","y_err")
        coord_summary <- c(coord,max_coord)
        coord_summary
        
}





# Find the coordinates of a fiducial point -------------------------------------

# Important: the center row and center column must both be 2. If not, re-center the pix_coord matrix.

# Fiducial 1:
fid1_coord <-fid.coord(fid1_dn,pix_coord,1)

# Fiducial 2:
fid2_coord <- fid.coord(fid2_dn,pix_coord,2)

# Fiducial 3:
fid3_coord <- fid.coord(fid3_dn,pix_coord,3)

# Fiducial 4:
fid4_coord <- fid.coord(fid4_dn,pix_coord,4)

fid_coord <- rbind(fid1_coord,fid2_coord,fid3_coord,fid4_coord)
fid_coord <- as.data.frame(fid_coord)


# Compute pixel resolution -----------------------------------------------------

# Compute the x-axis distance between fiducial marks in untis of pixels
fx_1_2 <- fid_coord$x[1] - fid_coord$x[2]
fx_1_3 <- fid_coord$x[1] - fid_coord$x[3]
fx_1_4 <- fid_coord$x[1] - fid_coord$x[4]
fx_2_3 <- fid_coord$x[2] - fid_coord$x[3]
fx_2_4 <- fid_coord$x[2] - fid_coord$x[4]
fx_3_4 <- fid_coord$x[3] - fid_coord$x[4]

# Compute the y-axis distance between fiducial marks in untis of pixels
fy_1_2 <- fid_coord$y[1] - fid_coord$y[2]
fy_1_3 <- fid_coord$y[1] - fid_coord$y[3]
fy_1_4 <- fid_coord$y[1] - fid_coord$y[4]
fy_2_3 <- fid_coord$y[2] - fid_coord$y[3]
fy_2_4 <- fid_coord$y[2] - fid_coord$y[4]
fy_3_4 <- fid_coord$y[3] - fid_coord$y[4]

# Gather the pixel distances into a df:
dist2fid_px <-
        data.frame(
                f = c("1_2", "1_3", "1_4", "2_3", "2_4", "3_4"), # Fiducial pair with computed distance
                fx = c(fx_1_2, fx_1_3, fx_1_4, fx_2_3, fx_2_4, fx_3_4), # Computed distance along x-axis
                fy = c(fy_1_2, fy_1_3, fy_1_4, fy_2_3, fy_2_4, fy_3_4) # Computed distance along y-axis
        )

dist2fid_px$dist <- sqrt(dist2fid_px$fx ^ 2 + dist2fid_px$fy ^ 2) # Compute 
dist2fid_px$dist_err <- dist_err
dist2fid_px$res <- dist2fid_mm$dist / dist2fid_px$dist
scan_res <- mean(dist2fid_px$res) # Scanning resolution in mm
scan_sd <- sd(dist2fid_px$res)





