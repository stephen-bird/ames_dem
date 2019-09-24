
# This is a Packrat project.

library(here)






# Enter data from the camera calibration report and a given image --------------

# Fid = fiducial
# Numbers refer to the numbering system in the report (fiducial 1, 2, 3, and 4...usually)
# Sub_x = x coordinate, sub_y = y coordinate
# Px = units of pixels, mm = units of millimeters

# Enter the number of pixels to the fiducial mark and the film coordinates for
# the same mark in mm. Use a pixel/image coordinate system for the pixels (0,0
# at top left). Note that the origin of the image and film coordinate systems do
# not need to match.
pix_coord <-
        data.frame(
                fid = c(1, 2, 3, 4), # Fiducial number
                px_x = c(798, 15947, 801, 15943), # x-coordinate in pixels
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

fid1_dn <- matrix(c(255,249,222,253,252,231,255,250,224, # Red
                     255,255,234,255,255,243,255,251,228, # Blue
                     255,250,223,255,255,234,255,250,220) # Green
                   ,nrow = 3)


fid2_dn <- matrix(c(255,248,212,255,252,230,255,255,233, # Red
                    255,241,204,255,252,238,248,246,231, # Blue
                    255,232,201,255,246,230,255,246,228) # Green
                  ,nrow = 3)

fid3_dn <- matrix(c(255,242,219,255,248,219,255,249,216, # Red
                     255,245,236,255,253,246,255,245,236, # Blue
                     255,241,235,255,244,241,255,243,239) # Green
                   ,nrow = 3)

fid4_dn <- matrix(c(255,242,224,255,251,225,255,238,200, # Red
                    255,250,228,255,241,221,255,243,201, # Blue
                    255,230,207,255,232,203,255,206,169) # Green
                  ,nrow = 3)


# Load functions ---------------------------------------------------------------

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
line.x <- function(lm1,lm2){
        slope_diff <- lm1$coefficients[2] - lm2$coefficients[2]
        intercept_diff <- lm2$coefficients[1] - lm1$coefficients[1]
        solved_x <- intercept_diff/slope_diff
        solved_y <- lm1$coefficients[2] * solved_x + lm1$coefficients[1]
        point <- c(solved_x,solved_y)
        names(point) <- c("x","y")
        point
}


# Function to find the coordinates of a fiducial point
fid.coord <- function(fid_dn,top_left){
        
        # Sum the DN's for each channel (e.g. red + green + blue) and save as a matrix:
        fid_mx <- subfid.mx(1,fid_dn) + subfid.mx(2,fid_dn) + subfid.mx(3,fid_dn)
        
        # Find the sub-pixel center for each row and column:
        row1_fid <- sub.px(fid_mx[1,1],fid_mx[1,2],fid_mx[1,3])
        row2_fid <- sub.px(fid_mx[2,1],fid_mx[2,2],fid_mx[2,3])
        row3_fid <- sub.px(fid_mx[3,1],fid_mx[3,2],fid_mx[3,3])
        col1_fid <- sub.px(fid_mx[1,1],fid_mx[2,1],fid_mx[3,1])
        col2_fid <- sub.px(fid_mx[1,2],fid_mx[2,2],fid_mx[3,2])
        col3_fid <- sub.px(fid_mx[1,3],fid_mx[2,3],fid_mx[3,3])
        
        # Use the sub-pixel column centers to plot a single horizontal line:
        horz_fid <- data.frame(x=c(1,2,3),y=c(col1_fid,col2_fid,col3_fid))
        horz_lm <- lm(y ~ x, data = horz_fid) # Fit a linear model
        
        # Use the sub-pixel row centers to plot a single vertical line:
        vert_fid <- data.frame(x=c(1,2,3),y=c(row1_fid,row2_fid,row3_fid))
        vert_lm <- lm(y ~ x, data = vert_fid) # Fit a linear model
        
        coord <- line.x(horz_lm,vert_lm) # Find the intersect between the two lines
        coord <- c(top_left[1] - coord[1], top_left[2] + coord[2]) # Translate to pixel coordinates (origin is in top left)
        error <- c(summary(horz_lm)$sigma,summary(vert_lm)$sigma) # Append the uncertainty estimates
        names(error) <- c("x_se","y_se")
        coord <- c(coord,error)
        coord
        
}



# Compute pixel resolution -----------------------------------------------------

# Compute the x-axis distance between fiducial marks in untis of pixels
fx_1_2 <- pix_coord$px_x[1] - pix_coord$px_x[2]
fx_1_3 <- pix_coord$px_x[1] - pix_coord$px_x[3]
fx_1_4 <- pix_coord$px_x[1] - pix_coord$px_x[4]
fx_2_3 <- pix_coord$px_x[2] - pix_coord$px_x[3]
fx_2_4 <- pix_coord$px_x[2] - pix_coord$px_x[4]
fx_3_4 <- pix_coord$px_x[3] - pix_coord$px_x[4]

# Compute the y-axis distance between fiducial marks in untis of pixels
fy_1_2 <- pix_coord$px_y[1] - pix_coord$px_y[2]
fy_1_3 <- pix_coord$px_y[1] - pix_coord$px_y[3]
fy_1_4 <- pix_coord$px_y[1] - pix_coord$px_y[4]
fy_2_3 <- pix_coord$px_y[2] - pix_coord$px_y[3]
fy_2_4 <- pix_coord$px_y[2] - pix_coord$px_y[4]
fy_3_4 <- pix_coord$px_y[3] - pix_coord$px_y[4]

# Gather the pixel distances into a df:
dist2fid_px <-
        data.frame(
                f = c("1_2", "1_3", "1_4", "2_3", "2_4", "3_4"), # Fiducial pair with computed distance
                fx = c(fx_1_2, fx_1_3, fx_1_4, fx_2_3, fx_2_4, fx_3_4), # Computed distance along x-axis
                fy = c(fy_1_2, fy_1_3, fy_1_4, fy_2_3, fy_2_4, fy_3_4) # Computed distance along y-axis
        )

dist2fid_px$dist <- sqrt(dist2fid_px$fx ^ 2 + dist2fid_px$fy ^ 2) # Compute 
dist2fid_px$dist_err <- dist_err
scan_res <- mean(dist2fid_mm$dist / dist2fid_px$dist) # Scanning resolution in mm

# pix_coord$pix_dif <- pix_coord$mm_x - pix_coord$mm_y
# pix_coord$mm_dif <- pix_coord$px_x - pix_coord$px_y


# Find the coordinates of a fiducial point -------------------------------------

# Fiducial 1:
top_left <- c(pix_coord$px_x[1],pix_coord$px_y[1])
fid1_coord <-fid.coord(fid1_dn,top_left)

# Fiducial 2:
top_left <- c(pix_coord$px_x[2],pix_coord$px_y[2])
fid2_coord <- fid.coord(fid2_dn,top_left)

# Fiducial 3:
top_left <- c(pix_coord$px_x[3],pix_coord$px_y[3])
fid3_coord <- fid.coord(fid3_dn,top_left)

# Fiducial 4:
top_left <- c(pix_coord$px_x[4],pix_coord$px_y[4])
fid4_coord <- fid.coord(fid4_dn,c(0,0))


























horz_lm$coefficients[2] - vert_lm$coefficients[2]
vert_lm$coefficients[1] - horz_lm$coefficients[1]







horz_lm$coefficients[1] + 



plot(vert_lm)

predict(vert_lm,newdata = data.frame(x=2,y=1))



r_val <- fid_val[1,]
r_val <- matrix(r_val,ncol=3)
r_val <- t(r_val)
g_val <- fid_val[2,]
g_val <- matrix(g_val,ncol=3)
g_val <- t(g_val)
b_val <- fid_val[3,]
b_val <- matrix(b_val,ncol=3)
b_val <- t(b_val)
r_val + g_val + b_val



sub.px(242,248,249)
sub.px(219,219,216)

sub.px(245,253,236)
sub.px(236,246,236)

sub.px(245,253,236)
sub.px(236,246,236)

