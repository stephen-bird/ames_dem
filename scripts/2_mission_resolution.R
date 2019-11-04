
# Projet Name: Frame Camera Calibration
# Script Name: 2_mission_resolution

# Stephen Bird, 2019-11-04

# This script calculates resolution of all images acquired as part of a single mission.


# Manual data entry ------------------------------------------------------------

# Enter images in the direction of the flight line:

image_no <- c("bcc04050_043","bcc04050_042","bcc04050_041")
image_dimensions <- data.frame(image_no, x_pixels = rep(16730,3), y_pixels = rep(16816,3)) # 


# Enter the DN's for each pixel in a 9x9 window centered at the pixel at the
# perceived center of the fiducial. Start in the top right position in the
# matrix (assuming we're viewing the image from the focal plane). If the image
# is in RGB, do all three channels before moving to the next cell. Filter the
# fiducial by Outliers, Despekle, and Gaussian Blur using the Fiji tool.


# Image bcc04050_043 -----------------------------------------------------------

fid1_pxy <- c(15946,15983)
fid1_dn <- matrix(c(84,55,49,97,63,56,96,60,53, # Row 1,
                    96,64,58,111,73,66,110,70,62, # Row 2
                    96,63,57,111,72,65,110,69,61), # Row 3
                  nrow = 3)

fid2_pxy <- c(798,842)
fid2_dn <- matrix(c(131,97,79,133,98,78,119,83,66, # Row 1,
                    145,112,90,147,113,90,131,97,75, # Row 2
                    141,111,88,143,112,88,129,96,74), # Row 3
                  nrow = 3)

fid3_pxy <- c(15943,838)
fid3_dn <- matrix(c(193.999,130,110,193,131,110,178,113,95, # Row 1, rgb,rgb,...
                    195.999,137,115,196,137,115,181,118,99, # Row 2
                    184,125,104,184,125,104,171,108,89) # Row 3
                  ,nrow = 3)

fid4_pxy <- c(802,15985)
fid4_dn <- matrix(c(185,135,110,191,145,118,183,135,109, # Row 1,
                    186,135,111,192,145,119,183,135,110, # Row 2
                    172,117,97,177,125,103,169,117,95), # Row 3
                  nrow = 3)

# Use the "mission.res" function (4) to calculate the image resolution:
bcc04050_043 <- mission.res(fid1_pxy,fid2_pxy,fid3_pxy,fid4_pxy,film_coord,film_fids,n)
bcc04050_043


# Image bcc04050_042 -----------------------------------------------------------

fid1_pxy <- c(15944,15986)
fid1_dn <- matrix(c(89,58,51,103,69,60,102,68,59, # Row 1,
                    98,63,55,113,74,65,111,73,64, # Row 2
                    93,57,50,107,67,59,105,66,58), # Row 3
                  nrow = 3)

fid2_pxy <- c(801,840)
fid2_dn <- matrix(c(127,90,79,139,102,90,135,97,86, # Row 1,
                    131,96,84,144,108,95,140,103,91, # Row 2
                    120,86,75,132,97,85,128,93,81), # Row 3
                  nrow = 3)

fid3_pxy <- c(15946,840)
fid3_dn <- matrix(c(174,112,88,181,121,95,175,114,88, # Row 1, rgb,rgb,...
                    186,128,102,193,138,110,187,131,103, # Row 2
                    183,126,103,191,137,111,185,129,104) # Row 3
                  ,nrow = 3)

fid4_pxy <- c(800,15982)
fid4_dn <- matrix(c(178,131,113,183,135,116,172,119,102, # Row 1,
                    187,140,119,192,143,122,180,127,108, # Row 2
                    180,130,108,185,133,111,174,118,98), # Row 3
                  nrow = 3)

# Use the "mission.res" function (4) to calculate the image resolution:
bcc04050_042 <- mission.res(fid1_pxy,fid2_pxy,fid3_pxy,fid4_pxy,film_coord,film_fids,n)
bcc04050_042


# Image bcc04050_041 -----------------------------------------------------------

fid1_pxy <- c(15943,15985)
fid1_dn <- matrix(c(106,74,65,111,77,67,101,68,59, # Row 1,
                    108,77,66,113,80,69,102,71,60, # Row 2
                    95,67,57,99,70,60,90,62,52), # Row 3
                  nrow = 3)

fid2_pxy <- c(801,836)
fid2_dn <- matrix(c(127,84,74,141,97,85,139,95,81, # Row 1,
                    138,94,83,154,107,95,151,104,91, # Row 2
                    133,88,78,148,101,89,146,97,86), # Row 3
                  nrow = 3)

fid3_pxy <- c(15947,839)
fid3_dn <- matrix(c(188,123,102,188,125,102,174,110,88, # Row 1, rgb,rgb,...
                    194,134,112,195,137,113,180,120,98, # Row 2
                    187,128,107,187,130,108,174,114,93) # Row 3
                  ,nrow = 3)

fid4_pxy <- c(798,15978)
fid4_dn <- matrix(c(181,123,104,185,129,108,174,118,96, # Row 1,
                    193,135,114,197,142,118,185,129,105, # Row 2
                    190,129,108,193,135,112,181,123,99), # Row 3
                  nrow = 3)

# Use the "mission.res" function (4) to calculate the image resolution:
bcc04050_041 <- mission.res(fid1_pxy,fid2_pxy,fid3_pxy,fid4_pxy,film_coord,film_fids,n)
bcc04050_041


# Derive the resolution of each image and then average the results -------------

mean_res <- c(bcc04050_041[[2]],bcc04050_042[[2]],bcc04050_043[[2]]) # Extract element 2 from each list for the mean
sd_res <- c(bcc04050_041[[3]],bcc04050_042[[3]],bcc04050_043[[3]]) # Extract element 3 from each list for the sd
# Use the weighted mean function (5) to calculate the weigthed mean of all calculated image resolutions based on Taylor
mission_res <- weighted.mean(mean_res,sd_res) # Calculate weighted mean in units of mm
print(mission_res)

