rm(list = ls())
library(tidyverse)
library(moderndive)
library(ISLR)
library(infer)
library(rgl)
library(ggplot2movies)
library(nycflights13)
library(plot3D)
library(plot3Drgl)
options(rgl.printRglwidget = TRUE)
# mvmeta for symmetric matrix vectorization (vechMat)
library(mvmeta)
library(data.table)
# from: https://search.r-project.org/CRAN/refmans/dplR/html/pass.filt.html
library(dplR) # per filtering
library(pracma) # per cross product
# from: https://search.r-project.org/CRAN/refmans/matlib/html/vectors3d.html
library(matlib)

toRad <- pi/180
toDeg <- 1/toRad

##################################################################################################
#### Funzioni base per matrici di rotazione
##################################################################################################
eucMatrix <- function(roll, pitch, yaw) {
  RES <- matrix(c(
    cos(yaw)*cos(pitch), cos(yaw)*sin(pitch)*sin(roll)-sin(yaw)*cos(roll), cos(yaw)*sin(pitch)*cos(roll)+sin(yaw)*sin(roll),
    sin(yaw)*cos(pitch), sin(yaw)*sin(pitch)*sin(roll)+cos(yaw)*cos(roll), sin(yaw)*sin(pitch)*cos(roll)-cos(yaw)*sin(roll),
    -sin(pitch), cos(pitch)*sin(roll), cos(pitch)*cos(roll)
  ), nrow=3, ncol=3, byrow = TRUE)
  return(RES)
}
eucMatrixRPY <- function(rpy) {
  roll <- rpy[1]
  pitch <- rpy[2]
  yaw =   rpy[3]
  return(eucMatrix(roll, pitch, yaw))
}
crossMatrix <- function(u) {
  RES <- matrix(c(
    0, -u[3], u[2],
    u[3], 0, -u[1],
    -u[2], u[1], 0
  ), nrow=3, ncol=3, byrow=TRUE)
  return(RES)
}
rotMatrix <- function(teta, u) {
  I <- diag(1,3,3)
  RES <- I*cos(teta) + sin(teta)*crossMatrix(u) + (1-cos(teta))*u%*%t(u)
  return(RES)
}
toIFMatrix <- function(rpy) {
  return(eucMatrix(rpy[1],rpy[2],rpy[3]))
}
toIF <- function(rpy, v) {
  return(toIFMatrix(rpy) %*% v)
}
toBFMatrix <- function(rpy) {
  return(t(toIFMatrix(rpy)))
}
toBF <- function(rpy, v) {
  return(toBFMatrix(rpy) %*% v)
}
fromMFFToIFMatrix <- function(declination) {
  return(eucMatrix(0, (pi/2+declination),0))
}
fromIFToMFFMatrix <- function(declination) {
  return(t(fromMFFToIFMatrix(declination)))  
}
fromIFToMFF <- function(declination, v) {
  return(fromIFToMFFMatrix(declination) %*% v)
}
fromMFFToIF <- function(declination, v) {
  return(fromMFFToIFMatrix(declination) %*% v)
}
makeDecRefInIF <- function(declination) {
  return(fromMFFToIF(declination, as.matrix(c(0,0,1))))
}
# RPY reference is vector c(0,0,1) = -g
calcRP <- function(v) {
  ipitch <- -asin(v[1])
  
  arg <- min(1,max(-1,v[3]/cos(ipitch)))
  iroll <- acos(arg)
  if(v[2] < 0) {
    iroll <- -iroll
  }
  
  if(iroll > pi) {
    iroll = iroll - 2*pi
  } else if(iroll < -pi) {
    iroll = iroll + 2*pi
  }
  iyaw <- -atan2(v[2], v[1])
  result <- as.matrix(c(iroll, ipitch, iyaw))
  return(result)
}
calcYaw <- function(v) {
  iyaw <- -atan2(v[2], v[1])
  result <- as.matrix(c(0, 0, iyaw))
  return(result)
}
toBFTest <- function(rpy, v, v_expected) {
  v_actual <- toBF(rpy, v)
  stopifnot(round(v_actual,7) == round(v_expected,7))
  return(TRUE)
}
toIFTest <- function(rpy, v, v_expected) {
  v_actual <- toIF(rpy, v)
  stopifnot(round(v_actual,7) == round(v_expected,7))
  return(TRUE)
}
calcRPTest <- function(rpy) {
  rpy[3]<-0
  v <- toBF(rpy, as.matrix(c(0,0,1)))
  new_rpy <- calcRP(v)
  new_rpy[3]<-0
  v1 <- toBF(new_rpy, as.matrix(c(0,0,1)))
  stopifnot(round(v - v1, 6) == c(0,0,0))
}
calcYawTest <- function(rpy) {
  rpy[1]<-0
  rpy[2]<-0
  v <- toBF(rpy, as.matrix(c(1,0,0)))
  new_rpy <- calcYaw(v)
  stopifnot(round(rpy - new_rpy, 10) == c(0,0,0))
}
toBFTest(as.matrix(c(0,0,90*toRad)), as.matrix(c(1,0,0)),as.matrix(c(0,-1,0)))
toIFTest(as.matrix(c(0,0,90*toRad)),as.matrix(c(0,-1,0)), as.matrix(c(1,0,0)))
toBFTest(as.matrix(c(0,90*toRad,0)), as.matrix(c(1,0,0)),as.matrix(c(0,0,1)))
toIFTest(as.matrix(c(0,90*toRad,0)),as.matrix(c(0,0,1)), as.matrix(c(1,0,0)))
toBFTest(as.matrix(c(0,90*toRad,90*toRad)), as.matrix(c(1,0,0)),as.matrix(c(0,-1,0)))
toIFTest(as.matrix(c(0,90*toRad,90*toRad)),as.matrix(c(0,-1,0)), as.matrix(c(1,0,0)))
toBFTest(as.matrix(c(90*toRad,0,0*toRad)), as.matrix(c(0,0,1)),as.matrix(c(0,1,0)))
toBFTest(as.matrix(c(10*toRad,0,0*toRad)), as.matrix(c(0,0,1)),as.matrix(c(0,0.1736482,0.9848078)))

calcRPTest(as.matrix(c(0,0,90*toRad)))
calcRPTest(as.matrix(c(10*toRad,20*toRad,90*toRad)))
calcRPTest(as.matrix(c(-10*toRad,-20*toRad,90*toRad)))
calcRPTest(as.matrix(c(100*toRad,0*toRad,0*toRad)))
calcRPTest(as.matrix(c(-100*toRad,0*toRad,0*toRad)))
calcYawTest(as.matrix(c(0,0,90*toRad)))
calcYawTest(as.matrix(c(10*toRad,20*toRad,90*toRad)))
calcYawTest(as.matrix(c(-10*toRad,-20*toRad,90*toRad)))
calcYawTest(as.matrix(c(100*toRad,0*toRad,90*toRad)))
calcYawTest(as.matrix(c(-100*toRad,.200*toRad,90*toRad)))

for(roll in -180:180) {
  for(pitch in -180:180) {
    calcRPTest(as.matrix(c(roll*toRad,pitch*toRad,0)))
  }
}

#######################################################################################
### Load Data
#######################################################################################
#imu.cal.data.gyro <- read.csv('/hd/eclipse-cpp-2020-12/eclipse/workspace/my_mpu9250_thing/examples/imu-cal-data-gyro-pitch-360.csv')
#imu.cal.data.gyro <- read.csv('/hd/eclipse-cpp-2020-12/eclipse/workspace/my_mpu9250_thing/examples/imu-cal-data-gyro-rpy-360.csv')
imu.cal.data.gyro <- read.csv('/hd/eclipse-cpp-2020-12/eclipse/workspace/my_mpu9250_thing/examples/imu-cal-data-gyro-pitch.csv')
imu.cal.data.gyro <- imu.cal.data.gyro[2:dim(imu.cal.data.gyro)[1],]

# plot original data
scatter3D(imu.cal.data.gyro$MX, imu.cal.data.gyro$MY, imu.cal.data.gyro$MZ, colvar = imu.cal.data.gyro$MZ, col = NULL, add = FALSE, ticktype = "detailed", scale = FALSE)
scatter3D(imu.cal.data.gyro$AX, imu.cal.data.gyro$AY, imu.cal.data.gyro$AZ, colvar = imu.cal.data.gyro$AZ, col = NULL, add = TRUE, ticktype = "detailed", scale = FALSE)
plotrgl()
rglwidget()

mag_plot_data <- function(mag_data, sphere_radius = -1, title = "") {
  scatter3D(mag_data[,1], mag_data[,2], mag_data[,3], colvar = mag_data[,3], col = NULL, add = FALSE, scale = FALSE, ticktype = "detailed", main = title)
  plotrgl()
  rglwidget()
}

############################################################################################################################
##### Remove Gyro Bias
##### Calculate Roll, Pitch, Yaw from Accelerometer & Magnetometer
############################################################################################################################
gyro_bias = c(mean(imu.cal.data.gyro$GX[1:10000]), mean(imu.cal.data.gyro$GY[1:10000]), mean(imu.cal.data.gyro$GZ[1:10000]))
imu.cal.data.gyro.cnt <- imu.cal.data.gyro %>% mutate(GX = GX - gyro_bias[1],GY = GY - gyro_bias[2],GZ = GZ - gyro_bias[3])

imu.cal.data.gyro_rp <- imu.cal.data.gyro.cnt %>% 
  mutate(NRMA = sqrt(AX^2+AY^2+AZ^2)) %>%
  mutate(NRMM = sqrt(MX^2+MY^2+MZ^2)) %>%
  mutate(MX = MX/NRMM, MY = MY/NRMM, MZ = MZ/NRMM) %>%
  mutate(AX = AX/NRMA, AY = AY/NRMA, AZ = AZ/NRMA) %>%
  mutate(PA = -asin(AX)*toDeg) %>%
  mutate(RA = acos(round(AZ/cos(PA*toRad),10))*toDeg) %>%
  mutate(RA = case_when(AY < 0 ~ -RA, TRUE ~ RA)) %>%
  mutate(RA = case_when(RA*toRad > pi ~ (RA*toRad-2*pi)*toDeg,RA*toRad < -pi ~ (RA*toRad+2*pi)*toDeg, TRUE ~ RA)) %>%
  mutate(
    ICMX = MX * cos(PA*toRad)   + MY*sin(PA*toRad)*sin(RA*toRad)  + MZ * sin(PA*toRad)*cos(RA*toRad),
    ICMY = 0                + MY*cos(RA*toRad)            - MZ*sin(RA*toRad),
    ICMZ = -MX*sin(PA*toRad)    + MY*cos(PA*toRad)*sin(RA*toRad)  + MZ*cos(PA*toRad)*cos(RA*toRad), 
    YM = atan2(-ICMY,ICMX)*toDeg,
    ICAX = AX * cos(PA*toRad)   + AY*sin(PA*toRad)*sin(RA*toRad)  + AZ * sin(PA*toRad)*cos(RA*toRad),
    ICAY = 0                + AY*cos(RA*toRad)            - AZ*sin(RA*toRad),
    ICAZ = -AX*sin(PA*toRad)    + AY*cos(PA*toRad)*sin(RA*toRad)  + AZ*cos(PA*toRad)*cos(RA*toRad), 
  ) %>% 
  mutate(
    IICMX = ICMX * cos(YM*toRad) - ICMY*sin(YM*toRad),
    IICMY = ICMX*sin(YM*toRad) + ICMY*cos(YM*toRad),
    IICMZ = ICMZ,
    IICAX = ICAX * cos(YM*toRad) - ICAY*sin(YM*toRad),
    IICAY = ICAX*sin(YM*toRad) + ICAY*cos(YM*toRad),
    IICAZ = ICAZ
  ) %>%
  mutate(DEC = acos(IICMX/sqrt(IICMX^2+IICMY^2+IICMZ^2))*toDeg) %>%
  rename(AROLL = RA, APITCH = PA, AYAW = YM, IMX = IICMX, IMY = IICMY, IMZ = ICMZ, IAX = IICAX, IAY = IICAY, IAZ = IICAZ) %>%
  mutate(DAROLL = AROLL - shift(AROLL, 1, fill = AROLL[1], type = "lag")) %>%
  mutate(DAPITCH = APITCH - shift(APITCH, 1, fill = APITCH[1], type = "lag")) %>%
  mutate(DAYAW = AYAW - shift(AYAW, 1, fill = AYAW[1], type = "lag")) %>%
  mutate(DAROLL = case_when(DAROLL > pi*toDeg ~ 2*pi*toDeg - DAROLL,  DAROLL < -pi*toDeg ~  2*pi*toDeg + DAROLL, TRUE ~ DAROLL)) %>%
  mutate(DPITCH = case_when(DAPITCH > pi*toDeg ~ 2*pi*toDeg - DAPITCH,  DAPITCH < -pi*toDeg ~  2*pi*toDeg + DAPITCH, TRUE ~ DAPITCH)) %>%
  mutate(DAYAW = case_when(DAYAW > pi*toDeg ~ 2*pi*toDeg - DAYAW,  DAYAW < -pi*toDeg ~  2*pi*toDeg + DAYAW, TRUE ~ DAYAW)) %>%
  select(TS, MX, MY, MZ, AX, AY, AZ, GX, GY, GZ, AROLL, APITCH, AYAW, DAROLL, DAPITCH, DAYAW, DEC, IMX, IMY, IMZ, IAX, IAY, IAZ, MV)
  
############################################################################################################################
##### Gyroscope data
############################################################################################################################
# gyro data are on 1000DPS => LSB = 32.768
gyro_data <- imu.cal.data.gyro_rp %>% 
  select(GX, GY, GZ, MV) %>%
  mutate(GX = GX/32.768, GY = GY/32.768, GZ = GZ/32.768)
# group data before MV == 1
grp = 0
for( i in 1:dim(gyro_data)[1]) {
  gyro_data$GRP[i]=grp
  if(gyro_data$MV[i] == 1) {
    grp = grp+1
  }
} 
gyro_data_means <- gyro_data %>% 
  group_by(GRP) %>%
  summarise(across(c(GX,GY,GZ, MV), list(MEAN = mean, SUM = sum))) %>%
  mutate(NUM_SMP = 1/MV_MEAN) %>%
  filter(MV_SUM == 1) %>%
  select(GX_MEAN, GY_MEAN, GZ_MEAN)

ts_data <- imu.cal.data.gyro_rp %>% filter(MV == 1) %>% select(TS)

imu.cal.data.gyro_rp_amg <- cbind(ts_data, imu.cal.data.gyro_rp %>% filter(MV == 1) %>% select(-GX, -GY, -GZ, -TS), gyro_data_means) %>%
  rename(GX = GX_MEAN, GY = GY_MEAN, GZ = GZ_MEAN) %>%
  mutate(DT = (TS - shift(TS, 1, fill = TS[1], type = "lag"))/1000) %>%
  mutate(DGROLL = GX*DT,DGPITCH = GY*DT, DGYAW = GZ*DT) %>%
  mutate(GROLL = cumsum(DGROLL) + mean(AROLL[1:3000])) %>%
  mutate(GPITCH = cumsum(DGPITCH) + mean(APITCH[1:3000])) %>%
  mutate(GYAW = cumsum(DGYAW) + mean(AYAW[1:3000])) %>%
  mutate(
    AXG= -sin(GPITCH*toRad),
    AYG= sin(GROLL*toRad)*cos(GPITCH*toRad),
    AZG= cos(GROLL*toRad)*cos(GPITCH*toRad),
  ) %>%
  mutate(
    IGMX = MX * cos(GPITCH*toRad)   + MY*sin(GPITCH*toRad)*sin(GROLL)  + MZ * sin(GPITCH*toRad)*cos(GROLL*toRad),
    IGMY = 0                + MY*cos(GROLL*toRad)            - MZ*sin(GROLL*toRad),
    IGMZ = -MX*sin(GPITCH*toRad)    + MY*cos(GPITCH*toRad)*sin(GROLL*toRad)  + MZ*cos(GPITCH*toRad)*cos(GROLL*toRad), 
  ) %>% 
  mutate(
    IIGMX = IGMX * cos(GYAW*toRad) - IGMY*sin(GYAW*toRad),
    IIGMY = IGMX*sin(GYAW*toRad) + IGMY*cos(GYAW*toRad),
    IIGMZ = IGMZ,
  ) %>%
  mutate(GDEC = acos(IIGMX/sqrt(IIGMX^2+IIGMY^2+IIGMZ^2))*toDeg) %>%
  select(TS, MX, MY, MZ, AX, AY, AZ, AXG, AYG, AZG, GX, GY, GZ, DT, AROLL, APITCH, AYAW, GROLL, GPITCH, GYAW, DAROLL, DGROLL, DAPITCH, DGPITCH, DAYAW, DGYAW, DEC, GDEC)
  
plot(imu.cal.data.gyro_rp_amg$AX, type="l")
lines(imu.cal.data.gyro_rp_amg$AXG, col="red")

ylim_min <- min(min(imu.cal.data.gyro_rp_amg$GROLL), min(imu.cal.data.gyro_rp_amg$AROLL)) - 0.5
ylim_max <- max(max(imu.cal.data.gyro_rp_amg$GROLL), max(imu.cal.data.gyro_rp_amg$AROLL)) + 0.5
plot(imu.cal.data.gyro_rp_amg$AROLL, type="l", main = "AROLL (Black) vs GROLL (Red)", ylab = "Roll", ylim = c(ylim_min,ylim_max))
lines(imu.cal.data.gyro_rp_amg$GROLL, col="red")

ylim_min <- min(min(imu.cal.data.gyro_rp_amg$GPITCH), min(imu.cal.data.gyro_rp_amg$APITCH)) - 0.5
ylim_max <- max(max(imu.cal.data.gyro_rp_amg$GPITCH), max(imu.cal.data.gyro_rp_amg$APITCH)) + 0.5
plot(imu.cal.data.gyro_rp_amg$APITCH, type="l", main = "APITCH (Black) vs GPITCH (Red)", ylab = "Pitch", ylim = c(ylim_min,ylim_max))
lines(imu.cal.data.gyro_rp_amg$GPITCH, col="red")

ylim_min <- min(min(imu.cal.data.gyro_rp_amg$GYAW), min(imu.cal.data.gyro_rp_amg$AYAW)) - 0.5
ylim_max <- max(max(imu.cal.data.gyro_rp_amg$GYAW), max(imu.cal.data.gyro_rp_amg$AYAW)) + 0.5
plot(imu.cal.data.gyro_rp_amg$AYAW, type="l", main = "AYAW (Black) vs GYAW (Red)", ylab = "Yaw", ylim = c(ylim_min,ylim_max))
lines(imu.cal.data.gyro_rp_amg$GYAW, col="red")

ylim_min <- min(imu.cal.data.gyro_rp_amg$MX, imu.cal.data.gyro_rp_amg$AX, imu.cal.data.gyro_rp_amg$AXG) - 0.5
ylim_max <- max(imu.cal.data.gyro_rp_amg$MX, imu.cal.data.gyro_rp_amg$AX, imu.cal.data.gyro_rp_amg$AXG) + 0.5
plot(imu.cal.data.gyro_rp_amg$MX, type="l", main = "MX (Black) vs AX (Red) vs AXG (Blue)", ylab = "MX,AX,AXG", ylim = c(ylim_min,ylim_max))
lines(imu.cal.data.gyro_rp_amg$AX, col="red")
lines(imu.cal.data.gyro_rp_amg$AXG, col="blue")

ylim_min <- min(imu.cal.data.gyro_rp_amg$MY, imu.cal.data.gyro_rp_amg$AY, imu.cal.data.gyro_rp_amg$AYG) - 0.5
ylim_max <- max(imu.cal.data.gyro_rp_amg$MY, imu.cal.data.gyro_rp_amg$AY, imu.cal.data.gyro_rp_amg$AYG) + 0.5
plot(imu.cal.data.gyro_rp_amg$MY, type="l", main = "MY (Black) vs AY (Red) vs AYG (Blue)", ylab = "MY,AY,AYG", ylim = c(ylim_min,ylim_max))
lines(imu.cal.data.gyro_rp_amg$AY, col="red")
lines(imu.cal.data.gyro_rp_amg$AYG, col="blue")

ylim_min <- min(imu.cal.data.gyro_rp_amg$MZ, imu.cal.data.gyro_rp_amg$AZ, imu.cal.data.gyro_rp_amg$AZG) - 0.5
ylim_max <- max(imu.cal.data.gyro_rp_amg$MZ, imu.cal.data.gyro_rp_amg$AZ, imu.cal.data.gyro_rp_amg$AZG) + 0.5
plot(imu.cal.data.gyro_rp_amg$MZ, type="l", main = "MZ (Black) vs -AZ (Red) vs -AZG (Blue)", ylab = "MZ,-AZ,-AZG", ylim = c(ylim_min,ylim_max))
lines(-imu.cal.data.gyro_rp_amg$AZ, col="red")
lines(-imu.cal.data.gyro_rp_amg$AZG, col="blue")

plot(imu.cal.data.gyro_rp_amg$DEC - imu.cal.data.gyro_rp_amg$GDEC, type="l", main = "Declination Error", ylab = "ADEC - GDEC")
plot(imu.cal.data.gyro_rp_amg$AX - imu.cal.data.gyro_rp_amg$AXG, type="l", main = "AX Error", ylab = "AX - AXG")
plot(imu.cal.data.gyro_rp_amg$AY - imu.cal.data.gyro_rp_amg$AYG, type="l", main = "AY Error", ylab = "AY - AYG")
plot(imu.cal.data.gyro_rp_amg$AZ - imu.cal.data.gyro_rp_amg$AZG, type="l", main = "AZ Error", ylab = "AZ - AZG")

plot(imu.cal.data.gyro_rp_amg$AX - imu.cal.data.gyro_rp_amg$AXG, type="l", main = "AX Error", ylab = "AX - AXG")
plot(imu.cal.data.gyro_rp_amg$AROLL, imu.cal.data.gyro_rp_amg$GROLL, xlab = "Roll from Accelerometer", ylab = "Roll from Gyroscope", main = "Roll: Accelerometer vs Gyro")
plot(imu.cal.data.gyro_rp_amg$APITCH, imu.cal.data.gyro_rp_amg$GPITCH, xlab = "Pitch from Accelerometer", ylab = "Pitch from Gyroscope", main = "Pitch: Accelerometer vs Gyro", xlim = c(-180,180), ylim = c(-180,180))
plot(imu.cal.data.gyro_rp_amg$AYAW, imu.cal.data.gyro_rp_amg$GYAW, xlab = "Yaw from Magnetometer", ylab = "Yaw from Gyroscope", main = "Yaw: Magnetometer vs Gyro")
plot(imu.cal.data.gyro_rp_amg$DEC, imu.cal.data.gyro_rp_amg$GDEC, xlab = "Declination from Accelerometer", ylab = "Declination from Gyroscope", main = "Declination: Accelerometer vs Gyro")
plot(imu.cal.data.gyro_rp_amg$DEC, main = "Declination from Accelerometer", type = "l")

################################################################################################
### Roll, Pitch, Yaw from Accelerometer + Magnetometer and from Gyroscope + Magnetometer
################################################################################################
# mag: vector from magnetometer
# ng: (negative g o -g) vector from accelerometer or gyroscope
calcRPByNG <- function(mag, ng) {
  # calc IF axis coordinates in BF
  y <- cross(ng, mag)
  y <- y/norm(y, "2")
  x <- cross(y, ng)
  x <- x/norm(x, "2")
  z <- cross(x,y)
  z <- z/norm(z, "2")

  # Applico Roll e Pitch e calcolo yaw
  rpy <- calcRP(z)
  rpy[3] <- 0
  rpy[3] <- calcRP(eucMatrixRPY(rpy) %*% mag)[3]
  return(rpy)
}
normalize_vector <- function(v) {
  return(v/norm(v, "2"))
}
calcRPYByData <- function(data) {
  for(i in 1:dim(data)[1]) {
    currRow <- data[i,]
    mag_acc_rpy <- calcRPByNG(normalize_vector(as.matrix(c(currRow$MX, currRow$MY, currRow$MZ))), normalize_vector(as.matrix(c(currRow$AX, currRow$AY, currRow$AZ))))
    mag_gyro_rpy <- calcRPByNG(normalize_vector(as.matrix(c(currRow$MX, currRow$MY, currRow$MZ))), normalize_vector(as.matrix(c(currRow$AXG, currRow$AYG, currRow$AZG))))
    data[i,"MAROLL"] <- mag_acc_rpy[1]*toDeg
    data[i,"MAPITCH"] <- mag_acc_rpy[2]*toDeg
    data[i,"MAYAW"] <- mag_acc_rpy[3]*toDeg
    data[i,"MGROLL"] <- mag_gyro_rpy[1]*toDeg
    data[i,"MGPITCH"] <- mag_gyro_rpy[2]*toDeg
    data[i,"MGYAW"] <- mag_gyro_rpy[3]*toDeg
  }
  return(data)
}

imu.cal.data.gyro_4rpy <- calcRPYByData(imu.cal.data.gyro_rp_amg) %>% 
  select(TS, DT, MX, MY, MZ, AX, AY, AZ, AXG, AYG, AZG, AROLL, GROLL, MAROLL, MGROLL, APITCH, GPITCH, MAPITCH, MGPITCH, AYAW, GYAW, MAYAW, MGYAW, DEC, GDEC)

ylim_min <- min(imu.cal.data.gyro_4rpy$AROLL) - 0.5
ylim_max <- max(imu.cal.data.gyro_4rpy$AROLL) + 0.5
plot(imu.cal.data.gyro_4rpy$AROLL, type="l", main = "AROLL (Black) vs MGROLL (Red) vs GROLL (blue) ", ylim = c(ylim_min,ylim_max))
lines(imu.cal.data.gyro_4rpy$MGROLL, col="red")
lines(imu.cal.data.gyro_4rpy$GROLL, col="blue")

ylim_min <- min(imu.cal.data.gyro_4rpy$APITCH, imu.cal.data.gyro_4rpy$MGPITCH,imu.cal.data.gyro_4rpy$GPITCH) - 0.5
ylim_max <- max(imu.cal.data.gyro_4rpy$APITCH, imu.cal.data.gyro_4rpy$MGPITCH,imu.cal.data.gyro_4rpy$GPITCH) + 0.5
plot(imu.cal.data.gyro_4rpy$APITCH, type="l", main = "APITCH (Black) vs MGPITCH (Red) vs GPITCH (blue) ", ylim = c(ylim_min,ylim_max))
lines(imu.cal.data.gyro_4rpy$MGPITCH, col="red")
lines(imu.cal.data.gyro_4rpy$GPITCH, col="blue")

ylim_min <- min(imu.cal.data.gyro_4rpy$AYAW, imu.cal.data.gyro_4rpy$MGYAW,imu.cal.data.gyro_4rpy$GYAW) - 0.5
ylim_max <- max(imu.cal.data.gyro_4rpy$AYAW, imu.cal.data.gyro_4rpy$MGYAW,imu.cal.data.gyro_4rpy$GYAW) + 0.5
plot(imu.cal.data.gyro_4rpy$AYAW, type="l", main = "AYAW (Black) vs MGYAW (Red) vs GYAW (blue) ", ylim = c(ylim_min,ylim_max))
lines(imu.cal.data.gyro_4rpy$MGYAW, col="red")
lines(imu.cal.data.gyro_4rpy$GYAW, col="blue")


################################################################################################
### Experimental
################################################################################################
# 
# With Inertial Frame Axis:
#  X oriented to North
#  Y oriented to West
#  Z oriented to Zenith
# this is the angle between vectors of magnetic field mr <- c(x,0,z) and X axis versor c(1,0,0)
declination = 58.21*toRad

# estimation of roll, pitch from accel and yaw from mag in body frame
rpy <- c(0*toRad,0*toRad,0*toRad)

# real rpy (I force some roll,pitch,yaw difference)
rpy_real = rpy + c(0*toRad,0*toRad,-20*toRad)

# declination coordinate reference in inertial frame
mr_ref_if <- makeDecRefInIF(declination) 

# expected meausure from Mag in body frame
mr_bf_expected <- toBF(rpy, mr_ref_if) 
# real measures in body frame
mr_bf_real <- toBF(rpy_real,mr_ref_if)

# magnetic field reference in ipothetic inertial frame (if_tilde)
mr_ref_if_tilde <- toIF(rpy, mr_bf_real)

# magnetic field reference in ipothetic magnetic field frame (mff_tilde)
mr_ref_mff_tilde <- fromIFToMFF(declination, mr_ref_if_tilde)
# Roll, Pitch Error Calculation in mff_tilde
irpy_err <- calcRP(mr_ref_mff_tilde)
#irpy_err[3] <- 0
#irpy_err[3] <- calcYaw(mr_ref_mff_tilde)

# magnetic field reference in magnetic field frame (mff) => c(0,0,1)
mr_ref_mff <- eucMatrix(irpy_err[1], irpy_err[2], irpy_err[3]) %*% mr_ref_mff_tilde

# magnetic field reference in inertial frame obtained by correction
mr_ref_if_recalculated <- fromMFFToIF(declination, mr_ref_mff)

print("Error in magnetic field frame after Roll,Pitch correction")
t(round(mr_ref_mff - as.matrix(c(0,0,1)),10))
print("Error in magnetic field frame after Yaw correction")
t(round(mr_ref_if - mr_ref_if_recalculated,10))

# from mag reference in Inertial Frame to mag reference in Body Frame with error correction
# must be equal to mrr_real (real mag reference)
IFCorrectionMatrix <- fromMFFToIFMatrix(declination) %*% eucMatrix(irpy_err[1], irpy_err[2], irpy_err[3]) %*% fromIFToMFFMatrix(declination) 
IF_NEW_RPY_MATRIX <- IFCorrectionMatrix %*% toIFMatrix(rpy)

# Correzione in inertial frame
# OK! this can work.
mr_ref_if_recalculated_1 <- IF_NEW_RPY_MATRIX %*% mr_bf_real
print("Error from mag reference in Inertial Frame before correction")
t(round(mr_ref_if - mr_ref_if_tilde,10))
print("Error from mag reference in Body Frame after correction")
t(round(mr_bf_real - t(IF_NEW_RPY_MATRIX) %*% mr_ref_if, 10))

print("Magnetic Field RPY before correction")
t(round(calcRP(mr_ref_if_tilde)*toDeg,7))
print("Magnetic Field RPY after correction")
t(round(calcRP(mr_ref_if_recalculated_1)*toDeg,7))

#################################################################################
#### Some Test
#################################################################################
# Check with gravity
newRPY <- calcRP(t(IF_NEW_RPY_MATRIX) %*% as.matrix(c(0,0,1)))
newRPY[3] <- 0
print("Real Roll, Pitch, Yaw")
t(newRPY*toDeg)

# this is correct
v <- toBF(rpy_real, fromMFFToIF(declination, as.matrix(c(0,0,1))))
calcRP(round(IF_NEW_RPY_MATRIX %*% v,10))*toDeg

calcRP(round(IF_NEW_RPY_MATRIX %*% mr_bf_real,10))*toDeg
calcRP(toBF(rpy_real, as.matrix(c(0,0,1))))*toDeg

ng_bf <- t(IF_NEW_RPY_MATRIX) %*% as.matrix(c(0,0,1))
ng_bf_angle_mr <- acos(t(ng_bf)%*%mr_bf_real)*toDeg
ng_bf_angle_mr
ng_bf_angle_mr - declination*toDeg

mr_ref_if_recalculated <- IF_NEW_RPY_MATRIX %*% mr_bf_real

ng_bf_real <- toBF(rpy_real, as.matrix(c(0,0,1)))
acos(t(mr_bf_real)%*%ng_bf_real)*toDeg
calcRP(cross(mr_bf_real, ng_bf_real))*toDeg
calcRP(fromIFToMFF(declination, toIF(rpy_real, cross(mr_bf_real, ng_bf_real))))*toDeg
acos(t(ng_bf)%*%ng_bf_real)*toDeg

cross(ng_bf, ng_bf_real)
calcRP(cross(cross(mr_bf_real, ng_bf_real), cross(ng_bf, ng_bf_real)))*toDeg
acos(t(cross(mr_bf_real, ng_bf_real))%*%cross(ng_bf, ng_bf_real))*toDeg

calcRP(cross(ng_bf, ng_bf_real))*toDeg
calcRP(fromIFToMFF(declination, toIF(rpy_real, cross(ng_bf, ng_bf_real))))*toDeg
calcRP(cross(ng_bf, ng_bf_real))*toDeg

calcRP(cross(cross(mr_bf_real, ng_bf_real), cross(ng_bf, mr_bf_real)))*toDeg

calcRP(toIF(rpy_real, cross(ng_bf, mr_bf_real)))*toDeg
calcRP(fromIFToMFF(declination, toIF(rpy_real, cross(ng_bf, mr_bf_real))))*toDeg
calcRP(toIF(rpy_real, ng_bf))*toDeg
calcRP(fromIFToMFF(declination, toIF(rpy_real, ng_bf)))*toDeg

###########################################################################
#### Visualization
###########################################################################
if_base <- diag(3)
rownames(if_base) <- c("X", "Y", "Z")

bfm <- toIFMatrix(rpy_real)
bf_base <- rbind(t(bfm%*%c(1,0,0)),t(bfm%*%c(0,1,0)),t(bfm%*%c(0,0,1)))
rownames(bf_base) <- c("X", "Y", "Z")

mfm <- fromMFFToIFMatrix(declination)
mff_base <- rbind(t(mfm%*%c(1,0,0)), t(mfm%*%c(0,1,0)),t(mfm%*%c(0,0,1)))
rownames(mff_base) <- c("X", "Y", "Z")

open3d()
vectors3d(if_base, color=c(rep("black",3)), lwd=2, draw = TRUE)
vectors3d(bf_base, color=c(rep("blue",3)), lwd=2)
vectors3d(mff_base, color=c(rep("green",3)), lwd=2)

v <- rbind(t(toIF(rpy_real, mr_bf_real)), t(toIF(rpy_real, ng_bf)), t(mr_ref_if_recalculated), t(toIF(rpy, mr_bf_real)))
rownames(v) <- c("mr1", "ng_bf", "mr2", "mr3")
vectors3d(v, color=c("red","darkgreen","goldenrod3", "deepskyblue"))
planes3d(0, 0, 1, 0, col="gray", alpha=0.2)
highlevel() 
rgl.bringtotop()

##############################################################
## appunti matrici rotazione
##############################################################
v_rpy <- calcRP(mr_bf_real)
v_rpy*toDeg
v_rot <- eucMatrixRPY(v_rpy) # matrice di rotazione: BF -> MFF
v_rot %*% mr_bf_real # BF -> MFF => equals to c(0,0,1)
t(v_rot) %*% c(0,0,1) # MFF->BF => equal to mr_bf_real
mr_bf_real
calcRP(t(v_rot) %*% c(0,0,1))*toDeg
calcRP(mr_bf_real)*toDeg
v1 <- fromMFFToIF(declination, as.matrix(c(0,0,1)))
calcRP(toBF(calcRP(mr_bf_real)-calcRP(v1), as.matrix(c(0,0,1))))*toDeg
v2 <- cross(as.matrix(c(0,1,0)), v1)
v3 <- cross(v1, v2)
new_base <- cbind(v2, v3, v1) # this matrix is equal to fromMFFToIFMatrix(declination)

# Costruisco una base che abbia mr_bf_real come asse Z
v1 <- mr_bf_real
v2 <- cross(as.matrix(c(0,1,0)), v1)
v2 <- v2/norm(v2, "2")
v3 <- cross(v1, v2)
v3 <- v3/norm(v3, "2")
new_base <- cbind(v2, v3, v1) # from MFFToBF Matrix
Rc <- fromMFFToIFMatrix(declination) %*% t(new_base) %*% toBFMatrix(rpy) # this is the correction matrix from IF tilde to IF
new_toif_matrix <- Rc %*% toIFMatrix(rpy)
new_toif_matrix %*% mr_bf_real # this is equal to mr_ref_if
new_tobf_matrix <- t(new_toif_matrix)
new_toif_matrix %*% ng_bf_real # what is this?
new_ng_if <- new_toif_matrix %*% ng_bf_real
new_ng_ierr <- calcRP(new_ng_if)
new_ng_if_base <- eucMatrixRPY(new_ng_ierr) %*% new_toif_matrix 
new_ng_if_base %*% ng_bf_real

# Visualizzo i vettori
if_base <- diag(3)
rownames(if_base) <- c("X", "Y", "Z")

bfm <- toIFMatrix(rpy_real)
bfmt <- t(bfm)
bf_base <- rbind(t(bfm%*%c(1,0,0)),t(bfm%*%c(0,1,0)),t(bfm%*%c(0,0,1)))
bf_base_t <- rbind(t(bfmt%*%c(1,0,0)),t(bfmt%*%c(0,1,0)),t(bfmt%*%c(0,0,1)))
rownames(bf_base) <- c("X", "Y", "Z")

bf_to_mff <- t(new_base)
bf_to_mmf_base <- rbind(t(bf_to_mff%*%c(1,0,0)),t(bf_to_mff%*%c(0,1,0)),t(bf_to_mff%*%c(0,0,1)))
rownames(bf_to_mmf_base) <- c("X", "Y", "Z")

bf_to_if <- new_toif_matrix
bf_to_if_base <- rbind(t(bf_to_if%*%c(1,0,0)),t(bf_to_if%*%c(0,1,0)),t(bf_to_if%*%c(0,0,1)))
rownames(bf_to_if_base) <- c("X", "Y", "Z")

open3d()
vectors3d(if_base, color=c(rep("black",3)), lwd=2, draw = TRUE)
vectors3d(bf_base, color=c(rep("blue",3)), lwd=2)
vectors3d(bf_to_mmf_base, color=c(rep("darkgoldenrod1",3)), lwd=2)
vectors3d(bf_to_if_base, color=c(rep("green",3)), lwd=2)
if_base_plan <- t(if_base%*%c(0,0,1))
bf_base_plan <- t(bf_base_t%*%c(0,0,1))
bf_to_if_plan <- t(bf_to_if%*%c(0,0,1))
bf_to_mmf_plan <- t(bf_to_mff%*%c(0,0,1))

ng_bf_plan <- t(toIF(rpy_real, ng_bf))
my_planes <- rbind(if_base_plan, bf_base_plan,bf_to_if_plan,bf_to_mmf_plan, ng_bf_plan, t(mr_ref_if))

v <- rbind(t(toIF(rpy_real, mr_bf_real)), t(toIF(rpy_real, ng_bf)), t(mr_ref_if_recalculated), t(toIF(rpy, mr_bf_real)))
rownames(v) <- c("mr1", "ng_bf", "mr2", "mr3")
vectors3d(v, color=c("red","darkgreen","goldenrod3", "deepskyblue"))

#planes3d(my_planes[1,],d= 0, col="gray", alpha=0.2)  # IF
#planes3d(my_planes[2,],d= 0, col="blue", alpha=0.2) # BF
#planes3d(my_planes[3,],d= 0, col="green", alpha=0.2) # 
planes3d(my_planes[4,],d= 0, col="darkgoldenrod1", alpha=0.2)
#planes3d(my_planes[5,],d= 0, col="darkgreen", alpha=0.4)
planes3d(my_planes[6,],d= 0, col="magenta", alpha=0.4)

highlevel() 
rgl.bringtotop()

mr_x <- cross(mr_bf_real, mr_bf_expected)
mr_x <- mr_x/norm(mr_x, "2")
mr_x_rpy <- calcRP(mr_x)*toDeg
mr_y <- cross(mr_bf_expected, mr_x)
mr_y <- mr_y/norm(mr_y, "2")
mr_z <- mr_bf_expected
mr_base <- cbind(mr_x, mr_y, mr_z) # MR_BF -> BF
t(mr_base) %*% mr_bf_expected # c(0,0,1)
calcRP(t(mr_base) %*% mr_bf_real)*toDeg # only roll (pitch must be 0, and yaw??)

# TODO: cercare formula adeguata per il caso generale
bf_rpy_real <- as.matrix(c(
  -calcRP(mr_bf_expected)[1]*toDeg+calcRP(mr_bf_real)[1]*toDeg,
  calcRP(mr_bf_expected)[2]*toDeg-calcRP(mr_bf_real)[2]*toDeg,
  calcRP(mr_bf_expected)[3]*toDeg+calcRP(mr_bf_real)[3]*toDeg
))
bf_rpy_real

# Nuova Rotazione (IF -> BF)
bf_tilde_z <- toBFMatrix(rpy)%*%t(Rc)%*%as.matrix(c(0,0,1))
bf_tilde_y <- toBFMatrix(rpy)%*%t(Rc)%*%as.matrix(c(0,1,0))
bf_tilde_x <- toBFMatrix(rpy)%*%t(Rc)%*%as.matrix(c(1,0,0))
calcRP(bf_tilde_z)*toDeg

# Caso 1: rotazione yaw su un asse nel piano <X,Y>
v_z <- eucMatrix(0,0,-30*toRad) %*% as.matrix(c(1,0,0))
v_x <- as.matrix(c(0,0,1))
v_y <- cross(v_z, v_x)
v_y <- v_y/norm(v_y, "2")
v_to_if <- cbind(v_x, v_y, v_z)
v_to_if %*% as.matrix(c(1,0,0)) # x in VF is Z in IF => result = c(0,0,1)
v_to_if %*% as.matrix(c(0,0,1)) # z in VF is a -30Deg yaw of x in IF => result = v_z
v <- v_to_if %*% eucMatrix(0,0,90*toRad) %*% as.matrix(c(1,0,0)) # yaw rotation of X in VF then convert to IF 
v_rpy <- calcRP(v)
v_rpy*toDeg

# Using
#https://math.stackexchange.com/questions/1560039/closed-formula-to-transform-roll-pitch-yaw-angles-into-axis-angle-representation
k <- as.matrix(c(cos(v_rpy[3]/2)*cos(v_rpy[2]/2)*sin(v_rpy[1]/2)-sin(v_rpy[3]/2)*sin(v_rpy[2]/2)*cos(v_rpy[1]/2),
                 cos(v_rpy[3]/2)*sin(v_rpy[2]/2)*cos(v_rpy[1]/2)+sin(v_rpy[3]/2)*cos(v_rpy[2]/2)*sin(v_rpy[1]/2),
                 sin(v_rpy[3]/2)*cos(v_rpy[2]/2)*cos(v_rpy[1]/2)-cos(v_rpy[3]/2)*sin(v_rpy[2]/2)*sin(v_rpy[1]/2)
)) # this is the rotation axis
acos(prod(cos(v_rpy/2))+prod(sin(v_rpy/2)))*toDeg*2 # this is the angle of rotation
k <- k/norm(k,"2")
k
# TODO: visualizzarlo nel diagramma rispetto al rpy_real del BF

####################################################################################################
### Ricalcolo Roll, Pitch, Yaw da mag + accel e da mag + gyro
####################################################################################################
# ng_bf_real e mr_bf_real derivano dalle rilevazioni effettuate con i sensori
# cross(ng_bf_real, mr_bf_real) normalizzato è il vettore y=c(0,1,0) in IF espresso in coordinate BF
y <- cross(ng_bf_real, mr_bf_real)
y <- y/norm(y, "2")
toIF(rpy_real, y) # c(0,1,0)
# cross(y, ng_bf_real) normalizzato è il vettore x=c(1,0,0) in IF espresso in coordinate BF
x <- cross(y, ng_bf_real)
x <- x/norm(x, "2")
toIF(rpy_real, x)
# cross(x, y) normalizzato è il vettore z=c(0,0,1) in IF espresso in coordinate BF
z <- cross(x,y)
z <- z/norm(z, "2")
toIF(rpy_real, z)

# Applico Roll e Pitch e calcolo yaw
new_rpy <- calcRP(z)
new_rpy[3] <- 0
v <- eucMatrixRPY(new_rpy) %*% mr_bf_real
new_rpy[3] <- calcRP(v)[3]
new_rpy*toDeg

# Declination
acos(t(z)%*%mr_bf_real)*toDeg - 90

# Verifico il mag reference in BF
# deve corrispondere a mr_bf_real
stopifnot(round(toBF(new_rpy, mr_ref_if) - mr_bf_real) == as.matrix(c(0,0,0)))
toBF(new_rpy, mr_ref_if)
mr_bf_real

###############################################################################
### Simplified Attitude Determination Algorithm
###############################################################################
SADA <- imu.cal.data.gyro_4rpy %>%
  mutate(TDEC = -asin(MX*AX+MY*AY+MZ*AZ)*toDeg) %>%
  mutate(MD = MX*AX+MY*AY+MZ*AZ, MN = sqrt(1-MD^2)) %>%
  mutate(B11 = 0.5*MN*MX, B13 = 0.5*(AX+MD*MX), B21 = 0.5*MN*MY, B23 = 0.5*(AY+MD*MY), B31 = 0.5*MN*MZ, B33 = 0.5*(AZ+MD*MZ)) %>%
  mutate(tau = B13+B31) %>%
  mutate(a1 = B11 - B33 -1, Y23 = tau/a1, a2 = B21^2/a1-B11-B33-1, p1 = a1, p2 = a2) %>% # TODO: Verificare p1, p2
  mutate(Y11 = -1/a1, Y12 = B21/a1, Y13 = tau/a1, Y21 = -B21/(a1*a2), Y22 = 1/a2, Y23 = (B23+B21*Y13)/a2, a3 = p1 -2 + tau^2/p1+Y23^2*p2) %>%
  mutate(Y31 = (tau+B21*Y23)/(a1*a3), Y32 = Y23/a3, Y33 = 1/a3) %>%
  mutate(a = B23*(Y11+Y12*(Y21+Y23*Y31) + Y13*Y31) - (B13 - B31)*(Y21+Y23*Y31) - Y31*B21) %>%
  mutate(b = B23*(Y12*(Y22+Y23*Y32) + Y13*Y32) - (B13 - B31)*(Y22+Y23*Y32) - Y31*B21) %>%
  mutate(c = B23*(Y13*Y33+Y12*Y23*Y33) - Y33*B21 - Y23*Y33*(B13 - B31)) %>%
  mutate(qn = 1/sqrt(a^2+b^2+c^2+1), a = a*qn, b = b*qn, c = c*qn) %>%
  mutate(sada_roll = -atan2(2*a*(-1)+2*b*c, 1-2*(a^2-b^2))*toDeg) %>%
  mutate(sada_pitch = atan2(2*b*(-1)-2*a*c, 1-2*(b^2+c^2))*toDeg) %>%
  mutate(sada_yaw = atan2(2*((-1)*c+a*b), 1-2*(b^2+c^2))*toDeg)


ylim_min <- min(SADA$AROLL, SADA$sada_roll, SADA$GROLL) - 0.5
ylim_max <- max(SADA$AROLL, SADA$sada_roll, SADA$GROLL) + 0.5
plot(SADA$AROLL, type="l", main = "AROLL (Black) vs SADA_ROLL (Red) vs GROLL (Blue)", ylim = c(ylim_min,ylim_max))
lines(SADA$sada_roll, col="red")
lines(SADA$GROLL, col="blue")

ylim_min <- min(SADA$APITCH, SADA$sada_pitch, SADA$GPITCH) - 0.5
ylim_max <- max(SADA$APITCH, SADA$sada_pitch, SADA$GPITCH) + 0.5
plot(SADA$APITCH, type="l", main = "APITCH (Black) vs SADA_PITCH (Red) vs GPITCH (Blue)", ylim = c(ylim_min,ylim_max))
lines(SADA$sada_pitch, col="red")
lines(SADA$GPITCH, col="blue")

ylim_min <- min(SADA$AYAW, SADA$sada_yaw, SADA$GYAW) - 0.5
ylim_max <- max(SADA$AYAW, SADA$sada_yaw, SADA$GYAW) + 0.5
plot(SADA$AYAW, type="l", main = "AYAW (Black) vs SADA_YAW (Red) vs GYAW (Blue)", ylim = c(ylim_min,ylim_max))
lines(SADA$sada_yaw, col="red")
lines(SADA$GYAW, col="blue")

# from: https://answers.unity.com/questions/416169/finding-pitchrollyaw-from-quaternions.html
# from: https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
#float Pitch = Mathf.Rad2Deg * Mathf.Atan2(2 * q.x * q.w - 2 * q.y * q.z, 1 - 2 * q.x * q.x - 2 * q.z * q.z);
#float Yaw = Mathf.Rad2Deg * Mathf.Atan2(2 * q.y * q.w - 2 * q.x * q.z, 1 - 2 * q.y * q.y - 2 * q.z * q.z);
#float Roll = Mathf.Rad2Deg * Mathf.Asin(2 * q.x * q.y + 2 * q.z * q.w);
sr <- SADA[5846,]
w = -1
x = sr$a
y = sr$b
z = sr$c
roll <- atan2(2*x*w+2*y*z, 1-2*(x^2-y^2))*toDeg
roll
#pitch <- (2*atan2(sqrt(1+2*(w*x-z*y)), sqrt(1-2*(w*x-z*y))) - pi/2)*toDeg
pitch <- atan2(2*y*w-2*x*z, 1-2*y*y-2*z*z)*toDeg
pitch
yaw <- atan2(2*(w*z+x*y), 1 - 2*(y^2+z^2))*toDeg
yaw
