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

imu.cal.data.gyro <- read.csv('/hd/eclipse-cpp-2020-12/eclipse/workspace/my_mpu9250_thing/examples/imu-cal-data-gyro-roll.csv')
imu.cal.data.gyro <- imu.cal.data.gyro[2:dim(imu.cal.data.gyro)[1],]

# plot original data
scatter3D(imu.cal.data.gyro$MX, imu.cal.data.gyro$MY, imu.cal.data.gyro$MZ, colvar = imu.cal.data.gyro$MZ, col = NULL, add = FALSE, ticktype = "detailed", scale = FALSE)
plotrgl()
rglwidget()

mag_plot_data <- function(mag_data, sphere_radius = -1, title = "") {
  scatter3D(mag_data[,1], mag_data[,2], mag_data[,3], colvar = mag_data[,3], col = NULL, add = FALSE, scale = FALSE, ticktype = "detailed", main = title)
  plotrgl()
  rglwidget()
}

############################################################################################################################
##### Check pitch, roll
############################################################################################################################
gyro_bias = c(mean(imu.cal.data.gyro$GX[1:10000]), mean(imu.cal.data.gyro$GY[1:10000]), mean(imu.cal.data.gyro$GZ[1:10000]))
imu.cal.data.gyro.cnt <- imu.cal.data.gyro %>% mutate(GX = GX - gyro_bias[1],GY = GY - gyro_bias[2],GZ = GZ - gyro_bias[3])

f = 180/pi
imu.cal.data.gyro_rp <- imu.cal.data.gyro.cnt %>% 
  mutate(NRMA = sqrt(AX^2+AY^2+AZ^2)) %>%
  mutate(NRMM = sqrt(MX^2+MY^2+MZ^2)) %>%
  mutate(MX = MX/NRMM, MY = MY/NRMM, MZ = MZ/NRMM) %>%
  mutate(AX = AX/NRMA, AY = AY/NRMA, AZ = AZ/NRMA) %>%
  mutate(PA = -asin(AX)*f) %>%
  mutate(RA = acos(round(AZ/cos(PA/f),10))*f) %>%
  mutate(RA = case_when(AY < 0 ~ -RA, TRUE ~ RA)) %>%
  mutate(RA = case_when(RA/f > pi ~ (RA/f-2*pi)*f,RA/f < -pi ~ (RA/f+2*pi)*f, TRUE ~ RA)) %>%
  mutate(
    ICMX = MX * cos(PA/f)   + MY*sin(PA/f)*sin(RA/f)  + MZ * sin(PA/f)*cos(RA/f),
    ICMY = 0                + MY*cos(RA/f)            - MZ*sin(RA/f),
    ICMZ = -MX*sin(PA/f)    + MY*cos(PA/f)*sin(RA/f)  + MZ*cos(PA/f)*cos(RA/f), 
    YM = atan2(-ICMY,ICMX)*f,
    ICAX = AX * cos(PA/f)   + AY*sin(PA/f)*sin(RA/f)  + AZ * sin(PA/f)*cos(RA/f),
    ICAY = 0                + AY*cos(RA/f)            - AZ*sin(RA/f),
    ICAZ = -AX*sin(PA/f)    + AY*cos(PA/f)*sin(RA/f)  + AZ*cos(PA/f)*cos(RA/f), 
  ) %>% 
  mutate(
    IICMX = ICMX * cos(YM/f) - ICMY*sin(YM/f),
    IICMY = ICMX*sin(YM/f) + ICMY*cos(YM/f),
    IICMZ = ICMZ,
    IICAX = ICAX * cos(YM/f) - ICAY*sin(YM/f),
    IICAY = ICAX*sin(YM/f) + ICAY*cos(YM/f),
    IICAZ = ICAZ
  ) %>%
  mutate(DEC = acos(IICMX/sqrt(IICMX^2+IICMY^2+IICMZ^2))*f) %>%
  rename(AROLL = RA, APITCH = PA, AYAW = YM, IMX = IICMX, IMY = IICMY, IMZ = ICMZ, IAX = IICAX, IAY = IICAY, IAZ = IICAZ) %>%
  mutate(DAROLL = AROLL - shift(AROLL, 1, fill = AROLL[1], type = "lag")) %>%
  mutate(DAPITCH = APITCH - shift(APITCH, 1, fill = APITCH[1], type = "lag")) %>%
  mutate(DAYAW = AYAW - shift(AYAW, 1, fill = AYAW[1], type = "lag")) %>%
  mutate(DAROLL = case_when(DAROLL > pi*f ~ 2*pi*f - DAROLL,  DAROLL < -pi*f ~  2*pi*f + DAROLL, TRUE ~ DAROLL)) %>%
  mutate(DPITCH = case_when(DAPITCH > pi*f ~ 2*pi*f - DAPITCH,  DAPITCH < -pi*f ~  2*pi*f + DAPITCH, TRUE ~ DAPITCH)) %>%
  mutate(DAYAW = case_when(DAYAW > pi*f ~ 2*pi*f - DAYAW,  DAYAW < -pi*f ~  2*pi*f + DAYAW, TRUE ~ DAYAW)) %>%
  select(TS, MX, MY, MZ, AX, AY, AZ, GX, GY, GZ, AROLL, APITCH, AYAW, DAROLL, DAPITCH, DAYAW, DEC, IMX, IMY, IMZ, IAX, IAY, IAZ, MV)
  

############################################################################################################################
##### Gyroscope calibration
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
    AXG= -sin(GPITCH/f),
    AYG= sin(GROLL/f)*cos(GPITCH/f),
    AZG= cos(GROLL/f)*cos(GPITCH/f),
  ) %>%
  mutate(
    IGMX = MX * cos(GPITCH/f)   + MY*sin(GPITCH/f)*sin(GROLL)  + MZ * sin(GPITCH/f)*cos(GROLL/f),
    IGMY = 0                + MY*cos(GROLL/f)            - MZ*sin(GROLL/f),
    IGMZ = -MX*sin(GPITCH/f)    + MY*cos(GPITCH/f)*sin(GROLL/f)  + MZ*cos(GPITCH/f)*cos(GROLL/f), 
  ) %>% 
  mutate(
    IIGMX = IGMX * cos(GYAW/f) - IGMY*sin(GYAW/f),
    IIGMY = IGMX*sin(GYAW/f) + IGMY*cos(GYAW/f),
    IIGMZ = IGMZ,
  ) %>%
  mutate(GDEC = acos(IIGMX)*f) %>%
  select(TS, MX, MY, MZ, AX, AY, AZ, AXG, AYG, AZG, GX, GY, GZ, DT, AROLL, APITCH, AYAW, GROLL, GPITCH, GYAW, DAROLL, DGROLL, DAPITCH, DGPITCH, DAYAW, DGYAW, DEC, GDEC)
  
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
plot(imu.cal.data.gyro_rp_amg$APITCH, imu.cal.data.gyro_rp_amg$GPITCH, xlab = "Pitch from Accelerometer", ylab = "Pitch from Gyroscope", main = "Pitch: Accelerometer vs Gyro", ylim = c(-60,80))
plot(imu.cal.data.gyro_rp_amg$AYAW, imu.cal.data.gyro_rp_amg$GYAW, xlab = "Yaw from Magnetometer", ylab = "Yaw from Gyroscope", main = "Yaw: Magnetometer vs Gyro")
plot(imu.cal.data.gyro_rp_amg$DEC, imu.cal.data.gyro_rp_amg$GDEC, xlab = "Declination from Accelerometer", ylab = "Declination from Gyroscope", main = "Declination: Accelerometer vs Gyro")
plot(imu.cal.data.gyro_rp_amg$DEC, main = "Declination from Accelerometer", type = "l")

##################################################################################################
#### Linear Regression: Pitch Accelerometer from Gyroscope
##################################################################################################
imu.cal.data.gyro_rp_amg_model_fit_mrx <- lm(APITCH ~ ., imu.cal.data.gyro_rp_amg %>% select(GPITCH, APITCH))
imu.cal.data.gyro_rp_amg_model_coeff_mrx <- imu.cal.data.gyro_rp_amg_model_fit_mrx$coefficients
P <- imu.cal.data.gyro_rp_amg %>% select(GPITCH, APITCH) %>% mutate(RAPITCH = imu.cal.data.gyro_rp_amg_model_coeff_mrx[1]+GPITCH*imu.cal.data.gyro_rp_amg_model_coeff_mrx[2], PITCH_ERR = APITCH - RAPITCH)

plot(imu.cal.data.gyro_rp_amg$APITCH, imu.cal.data.gyro_rp_amg$GPITCH, xlab = "Pitch from Accelerometer", ylab = "Pitch from Gyroscope", main = "Pitch: Linear Regression Accelerometer vs Gyro", ylim = c(-60,80))
lines(P$APITCH,P$RAPITCH, col="red")

ylim_min <- min(P$GPITCH, P$APITCH, P$RAPITCH) - 0.5
ylim_max <- max(P$GPITCH, P$APITCH, P$RAPITCH) + 0.5
plot(P$APITCH, type="l", main = "APITCH (Black) vs GPITCH (Red) vs RAPITCH (Blue)", ylab = "Pitch", ylim = c(ylim_min,ylim_max))
lines(P$GPITCH, col="red")
lines(P$RAPITCH, col="blue")

imu.cal.data.gyro_rp_amg_model_fit_mrx <- lm(AROLL ~ ., imu.cal.data.gyro_rp_amg %>% select(GROLL, AROLL))
imu.cal.data.gyro_rp_amg_model_coeff_mrx <- imu.cal.data.gyro_rp_amg_model_fit_mrx$coefficients
P <- imu.cal.data.gyro_rp_amg %>% select(GROLL, AROLL) %>% mutate(RAROLL = imu.cal.data.gyro_rp_amg_model_coeff_mrx[1]+GROLL*imu.cal.data.gyro_rp_amg_model_coeff_mrx[2], PITCH_ERR = AROLL - RAROLL)

plot(imu.cal.data.gyro_rp_amg$AROLL, imu.cal.data.gyro_rp_amg$GROLL, xlab = "Roll from Accelerometer", ylab = "Roll from Gyroscope", main = "Roll: Linear Regression Accelerometer vs Gyro", ylim = c(-60,80))
##lines(P$AROLL,P$RAROLL, col="red")
lines(seq(-100,100),imu.cal.data.gyro_rp_amg_model_coeff_mrx[1] + seq(-100,100)*imu.cal.data.gyro_rp_amg_model_coeff_mrx[2], col="red")

ylim_min <- min(P$GROLL, P$AROLL, P$RAROLL) - 0.5
ylim_max <- max(P$GROLL, P$AROLL, P$RAROLL) + 0.5
plot(P$AROLL, type="l", main = "AROLL (Black) vs GROLL (Red) vs RAROLL (Blue)", ylab = "Roll", ylim = c(ylim_min,ylim_max))
lines(P$GROLL, col="red")
lines(P$RAROLL, col="blue")

imu.cal.data.gyro_rp_amg_model_fit_mrx <- lm(AYAW ~ ., imu.cal.data.gyro_rp_amg %>% select(GYAW, AYAW))
imu.cal.data.gyro_rp_amg_model_coeff_mrx <- imu.cal.data.gyro_rp_amg_model_fit_mrx$coefficients
P <- imu.cal.data.gyro_rp_amg %>% select(GYAW, AYAW) %>% mutate(RAYAW = imu.cal.data.gyro_rp_amg_model_coeff_mrx[1]+GYAW*imu.cal.data.gyro_rp_amg_model_coeff_mrx[2], PITCH_ERR = AYAW - RAYAW)

plot(imu.cal.data.gyro_rp_amg$AYAW, imu.cal.data.gyro_rp_amg$GYAW, xlab = "YAW from Accelerometer", ylab = "YAW from Gyroscope", main = "YAW: Linear Regression Accelerometer vs Gyro", ylim = c(-60,80))
lines(P$AYAW,P$RAYAW, col="red")

ylim_min <- min(P$GYAW, P$AYAW, P$RAYAW) - 0.5
ylim_max <- max(P$GYAW, P$AYAW, P$RAYAW) + 0.5
plot(P$AYAW, type="l", main = "AYAW (Black) vs GYAW (Red) vs RAYAW (Blue)", ylab = "YAW", ylim = c(ylim_min,ylim_max))
lines(P$GYAW, col="red")
lines(P$RAYAW, col="blue")

##################################################################################################
#### Qualche relazione sui dati Mag, Accel, Gyro
##################################################################################################
# TODO: Calcolare Roll,Pitch,Yaw da magnetometro
# Provare a sfruttare la declinazione teorica (58.21deg) come reference

# In assenza di accelerazioni, u è l'asse di rotazione
u <- pracma::cross(c(imu.cal.data.gyro_rp_amg$AX[5001], imu.cal.data.gyro_rp_amg$AY[5001], imu.cal.data.gyro_rp_amg$AZ[5001]),c(imu.cal.data.gyro_rp_amg$MX[5001], imu.cal.data.gyro_rp_amg$MY[5001], imu.cal.data.gyro_rp_amg$MZ[5001]))
u <- u/norm(u, "2")
u

toRad <- pi/180
toDeg <- 1/toRad
eucMatrix <- function(roll, pitch, yaw) {
  RES <- matrix(c(
           cos(yaw)*cos(pitch), cos(yaw)*sin(pitch)*sin(roll)-sin(yaw)*cos(roll), cos(yaw)*sin(pitch)*cos(roll)+sin(yaw)*sin(roll),
           sin(yaw)*cos(pitch), sin(yaw)*sin(pitch)*sin(roll)+cos(yaw)*cos(roll), sin(yaw)*sin(pitch)*cos(roll)-cos(yaw)*sin(roll),
           -sin(pitch), cos(pitch)*sin(roll), cos(pitch)*cos(roll)
           ), nrow=3, ncol=3, byrow = TRUE)
  return(RES)
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
  return(eucMatrix(0,(pi/2+declination),0))
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
  return(round(fromMFFToIF(declination, matrix(c(0,0,1))),10))
}
# RPY reference is vector c(0,0,1) = -g
calcRPY <- function(v) {
  ipitch <- -asin(v[1])
  arg <- v[3]/cos(ipitch)
  if(abs(arg) <= 1) {
    iroll <- acos(arg)
  } else {
    iroll <- 0
  }
  if(v[2] < 0) {
    iroll = -iroll
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

toBFTest(as.matrix(c(0,0,90*toRad)), as.matrix(c(1,0,0)),as.matrix(c(0,-1,0)))
toIFTest(as.matrix(c(0,0,90*toRad)),as.matrix(c(0,-1,0)), as.matrix(c(1,0,0)))
toBFTest(as.matrix(c(0,90*toRad,0)), as.matrix(c(1,0,0)),as.matrix(c(0,0,1)))
toIFTest(as.matrix(c(0,90*toRad,0)),as.matrix(c(0,0,1)), as.matrix(c(1,0,0)))
toBFTest(as.matrix(c(0,90*toRad,90*toRad)), as.matrix(c(1,0,0)),as.matrix(c(0,-1,0)))
toIFTest(as.matrix(c(0,90*toRad,90*toRad)),as.matrix(c(0,-1,0)), as.matrix(c(1,0,0)))
toBFTest(as.matrix(c(90*toRad,0,0*toRad)), as.matrix(c(0,0,1)),as.matrix(c(0,1,0)))
toBFTest(as.matrix(c(10*toRad,0,0*toRad)), as.matrix(c(0,0,1)),as.matrix(c(0,0.1736482,0.9848078)))

# 
# With Inertial Frame Axis:
#  X oriented to North
#  Y oriented to West
#  Z oriented to Zenith
# this is the angle between vectors of magnetic field mr <- c(x,0,z) and X axis versor c(1,0,0)
declination = 58.21*toRad

# estimation of roll, pitch from accel and yaw from mag in body frame
rpy <- c(37*toRad,0*toRad,0*toRad)

# real rpy (I force some roll,pitch,yaw difference)
rpy_real = rpy + c(6*toRad,0*toRad,0*toRad)

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
irpy_err <- calcRPY(mr_ref_mff_tilde)
#irpy_err[3] <- 0

# magnetic field reference in magnetic field frame (mff)
mr_ref_mff <- eucMatrix(irpy_err[1], irpy_err[2], irpy_err[3]) %*% mr_ref_mff_tilde

# yaw calculation in Inertial Frame
mr_ref_if_to_correct_yaw <- fromMFFToIF(declination, mr_ref_mff)
iyaw_err <- calcRPY(mr_ref_if_to_correct_yaw)[3]

# magnetic field reference in inertial frame obtained by correction
mr_ref_if_recalculated <- eucMatrix(0,0,iyaw_err)%*% mr_ref_if_to_correct_yaw

print("Error in magnetic field frame after Roll,Pitch correction")
t(round(mr_ref_mff - as.matrix(c(0,0,1)),10))
print("Error in magnetic field frame after Yaw correction")
t(round(mr_ref_if - mr_ref_if_recalculated,10))

# from mag reference in Inertial Frame to mag reference in Body Frame with error correction
# must be equal to mrr_real (real mag reference)
IFCorrectionMatrix <- eucMatrix(0, 0, iyaw_err) %*% fromMFFToIFMatrix(declination) %*% eucMatrix(irpy_err[1], irpy_err[2], 0) %*% fromIFToMFFMatrix(declination) 
IF_NEW_RPY_MATRIX <- IFCorrectionMatrix %*% toIFMatrix(rpy)

# Correzione in inertial frame
# OK! this can work.
mr_ref_if_recalculated_1 <- IF_NEW_RPY_MATRIX %*% mr_bf_real
print("Error from mag reference in Inertial Frame before correction")
t(round(mr_ref_if - mr_ref_if_tilde,10))
print("Error from mag reference in Body Frame after correction")
t(round(mr_bf_real - t(IF_NEW_RPY_MATRIX) %*% mr_ref_if_recalculated_1, 10))

print("Magnetic Field RPY before correction")
t(round(calcRPY(mr_ref_if_tilde)*toDeg,7))
print("Magnetic Field RPY after correction")
t(round(calcRPY(mr_ref_if_recalculated_1)*toDeg,7))

#################################################################################
#### Some Test
# FIXME: with pitch = 0, there is some problem with mr_ref_mff_tilde
#> mr_ref_mff_tilde
#[,1]
#[1,] -0.002452976
#[2,] -0.088847589
#[3,]  0.996042212
#> -asin(-0.002452976)*toDeg
#[1] 0.1405453  <----------- must be 0!
#################################################################################
# Check with gravity
newRPY <- calcRPY(t(IF_NEW_RPY_MATRIX) %*% as.matrix(c(0,0,1)))
newRPY[3] <- iyaw_err + rpy[3]
print("Real Roll, Pitch, Yaw")
t(newRPY*toDeg)

