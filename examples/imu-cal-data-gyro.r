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
library(RSpincalc)

toRad <- pi/180
toDeg <- 1/toRad
accel_reference <- as.matrix(c(0,0,-1))
declination = -54.83*toRad*accel_reference[3]

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

calcRP <- function(v) {
  ipitch <- -asin(v[1])*accel_reference[3]
  
  arg <- min(1,max(-1,v[3]/cos(ipitch)))*accel_reference[3]
  iroll <- acos(arg)
  if(v[2]*accel_reference[3] < 0) {
    iroll <- -iroll
  }
  
  if(iroll > pi) {
    iroll = iroll - 2*pi
  } else if(iroll < -pi) {
    iroll = iroll + 2*pi
  }
  iyaw <- atan2(v[2], v[1])*accel_reference[3]
  result <- as.matrix(c(iroll, ipitch, iyaw))
  return(result)
}
calcYaw <- function(v) {
  iyaw <- atan2(v[2], v[1])*accel_reference[3]
  result <- as.matrix(c(0, 0, iyaw))
  return(result)
}

# m: magnetometer in body frame
# a: accelerometer in body frame
# mr: magnetometer reference in inertial frame 
# ar: accelerometer reference in inertial frame
sada_B_matrix <- function(m, a, mr, ar) {
  B <- 0.5*(a%*%t(ar) + m%*%t(mr))
  return(B)    
}
sada_z_vector <- function(B) {
  return(as.matrix(c(B[2,3]-B[3,2],B[3,1]-B[1,3],B[1,2]-B[2,1])))
}
sada_K_matrix <- function(B, z) {
  tr_B <- as.double(tr(B))
  K <- rbind(cbind(B + t(B)-tr_B*diag(1,3,3),z),cbind(t(z), tr_B))
  return(K)
}
sada_S_matrix <- function(K) {
  return (K - diag(1,4,4))
}

sada_L_parameters <- function(B) {
  a1 <- B[1,1]-B[3,3]-1
  a2 <- -B[2,1]^2/a1-B[1,1]-B[3,3]-1
  a3 <- -(B[2,3]-B[2,1]*(B[1,3]+B[3,1])/a1)^2/a2-(B[1,3]+B[3,1])^2/a1+B[3,3]-B[1,1]-1
  return(as.matrix(c(a1, a2, a3)))
}

sada_L1_matrix <- function(S, pv) {
  H1 <- diag(1,4,4)
  H1[4,1] <- -S[4,1]
  H2 <- diag(1,4,4)
  H2[3,1] <- -S[3,1]
  H3 <- diag(1,4,4)
  H3[2,1] <- -S[2,1]
  H4 <- diag(1,4,4)
  H4[1,1] <- 1/pv
  L1 <- H1 %*% H2 %*% H3 %*% H4
  return(L1)
}
sada_L2_matrix <- function(S, pv) {
  H1 <- diag(1,4,4)
  H1[4,2] <- -S[4,2]
  H2 <- diag(1,4,4)
  H2[3,2] <- -S[3,2]
  H3 <- diag(1,4,4)
  H3[1,2] <- -S[1,2]
  H4 <- diag(1,4,4)
  H4[2,2] <- 1/pv
  L2 <- H1 %*% H2 %*% H3 %*% H4
  return(L2)
}
sada_L3_matrix <- function(S, pv) {
  H1 <- diag(1,4,4)
  H1[4,3] <- -S[4,3]
  H2 <- diag(1,4,4)
  H2[2,3] <- -S[2,3]
  H3 <- diag(1,4,4)
  H3[1,3] <- -S[1,3]
  H4 <- diag(1,4,4)
  H4[3,3] <- 1/pv
  L3 <- H1 %*% H2 %*% H3 %*% H4
  return(L3)
}

sada_T_matrix <- function(K,B) {
  S <- sada_S_matrix(K)
  pars <- sada_L_parameters(B)
  a1 <- pars[1]
  a2 <- pars[2]
  a3 <- pars[3]

  L1 <- sada_L1_matrix(S,a1)
  S1 <- L1 %*% S
  
  L2 <- sada_L2_matrix(S1,a2)
  S2 <- L2 %*% S1
  
  L3 <- sada_L3_matrix(S2,a3)
  S3 <- L3 %*% S2
  
  TR <- L3 %*% L2 %*% L1 %*% S
  return(TR)
}

sada_quaternion <- function(B) {
  pars <- sada_L_parameters(B)
  tau <- B[1,3]+B[3,1]
  Y11 = 1/pars[1]
  Y22 = 1/pars[2]
  Y33 = 1/pars[3]
  Y12 = -B[2,1]/pars[1]
  Y13 = -tau/pars[1]
  Y21 = -B[2,1]/(pars[1]*pars[2])
  Y23 = -(B[2,1]*(-tau/pars[1])+B[2,3])/pars[2]
  Y31 = -tau/(pars[1]*pars[3]) - tau*B[2,1]^2/(pars[1]^2*pars[2]*pars[3])+B[2,3]*B[2,1]/(pars[1]*pars[2]*pars[3])
  Y32 = (tau*B[2,1]/pars[1]-B[2,3])/(pars[2]*pars[3])
  a <- B[2,3]*(Y11+Y12*(Y21+Y23*Y31)+Y13*Y31)-(B[1,3]-B[3,1])*(Y21+Y23*Y31)-Y31*B[2,1]
  b <- B[2,3]*(Y12*(Y22+Y23*Y32)+Y13*Y32)-(B[1,3]-B[3,1])*(Y22+Y23*Y32)-Y32*B[2,1]
  c <- B[2,3]*(Y13*Y33+Y12*Y23*Y33)-Y33*B[2,1]-Y23*Y33*(B[1,3]-B[3,1])
  n <- sqrt(a^2+b^2+c^2+1)
  q <- c(-1, a, b, c)/n
  return(q)
}

## Test 1
#mr <- as.matrix(c(cos(declination), 0, -sin(declination)))
#ar <- as.matrix(c(0,0,1))
#m <- as.matrix(c(0.96361408, 0.25812689, -0.06941483))
#a <- as.matrix(c(-0.77144974, -0.04735083, 0.63452597))

# North, East, Down
mr <- as.matrix(c(cos(declination), 0, sin(declination)))
ar <- as.matrix(c(0,0,-1))
m <- as.matrix(c(0.96361408, -0.25812689, 0.06941483))
a <- as.matrix(c(-0.77144974, 0.04735083, -0.63452597))
MD <- -t(m)%*%a
MN <- sqrt(1-(MD)^2)
#mr <- as.matrix(c(MN, 0, MD))
B <- sada_B_matrix(m, a, mr, ar)
z <- sada_z_vector(B)
K <- sada_K_matrix(B, z)
S <- sada_S_matrix(K)
pars <- sada_L_parameters(B)
L1 <- sada_L1_matrix(S, pars[1])
S1 <- L1 %*%S
S1
#L2 <- sada_L2_matrix(S1, -1.59635357)
L2 <- sada_L2_matrix(S1, pars[2])
S2 <- L2 %*%S1
S2
#L3 <- sada_L3_matrix(S2, -0.28029031)
L3 <- sada_L3_matrix(S2, pars[3])
S3 <- L3 %*%S2
S3
TR <- L3 %*% L2 %*% L1 %*% S
TR
x <- TR[1,4]
y <- TR[2,4]
z <- TR[3,4]
n <- sqrt(x^2+y^2+z^2+1)
x <- x/n
y <- y/n
z <- z/n
q1 <- t(as.matrix(c(1/n, x, y, z)))
Q2EA(q1,EulerOrder = "xyz")*toDeg # roll, pitch, yaw hanno segno inverso
Q2EA(Qconj(q1),EulerOrder = "zyx")*toDeg

## Test 2 North-East-Down & Quaternions
mr <- as.matrix(c(cos(declination), 0, sin(declination)))
ar <- as.matrix(c(0,0,-1))
rm <- toBFMatrix(as.matrix(c(21.4*toRad,14.38*toRad, 25.01*toRad)))
q <- DCM2Q(t(rm))
m <- rm %*% mr
a <- rm %*% ar
B <- sada_B_matrix(m, a, mr, ar)
z <- sada_z_vector(B)
K <- sada_K_matrix(B, z)
S <- sada_S_matrix(K)
pars <- sada_L_parameters(B)
L1 <- sada_L1_matrix(S, pars[1])
S1 <- L1 %*%S
S1
L2 <- sada_L2_matrix(S1, pars[2]) # Calcolare a2 correttamente (in sada_L_parameters)
S2 <- L2 %*%S1
S2
L3 <- sada_L3_matrix(S2, pars[3]) # Calcolare a3 correttamente (n sada_L_parameters)
S3 <- L3 %*%S2
S3
TR <- L3 %*% L2 %*% L1 %*% S
TR
q
x <- TR[1,4]
y <- TR[2,4]
z <- TR[3,4]
n <- sqrt(x^2+y^2+z^2+1)
x <- x/n
y <- y/n
z <- z/n
q1 <- t(as.matrix(c(1/n, x, y, z))) # x, y, z sono negati (oppure negare w)
Q2EA(q1,EulerOrder = "xyz")*toDeg # roll, pitch, yaw hanno segno inverso
Q2EA(Qconj(q1),EulerOrder = "zyx")*toDeg
Q2EA(q,EulerOrder = "xyz")*toDeg
Q2EA(Qconj(q),EulerOrder = "zyx")*toDeg # Il coniugato ha EulerOrder inverso
q1
q
stopifnot(round(q - q1, 8) == 0)

################################################################################
### Helpers for tests
################################################################################
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
  v <- toBF(rpy, accel_reference)
  new_rpy <- calcRP(v)
  new_rpy[3]<-0
  v1 <- toBF(new_rpy, accel_reference)
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
#imu.cal.data.gyro <- read.csv('/hd/eclipse-cpp-2020-12/eclipse/workspace/my_mpu9250_thing/examples/imu-cal-data-gyro-yaw.csv')

#imu.cal.data.gyro <- imu.cal.data.gyro[2:dim(imu.cal.data.gyro)[1],] %>% mutate(MY = -MY, MZ = -MZ, AY <- -AY, AZ = -AZ, GX = -GX, GY = -GY, GZ = -GZ)

# I dati forniti sono allineati agli assi dell'accelerometro.
# Li converto in NED
imu.cal.data.gyro <- imu.cal.data.gyro[2:dim(imu.cal.data.gyro)[1],] %>% mutate(MY = -MY, MZ = -MZ,AY = -AY, AZ = -AZ,GY = -GY, GZ = -GZ)

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
  mutate(PA = -asin(AX)*accel_reference[3]*toDeg) %>%
  mutate(RA = acos(round(AZ/cos(PA*toRad)*accel_reference[3],10))*toDeg) %>%
  mutate(RA = case_when(AY*accel_reference[3] < 0 ~ -RA, TRUE ~ RA)) %>%
  mutate(RA = case_when(RA*toRad > pi ~ (RA*toRad-2*pi)*toDeg,RA*toRad < -pi ~ (RA*toRad+2*pi)*toDeg, TRUE ~ RA)) %>%
  mutate(
    ICMX = MX * cos(PA*toRad)   + MY*sin(PA*toRad)*sin(RA*toRad)  + MZ * sin(PA*toRad)*cos(RA*toRad),
    ICMY = 0                + MY*cos(RA*toRad)            - MZ*sin(RA*toRad),
    ICMZ = -MX*sin(PA*toRad)    + MY*cos(PA*toRad)*sin(RA*toRad)  + MZ*cos(PA*toRad)*cos(RA*toRad), 
    YM = atan2(ICMY,ICMX)*accel_reference[3]*toDeg,
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
  select(TS, MX, MY, MZ, AX, AY, AZ, GX, GY, GZ, AROLL, APITCH, AYAW, DEC, IMX, IMY, IMZ, IAX, IAY, IAZ, MV)
  
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

df_gyro <- imu.cal.data.gyro_rp %>%
  mutate(DT = (TS - shift(TS, 1, fill = TS[1], type = "lag"))/1000) %>%
  mutate(GX = GX/32.768, GY = GY/32.768, GZ = GZ/32.768) %>%
  mutate(DGROLL = GX*DT,DGPITCH = GY*DT, DGYAW = GZ*DT) %>%
  mutate(GROLL = cumsum(DGROLL) + mean(AROLL[1:3000])) %>%
  mutate(GPITCH = cumsum(DGPITCH) + mean(APITCH[1:3000])) %>%
  mutate(GYAW = cumsum(DGYAW) + mean(AYAW[1:3000]))
  
ylim_min <- min(df_gyro$GROLL,df_gyro$AROLL) - 0.5
ylim_max <- max(df_gyro$GROLL,df_gyro$AROLL) + 0.5
plot(df_gyro$AROLL, type="l", main = "AROLL (Black) vs GROLL (Red)", ylab = "Roll", ylim = c(ylim_min,ylim_max))
lines(df_gyro$GROLL, col="red", )

ylim_min <- min(df_gyro$GPITCH,df_gyro$APITCH) - 0.5
ylim_max <- max(df_gyro$GPITCH,df_gyro$APITCH) + 0.5
plot(df_gyro$APITCH, type="l", main = "APITCH (Black) vs GPITCH (Red)", ylab = "Pitch", ylim = c(ylim_min,ylim_max))
lines(df_gyro$GPITCH, col="red", )

ylim_min <- min(df_gyro$GYAW,df_gyro$AYAW) - 0.5
ylim_max <- max(df_gyro$GYAW,df_gyro$AYAW) + 0.5
plot(df_gyro$AYAW, type="l", main = "AYAW (Black) vs GYAW (Red)", ylab = "Yaw", ylim = c(ylim_min,ylim_max))
lines(df_gyro$GYAW, col="red", )

imu.cal.data.gyro_rp_amg <- cbind(ts_data, imu.cal.data.gyro_rp %>% filter(MV == 1) %>% select(-GX, -GY, -GZ, -TS), gyro_data_means) %>%
  rename(GX = GX_MEAN, GY = GY_MEAN, GZ = GZ_MEAN) %>%
  mutate(DT = (TS - shift(TS, 1, fill = TS[1], type = "lag"))/1000) %>%
  mutate(DMX = (MX - shift(MX, 1, fill = MX[1], type = "lag"))) %>%
  mutate(DMY = (MY - shift(MY, 1, fill = MY[1], type = "lag"))) %>%
  mutate(DMZ = (MZ - shift(MZ, 1, fill = MZ[1], type = "lag")))
  

###############################################
### Simplified Attitude Determination Algorithm
###############################################
sada_update_df <- function(df) {
  for(row_idx in 1:dim(df)[1]) {
    curr_data <- df[row_idx,]
    m <- as.matrix(c(curr_data$MX, curr_data$MY,curr_data$MZ))
    a <- as.matrix(c(curr_data$AX, curr_data$AY,curr_data$AZ))
    m <- m/norm(m, "2")
    a <- a/norm(a, "2")
    ar <- as.matrix(c(0,0,-1))
    MD = -t(m)%*%a # negative becose accelerometer point to '-g'
    MN = sqrt(1-MD^2)
    mr <- as.matrix(c(MN, 0, MD))
    #mr <- as.matrix(c(cos(declination), 0, sin(declination)))
    
    B <- sada_B_matrix(m, a, mr, ar)
    q <- sada_quaternion(B)
    rpy <- (Q2EA.Xiao(q,EulerOrder = "xyz")*toDeg)
    df$sada_roll[row_idx] <- rpy[1]
    df$sada_pitch[row_idx] <- rpy[2]
    df$sada_yaw[row_idx] <- rpy[3]
    df$w[row_idx] <- q[1]
    df$a[row_idx] <- q[2]
    df$b[row_idx] <- q[3]
    df$c[row_idx] <- q[4]
    df$MD[row_idx] <- MD
    df$MN[row_idx] <- MN
    ng <- Q2DCM(t(q))%*%(accel_reference)
    df$sada_ngx[row_idx] <- ng[1]
    df$sada_ngy[row_idx] <- ng[2]
    df$sada_ngz[row_idx] <- ng[3]
  }
  return(df)
}

########################################
### Gyroscope Accelerometer Fusion (GAF)
########################################
gyroacc_fusion <- function(df, gamma) {
  qprev = NULL
  prev_data = NULL
  for(row_idx in 1:dim(df)[1]) {
    curr_data <- df[row_idx,]
    if(is.null(qprev)) {
      qprev <- as.matrix(c(curr_data$w, curr_data$a, curr_data$b, curr_data$c))
      prev_data <- curr_data
      q <- qprev
      ng <- Q2DCM(t(q))%*%(accel_reference)
      df$ngx[row_idx] <- ng[1]
      df$ngy[row_idx] <- ng[2]
      df$ngz[row_idx] <- ng[3]
      rpy <- Q2EA.Xiao(t(q), EulerOrder = "xyz")
      df$ng_roll[row_idx] = rpy[1]*toDeg
      df$ng_pitch[row_idx] = rpy[2]*toDeg
      df$ng_yaw[row_idx] = rpy[3]*toDeg
      next
    }    
    omega <- matrix(c(
      0,             -curr_data$GX,  -curr_data$GY, -curr_data$GZ,
      curr_data$GX,   0,             -curr_data$GZ,  -curr_data$GY,
      curr_data$GY,  curr_data$GZ,  0,             curr_data$GX,
      curr_data$GZ,  curr_data$GY,   -curr_data$GX, 0
    ), nrow = 4, ncol = 4, byrow = T)
    omega <- omega*toRad
    ds <- -as.matrix(c(
      curr_data$sada_ngx,
      curr_data$sada_ngy,
      curr_data$sada_ngz
    ))
    W <- matrix(c(
      ds[3],    -ds[2],  -ds[1],    0,
      ds[2],    -ds[3], 0,         ds[1],
      -ds[1],   0,      -ds[3],    -ds[2],
      0,        ds[1],  ds[2],     ds[3]
    ), nrow = 4, ncol = 4, byrow = T)
    sada_quat <- as.matrix(c(curr_data$w, curr_data$a, curr_data$b, curr_data$c))
    dt <- as.double(curr_data$DT)
    ####################################################################################
    # This is the original formula, but does not works. Yaw from accelerometer is zero.
    # q <- (1-gamma)*(dt/2*omega + diag(1, 4))%*%qprev + gamma*(W/2 + diag(1,4))%*%qprev
    # I use sada quaternion for gyroscope fusion
    ####################################################################################
    q <- (1-gamma)*(dt/2*omega + diag(1, 4))%*%qprev + gamma*sada_quat
    q <- q/norm(q, "2")
    qprev <- q
    prev_data <- curr_data
    ng <- Q2DCM(t(q))%*%(accel_reference)
    df$ngx[row_idx] <- ng[1]
    df$ngy[row_idx] <- ng[2]
    df$ngz[row_idx] <- ng[3]
    rpy <- Q2EA.Xiao(t(q), EulerOrder = "xyz")
    df$ng_roll[row_idx] = rpy[1]*toDeg
    df$ng_pitch[row_idx] = rpy[2]*toDeg
    df$ng_yaw[row_idx] = rpy[3]*toDeg
  }  
  return(df)
}

generate_df <- function() {
  mr <- as.matrix(c(cos(declination), 0, sin(declination)))
  ar <- as.matrix(c(0,0,-1))
  mx <- c()
  my <- c()
  mz <- c()
  ax <- c()
  ay <- c()
  az <- c()
  aroll <- c()
  apitch <- c()
  ayaw <- c()
  for(roll_deg in seq(-179,179,3)) {
    for(pitch_deg in seq(-89,89,3)) {
      ## Test 2 North-East-Down & Quaternions
      rm <- toBFMatrix(as.matrix(c(roll_deg*toRad,pitch_deg*toRad, 25.01*toRad)))
      m <- rm %*% mr
      a <- rm %*% ar
      mx = append(mx, m[1])
      my = append(my, m[2])
      mz = append(mz, m[3])
      ax = append(ax, a[1])
      ay = append(ay, a[2])
      az = append(az, a[3])
      aroll = append(aroll, roll_deg)
      apitch = append(apitch, pitch_deg)
      ayaw = append(ayaw, 25.01)
    }
  }
  df <- data.frame(mx, my, mz, ax, ay, az, aroll, apitch, ayaw)
  names(df) <- c('MX', 'MY', 'MZ', 'AX', 'AY', 'AZ', 'AROLL', 'APITCH', 'AYAW')
  return(df)
}
imu.cal.data.gyro_sada <- gyroacc_fusion(sada_update_df(imu.cal.data.gyro_rp_amg), 0.01)

SADA <- imu.cal.data.gyro_sada %>%
  mutate(TDEC = -asin(MX*AX+MY*AY+MZ*AZ)*toDeg) %>%
  mutate(DGROLL = GX*DT,DGPITCH = GY*DT, DGYAW = GZ*DT) %>%
  mutate(GROLL = cumsum(DGROLL) + mean(sada_roll[1:3000])) %>%
  mutate(GPITCH = cumsum(DGPITCH) + mean(sada_pitch[1:3000])) %>%
  mutate(GYAW = cumsum(DGYAW) + mean(sada_yaw[1:3000])) %>%
  select(TS, DT, MX, MY, MZ, AX, AY, AZ, GX, GY, GZ,ngx, ngy, ngz,sada_ngx, sada_ngy, sada_ngz, AROLL, APITCH, AYAW, GROLL, GPITCH, GYAW, sada_roll, sada_pitch, sada_yaw, ng_roll, ng_pitch, ng_yaw, DEC, TDEC, MN, MD, w, a, b, c)

ylim_min <- min(SADA$ng_roll, SADA$sada_roll, SADA$GROLL,SADA$AROLL) - 0.5
ylim_max <- max(SADA$ng_roll, SADA$sada_roll, SADA$GROLL,SADA$AROLL) + 0.5
plot(SADA$AROLL, type="l", main = "AROLL (Black) vs GROLL (Red) vs SADA_ROLL (blue)", ylab = "Roll", ylim = c(ylim_min,ylim_max))
lines(SADA$sada_roll, col="blue", )
lines(SADA$GROLL, col="red")
lines(SADA$ng_roll, col="green")

ylim_min <- min(SADA$ng_pitch, SADA$sada_pitch, SADA$GPITCH,SADA$APITCH) - 0.5
ylim_max <- max(SADA$ng_pitch, SADA$sada_pitch, SADA$GPITCH,SADA$APITCH) + 0.5
plot(SADA$APITCH, type="l", main = "APITCH (Black) vs GPITCH (Red) vs SADA_PITCH (blue)", ylab = "Pitch", ylim = c(ylim_min,ylim_max))
lines(SADA$sada_pitch, col="blue", )
lines(SADA$GPITCH, col="red")
lines(SADA$ng_pitch, col="green")

ylim_min <- min(SADA$ng_yaw, SADA$sada_yaw, SADA$GYAW,SADA$AYAW) - 0.5
ylim_max <- max(SADA$ng_yaw, SADA$sada_yaw, SADA$GYAW,SADA$AYAW) + 0.5
plot(SADA$AYAW, type="l", main = "AYAW (Black) vs GYAW (Red) vs SADA_YAW (blue)", ylab = "Yaw", ylim = c(ylim_min,ylim_max))
lines(SADA$sada_yaw, col="blue", )
lines(SADA$GYAW, col="red")
lines(SADA$ng_yaw, col="green")

ylim_min <- min(SADA$AX, SADA$sada_ngx) - 0.5
ylim_max <- max(SADA$AX, SADA$sada_ngx) + 0.5
plot(SADA$AX, type="l", main = "AX (Black) vs NGX (Red)", ylab = "AX vs. NGX", ylim = c(ylim_min,ylim_max))
lines(SADA$ngx, col="red", )

ylim_min <- min(SADA$AY, SADA$sada_ngy) - 0.5
ylim_max <- max(SADA$AY, SADA$sada_ngy) + 0.5
plot(SADA$AY, type="l", main = "AY (Black) vs NGY (Red)", ylab = "AY vs. NGY", ylim = c(ylim_min,ylim_max))
lines(SADA$ngy, col="red", )

ylim_min <- min(SADA$AZ, SADA$sada_ngz) - 0.5
ylim_max <- max(SADA$AZ, SADA$sada_ngz) + 0.5
plot(SADA$AZ, type="l", main = "AZ (Black) vs NGZ (Red)", ylab = "AZ vs. NGZ", ylim = c(ylim_min,ylim_max))
lines(SADA$ngz, col="red", )

mdec<- mean(SADA$TDEC)
vardec <- var(SADA$TDEC)
tdec_err <- (SADA$TDEC - mdec)
plot(tdec_err, type="l", main = paste("Declination Error (mean = ", round(mdec, 2), ")") , ylab = "TDEC Error")

plot(SADA$sada_roll, SADA$GROLL, xlab = "Roll from SADA", ylab = "Roll from Gyroscope", main = "Roll: SADA vs Gyro")
plot(SADA$sada_pitch, SADA$GPITCH, xlab = "Pitch from SADA", ylab = "Pitch from Gyroscope", main = "Pitch: SADA vs Gyro", xlim = c(-180,180), ylim = c(-180,180))
plot(SADA$sada_yaw, SADA$GYAW, xlab = "Yaw from SADA", ylab = "Yaw from Gyroscope", main = "Yaw: SADA vs Gyro")
plot(SADA$DEC, SADA$TDEC, xlab = "Declination from Accelerometer", ylab = "Declination from SADA", main = "Declination: Accelerometer vs SADA")

################################################################################################
### Experimental
################################################################################################
# 
# With Inertial Frame Axis:
#  X oriented to North
#  Y oriented to West
#  Z oriented to Zenith
# this is the angle between vectors of magnetic field mr <- c(x,0,z) and X axis versor c(1,0,0)

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

###########################################################################
#### Visualization
###########################################################################
if_base <- diag(3)
rownames(if_base) <- c("X", "Y", "Z")

bfm <- toIFMatrix(rpy_real)
bf_base <- rbind(t(bfm%*%c(1,0,0)),t(bfm%*%c(0,1,0)),t(bfm%*%c(0,0,1)))
rownames(bf_base) <- c("X", "Y", "Z")

open3d()
vectors3d(if_base, color=c(rep("black",3)), lwd=2, draw = TRUE)
vectors3d(bf_base, color=c(rep("blue",3)), lwd=2)

planes3d(0, 0, 1, 0, col="gray", alpha=0.2)
highlevel() 
rgl.bringtotop()

##############################################################
## appunti Visualizzazione Vettori
##############################################################
sr <- SADA[15000,]
w = -1
x = sr$a
y = sr$b
z = sr$c

#roll <- -atan2(2*x*(-1)+2*y*z, 1-2*(x^2-y^2))*toDeg
roll <- sr$sada_roll
roll
#pitch <- atan2(2*y*(-1)-2*x*z, 1-2*(y^2+z^2))*toDeg
pitch <- sr$sada_pitch
pitch
#yaw <- atan2(2*((-1)*z+x*y), 1-2*(y^2+z^2))*toDeg
yaw <- sr$sada_yaw
yaw
sada_rpy <- as.matrix(c(roll*toRad, pitch*toRad, yaw*toRad))
vq <- as.matrix(c(x, y, z))
vm <- as.matrix(c(sr$MX,sr$MY,sr$MZ))
va <- as.matrix(c(sr$AX,sr$AY,sr$AZ))


# Visualizzo i vettori
if_base <- diag(3)
rownames(if_base) <- c("X", "Y", "Z")

rm <- toBFMatrix(sada_rpy)
bf_base <- rbind(t(rm%*%c(1,0,0)),t(rm%*%c(0,1,0)),t(rm%*%c(0,0,1)))
rownames(bf_base) <- c("X", "Y", "Z")

open3d()
vectors3d(if_base, color=c(rep("black",3)), lwd=2, draw = TRUE)
vectors3d(bf_base, color=c(rep("blue",3)), lwd=2)
if_base_plan <- t(if_base%*%c(0,0,1))
bf_base_plan <- t(bf_base%*%c(0,0,1))

my_planes <- rbind(if_base_plan, bf_base_plan)

v <- rbind(t(va), t(toIF(sada_rpy, vm)), t(vq))
rownames(v) <- c("va", "vm", "vq")
vectors3d(v, color=c("red","blue","goldenrod3"))

planes3d(my_planes[1,],d= 0, col="gray", alpha=0.2)  # IF
planes3d(my_planes[2,],d= 0, col="blue", alpha=0.2) # BF

highlevel() 
rgl.bringtotop()

# from: https://answers.unity.com/questions/416169/finding-pitchrollyaw-from-quaternions.html
# from: https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
#float Pitch = Mathf.Rad2Deg * Mathf.Atan2(2 * q.x * q.w - 2 * q.y * q.z, 1 - 2 * q.x * q.x - 2 * q.z * q.z);
#float Yaw = Mathf.Rad2Deg * Mathf.Atan2(2 * q.y * q.w - 2 * q.x * q.z, 1 - 2 * q.y * q.y - 2 * q.z * q.z);
#float Roll = Mathf.Rad2Deg * Mathf.Asin(2 * q.x * q.y + 2 * q.z * q.w);

