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
imu.cal.data.gyro <- read.csv('/hd/eclipse-cpp-2020-12/eclipse/workspace/my_mpu9250_thing/examples/imu-cal-data-gyro.csv')
  
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
  mutate(RA = acos(min(1,max(-1, AZ/cos(PA/f))))*f) %>%
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
  mutate(GROLL = cumsum(DGROLL) + AROLL[1]) %>%
  mutate(GPITCH = cumsum(DGPITCH) + APITCH[1]) %>%
  mutate(GYAW = cumsum(DGYAW) + AYAW[1]) %>%
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
  

plot(imu.cal.data.gyro_rp_amg$AX, type="l")
lines(imu.cal.data.gyro_rp_amg$AXG, col="red")

plot(imu.cal.data.gyro_rp_amg$AY, type="l")
lines(imu.cal.data.gyro_rp_amg$AYG, col="red")

plot(imu.cal.data.gyro_rp_amg$AZ, type="l")
lines(imu.cal.data.gyro_rp_amg$AZG, col="red")

plot(imu.cal.data.gyro_rp_amg$AROLL, imu.cal.data.gyro_rp_amg$GROLL)
plot(imu.cal.data.gyro_rp_amg$APITCH, imu.cal.data.gyro_rp_amg$GPITCH)
plot(imu.cal.data.gyro_rp_amg$AYAW, imu.cal.data.gyro_rp_amg$GYAW)

plot(imu.cal.data.gyro_rp_amg$DEC, type="l")
lines(imu.cal.data.gyro_rp_amg$GDEC, col="red")


