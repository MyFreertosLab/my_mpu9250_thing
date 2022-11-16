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
imu.cal.data <- read.csv('/hd/eclipse-cpp-2020-12/eclipse/workspace/my_mpu9250_thing/imu-cal-data.csv') %>% filter(MV == 1) %>% mutate(RADM = sqrt(MX^2+MY^2+MZ^2)) %>% filter(RADM <= mean(RADM) + 0.1) %>% filter(RADM >= mean(RADM) - 0.1)  %>% mutate(RADA = sqrt(AX^2+AY^2+AZ^2)) %>% filter(RADA <= mean(RADA) + 0.3) %>% filter(RADA >= mean(RADA) - 0.3) %>% mutate(ANG = acos((MX*AX+MY*AY+MZ*AZ)/(sqrt(MX^2+MY^2+MZ^2)*(sqrt(AX^2+AY^2+AZ^2))))/pi*180)

mean(imu.cal.data$RADM)
sd(imu.cal.data$RADM)
median(imu.cal.data$RADM)

mean(imu.cal.data$RADA)
sd(imu.cal.data$RADA)
median(imu.cal.data$RADA)

mean(imu.cal.data$ANG)
sd(imu.cal.data$ANG)
median(imu.cal.data$ANG)

scatter3D(imu.cal.data$MX, imu.cal.data$MY, imu.cal.data$MZ, colvar = imu.cal.data$MZ, col = NULL, add = FALSE)
plotrgl()
rglwidget()

scatter3D(imu.cal.data$AX, imu.cal.data$AY, imu.cal.data$AZ, colvar = imu.cal.data$AZ, col = NULL, add = FALSE)
plotrgl()
rglwidget()

# aggiungere yaw alcolato da MX.MY.MZ
# calcolare rotazione assi da MX, MY, MZ (roll, pitch, yaw)
# confrontarlo con rotazione assi da AX, AY (roll, pitch)
imu.cal.data.angle <- imu.cal.data %>% mutate(ROLLA = asin(AY/RADA)/pi*180) %>% mutate(PITCHA = asin(-AX/RADA)/pi*180)
