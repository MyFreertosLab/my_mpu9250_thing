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

imu.data.body <- imu.data.mag %>% select(MY, MX, MZ, MV) %>% filter(MV == 1) %>% select(MY, MX, MZ)

centerY = mean(imu.data.body$MX)
centerX = mean(imu.data.body$MY)
centerZ = -mean(imu.data.body$MZ)

imu.data.body <- imu.data.body %>% rename(MX = MY, MY = MX) %>% mutate(MZ = -MZ) %>% mutate(MX = (MX -centerX), MY = (MY - centerY), MZ = (MZ - centerZ))

#acc_mag_cos <-(imu.data.body$AX*imu.data.body$MX + imu.data.body$AY*imu.data.body$MY + imu.data.body$AZ*imu.data.body$MZ)/(sqrt(imu.data.body$AX**2 + imu.data.body$AY**2 + imu.data.body$AZ**2)*sqrt(imu.data.body$MX**2 + imu.data.body$MY**2 + imu.data.body$MZ**2))
#T = c(1:1:length(acc_mag_cos))
#acc_mag_angle <- 90 - acos(acc_mag_cos)/(2*pi)*360
#plot(T,acc_mag_angle)
#mean(acc_mag_angle)
#var(acc_mag_angle)
#sd(acc_mag_angle)


scatter3D(imu.data.body$MX, imu.data.body$MY, imu.data.body$MZ, colvar = imu.data.body$MZ, col = NULL, add = FALSE)
plotrgl()
rglwidget()

imu.data.mag.poly <- imu.data.body %>% select(MX, MY, MZ) %>% mutate(One = 1, MX2 = MX**2, MY2 = MY**2, MZ2=MZ**2, MXY=MX*MY, MXZ=MX*MZ, MYZ=MY*MZ)
fit <- lm(One ~ . - 1, imu.data.mag.poly)
coeff.hat <- coef(fit)
P <- as.matrix(imu.data.mag.poly %>% select(-One))

open3d()
plot3d(ellipse3d(fit, level = 0.90), col = "blue", alpha = 0.5, aspect = TRUE)

a=coeff.hat["MX2"]
b=coeff.hat["MXY"]/2
c=coeff.hat["MXZ"]/2
f=coeff.hat["MYZ"]/2
d=coeff.hat["MY2"]
e=coeff.hat["MZ2"]
A = matrix(c(a,b,c,b,d,f,c,f,e), nrow = 3, ncol=3)
B = matrix(c(coeff.hat["MX"],coeff.hat["MY"],coeff.hat["MZ"]))
X=as.matrix(imu.data.mag.poly %>% select(MX, MY, MZ))

# Calcola espressione algebrica
matrix(colSums(coeff.hat * t(imu.data.mag.poly %>% select(MX, MY, MZ, MX2, MY2, MZ2, MXY, MXZ, MYZ))))[1]
# Stesso calcolo in algebra lineare
t(matrix(X[1,1:3])) %*% A %*% matrix(X[1,1:3]) + t(B) %*% matrix(X[1,1:3])

# Provo il calcolo con tutti gli items
matrix(X, nrow=dim(X)[1], ncol = dim(X)[2])
# TODO: continuare
colSums(t(matrix(X,nrow = dim(X)[1], ncol = dim(X)[2]) * t(matrix(A %*% t(matrix(X, nrow = dim(X)[1], ncol = dim(X)[2])), nrow = dim(X)[2], ncol = dim(X)[1])))) + t(t(B) %*% t(matrix(X, dim(X)[1], dim(X)[2])))

C=-0.5*solve(A)%*%B
AE=1/(t(C)%*%A%*%C+1)[1,1]*A

# from Complete TriAxis Magnetometer Calibration
# Attenzione: 
#   per la definizione del symmetric kronecker vedere: https://it.mathworks.com/matlabcentral/answers/361452-symmetric-kronecker-product-in-matlab
#   0.5U%*%(kronecker(A,B) + kronecker(B,A))%*%t(U)
#   dove U è definita esplicitamente:
#   1 0 0 0 0 0 0 0 0
#   0 1/sqrt(2) 0 1/sqrt(2) 0 0 0 0 0 
#   0 0 1/sqrt(2) 0 0 0 1/sqrt(2) 0 0
#   0 0 0 0 1 0 0 0 0 
#   0 0 0 0 0 1/sqrt(2) 0 1/sqrt(2) 0
#   0 0 0 0 0 0 0 0 1
# Non sembra fornire il risultato corretto.
U=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 1/sqrt(2), 0, 1/sqrt(2), 0, 0, 0, 0, 0,
           0, 0, 1/sqrt(2), 0, 0, 0, 1/sqrt(2), 0, 0,
           0, 0, 0, 0, 1, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 1/sqrt(2), 0, 1/sqrt(2), 0,
           0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=6, ncol=9)

# Symmetric Kronecker Product
U=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 1, 0, 1, 0, 0, 0, 0, 0,
           0, 0, 1, 0, 0, 0, 1, 0, 0,
           0, 0, 0, 0, 1, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 1, 0, 1, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=6, ncol=9, byrow = TRUE)
U %*% kronecker(matrix(X[1,1:3]),matrix(X[1,1:3]))
xtAx=t(xk) %*% matrix(vechMat(A))
y=c(xk, X[1,1:3], 1)
psiOls=matrix(y) %*% t(matrix(y))

# T.B.D
Q=AE
k=t(C) %*% Q %*% C -1
u=-2*t(AE) %*% C
h=matrix(X[1,1:3])
t(h) %*% Q %*% (h) + t(u) %*% h + k


