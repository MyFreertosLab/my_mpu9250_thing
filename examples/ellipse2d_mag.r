library(tidyverse)
library(moderndive)
library(skimr)
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
library(quadprog)
library(expm)

# axis correction
imu.data.body.original <- imu.data.mag %>% select(MY, MX, MZ, MV) %>% filter(MV == 1) %>% select(MY, MX, MZ) %>% rename(MX = MY, MY = MX) %>% mutate(MZ = -MZ)
imu.data.body <- imu.data.body.original

# plot original data
scatter3D(imu.data.body$MX, imu.data.body$MY, imu.data.body$MZ, colvar = imu.data.body$MZ, col = NULL, add = FALSE)
plotrgl()
rglwidget()

# create model
imu.data.mag.poly <- imu.data.body %>% select(MX, MY, MZ) %>% mutate(One = -1, MX2 = MX**2, MY2 = MY**2, MZ2=MZ**2, MXY=MX*MY, MXZ=MX*MZ, MYZ=MY*MZ)
mag_model <- lm(One ~ . - 1, imu.data.mag.poly)
mag_model_coeff.hat <- coef(mag_model)
mag_model_data <- as.matrix(imu.data.body, ncol=3)

a=mag_model_coeff.hat["MX2"]
b=mag_model_coeff.hat["MXY"]/2
c=mag_model_coeff.hat["MXZ"]/2
e=mag_model_coeff.hat["MYZ"]/2
d=mag_model_coeff.hat["MY2"]
f=mag_model_coeff.hat["MZ2"]

mag_model_Q = xpndMat(c(a,b,c,d,e,f))
mag_model_u = matrix(c(mag_model_coeff.hat["MX"],mag_model_coeff.hat["MY"],mag_model_coeff.hat["MZ"]))
mag_model_k = 1
X=as.matrix(imu.data.mag.poly %>% select(MX, MY, MZ))
eigen(mag_model_Q)
D=matrix(c(eigen(mag_model_Q)$values[1], 0,0, 0,eigen(mag_model_Q)$values[2],0,0,0,eigen(mag_model_Q)$values[3]), nrow=3, ncol=3)
NV=eigen(mag_model_Q)$vectors/norm(eigen(mag_model_Q)$vectors, "F")
V=eigen(mag_model_Q)$vectors
TA=V%*%D%*%t(V)
stopifnot(abs(mag_model_Q - TA) < 1e-15)
TD=t(V)%*%mag_model_Q%*%V
stopifnot(abs(D - TD) < 1e-15)

# calc center

# from: https://math.stackexchange.com/questions/1217796/compute-center-axes-and-rotation-from-equation-of-ellipse
# Translation paragraph
# si risolve: D[1,1]*X^2 + D[2,2]*Y^2 + D[3,3]*Z^2+LB[1]*X+LB[2]*Y+LB[3]*Z+1=0
#             Le soluzioni forniscono il centro
# NOTA: si è posto F = 1, i coefficienti sono normalizzati rispetto al valore reale di F
#       Trovare un modo per stimare F
LB=t(mag_model_u)%*%V
Dmat <- 2*D
dvec <- -LB
Amat <- matrix(0,3,3)
bvec <- c(0,0,0)
res=solve.QP(Dmat,dvec, Amat)
S=res$unconstrained.solution
S%*%D%*%matrix(S) + LB%*%S
C=V%*%S
# center
C

# Translation
# plot before translation
open3d()
plot3d(imu.data.body, col = "blue")
spheres3d(c(0,0,0), radius = 227.9377, col = "red", alpha = 0.4)
# apply translation
imu.data.body <- imu.data.body %>% mutate(MX = MX - C[1],MY = MY - C[2],MZ = MZ - C[3])
imu.data.body.original <- imu.data.body.original %>% mutate(MX = MX - C[1],MY = MY - C[2],MZ = MZ - C[3])

# plot after translation
open3d()
plot3d(imu.data.body, col = "blue")
spheres3d(c(0,0,0), radius = 227.9377, col = "red", alpha = 0.4)

# Rotation
imu.data.body.rotated <- imu.data.body %>% mutate(multiplication= as.matrix(imu.data.body[,]) %*% V) %>% select(multiplication) %>% mutate(MX=multiplication[,1], MY=multiplication[,2], MZ=multiplication[,3]) %>% select(MX,MY,MZ)
open3d()
plot3d(imu.data.body.rotated, col = "blue")
spheres3d(c(0,0,0), radius = 227.9377, col = "red", alpha = 0.4)

scatter3D(imu.data.body.rotated$MX, imu.data.body.rotated$MY, imu.data.body.rotated$MZ, colvar = imu.data.body.rotated$MZ, col = NULL, add = FALSE)
plotrgl()
rglwidget()

# calc yaw pitch roll
pitch=asin(-NV[3,1])
roll=asin(NV[3,2]/cos(pitch))
yaw=asin(NV[2,1]/cos(pitch))
pitch/(2*pi)*360
roll/(2*pi)*360
yaw/(2*pi)*360

# calc axis lengths (from X^2/a^2+Y^2/b^2+Z^2/c^2=1) and scale factors
imu.data.mag.poly.rotated <- imu.data.body.rotated %>% select(MX, MY, MZ) %>% mutate(One = 1, MX2 = MX**2, MY2 = MY**2, MZ2=MZ**2) %>% select(One, MX2, MY2, MZ2)
fit_rotated <- lm(One ~ . - 1, imu.data.mag.poly.rotated)
coeff_rotated.hat <- coef(fit_rotated)
ellipsoid_axis = matrix(0,3,3)
c1=sqrt(1/coeff_rotated.hat[1])
c2=sqrt(1/coeff_rotated.hat[2])
c3=sqrt(1/coeff_rotated.hat[3])

diag(ellipsoid_axis) <- (1/((c1*c2*c3)^(2/3)))*c(c2*c3,c1*c3,c1*c2)
scale_factors_rotated = matrix(0,3,3)
sphere_radius = sqrt(t(matrix(diag(ellipsoid_axis))) %*% matrix(diag(ellipsoid_axis)))
diag(scale_factors_rotated) <- diag(ellipsoid_axis)/as.double(sphere_radius)
# plot sphere
imu.data.mag.sphere <- imu.data.body.rotated %>% mutate(multiplication= as.matrix(imu.data.body.rotated[,]) %*% scale_factors_rotated) %>% select(multiplication) %>% mutate(MX=multiplication[,1], MY=multiplication[,2], MZ=multiplication[,3]) %>% select(MX,MY,MZ)
open3d()
plot3d(imu.data.mag.sphere, col = "blue")
spheres3d(c(0,0,0), radius = sphere_radius, col = "red", alpha = 0.4)

# return to original space (from rotated)
scale_factors_1 = V %*% scale_factors_rotated %*% t(V)
#imu.data.mag.sphere <- imu.data.body %>% mutate(multiplication=  as.matrix(imu.data.body[,]) %*% scale_factors) %>% select(multiplication) %>% mutate(MX=multiplication[,1], MY=multiplication[,2], MZ=multiplication[,3]) %>% select(MX,MY,MZ)
imu.data.mag.sphere <- imu.data.mag.sphere %>% mutate(multiplication= as.matrix(imu.data.mag.sphere[,]) %*% t(V)) %>% select(multiplication) %>% mutate(MX=multiplication[,1], MY=multiplication[,2], MZ=multiplication[,3]) %>% select(MX,MY,MZ)
open3d()
plot3d(imu.data.mag.sphere, col = "blue")
spheres3d(c(0,0,0), radius = sphere_radius, col = "yellow", alpha = 0.4)

# FIXME: Why I need this?
# recalc model from centered data
imu.data.mag.poly <- imu.data.mag.sphere %>% select(MX, MY, MZ) %>% mutate(One = 1, MX2 = MX**2, MY2 = MY**2, MZ2=MZ**2) %>% select(One, MX2, MY2, MZ2)
fit_spheric <- lm(One ~ . - 1, imu.data.mag.poly)
coeff_spheric.hat <- coef(fit_spheric)

ellipsoid_axis = matrix(0,3,3)
c1=sqrt(1/coeff_spheric.hat[1])
c2=sqrt(1/coeff_spheric.hat[2])
c3=sqrt(1/coeff_spheric.hat[3])

diag(ellipsoid_axis) <- (1/((c1*c2*c3)^(2/3)))*c(c2*c3,c1*c3,c1*c2)
scale_factors_2 = matrix(0,3,3)
sphere_radius = sqrt(t(matrix(diag(ellipsoid_axis))) %*% matrix(diag(ellipsoid_axis)))
diag(scale_factors_2) <- diag(ellipsoid_axis)/as.double(sphere_radius)

# plot sphere
imu.data.mag.sphere <- imu.data.mag.sphere  %>% mutate(multiplication= as.matrix(imu.data.mag.sphere [,]) %*% scale_factors_2) %>% select(multiplication) %>% mutate(MX=multiplication[,1], MY=multiplication[,2], MZ=multiplication[,3]) %>% select(MX,MY,MZ)
open3d()
plot3d(imu.data.mag.sphere, col = "blue")
spheres3d(c(0,0,0), radius = sphere_radius, col = "red", alpha = 0.4)

scatter3D(imu.data.mag.sphere$MX, imu.data.mag.sphere$MY, imu.data.mag.sphere$MZ, colvar = imu.data.mag.sphere$MZ, col = NULL, add = FALSE)
plotrgl()
rglwidget()

# Check that it's a sphere
imu.data.mag.poly <- imu.data.mag.sphere %>% select(MX, MY, MZ) %>% mutate(One = 1, MX2 = MX**2, MY2 = MY**2, MZ2=MZ**2) %>% select(One, MX2, MY2, MZ2)
fit_spheric <- lm(One ~ . - 1, imu.data.mag.poly)
coeff_spheric2.hat <- coef(fit_spheric)
sphere_axis=sqrt(1/coeff_spheric2.hat )
stopifnot((sphere_axis[1] -  sphere_axis[2] < 1e-10 ))
stopifnot((sphere_axis[1] -  sphere_axis[3] < 1e-10 ))
scale_factors_1
scale_factors_2
scale_factors_final=scale_factors_2%*%scale_factors_1
sphere_radius = sqrt(t(scale_factors_final %*%as.matrix(sphere_axis)) %*% as.matrix(sphere_axis))
imu.data.mag.sphere.final <- imu.data.body  %>% mutate(multiplication= as.matrix(imu.data.body [,]) %*% scale_factors_final) %>% select(multiplication) %>% mutate(MX=multiplication[,1], MY=multiplication[,2], MZ=multiplication[,3]) %>% select(MX,MY,MZ)
open3d()
plot3d(imu.data.mag.sphere.final, col = "blue")
spheres3d(c(0,0,0), radius = sphere_radius, col = "red", alpha = 0.4)

scatter3D(imu.data.mag.sphere.final$MX, imu.data.mag.sphere.final$MY, imu.data.mag.sphere.final$MZ, colvar = imu.data.mag.sphere.final$MZ, col = NULL, add = FALSE)
plotrgl()
rglwidget()

#########################################################################################################################################
# from Complete Triaxis Magnetometer Caibration (is another implementation of the same problem)
#########################################################################################################################################
#center
mag_model_b = -1/2*solve(mag_model_Q)%*%mag_model_u
#matrix square root of Q
appo = mag_model_Q %*% t(mag_model_Q)
appo_eigen = eigen(appo)
mag_model_V = appo_eigen$vectors
mag_model_D = diag(eigen(mag_model_Q)$values)
mag_model_magnetic_norm = t(mag_model_b)%*%mag_model_Q%*%mag_model_b - mag_model_k 
mag_model_alpha = -as.double(4*mag_model_magnetic_norm/(4*mag_model_k - t((t(mag_model_V)%*%mag_model_u))%*%solve(mag_model_D)%*%(t(mag_model_V)%*%mag_model_u)))
mag_model_B = mag_model_V %*% sqrt(mag_model_alpha*mag_model_D) %*% t(mag_model_V) # the same of sqrtm(mag_model_Q)
mag_model_inv_A = mag_model_B 
mag_model_A = solve(mag_model_inv_A)
mag_model_data_transformed <- mag_model_data

#######################################################################
### T.B.D. Trovare qui i fattori di scala per normalizzare gli assi
### mag_model_magnetic_norm risuta: 3.076837. 
###    sqrt(mag_model_magnetic_norm) corrisponde al raggio medio dell'ellissoide
### devo normalizzare?
#######################################################################

##### T.B.D: L'algoritmo richiede un secondo passaggio per calcolare scale_factors_3 ed applicarlo alla trasformazione
#####        affinché sia realmente sferico
for(i in 1:dim(mag_model_data)[1]) {
  mag_model_data_transformed[i,] = (mag_model_inv_A %*% t(mag_model_data[i,] - t(mag_model_b)))
}

open3d()
points3d(mag_model_data_transformed, col = "blue")
spheres3d(c(0,0,0), radius = as.double(sqrt(mag_model_magnetic_norm)), col = "red", alpha = 0.4)
mag_model_data_transformed <- as.data.frame(mag_model_data_transformed)
scatter3D(mag_model_data_transformed$MX, mag_model_data_transformed$MY, mag_model_data_transformed$MZ, colvar = mag_model_data_transformed$MZ, col = NULL, add = FALSE)
plotrgl()
rglwidget()


# Check that it's a sphere
imu.data.mag.poly <- as.data.frame(mag_model_data_transformed) %>% select(MX, MY, MZ) %>% mutate(One = 1, MX2 = MX**2, MY2 = MY**2, MZ2=MZ**2) %>% select(One, MX2, MY2, MZ2)
fit_spheric <- lm(One ~ . - 1, imu.data.mag.poly)
coeff_spheric2.hat <- coef(fit_spheric)
sphere_axis=sqrt(1/coeff_spheric2.hat )
sphere_axis

## T.B.D: after the first estimation I need to calculate the scale_factors. Why? 
## T.B.D: the radius of sphere
ellipsoid_axis = matrix(0,3,3)
c1=sqrt(1/coeff_spheric2.hat[1])
c2=sqrt(1/coeff_spheric2.hat[2])
c3=sqrt(1/coeff_spheric2.hat[3])

diag(ellipsoid_axis) <- (1/((c1*c2*c3)^(2/3)))*c(c2*c3,c1*c3,c1*c2)
sphere_radius = sqrt(t(matrix(diag(ellipsoid_axis))) %*% matrix(diag(ellipsoid_axis)))
scale_factors_3 = matrix(0,3,3)
diag(scale_factors_3) <- diag(ellipsoid_axis)/as.double(sphere_radius)
scale_factors_3
#sphere_radius
mean(radius_data)
for(i in 1:dim(mag_model_data)[1]) {
  mag_model_data_transformed[i,] = scale_factors_3 %*% (mag_model_inv_A %*% t(mag_model_data[i,] - t(mag_model_b)))
}
open3d()
plot3d(mag_model_data_transformed, col = "blue")
spheres3d(c(0,0,0), radius = 1, col = "red", alpha = 0.4)
mean(mag_model_data_transformed$MX)
mean(mag_model_data_transformed$MY)
mean(mag_model_data_transformed$MZ)
