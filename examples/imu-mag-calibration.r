library(tidyverse)
library(moderndive)
library(skimr)
library(ISLR)
library(infer)
library(rgl)
library(plot3D)
library(plot3Drgl)
options(rgl.printRglwidget = TRUE)
# mvmeta for symmetric matrix vectorization (vechMat)
library(mvmeta)
mag_data_correction <- function(mag_data) {
  # axis correction
  result <- mag_data %>% select(MY, MX, MZ, MV) %>% filter(MV == 1) %>% select(MY, MX, MZ) %>% rename(MX = MY, MY = MX) %>% mutate(MZ = -MZ)
  return(result)
}

mag_create_model <- function(mag_data) {
  # create model
  # t(x)%*%Q%*x+t(u)%*%x + k = 0  where k = 1
  imu.data.mag.poly <- mag_data %>% select(MX, MY, MZ) %>% mutate(One = -1, MX2 = MX**2, MY2 = MY**2, MZ2=MZ**2, MXY=MX*MY, MXZ=MX*MZ, MYZ=MY*MZ)
  mag_model <- lm(One ~ . - 1, imu.data.mag.poly)
  mag_model_coeff.hat <- coef(mag_model)

  a=mag_model_coeff.hat["MX2"]
  b=mag_model_coeff.hat["MXY"]/2
  c=mag_model_coeff.hat["MXZ"]/2
  e=mag_model_coeff.hat["MYZ"]/2
  d=mag_model_coeff.hat["MY2"]
  f=mag_model_coeff.hat["MZ2"]

  mag_model_Q = xpndMat(c(a,b,c,d,e,f))
  mag_model_u = matrix(c(mag_model_coeff.hat["MX"],mag_model_coeff.hat["MY"],mag_model_coeff.hat["MZ"]))
  mag_model_k = 1
  mag_model_b = -1/2*solve(mag_model_Q)%*%mag_model_u

  # Check result
  eigen(mag_model_Q)
  D=diag(eigen(mag_model_Q)$values)
  V=eigen(mag_model_Q)$vectors
  TA=V%*%D%*%t(V)
  stopifnot(abs(mag_model_Q - TA) < 1e-15)
  TD=t(V)%*%mag_model_Q%*%V
  stopifnot(abs(D - TD) < 1e-15)

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
  mag_data_source <- mag_data
  mag_data_target <- mag_data
  
  # apply model to data
  for(i in 1:dim(mag_model_data)[1]) {
    mag_data_target[i,] = (mag_model_inv_A %*% t(mag_data_source[i,] - t(mag_model_b)))
  }
  mag_data_target <- as.data.frame(mag_data_target)
  
#  open3d()
#  points3d(mag_data_target, col = "blue")
#  spheres3d(c(0,0,0), radius = as.double(sqrt(mag_model_magnetic_norm)), col = "red", alpha = 0.4)
#  scatter3D(mag_data_target$MX, mag_data_target$MY, mag_data_target$MZ, colvar = mag_data_target$MZ, col = NULL, add = FALSE)
#  plotrgl()
#  rglwidget()

  
  # Check that it's a sphere
  imu.data.mag.poly <- mag_data_target %>% select(MX, MY, MZ) %>% mutate(One = 1, MX2 = MX**2, MY2 = MY**2, MZ2=MZ**2) %>% select(One, MX2, MY2, MZ2)
  fit_spheric <- lm(One ~ . - 1, imu.data.mag.poly)
  coeff_spheric2.hat <- coef(fit_spheric)
  sphere_axis=sqrt(1/coeff_spheric2.hat )

  ellipsoid_axis = matrix(0,3,3)
  c1=sqrt(1/coeff_spheric2.hat[1])
  c2=sqrt(1/coeff_spheric2.hat[2])
  c3=sqrt(1/coeff_spheric2.hat[3])
  
  diag(ellipsoid_axis) <- (1/((c1*c2*c3)^(2/3)))*c(c2*c3,c1*c3,c1*c2)
  sphere_radius = sqrt(t(matrix(diag(ellipsoid_axis))) %*% matrix(diag(ellipsoid_axis)))
  scale_factors_3 = matrix(0,3,3)
  diag(scale_factors_3) <- diag(ellipsoid_axis)/as.double(sphere_radius)

  
  result <- {}
  result$Q <- mag_model_Q
  result$b <- mag_model_b
  result$k <- mag_model_k
  result$V <- mag_model_V
  result$D <- mag_model_D
  result$B <- mag_model_B
  result$Hm2 <- mag_model_magnetic_norm
  result$alpha <- mag_model_alpha
  result$invA <- mag_model_inv_A
  result$scale_factors <- scale_factors_3
  result$model <- mag_model
  result$data_source <- mag_data
  result$last_coeff <- coeff_spheric2.hat
  return(result)
}

mag_apply_model <- function(mag_model) {
  mag_data <- mag_model$data_source
  mag_model_inv_A <- mag_model$invA
  mag_model_b <- mag_model$b
  mag_data_target <- mag_data
  mag_model_scale_factors <- mag_model$scale_factors
  # apply model to data
  for(i in 1:dim(mag_data)[1]) {
    x <- t(mag_data[i,] - t(mag_model_b))
    
    mag_data_target[i,] = mag_model_scale_factors %*% (mag_model_inv_A %*% x)
  }
  return(as.data.frame(mag_data_target))
}

mag_plot_data <- function(mag_data, sphere_radius = -1) {
  scatter3D(mag_data$MX, mag_data$MY, mag_data$MZ, colvar = mag_data$MZ, col = NULL, add = FALSE)
  plotrgl()
  rglwidget()
}

mag_model <- mag_create_model(mag_data_correction(imu.data.mag))
mag_data_target <-mag_apply_model(mag_model)
mag_plot_data(mag_model$data_source)

open3d()
plot3d(mag_model$data_source, col = "blue")
spheres3d(c(0,0,0), radius = 227.9377, col = "red", alpha = 0.4)

mag_plot_data(mag_data_target, sphere_radius=1)

open3d()
plot3d(mag_data_target, col = "blue")
spheres3d(c(0,0,0), radius = 1, col = "red", alpha = 0.4)

# Check that it's a sphere
imu.data.mag.poly <- mag_data_target %>% select(MX, MY, MZ) %>% mutate(One = 1, MX2 = MX**2, MY2 = MY**2, MZ2=MZ**2) %>% select(One, MX2, MY2, MZ2)
fit_spheric <- lm(One ~ . - 1, imu.data.mag.poly)
coeff_spheric2.hat <- coef(fit_spheric)
sphere_axis=sqrt(1/coeff_spheric2.hat )
stopifnot(sphere_axis[1] != sphere_axis[2])
stopifnot(sphere_axis[1] != sphere_axis[3])

