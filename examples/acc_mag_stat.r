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
imu.data.original <- read.csv('/hd/eclipse-cpp-2020-12/eclipse/workspace/my_mpu9250_thing/examples/imu-data-mag.csv')
# Axis North-East-Down
imu.data.NED <- imu.data.original %>% mutate(MX2 = MY, MY = -MX, MX = MX2, AY = -AY, AZ = -AZ, GY = -GY, GZ = -GZ) %>% select(-MX2)

imu.data.mag <- imu.data.NED %>% select(MX, MY, MZ, MV) %>% filter(MV == 1) %>% select(MX, MY, MZ)

#######################################
######## Filter Outlier of mag ########
#######################################
imu.data.mag.diff <- imu.data.mag %>%
  mutate(DX = (MX - shift(MX, 1, fill = MX[1], type = "lag"))) %>%
  mutate(DY = (MY - shift(MY, 1, fill = MY[1], type = "lag"))) %>%
  mutate(DZ = (MZ - shift(MZ, 1, fill = MZ[1], type = "lag")))
qx <- as.array(quantile(imu.data.mag.diff$DX))
qy <- as.array(quantile(imu.data.mag.diff$DY))
qz <- as.array(quantile(imu.data.mag.diff$DZ))
minx = 2*(2.5*qx[2]-1.5*qx[4])
maxx = 2*(2.5*qx[4]-1.5*qx[2])
miny = 2*(2.5*qy[2]-1.5*qy[4])
maxy = 2*(2.5*qy[4]-1.5*qy[2])
minz = 2*(2.5*qz[2]-1.5*qz[4])
maxz = 2*(2.5*qz[4]-1.5*qz[2])
imu.data.mag.diff.1 <- imu.data.mag.diff %>%
  mutate(EX = case_when(between(DX, minx,maxx) ~ 0, TRUE ~ 1)) %>%
  mutate(EX = (EX - shift(EX, 1, fill = 0, type = "lag"))) %>%
  mutate(EX = case_when(EX < 0 ~ 0, TRUE ~ EX)) %>%
  mutate(EY = case_when(between(DY, miny,maxy) ~ 0, TRUE ~ 1)) %>%
  mutate(EY = (EY - shift(EY, 1, fill = 0, type = "lag"))) %>%
  mutate(EY = case_when(EY < 0 ~ 0, TRUE ~ EY)) %>%
  mutate(EZ = case_when(between(DZ, minz,maxz) ~ 0, TRUE ~ 1)) %>%
  mutate(EZ = (EZ - shift(EZ, 1, fill = 0, type = "lag"))) %>%
  mutate(EZ = case_when(EZ < 0 ~ 0, TRUE ~ EZ))

imu.data.mag.diff.1 %>% filter(EX == 1)
imu.data.mag.diff.1 %>% filter(EY == 1)
imu.data.mag.diff.1 %>% filter(EZ == 1)

# Force previous measurement for outliers
imu.data.mag <- imu.data.mag.diff.1 %>%
  mutate(MX = case_when(EX == 0 ~ MX, TRUE ~ shift(MX, 1, fill = 0, type = "lag"))) %>%
  mutate(MY = case_when(EY == 0 ~ MY, TRUE ~ shift(MY, 1, fill = 0, type = "lag"))) %>%
  mutate(MZ = case_when(EZ == 0 ~ MZ, TRUE ~ shift(MZ, 1, fill = 0, type = "lag"))) %>%
  select(MX, MY, MZ)  

# plot original data
scatter3D(imu.data.mag$MX, imu.data.mag$MY, imu.data.mag$MZ, colvar = imu.data.mag$MZ, col = NULL, add = FALSE, ticktype = "detailed", scale = FALSE)
plotrgl()
rglwidget()

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
# Applicandola ad una matrice con diagonale = x,
# il risultato desiderato si ottiene con la seguente
# senza manipolazioni di x
make_3d_matrix_U <- function() {
  U=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 1, 0, 1, 0, 0, 0, 0, 0,
             0, 0, 1, 0, 0, 0, 1, 0, 0,
             0, 0, 0, 0, 1, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 1, 0, 1, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=6, ncol=9, byrow = TRUE)
  return(U)
}

symm_kronecker_x <- function(x,U) {
  n = sqrt(dim(U)[2])
  K <- matrix(0,n,n)
  diag(K) <- x
  as.matrix(diag((U %*% kronecker(K,K)%*%t(U))))
}

########################################################################
## Algorithms for adjusted least squares estimation
########################################################################
###############################
## Step 1): Form Tensor T 5xnxm
###############################
t0 <- function(x, variance) {
   1
}
t1 <- function(x, variance) {
  x
}
t2 <- function(x, variance) {
  x^2 - variance
}
t3 <- function(x, variance) {
  x^3 - 3*x*variance
}
t4 <- function(x, variance) {
  x^4 - 6*x^2*variance + 3*variance^2
}
T <- function(k, i, l, X, variance) {
    stopifnot(k >= 0 && k <= 4)
    res=switch(as.integer(k+1), t0(X[l,i], variance), t1(X[l,i], variance),t2(X[l,i], variance),t3(X[l,i], variance),t4(X[l,i], variance))
    as.double(res)
}
###############################
## Step 2): Define index matrix
###############################
make_index_matrix <- function(n) {
  M <- matrix(c(
    1, 1,
    1, 2,
    1, 3,
    2, 2,
    2, 3,
    3, 3,
    1, 0,
    2, 0,
    3, 0,
    0,0), nrow=10, ncol=2, byrow = TRUE
  )
  return(M)
}
#####################################
## Step 3): Form Tensor R 
#####################################
is_equal <- function(v1, v2) {
  as.integer(v1 == v2)
}
M <- make_index_matrix(3)
R <- function(p,q,i, M) {
  is_equal(M[p,1], i) + is_equal(M[p,2], i)  + is_equal(M[q,1], i)  + is_equal(M[q,2], i) 
}
#####################################
## Step 4): Compute ni_als
#####################################
ni_als <- function(p,q,X,variance,M) {
   m = dim(X)[1]
   n = dim(X)[2]
   s = 0
   for(l in 1:m) {
     f = 1
     for(i in 1:n) {
       f = f*T(R(p,q,i,M),i,l,X,variance)
     }
     s = s + f
   }
   s
}
#####################################
## Step 5): Define index off-diagonal
#####################################
index_off_diagonal <- function(n) {
  # only for n=3
  return(c(2,3,5))
}

#####################################
## Step 6): form matrix psi_als
#####################################
make_psi_als <- function(X,variance,M) {
  n = dim(X)[2]
  n1 = dim(M)[1]
  D = index_off_diagonal(n)
  psi_als = matrix(nrow = n1, ncol = n1)
  for(p in 1:n1) {
    for(q in p:n1) {
      if(is.element(p,D) && is.element(q,D)) {
        psi_als[p,q] = 4*ni_als(p,q,X,variance, M)
      } else if(!is.element(p,D) && !is.element(q,D)) {
        psi_als[p,q] = ni_als(p,q,X,variance, M)
      } else {
        psi_als[p,q] = 2*ni_als(p,q,X,variance, M)
      }
      psi_als[q,p] = psi_als[p,q]
    }
  }
  psi_als
}
################################################
## Step 7)8): Find eigenvector of min eigenvalue
################################################
calc_eigenvector_psi_als <- function(psi_als) {
  eigen_als = eigen(psi_als)
  eigen_values = eigen_als$values
  eigen_vectors = eigen_als$vectors 
  idx = which.min(abs(eigen_values))
  beta = matrix(eigen_vectors[, idx])
  result = beta/norm(beta, "2")
  return(result)
}
################################################
## Step 9): Estimate A, b, d
################################################
estimate_A <- function(beta, n) {
  xpndMat(beta[1:(n*(n+1)/2)])
}
estimate_b <- function(beta, n) {
  beta[(n*(n+1)/2+1):(dim(beta)[1]-1)]
}
estimate_d <- function(beta, n) {
  beta[dim(beta)[1]]
}
estimate_c <- function(A,b,n) {
  c=-0.5*solve(A)%*%b
}
estimate_Ae <- function(A,c,d) {
   as.double((1/(t(c) %*% A %*% c - d)))*A
}
to_definite_positive_factor <- function(A) {
  val <- eigen(A)$values
  result = 0
  if(prod(val >= 0) == 1) {
    result = 1
  } else if(prod(val < 0) == 1) {
    result = -1
  }
  return (result)
}
estimate_all <- function(psi_als, n) {
  beta = calc_eigenvector_psi_als(psi_als)
  A=estimate_A(beta, n)
  # force matrix to definite positive if it's negative
  factor <- to_definite_positive_factor(A)
  stopifnot(factor != 0)
  A = factor*A
  
  b=factor*estimate_b(beta, n)
  d=factor*estimate_d(beta, n)
  c=estimate_c(A, b, n)
  Ae=estimate_Ae(A,c,d)

  als = {}
  als$A=A
  als$b=b
  als$d=d
  als$c=c
  als$Ae=Ae
  als
}

mag_estimate <- function(als, mag_data) {
  mag_model_Q = als$A
  mag_model_u = als$b
  mag_model_k = als$d
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
  mag_model_B = mag_model_V %*% sqrt(mag_model_alpha*mag_model_D) %*% t(mag_model_V) # the same of sqrt(mag_model_Q)
  mag_model_inv_A = mag_model_B 
  mag_model_A = solve(mag_model_inv_A)
  f <- 1/sqrt(mag_model_magnetic_norm)
  mag_model_scale_factors <- matrix(0,3,3)
  diag(mag_model_scale_factors) <- c(f,f,f)
  
  result <- {}
  result$Q <- mag_model_Q
  result$b <- mag_model_b
  result$k <- mag_model_k
  result$offset <- mag_model_b
  result$V <- mag_model_V
  result$D <- mag_model_D
  result$B <- mag_model_B
  result$Hm2 <- mag_model_magnetic_norm
  result$alpha <- mag_model_alpha
  result$invA <- mag_model_inv_A
  result$data_source <- mag_data
  result$scale_factors <- mag_model_scale_factors
  result$als <- als
  return(result)
}

mag_apply_estimator <- function(mag_model) {
  mag_data <- mag_model$data_source
  mag_model_inv_A <- mag_model$invA
  mag_model_b <- mag_model$b
  mag_data_target <- mag_data
  # apply model to data
  for(i in 1:dim(mag_data)[1]) {
    x <- t(mag_data[i,] - t(mag_model_b))
    mag_data_target[i,] = mag_model$scale_factors %*% (mag_model_inv_A %*% x)
  }
  return(as.data.frame(mag_data_target))
}

mag_plot_data <- function(mag_data, sphere_radius = -1, title = "") {
  scatter3D(mag_data[,1], mag_data[,2], mag_data[,3], colvar = mag_data[,3], col = NULL, add = FALSE, scale = FALSE, ticktype = "detailed", main = title)
  plotrgl()
  rglwidget()
}

######################################################################################
# Magnetometer Estimation
######################################################################################
psi_als_mag = make_psi_als(as.matrix(imu.data.mag), 5.618457^2,make_index_matrix(3))
als_mag=estimate_all(psi_als_mag, 3)
imu.mag.estimator <- mag_estimate(als_mag, imu.data.mag)
imu.data.mag.estimated <- mag_apply_estimator(imu.mag.estimator)
# remove outlier
imu.data.mag.estimated <- imu.data.mag.estimated %>% mutate(RAD = sqrt(MX**2+MY**2+MZ**2)) 

# Plotting Magnetometer
mag_plot_data(imu.mag.estimator$data_source)
mag_plot_data(imu.data.mag.estimated)
sphere_radius_mag = mean(imu.data.mag.estimated$RAD)

scatter3D(imu.data.mag$MX-als_mag$c[1], imu.data.mag$MY-als_mag$c[2], imu.data.mag$MZ-als_mag$c[3], col = "green", add = FALSE, scale=FALSE)
scatter3D(imu.data.mag.estimated$MX*250, imu.data.mag.estimated$MY*250, imu.data.mag.estimated$MZ*250, col = "red", add = TRUE, scale=FALSE)
plotrgl()
rglwidget()

imu.data.mag.poly <- imu.data.mag.estimated %>% select(MX, MY, MZ) %>% mutate(One = 1, MX2 = MX**2, MY2 = MY**2, MZ2=MZ**2) %>% select(One, MX2, MY2, MZ2)
fit_spheric_mag <- lm(One ~ . - 1, imu.data.mag.poly)
coeff_spheric_mag.hat <- coef(fit_spheric_mag)
sphere_axis=sqrt(1/coeff_spheric_mag.hat )

print(c("mag sphere radius estimated: ", sphere_radius_mag), quote = FALSE)
print(c("mag sphere radius linear model: ", sphere_axis), quote = FALSE)


imu.data.acc <- imu.data.NED %>% filter(MV == 1) %>% select(AX, AY, AZ)

######################################################################################
# Accelerometer Estimation
######################################################################################
psi_als_acc = make_psi_als(as.matrix(imu.data.acc), 86.48833^2,make_index_matrix(3))
als_acc=estimate_all(psi_als_acc, 3)
imu.acc.estimator <- mag_estimate(als_acc, imu.data.acc)
imu.data.acc.estimated <- mag_apply_estimator(imu.acc.estimator)

# Plotting Accelerometer
mag_plot_data(imu.acc.estimator$data_source)
mag_plot_data(imu.data.acc.estimated)
sphere_radius_acc = sqrt(mean(imu.data.acc.estimated$AX^2+imu.data.acc.estimated$AY^2+imu.data.acc.estimated$AZ^2))

scatter3D(imu.data.acc$AX, imu.data.acc$AY, imu.data.acc$AZ, col = "blue", add = FALSE, ticktype = "detailed", scale=FALSE)
scatter3D(imu.data.acc$AX-als_acc$c[1], imu.data.acc$AY-als_acc$c[2], imu.data.acc$AZ-als_acc$c[3], col = "green", add = TRUE, scale=FALSE)
scatter3D(imu.data.acc.estimated$AX*4143, imu.data.acc.estimated$AY*4143, imu.data.acc.estimated$AZ*4143, col = "red", add = TRUE, scale=FALSE)

plotrgl()
rglwidget()

# Check that it's a sphere
imu.data.acc.poly <- imu.data.acc.estimated %>% select(AX, AY, AZ) %>% mutate(One = 1, AX2 = AX**2, AY2 = AY**2, AZ2=AZ**2) %>% select(One, AX2, AY2, AZ2)
fit_spheric_acc <- lm(One ~ . - 1, imu.data.acc.poly)
coeff_spheric_acc.hat <- coef(fit_spheric_acc)
sphere_axis=sqrt(1/coeff_spheric_acc.hat )

print(c("acc sphere radius estimated: ", sphere_radius_acc), quote = FALSE)
print(c("acc sphere radius linear model: ", sphere_axis), quote = FALSE)

scale_factors_3 = matrix(0,3,3)
diag(scale_factors_3) <- 1/sphere_axis

###########################################################################
### Gyroscope Bias
###########################################################################
gyro_bias <- as.matrix(c(
  mean(imu.data.NED$GX),
  mean(imu.data.NED$GY),
  mean(imu.data.NED$GZ)
))

###########################################################################
### Results
###########################################################################
print("Results: ", quote = F)
print("  Magnetometer Matrix: ")
print(imu.mag.estimator$scale_factors %*% imu.mag.estimator$invA)
print("  Magnetometer Center: ")
print(imu.mag.estimator$offset)
print("  Accelerometer Matrix: ")
print(imu.acc.estimator$scale_factors %*% imu.acc.estimator$invA)
print("  Accelerometer Center: ")
print(imu.acc.estimator$offset)
print("  Gyroscope Bias: ")
print(gyro_bias)
