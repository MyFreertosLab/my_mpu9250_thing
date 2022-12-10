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
imu.data.mag <- read.csv('/hd/eclipse-cpp-2020-12/eclipse/workspace/my_mpu9250_thing/examples/imu-data-mag.csv')
  
imu.data.body <- imu.data.mag %>% select(MY, MX, MZ, MV) %>% filter(MV == 1) %>% select(MY, MX, MZ)
imu.data.body <- imu.data.body %>% rename(MX = MY, MY = MX) %>% mutate(MZ = -MZ)

# plot original data
scatter3D(imu.data.body$MX, imu.data.body$MY, imu.data.body$MZ, colvar = imu.data.body$MZ, col = NULL, add = FALSE, ticktype = "detailed", scale = FALSE)
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
  mag_model_B = mag_model_V %*% sqrt(mag_model_alpha*mag_model_D) %*% t(mag_model_V) # the same of sqrtm(mag_model_Q)
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

# estimation Magnetometer
#psi_als_mag = make_psi_als(as.matrix(imu.data.body), 0.1720985^2,make_index_matrix(3))
#psi_als_mag = make_psi_als(as.matrix(imu.data.body), 4.213842^2,make_index_matrix(3))
psi_als_mag = make_psi_als(as.matrix(imu.data.body), 5.618457^2,make_index_matrix(3))
als_mag=estimate_all(psi_als_mag, 3)
prova <- mag_estimate(als_mag, imu.data.body)
prova_data <- mag_apply_estimator(prova)
# remove outlier
#prova_data <- prova_data %>% mutate(RAD = sqrt(MX**2+MY**2+MZ**2)) %>% filter(RAD < (mean(RAD) + 0.2)) %>% filter(RAD > (mean(RAD) - 0.2))
prova_data <- prova_data %>% mutate(RAD = sqrt(MX**2+MY**2+MZ**2)) 

# Plotting Magnetometer
mag_plot_data(prova$data_source)
mag_plot_data(prova_data)
sphere_radius_mag = mean(prova_data$RAD)

scatter3D(imu.data.body$MX-als_mag$c[1], imu.data.body$MY-als_mag$c[2], imu.data.body$MZ-als_mag$c[3], col = "green", add = FALSE, scale=FALSE)
scatter3D(prova_data$MX*250, prova_data$MY*250, prova_data$MZ*250, col = "red", add = TRUE, scale=FALSE)
plotrgl()
rglwidget()

print(c("sphere radius: ", sphere_radius_mag), quote = FALSE)

# calcola rotazione yaw, pitch, roll
V <- eigen(prova$invA)$vectors
v1 <- V[,1]
v2 <- V[,2]
v3 <- V[,3]
pitch <- asin(v1[3])
yaw <- acos(v1[1]/cos(pitch))
roll <- asin(v2[3]/cos(pitch))
print(c("yaw, pitch roll: ", yaw/pi*180, pitch/pi*180, roll/pi*180), quote = FALSE)

#imu.data.acc <- read.csv('/hd/eclipse-cpp-2020-12/eclipse/workspace/my_mpu9250_thing/examples/imu-data-acc.csv')
#imu.data.acc <- imu.data.acc %>% select(AX, AY, AZ) 
imu.data.acc <- imu.data.mag %>% filter(MV == 1) %>% select(AX, AY, AZ)

# estimation Accelerometer
##psi_als_acc = make_psi_als(as.matrix(imu.data.acc), 0.07497^2,make_index_matrix(3))
psi_als_acc = make_psi_als(as.matrix(imu.data.acc), 86.48833^2,make_index_matrix(3))
als_acc=estimate_all(psi_als_acc, 3)
prova_acc <- mag_estimate(als_acc, imu.data.acc)
prova_data_acc <- mag_apply_estimator(prova_acc)

# Plotting Accelerometer
mag_plot_data(prova_acc$data_source)
mag_plot_data(prova_data_acc)
sphere_radius_acc = sqrt(mean(prova_data_acc$AX^2+prova_data_acc$AY^2+prova_data_acc$AZ^2))

scatter3D(imu.data.acc$AX, imu.data.acc$AY, imu.data.acc$AZ, col = "blue", add = FALSE, ticktype = "detailed", scale=FALSE)
scatter3D(imu.data.acc$AX-als_acc$c[1], imu.data.acc$AY-als_acc$c[2], imu.data.acc$AZ-als_acc$c[3], col = "green", add = TRUE, scale=FALSE)
scatter3D(prova_data_acc$AX*4143, prova_data_acc$AY*4143, prova_data_acc$AZ*4143, col = "red", add = TRUE, scale=FALSE)

plotrgl()
rglwidget()

# Check that it's a sphere
imu.data.acc.poly <- prova_data_acc %>% select(AX, AY, AZ) %>% mutate(One = 1, AX2 = AX**2, AY2 = AY**2, AZ2=AZ**2) %>% select(One, AX2, AY2, AZ2)
fit_spheric_acc <- lm(One ~ . - 1, imu.data.acc.poly)
coeff_spheric_acc.hat <- coef(fit_spheric_acc)
sphere_axis=sqrt(1/coeff_spheric_acc.hat )

scale_factors_3 = matrix(0,3,3)
diag(scale_factors_3) <- 1/sphere_axis

############################################################################################################################
##### Check pitch, roll
############################################################################################################################
f = 180/pi
check_data_mag <- cbind(
  imu.data.body %>% mutate(MX = MX - prova$offset[1], MY = MY - prova$offset[2], MZ = MZ - prova$offset[3]) %>% select(MX, MY, MZ), 
  prova_data %>% rename(CMX = MX, CMY = MY, CMZ = MZ) %>% select(CMX, CMY, CMZ),
  imu.data.acc %>% mutate(AX = AX - prova_acc$offset[1], AY = AY - prova_acc$offset[2], AZ = AZ - prova_acc$offset[3]) %>% select(AX, AY, AZ), 
  prova_data_acc %>% rename(CAX = AX, CAY = AY, CAZ = AZ) %>% select(CAX, CAY, CAZ)
)
check_data_mag_rp <- check_data_mag %>% 
#  filter(abs(AX) > 0) %>%
#  filter(abs(AY) > 0) %>%
#  filter(abs(AZ) > 0) %>%
  mutate(NRMA = sqrt(AX^2+AY^2+AZ^2)) %>%
  mutate(CNRMA = sqrt(CAX^2+CAY^2+CAZ^2)) %>%
  mutate(NRMM = sqrt(MX^2+MY^2+MZ^2)) %>%
  mutate(CNRMM = sqrt(CMX^2+CMY^2+CMZ^2)) %>%
  mutate(CMX = CMX/CNRMM, CMY = CMY/CNRMM, CMZ = CMZ/CNRMM) %>%
  mutate(CAX = CAX/CNRMA, CAY = CAY/CNRMA, CAZ = CAZ/CNRMA) %>%
  mutate(PA = -asin(AX/NRMA)*f) %>%
  mutate(RA = acos(AZ/NRMA/cos(PA/f))*f) %>%
  mutate(RA = case_when(AY < 0 ~ -RA, TRUE ~ RA)) %>%
  mutate(RA = case_when(RA/f > pi ~ (RA/f-2*pi)*f,RA/f < -pi ~ (RA/f+2*pi)*f, TRUE ~ RA)) %>%
  mutate(CPA = -asin(CAX)*f) %>%
  mutate(CRA = acos(CAZ/cos(CPA/f))*f) %>%
  mutate(CRA = case_when(CAY < 0 ~ -CRA, TRUE ~ CRA)) %>%
  mutate(CRA = case_when(CRA/f > pi ~ (CRA/f-2*pi)*f,CRA/f < -pi ~ (CRA/f+2*pi)*f, TRUE ~ CRA)) %>%
  mutate(ERA = RA - CRA, EPA = PA - CPA) %>%
  ## FIXME: vedere meglio come coreggere i casi con segno inverso (es. 179 e -179 => ERA = +2 o -2)
  mutate(ERA = case_when(abs(ERA)/f > pi ~ RA + CRA,TRUE ~ ERA)) %>% 
  mutate(
    ICMX = CMX * cos(CPA/f) + CMY*sin(CPA/f)*sin(CRA/f) + CMZ * sin(CPA/f)*cos(CRA/f),
    ICMY = 0                + CMY*cos(CRA/f)            - CMZ*sin(CRA/f),
    ICMZ = -CMX*sin(CPA/f)  + CMY*cos(CPA/f)*sin(CRA/f) + CMZ*cos(CPA/f)*cos(CRA/f), 
    CYM = atan2(-ICMY,ICMX)*f,
    ICAX = CAX * cos(CPA/f) + CAY*sin(CPA/f)*sin(CRA/f) + CAZ * sin(CPA/f)*cos(CRA/f),
    ICAY = 0                + CAY*cos(CRA/f)            - CAZ*sin(CRA/f),
    ICAZ = -CAX*sin(CPA/f)  + CAY*cos(CPA/f)*sin(CRA/f) + CAZ*cos(CPA/f)*cos(CRA/f), 
  ) %>% 
  mutate(
    IICMX = ICMX * cos(CYM/f) - ICMY*sin(CYM/f),
    IICMY = ICMX*sin(CYM/f) + ICMY*cos(CYM/f),
    IICMZ = ICMZ,
    IICAX = ICAX * cos(CYM/f) - ICAY*sin(CYM/f),
    IICAY = ICAX*sin(CYM/f) + ICAY*cos(CYM/f),
    IICAZ = ICAZ
  ) %>%
  mutate(DEC = acos(IICMX/sqrt(IICMX^2+IICMY^2+IICMZ^2))*f)

plot(check_data_mag_rp$CRA, check_data_mag_rp$DEC)
plot(check_data_mag_rp$CPA, check_data_mag_rp$DEC)
plot(check_data_mag_rp$CYM, check_data_mag_rp$DEC)

############################################################################################################################
##### Prove filtro after calibration
############################################################################################################################
# angolo verticale magnetico = 58.21deg
# dati già forniti normalizzati ed in inertial frame
cal_data_rpy <- check_data_mag_rp %>% select(CPA, CRA, CYM, CMX, CMY, CMZ, CAX, CAY, CAZ, IICMX, IICMY, IICMZ, IICAX, IICAY, IICAZ) %>% rename(YM=CYM, PA=CPA, RA=CRA, MX=CMX, MY=CMY, MZ=CMZ, AX=CAX,AY=CAY, AZ=CAZ, IIMX=IICMX,IIMY=IICMY, IIMZ=IICMZ, IIAX=IICAX, IIAY=IICAY, IIAZ=IICAZ)
mrad <- 58.21/180*pi
mrx <- cos(mrad)
mrz <- -sin(mrad)
cal_data_rpy_model <- cal_data_rpy %>% mutate(MRX = mrx, MRY = 0, MRZ = mrz, ARX = 0, ARY = 0, ARZ = 1)
cal_data_rpy_model_fit_mrx <- lm(MRX ~ . - 1, cal_data_rpy_model %>% select(MRX, IIMX, IIMY, IIMZ, IIAX, IIAY, IIAZ))
cal_data_rpy_model_coeff_mrx <- cal_data_rpy_model_fit_mrx$coefficients

cal_data_rpy_model_fit_mry <- lm(MRY ~ . - 1, cal_data_rpy_model %>% select(MRY, IIMX, IIMY, IIMZ, IIAX, IIAY, IIAZ))
cal_data_rpy_model_coeff_mry <- cal_data_rpy_model_fit_mry$coefficients

cal_data_rpy_model_fit_mrz <- lm(MRZ ~ . - 1, cal_data_rpy_model %>% select(MRZ, IIMX, IIMY, IIMZ, IIAX, IIAY, IIAZ))
cal_data_rpy_model_coeff_mrz <- cal_data_rpy_model_fit_mrz$coefficients

cal_data_rpy_model_fit_arx <- lm(ARX ~ . - 1, cal_data_rpy_model %>% select(ARX, IIMX, IIMY, IIMZ, IIAX, IIAY, IIAZ))
cal_data_rpy_model_coeff_arx <- cal_data_rpy_model_fit_arx$coefficients

cal_data_rpy_model_fit_ary <- lm(ARY ~ . - 1, cal_data_rpy_model %>% select(ARY, IIMX, IIMY, IIMZ, IIAX, IIAY, IIAZ))
cal_data_rpy_model_coeff_ary <- cal_data_rpy_model_fit_ary$coefficients

cal_data_rpy_model_fit_arz <- lm(ARZ ~ . - 1, cal_data_rpy_model %>% select(ARZ, IIMX, IIMY, IIMZ, IIAX, IIAY, IIAZ))
cal_data_rpy_model_coeff_arz <- cal_data_rpy_model_fit_arz$coefficients

PMrx <- t(t(as.matrix(cal_data_rpy_model_coeff_mrx)) %*% t(as.matrix(cal_data_rpy_model %>% select(IIMX, IIMY, IIMZ, IIAX, IIAY, IIAZ))))
PMry <- t(t(as.matrix(cal_data_rpy_model_coeff_mry)) %*% t(as.matrix(cal_data_rpy_model %>% select(IIMX, IIMY, IIMZ, IIAX, IIAY, IIAZ))))
PMrz <- t(t(as.matrix(cal_data_rpy_model_coeff_mrz)) %*% t(as.matrix(cal_data_rpy_model %>% select(IIMX, IIMY, IIMZ, IIAX, IIAY, IIAZ))))

PArx <- t(t(as.matrix(cal_data_rpy_model_coeff_arx)) %*% t(as.matrix(cal_data_rpy_model %>% select(IIMX, IIMY, IIMZ, IIAX, IIAY, IIAZ))))
PAry <- t(t(as.matrix(cal_data_rpy_model_coeff_ary)) %*% t(as.matrix(cal_data_rpy_model %>% select(IIMX, IIMY, IIMZ, IIAX, IIAY, IIAZ))))
PArz <- t(t(as.matrix(cal_data_rpy_model_coeff_arz)) %*% t(as.matrix(cal_data_rpy_model %>% select(IIMX, IIMY, IIMZ, IIAX, IIAY, IIAZ))))

cal_data_rpy_corr <- cbind(cal_data_rpy, PMrx,PMry,PMrz,PArx,PAry,PArz)

plot(cal_data_rpy_corr$RA, acos(cal_data_rpy_corr$PMrx)*f, ylim = c(55,65))
plot(cal_data_rpy_corr$PA, acos(cal_data_rpy_corr$PMrx)*f, ylim = c(55,65))
plot(cal_data_rpy_corr$YM, acos(cal_data_rpy_corr$PMrx)*f, ylim = c(55,65))

# Rilevazioni corrette ruotate in body frame
cal_data_rpy <- cal_data_rpy_corr %>% 
  mutate(MDEC = acos(IIMX)*f) %>%
  mutate(CMDEC = acos(PMrx)*f) %>%
  mutate(EMDEC = CMDEC - MDEC) %>%
  mutate(
    YPMrx = cos(YM/f)*PMrx + sin(YM/f)*PMry,
    YPMry = -sin(YM/f)*PMrx + cos(YM/f)*PMry,
    YPMrz = PMrz,
    YPArx = cos(YM/f)*PArx + sin(YM/f)*PAry,
    YPAry = -sin(YM/f)*PArx + cos(YM/f)*PAry,
    YPArz = PArz
  ) %>%
  mutate(
   CMX = cos(PA/f)*YPMrx + 0                         - sin(PA/f)*YPMrz, 
   CMY = sin(PA/f)*sin(RA/f)*YPMrx + cos(RA/f)*YPMry + cos(PA/f)*sin(RA/f)*YPMrz,
   CMZ = sin(PA/f)*cos(RA/f)*YPMrx - sin(RA/f)*YPMry + cos(PA/f)*cos(RA/f)*YPMrz,
   CAX = cos(PA/f)*YPArx + 0                         - sin(PA/f)*YPArz, 
   CAY = sin(PA/f)*sin(RA/f)*YPArx + cos(RA/f)*YPAry + cos(PA/f)*sin(RA/f)*YPArz,
   CAZ = sin(PA/f)*cos(RA/f)*YPArx - sin(RA/f)*YPAry + cos(PA/f)*cos(RA/f)*YPArz,
  )

IMF <- matrix(c(
  cal_data_rpy_model_coeff_mrx, 
  cal_data_rpy_model_coeff_mry, 
  cal_data_rpy_model_coeff_mrz,
  cal_data_rpy_model_coeff_arx, 
  cal_data_rpy_model_coeff_ary, 
  cal_data_rpy_model_coeff_arz
), nrow = 6, ncol = 6, byrow = TRUE
)


# Ricalcolo la matrice dei coefficienti in body frame
cal_data_rpy_model_b_fit_cmx <- lm(CMX ~ . - 1, cal_data_rpy %>% select(CMX, MX, MY, MZ, AX, AY, AZ))
cal_data_rpy_model_b_coeff_cmx <- cal_data_rpy_model_b_fit_cmx$coefficients
cal_data_rpy_model_b_fit_cmy <- lm(CMY ~ . - 1, cal_data_rpy %>% select(CMY, MX, MY, MZ, AX, AY, AZ))
cal_data_rpy_model_b_coeff_cmy <- cal_data_rpy_model_b_fit_cmy$coefficients
cal_data_rpy_model_b_fit_cmz <- lm(CMZ ~ . - 1, cal_data_rpy %>% select(CMZ, MX, MY, MZ, AX, AY, AZ))
cal_data_rpy_model_b_coeff_cmz <- cal_data_rpy_model_b_fit_cmz$coefficients
cal_data_rpy_model_b_fit_cax <- lm(CAX ~ . - 1, cal_data_rpy %>% select(CAX, MX, MY, MZ, AX, AY, AZ))
cal_data_rpy_model_b_coeff_cax <- cal_data_rpy_model_b_fit_cax$coefficients
cal_data_rpy_model_b_fit_cay <- lm(CAY ~ . - 1, cal_data_rpy %>% select(CAY, MX, MY, MZ, AX, AY, AZ))
cal_data_rpy_model_b_coeff_cay <- cal_data_rpy_model_b_fit_cay$coefficients
cal_data_rpy_model_b_fit_caz <- lm(CAZ ~ . - 1, cal_data_rpy %>% select(CAZ, MX, MY, MZ, AX, AY, AZ))
cal_data_rpy_model_b_coeff_caz <- cal_data_rpy_model_b_fit_caz$coefficients

MF <- matrix(c(
  cal_data_rpy_model_b_coeff_cmx, 
  cal_data_rpy_model_b_coeff_cmy, 
  cal_data_rpy_model_b_coeff_cmz,
  cal_data_rpy_model_b_coeff_cax, 
  cal_data_rpy_model_b_coeff_cay, 
  cal_data_rpy_model_b_coeff_caz
  ), nrow = 6, ncol = 6, byrow = TRUE
)

mag_apply_filter <- function(mag_acc_data, mf) {
  # apply model to data
  mag_acc_data_result <- mag_acc_data
  for(i in 1:dim(mag_acc_data)[1]) {
    mag_acc_data_result[i,] = (mf %*% t(mag_acc_data[i,]))
  }
  return(as.data.frame(mag_acc_data_result))
}
calibrated_data_filtered <- mag_apply_filter(cal_data_rpy %>% select(MX, MY, MZ, AX, AY, AZ), MF)

#mag_plot_data(cal_data_rpy %>% select(MX, MY, MZ), title = "mag original")
#mag_plot_data(calibrated_data_filtered %>% select(MX, MY, MZ) , title = "mag filtered")
#mag_plot_data(cal_data_rpy %>% select(AX, AY, AZ), title = "acc original")
#mag_plot_data(calibrated_data_filtered %>% select(AX, AY, AZ), title = "acc filtered" )


############################################################################################################################
##### Gyroscope calibration
############################################################################################################################

inertial_frame_data <- calibrated_data_filtered %>% 
  select(MX, MY, MZ, AX, AY, AZ) %>% 
  mutate(PITCH = -asin(AX)*f) %>%
  mutate(ROLL = acos(AZ/cos(PITCH/f))*f) %>%
  mutate(ROLL = case_when(AY < 0 ~ -ROLL, TRUE ~ ROLL)) %>%
  mutate(ROLL = case_when(ROLL/f > pi ~ (ROLL/f-2*pi)*f,ROLL/f < -pi ~ (ROLL/f+2*pi)*f, TRUE ~ ROLL)) %>%
  mutate(
    # calculate Yaw
    ICMX = MX * cos(PITCH/f) + MY*sin(PITCH/f)*sin(ROLL/f) + MZ * sin(PITCH/f)*cos(ROLL/f),
    ICMY = 0                + MY*cos(ROLL/f)            - MZ*sin(ROLL/f),
    YAW = atan2(-ICMY,ICMX)*f
  ) %>% 
  select(-ICMX, -ICMY) %>%
  mutate(DROLL = ROLL - shift(ROLL, 1, fill = ROLL[1], type = "lag")) %>%
  mutate(DPITCH = PITCH - shift(PITCH, 1, fill = PITCH[1], type = "lag")) %>%
  mutate(DYAW = YAW - shift(YAW, 1, fill = YAW[1], type = "lag")) %>%
  mutate(DROLL = case_when(DROLL > pi*f ~ 2*pi*f - DROLL,  DROLL < -pi*f ~  2*pi*f + DROLL, TRUE ~ DROLL)) %>%
  mutate(DPITCH = case_when(DPITCH > pi*f ~ 2*pi*f - DPITCH,  DPITCH < -pi*f ~  2*pi*f + DPITCH, TRUE ~ DPITCH)) %>%
  mutate(DYAW = case_when(DYAW > pi*f ~ 2*pi*f - DYAW,  DYAW < -pi*f ~  2*pi*f + DYAW, TRUE ~ DYAW)) %>%
  select(MX, MY, MZ, AX, AY, AZ,ROLL, PITCH, YAW,DROLL, DPITCH, DYAW)

plot(inertial_frame_data$ROLL, inertial_frame_data$DROLL)
plot(inertial_frame_data$DROLL, inertial_frame_data$DPITCH)
plot(inertial_frame_data$DROLL, inertial_frame_data$DYAW)
plot(inertial_frame_data$DPITCH, inertial_frame_data$DYAW)

corr_dy_model <- lm(DYAW ~ ., inertial_frame_data %>% select(DROLL, DPITCH, DYAW))
summary(corr_dy_model)

n = dim(imu.data.mag)[1]
x = seq(1,n)
plot(x, imu.data.mag$GX)
plot(x, imu.data.mag$GY)
plot(x, imu.data.mag$GZ)

gyro_bias = c(mean(imu.data.mag$GX[1:10000]), mean(imu.data.mag$GY[1:10000]), mean(imu.data.mag$GZ[1:10000]))
# gyro data are on 1000DPS => LSB = 32.768
gyro_data <- imu.data.mag %>% 
  select(GX, GY, GZ, MV) %>%
  mutate(GX = (GX - gyro_bias[1])/32.768, GY = (GY - gyro_bias[2])/32.768, GZ = (GZ - gyro_bias[3])/32.768)
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

ts_data <- imu.data.mag %>% filter(MV == 1) %>% select(TS)
corrected_data_am <- check_data_mag_rp %>% 
  select(CMX, CMY, CMZ, CAX, CAY, CAZ, CRA, CPA, CYM) %>% 
  rename(MX = CMX, MY = CMY, MZ = CMZ, AX = CAX, AY = CAY, AZ = CAZ, AROLL = CRA, APITCH = CPA, AYAW = CYM) %>%
  mutate(DAROLL = AROLL - shift(AROLL, 1, fill = AROLL[1], type = "lag")) %>%
  mutate(DAPITCH = APITCH - shift(APITCH, 1, fill = APITCH[1], type = "lag")) %>%
  mutate(DAYAW = AYAW - shift(AYAW, 1, fill = AYAW[1], type = "lag")) %>%
  mutate(DAROLL = case_when(DAROLL > pi*f ~ 2*pi*f - DAROLL,  DAROLL < -pi*f ~  2*pi*f + DAROLL, TRUE ~ DAROLL)) %>%
  mutate(DAPITCH = case_when(DAPITCH > pi*f ~ 2*pi*f - DAPITCH,  DAPITCH < -pi*f ~  2*pi*f + DAPITCH, TRUE ~ DAPITCH)) %>%
  mutate(DAYAW = case_when(DAYAW > pi*f ~ 2*pi*f - DAYAW,  DAYAW < -pi*f ~  2*pi*f + DAYAW, TRUE ~ DAYAW))
  
corrected_data_amg <- cbind(ts_data, corrected_data_am, gyro_data_means) %>%
  rename(GX = GX_MEAN, GY = GY_MEAN, GZ = GZ_MEAN) %>%
  mutate(DT = (TS - shift(TS, 1, fill = TS[1], type = "lag"))/1000) %>%
  mutate(DGROLL = GX*DT,DGPITCH = GY*DT, DGYAW = GZ*DT) %>%
  mutate(
    RDGROLL  = DGROLL*cos(APITCH/f)*cos(AYAW/f)                                           + DGPITCH*cos(APITCH/f)*sin(AYAW/f)                                           - DGYAW*sin(APITCH/f),
    RDGPITCH = DGROLL*(sin(AROLL/f)*sin(APITCH/f)*cos(AYAW/f) - cos(AROLL/f)*sin(AYAW/f)) + DGPITCH*(sin(AROLL/f)*sin(APITCH/f)*sin(AYAW/f) + cos(AROLL/f)*cos(AYAW/f)) + DGYAW*(sin(AROLL/f)*cos(APITCH/f)),
    RDGYAW   = DGROLL*(cos(AROLL/f)*sin(APITCH/f)*cos(AYAW/f) + sin(AROLL/f)*sin(AYAW/f)) + DGPITCH*(cos(AROLL/f)*sin(APITCH/f)*sin(AYAW/f) - sin(AROLL/f)*cos(AYAW/f)) + DGYAW*(cos(AROLL/f)*cos(APITCH/f))
  ) %>%
  select(TS, MX, MY, MZ, AX, AY, AZ, GX, GY, GZ, AROLL, APITCH, AYAW, DAROLL, DGROLL, DAPITCH, DGPITCH, DAYAW, DGYAW, DT, RDGROLL, RDGPITCH, RDGYAW) %>%
  mutate(GROLL = cumsum(DGROLL) + AROLL[1]) %>%
  mutate(GPITCH = cumsum(DGPITCH) + APITCH[1]) %>%
  mutate(GYAW = cumsum(DGYAW) + AYAW[1]) %>%
  mutate(
    AXG= -sin(GPITCH/f),
    AYG= sin(GROLL/f)*cos(GPITCH/f),
    AZG= cos(GROLL/f)*cos(GPITCH/f),
  )

plot(corrected_data_amg$AX, type="l")
lines(corrected_data_amg$AXG, col="red")

plot(corrected_data_amg$AY, type="l")
lines(corrected_data_amg$AYG, col="red")

plot(corrected_data_amg$AZ, type="l")
lines(corrected_data_amg$AZG, col="red")

plot(corrected_data_amg$AROLL, corrected_data_amg$GROLL)
plot(corrected_data_amg$APITCH, corrected_data_amg$GPITCH)
plot(corrected_data_amg$AYAW, corrected_data_amg$GYAW)


###########################################################################################################################
# Qualcosa non va in corrected_data_amg per GROLL/GPITCH/GYAW
# Provo a calcolare riga per riga e confronto
# TODO: Eliminare al termine delle prove
###########################################################################################################################
check_corrected_data_gyro <- imu.data.mag %>% 
  select(GX, GY, GZ, MV, TS) %>%
  mutate(DT = (TS - shift(TS, 1, fill = TS[1], type = "lag"))/1000) %>%
  mutate(GX = (GX - gyro_bias[1])*DT/32.768, GY = (GY - gyro_bias[2])*DT/32.768, GZ = (GZ - gyro_bias[3])*DT/32.768) %>%
  mutate(GROLL = cumsum(GX) + corrected_data_amg$AROLL[1]) %>%
  mutate(GPITCH = cumsum(GY) + corrected_data_amg$APITCH[1]) %>%
  mutate(GYAW = cumsum(GZ) + corrected_data_amg$AYAW[1]) %>%
  filter(MV == 1)


###########################################################################################################################
# Dati già calibrati con unico movimento rotatorio
###########################################################################################################################
D <- read.csv('/hd/eclipse-cpp-2020-12/eclipse/workspace/my_mpu9250_thing/examples/imu-cal-data-pitch.csv')

# plot original data
scatter3D(D$MX, D$MY, D$MZ, colvar = D$MZ, col = NULL, add = FALSE, ticktype = "detailed", scale = FALSE)
plotrgl()
rglwidget()
