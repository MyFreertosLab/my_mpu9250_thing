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
imu.data.mag <- read.csv('/hd/eclipse-cpp-2020-12/eclipse/workspace/my_mpu9250_thing/examples/imu-data-mag.csv')
  
imu.data.body <- imu.data.mag %>% select(MY, MX, MZ, MV) %>% filter(MV == 1) %>% select(MY, MX, MZ)
imu.data.body <- imu.data.body %>% rename(MX = MY, MY = MX) %>% mutate(MZ = -MZ)

# plot original data
scatter3D(imu.data.body$MX, imu.data.body$MY, imu.data.body$MZ, colvar = imu.data.body$MZ, col = NULL, add = FALSE)
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
estimate_all <- function(psi_als, n) {
  beta = calc_eigenvector_psi_als(psi_als)
  A=estimate_A(beta, n)
  b=estimate_b(beta, n)
  d=estimate_d(beta, n)
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

mag_plot_data <- function(mag_data, sphere_radius = -1) {
  scatter3D(mag_data$MX, mag_data$MY, mag_data$MZ, colvar = mag_data$MZ, col = NULL, add = FALSE)
  plotrgl()
  rglwidget()
}

# estimation
psi_als = make_psi_als(as.matrix(imu.data.body), 0.1720985^2,make_index_matrix(3))
als=estimate_all(psi_als, 3)
prova <- mag_estimate(als, imu.data.body)
prova_data <- mag_apply_estimator(prova)

# Plotting
mag_plot_data(prova$data_source)
mag_plot_data(prova_data)
sphere_radius = sqrt(mean(prova_data$MX^2+prova_data$MY^2+prova_data$MZ^2))

scatter3D(imu.data.body$MX-als$c[1], imu.data.body$MY-als$c[2], imu.data.body$MZ-als$c[3], col = "green", add = FALSE)
scatter3D(prova_data$MX*250, prova_data$MY*250, prova_data$MZ*250, col = "red", add = TRUE)
plotrgl()
rglwidget()

print(c("sphere radius: ", sphere_radius), quote = FALSE)

# calcola rotazione yaw, pitch, roll
V <- eigen(prova$invA)$vectors
v1 <- V[,1]
v2 <- V[,2]
v3 <- V[,3]
pitch <- asin(-v1[3])
yaw <- acos(v1[1]/cos(pitch))
roll <- asin(v2[3]/cos(pitch))
print(c("yaw, pitch roll: ", yaw/pi*180, pitch/pi*180, roll/pi*180), quote = FALSE)
