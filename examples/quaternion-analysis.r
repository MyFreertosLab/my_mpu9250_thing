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

###############################################
### Quaternions
###############################################
quat_def <- function(axis, angle) {
  result = as.matrix(c(0,0,0,0))
  v <- axis/norm(axis, "2")
  result[1] = cos(angle/2)
  result[2:4] = v*sin(angle/2)
  return(result)
}
quat_sum <- function(q1, q2) {
  return(q1+q2)
}
quat_prod <- function(q, p) {
  result = as.matrix(c(0,0,0,0))
  result[1] = q[1]*p[1]-q[2]*p[2]-q[3]*p[3]-q[4]*p[4]
  result[2] = q[1]*p[2]+q[2]*p[1]+q[3]*p[4]-q[4]*p[3]
  result[3] = q[1]*p[3]-q[2]*p[4]+q[3]*p[1]+q[4]*p[2]
  result[4] = q[1]*p[4]+q[2]*p[3]-q[3]*p[2]+q[4]*p[1]
  return(result)
}
quat_conj <- function(q) {
  result = -q
  result[1] = q[1]
  return(result)
}
quat_inv <- function(q) {
  n = norm(q, "2")^2
  if(n > 0) {
    result = quat_conj(q)/n
  } else {
    result = as.matrix(c(0,0,0,0))
  }
  return(result)
}
quat_omegadt <- function(omega, dt) {
  axis = omega/norm(omega, "2")
  angle = dt*norm(omega, "2")
  result <- quat_def(axis, angle)
  return(result)
}

quat_rot <- function(q, p) {
  result <- quat_prod(q, quat_prod(p, quat_inv(q)))
  return(result)
} 

quat_toquat <- function(u) {
  result = as.matrix(c(0,0,0,0))
  result[2:4] = u
  return(result)
}
quat_ng <- function(u, gamma = 1) {
  v <- u/norm(u,"2")
  r <- as.matrix(c(0,0,-1))
  axis <- cross(v, r)
  angle <- acos(-v[3])
  result <- quat_def(axis, gamma*angle)
  return(result)
}

########################
## Queternion Tests
########################
p <- as.matrix(c(1,1,1))
q <- quat_ng(p)
qj <- quat_conj(q)
stopifnot(q[1] == qj[1], q[2:4] == -qj[2:4])
stopifnot(round(quat_prod(q,quat_inv(q)),14) == round(as.matrix(c(1,0,0,0)),14))
stopifnot(round(quat_rot(q,quat_toquat(p)), 14) == round(as.matrix(c(0,0,0,-norm(p,"2"))), 14))

###############################################
### Simplified Attitude Determination Algorithm
###############################################
sada_update_df <- function(df) {
  for(row_idx in 1:dim(df)[1]) {
    curr_data <- df[row_idx,]
    with_mag = as.integer(curr_data$MV)
    if(with_mag == 1) {
      m <- as.matrix(c(curr_data$MX, curr_data$MY,curr_data$MZ))
      a <- as.matrix(c(curr_data$AX, curr_data$AY,curr_data$AZ))
      m <- m/norm(m, "2")
      a <- a/norm(a, "2")
      ar <- as.matrix(accel_reference)
      MD = -t(m)%*%a # negative because accelerometer point to '-g'
      MN = sqrt(1-MD^2)
      mr <- as.matrix(c(MN, 0, MD))
      #mr <- as.matrix(c(cos(declination), 0, sin(declination)))
      
      B <- sada_B_matrix(m, a, mr, ar)
      q <- sada_quaternion(B)
      rpy <- (Q2EA.Xiao(q,EulerOrder = "xyz")*toDeg)
      df$sada_roll[row_idx] = rpy[1]
      df$sada_pitch[row_idx] = rpy[2]
      df$sada_yaw[row_idx] = rpy[3]
      df$w[row_idx] = q[1]
      df$a[row_idx] = q[2]
      df$b[row_idx] = q[3]
      df$c[row_idx] = q[4]
      df$MD[row_idx] = MD
      df$MN[row_idx] = MN
      ng <- Q2DCM(t(q))%*%(accel_reference)
      df$sada_ngx[row_idx] = ng[1]
      df$sada_ngy[row_idx] = ng[2]
      df$sada_ngz[row_idx] = ng[3]
    } else {
      df$sada_ngx[row_idx] = NA
      df$sada_ngy[row_idx] = NA
      df$sada_ngz[row_idx] = NA
      df$sada_roll[row_idx] = NA
      df$sada_pitch[row_idx] = NA
      df$sada_yaw[row_idx] = NA
      df$w[row_idx] = NA
      df$a[row_idx] = NA
      df$b[row_idx] = NA
      df$c[row_idx] = NA
      df$MD[row_idx] = NA
      df$MN[row_idx] = NA
    }
  }
  return(df)
}

########################################
### Gyroscope Accelerometer Fusion (GAF)
########################################
gyroacc_fusion <- function(df, gamma) {
  qprev = NULL
  prev_data = NULL
  cum_dt = 0
  cum_omega = as.matrix(c(0,0,0))
  cum_items = 0
  for(row_idx in 1:dim(df)[1]) {
    curr_data <- df[row_idx,]
    cum_dt = cum_dt+as.double(curr_data$DT)
    cum_items = cum_items + 1
    o = as.matrix(c(
      curr_data$GX*toRad,
      curr_data$GY*toRad,
      curr_data$GZ*toRad
    ))
    cum_omega = cum_omega+o
    if(curr_data$MV == 1) {
      if(is.null(qprev)) {
        qprev <- as.matrix(c(curr_data$w, curr_data$a, curr_data$b, curr_data$c))
        prev_data <- curr_data
        q <- qprev
        ng <- as.matrix(quat_rot(quat_inv(q), quat_toquat(accel_reference))[2:4])
        df$qfw[row_idx] = q[1]
        df$qfx[row_idx] = q[2]
        df$qfy[row_idx] = q[3]
        df$qfz[row_idx] = q[4]
        df$ngx[row_idx] = ng[1]
        df$ngy[row_idx] = ng[2]
        df$ngz[row_idx] = ng[3]
        rpy <- Q2EA.Xiao(t(q), EulerOrder = "xyz")
        df$ng_roll[row_idx] = rpy[1]*toDeg
        df$ng_pitch[row_idx] = rpy[2]*toDeg
        df$ng_yaw[row_idx] = rpy[3]*toDeg
        cum_dt = 0
        cum_omega = as.matrix(c(0,0,0))
        cum_items = 0
        next
      }
      omega_mean = cum_omega/cum_items
      cum_omega = as.matrix(c(0,0,0))
      cum_items = 0
      omega <- matrix(c(
        0,             -omega_mean[1],  -omega_mean[2], -omega_mean[3],
        omega_mean[1],   0,             -omega_mean[3],  -omega_mean[2],
        omega_mean[2],  omega_mean[3],  0,             omega_mean[1],
        omega_mean[3],  omega_mean[2],   -omega_mean[1], 0
      ), nrow = 4, ncol = 4, byrow = T)
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
      dt <- cum_dt
      cum_dt = 0
      ####################################################################################
      # This is the original formula, but does not works. Yaw from accelerometer is zero.
      # q <- (1-gamma)*(dt/2*omega + diag(1, 4))%*%qprev + gamma*(W/2 - diag(1,4))%*%qprev
      # I use sada quaternion for gyroscope fusion
      ####################################################################################
      #q <- (1-gamma)*(quat_prod(qprev, 0.5*dt*quat_toquat(omega_mean))) + gamma*sada_quat
      q <- (1-gamma)*(dt/2*omega + diag(1, 4))%*%qprev + gamma*sada_quat
      q <- q/norm(q, "2")
      prev_data <- curr_data
      ng <- as.matrix(quat_rot(quat_inv(q), quat_toquat(accel_reference))[2:4])
      df$qfw[row_idx] = q[1]
      df$qfx[row_idx] = q[2]
      df$qfy[row_idx] = q[3]
      df$qfz[row_idx] = q[4]
      df$ngx[row_idx] = ng[1]
      df$ngy[row_idx] = ng[2]
      df$ngz[row_idx] = ng[3]
      rpy <- tryCatch(
        {
          t(as.matrix(Q2EA.Xiao(t(q), EulerOrder = "xyz")))
        },
        error=function(cond){
          t(as.matrix(Q2EA.Xiao(t(qprev), EulerOrder = "xyz")))
        }
      )
      
      df$ng_roll[row_idx] = rpy[1]*toDeg
      df$ng_pitch[row_idx] = rpy[2]*toDeg
      df$ng_yaw[row_idx] = rpy[3]*toDeg
      qprev <- q
    } else {
      df$qfw[row_idx] = NA
      df$qfx[row_idx] = NA
      df$qfy[row_idx] = NA
      df$qfz[row_idx] = NA
      df$ngx[row_idx] = NA
      df$ngy[row_idx] = NA
      df$ngz[row_idx] = NA
      df$ng_roll[row_idx] = NA
      df$ng_pitch[row_idx] = NA
      df$ng_yaw[row_idx] = NA
    }
  }  
  return(df)
}

gyroacc_fusion_2 <- function(df, gamma) {
  qprev = NULL
  for(row_idx in 1:dim(df)[1]) {
    curr_data <- df[row_idx,]
    accel = as.matrix(c(
      curr_data$AX,
      curr_data$AY,
      curr_data$AZ
    ))
    omega = as.matrix(c(
      curr_data$GX*toRad,
      curr_data$GY*toRad,
      curr_data$GZ*toRad
    ))
    if(is.null(qprev)) {
      qprev <- quat_ng(accel)
      prev_data <- curr_data
      q <- qprev
      ng <- as.matrix(quat_rot(quat_inv(q), quat_toquat(accel_reference))[2:4])
      df$ng2x[row_idx] = ng[1]
      df$ng2y[row_idx] = ng[2]
      df$ng2z[row_idx] = ng[3]
      rpy <- Q2EA.Xiao(t(q), EulerOrder = "xyz")
      df$ng2_roll[row_idx] = rpy[1]*toDeg
      df$ng2_pitch[row_idx] = rpy[2]*toDeg
      df$ng2_yaw[row_idx] = rpy[3]*toDeg
      next
    }
    if(as.double(curr_data$MV) == 1) {
      q <- as.matrix(c(
        curr_data$qfw,
        curr_data$qfx,
        curr_data$qfy,
        curr_data$qfz
      ))
      df$ng2x[row_idx] = curr_data$ngx
      df$ng2y[row_idx] = curr_data$ngy
      df$ng2z[row_idx] = curr_data$ngz
      df$ng2_roll[row_idx] = curr_data$ng_roll
      df$ng2_pitch[row_idx] = curr_data$ng_pitch
      df$ng2_yaw[row_idx] = curr_data$ng_yaw
      qprev <- q
    } else {
      dt <- as.double(curr_data$DT)
      dq <- quat_omegadt(omega, dt) # gyro rotation
      q <- quat_prod(qprev, dq) # apply gyro rotation
      qa_body <- quat_toquat(accel)
      qa_if <- quat_rot(q, qa_body)
      qc <- quat_prod(quat_ng(as.matrix(qa_if[2:4]), gamma = gamma), q) # correction
      
      rpy <- tryCatch(
        {
          t(as.matrix(Q2EA.Xiao(t(qc), EulerOrder = "xyz")))
        },
        error=function(cond){
          t(as.matrix(Q2EA.Xiao(t(qprev), EulerOrder = "xyz")))
        }
      )
      qa_body_corrected <- quat_rot(quat_inv(qc), quat_toquat(accel_reference)) # recalc ax, ay, az from ng
      df$ng2x[row_idx] = qa_body_corrected[2]
      df$ng2y[row_idx] = qa_body_corrected[3]
      df$ng2z[row_idx] = qa_body_corrected[4]
      df$ng2_roll[row_idx] = rpy[1]*toDeg
      df$ng2_pitch[row_idx] = rpy[2]*toDeg
      df$ng2_yaw[row_idx] = rpy[3]*toDeg
      df$ng2_roll[row_idx] = NA
      df$ng2_pitch[row_idx] = NA
      df$ng2_yaw[row_idx] = NA
      qprev <- qc
    }
  }
  return(df)
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
gyro_bias = c(mean(imu.cal.data.gyro$GX[1:10000]), mean(-imu.cal.data.gyro$GY[1:10000]), mean(-imu.cal.data.gyro$GZ[1:10000]))
imu.cal.data.gyro.ned <- imu.cal.data.gyro %>% 
  mutate(MY = -MY, MZ = -MZ,AY = -AY, AZ = -AZ,GY = -GY, GZ = -GZ) %>%
  mutate(GX = GX - gyro_bias[1],GY = GY - gyro_bias[2],GZ = GZ - gyro_bias[3]) %>%
  mutate(GX = GX/32.768, GY = GY/32.768, GZ = GZ/32.768) %>%
  mutate(DT = (TS - shift(TS, 1, fill = TS[1], type = "lag"))/1000) %>%
  mutate(DGROLL = GX*DT,DGPITCH = GY*DT, DGYAW = GZ*DT) %>%
  mutate(GROLL = cumsum(DGROLL)) %>%
  mutate(GPITCH = cumsum(DGPITCH)) %>%
  mutate(GYAW = cumsum(DGYAW)) %>%
  mutate(NRMG = sqrt(GX^2+GY^2+GZ^2)) %>%
  mutate(NRMA = sqrt(AX^2+AY^2+AZ^2)) %>%
  mutate(NRMM = sqrt(MX^2+MY^2+MZ^2)) %>%
  mutate(MX = MX/NRMM, MY = MY/NRMM, MZ = MZ/NRMM) %>%
  mutate(AX = AX/NRMA, AY = AY/NRMA, AZ = AZ/NRMA) %>%
  select(-DGROLL, -DGPITCH, -DGYAW)
  
mag_plot_data <- function(mag_data, sphere_radius = -1, title = "", add = FALSE, widget = TRUE) {
  scatter3D(mag_data[,1], mag_data[,2], mag_data[,3], colvar = mag_data[,3], col = NULL, add = add, scale = FALSE, ticktype = "detailed", main = title)
  if(widget) {
    plotrgl()
    rglwidget()
  }
}
  
# plot original data
mag_plot_data(imu.cal.data.gyro.ned %>% select(MX, MY, MZ), widget = FALSE)
mag_plot_data(imu.cal.data.gyro.ned %>% select(AX, AY, AZ), widget = TRUE, add = TRUE)
imu.cal.data.gyro_sada <- sada_update_df(imu.cal.data.gyro.ned)

rpy <- as.matrix(imu.cal.data.gyro_sada %>% 
  filter(MV == 1) %>%
  select(sada_roll, sada_pitch, sada_yaw), ncol = 3, byrow = TRUE)

imu.cal.data.gyro_sada_gyro <-  imu.cal.data.gyro_sada %>%
  mutate(TDEC = -asin(MX*AX+MY*AY+MZ*AZ)*toDeg) %>%
  mutate(GROLL = GROLL + mean(rpy[1:3000,1])) %>%
  mutate(GPITCH = GPITCH + mean(rpy[1:3000,2])) %>%
  mutate(GYAW = GYAW + mean(rpy[1:3000,3])) %>%
  select(MV, TS, DT, MX, MY, MZ, AX, AY, AZ, GX, GY, GZ,sada_ngx, sada_ngy, sada_ngz, GROLL, GPITCH, GYAW, sada_roll, sada_pitch, sada_yaw, TDEC, MN, MD, w, a, b, c)

SADA <- gyroacc_fusion(imu.cal.data.gyro_sada_gyro, 0.0)
df <- SADA %>% filter(MV == 1)

df <- gyroacc_fusion_2(SADA, 0.0) %>% filter(MV == 1)
ylim_min <- min(df$sada_roll, df$GROLL) - 0.5
ylim_max <- max(df$sada_roll, df$GROLL) + 0.5
plot(df$sada_roll, type="l", main = "SADA_ROLL (Black) vs GROLL (Red) vs Fusion (Blue)", ylab = "Roll", ylim = c(ylim_min,ylim_max))
lines(df$GROLL, col="red")
lines(df$ng_roll, col="blue")
lines(df$ng2_roll, col="green")
lines(df$MZ*6, col="red", )

ylim_min <- min(df$sada_pitch, df$ng_pitch, df$sada_pitch) - 0.5
ylim_max <- max(df$sada_pitch, df$ng_pitch, df$sada_pitch) + 0.5
plot(df$sada_pitch, type="l", main = "SADA_PITCH (Black) vs GPITCH (Red) vs Fusion (Blue)", ylab = "Pitch", ylim = c(ylim_min,ylim_max))
lines(df$GPITCH, col="red")
lines(df$ng_pitch, col="blue")
lines(df$ng2_pitch, col="green")

ylim_min <- min(df$sada_yaw, df$ng_yaw, df$sada_yaw) - 0.5
ylim_max <- max(df$sada_yaw, df$ng_yaw, df$sada_yaw) + 0.5
plot(df$sada_yaw, type="l", main = "SADA_YAW (Black) vs GYAW (Red) vs Fusion (Blue)", ylab = "Yaw", ylim = c(ylim_min,ylim_max))
lines(df$GYAW, col="red")
lines(df$ng_yaw, col="blue")
lines(df$ng2_yaw, col="green")
lines(df$MZ*100, col="blue", )

ylim_min <- min(df$AX, df$sada_ngx) - 0.5
ylim_max <- max(df$AX, df$sada_ngx) + 0.5
plot(df$AX, type="l", main = "AX (Black) vs NGX (Red)", ylab = "AX vs. NGX vs NG2X", ylim = c(ylim_min,ylim_max))
lines(df$ngx, col="red", )
lines(df$ng2x, col="green", )
lines(df$MZ, col="blue", )

ylim_min <- min(df$AY, df$sada_ngy) - 0.5
ylim_max <- max(df$AY, df$sada_ngy) + 0.5
plot(df$AY, type="l", main = "AY (Black) vs NGY (Red)", ylab = "AY vs. NGY", ylim = c(ylim_min,ylim_max))
lines(df$ngy, col="red", )
lines(df$ng2y, col="green", )
lines(df$MZ, col="blue", )

ylim_min <- min(df$AZ, df$sada_ngz) - 0.5
ylim_max <- max(df$AZ, df$sada_ngz) + 0.5
plot(df$AZ, type="l", main = "AZ (Black) vs NGZ (Red)", ylab = "AZ vs. NGZ", ylim = c(ylim_min,ylim_max))
lines(df$ngz, col="red", )
lines(df$ng2z, col="green", )
lines(df$MZ, col="blue", )

mdec<- mean(df$TDEC)
vardec <- var(df$TDEC)
tdec_err <- (df$TDEC - mdec)
plot(tdec_err, type="l", main = paste("Declination Error (mean = ", round(mdec, 2), ")") , ylab = "TDEC Error")

plot(df$sada_roll, df$GROLL, xlab = "Roll from SADA", ylab = "Roll from Gyroscope", main = "Roll: SADA vs Gyro", xlim = c(-180,180), ylim = c(-180,180))
plot(df$sada_pitch, df$GPITCH, xlab = "Pitch from SADA", ylab = "Pitch from Gyroscope", main = "Pitch: SADA vs Gyro", xlim = c(-180,180), ylim = c(-180,180))
plot(df$sada_yaw, df$GYAW, xlab = "Yaw from SADA", ylab = "Yaw from Gyroscope", main = "Yaw: SADA vs Gyro", xlim = c(-180,180), ylim = c(-180,180))

ylim_min <- min(df$MX,df$ngx,df$AX) - 0.5
ylim_max <- max(df$MX,df$ngx,df$AX) + 0.5
plot(df$MX, type="l", main = "MX(black),NGX(red),AX(blue),SADA_NGX(green)",ylab = "", ylim = c(ylim_min,ylim_max))
lines(df$ngx, col="red")
lines(df$AX, col="blue")
lines(df$sada_ngx, col="green")

ylim_min <- min(df$MY,df$ngy,df$AY) - 0.5
ylim_max <- max(df$MY,df$ngy,df$AY) + 0.5
plot(df$MY, type="l", main = "MY(black),NGY(red),AY(blue),SADA_NGY(green)",ylab = "", ylim = c(ylim_min,ylim_max))
lines(df$ngy, col="red")
lines(df$AY, col="blue")
lines(df$sada_ngy, col="green")

ylim_min <- min(df$MZ,df$ngz,df$AZ) - 0.5
ylim_max <- max(df$MZ,df$ngz,df$AZ) + 0.5
plot(df$MZ, type="l", main = "MZ(black),NGZ(red),AZ(blue),SADA_NGZ(green)",ylab = "", ylim = c(ylim_min,ylim_max))
lines(df$ngz, col="red")
lines(df$AZ, col="blue")
lines(df$sada_ngz, col="green")

ylim_min <- min(df$sada_ngz-df$ngz) - 0.5
ylim_max <- max(df$sada_ngz-df$ngz) + 0.5
plot(df$sada_ngz-df$ngz, type="l", main = "SADA_NGZ - NGZ",ylab = "", ylim = c(ylim_min,ylim_max))
df1 <- df %>% mutate(NGZ_DIFF = sada_ngz - ngz )
# df at row 8462 (SADA 45077), sada_quat ~= quat_conj(q) and theta it's 180Deg

