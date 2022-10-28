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

# axis correction
imu.data.body.original <- imu.data.mag %>% select(MY, MX, MZ, MV) %>% filter(MV == 1) %>% select(MY, MX, MZ) %>% rename(MX = MY, MY = MX) %>% mutate(MZ = -MZ)
imu.data.body <- imu.data.body.original

# plot original data
scatter3D(imu.data.body$MX, imu.data.body$MY, imu.data.body$MZ, colvar = imu.data.body$MZ, col = NULL, add = FALSE)
plotrgl()
rglwidget()

# create model
imu.data.mag.poly <- imu.data.body %>% select(MX, MY, MZ) %>% mutate(One = -1, MX2 = MX**2, MY2 = MY**2, MZ2=MZ**2, MXY=MX*MY, MXZ=MX*MZ, MYZ=MY*MZ)
fit <- lm(One ~ . - 1, imu.data.mag.poly)
coeff.hat <- coef(fit)
P <- as.matrix(imu.data.mag.poly %>% select(-One))

a=coeff.hat["MX2"]
b=coeff.hat["MXY"]/2
c=coeff.hat["MXZ"]/2
e=coeff.hat["MYZ"]/2
d=coeff.hat["MY2"]
f=coeff.hat["MZ2"]

A = xpndMat(c(a,b,c,d,e,f))
B = matrix(c(coeff.hat["MX"],coeff.hat["MY"],coeff.hat["MZ"]))
X=as.matrix(imu.data.mag.poly %>% select(MX, MY, MZ))
eigen(A)
D=matrix(c(eigen(A)$values[1], 0,0, 0,eigen(A)$values[2],0,0,0,eigen(A)$values[3]), nrow=3, ncol=3)
NV=eigen(A)$vectors/norm(eigen(A)$vectors)
V=eigen(A)$vectors
TA=V%*%D%*%t(V)
stopifnot(abs(A - TA) < 1e-15)
TD=t(V)%*%A%*%V
stopifnot(abs(D - TD) < 1e-15)

# calc center

# from: https://math.stackexchange.com/questions/1217796/compute-center-axes-and-rotation-from-equation-of-ellipse
# Translation paragraph
# si risolve: D[1,1]*X^2 + D[2,2]*Y^2 + D[3,3]*Z^2+LB[1]*X+LB[2]*Y+LB[3]*Z+1=0
#             Le soluzioni forniscono il centro
# NOTA: si è posto F = 1, i coefficienti sono normalizzati rispetto al valore reale di F
#       Trovare un modo per stimare F
LB=t(B)%*%V
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
spheres3d(c(0,0,0), radius = 229.1845, col = "red", alpha = 0.4)
# apply translation
imu.data.body <- imu.data.body %>% mutate(MX = MX - C[1],MY = MY - C[2],MZ = MZ - C[3])
imu.data.body.original <- imu.data.body.original %>% mutate(MX = MX - C[1],MY = MY - C[2],MZ = MZ - C[3])

# plot after translation
open3d()
plot3d(imu.data.body, col = "blue")
spheres3d(c(0,0,0), radius = 229.1845, col = "red", alpha = 0.4)

# Rotation
imu.data.body.rotated <- imu.data.body %>% mutate(multiplication= as.matrix(imu.data.body[,]) %*% V) %>% select(multiplication) %>% mutate(MX=multiplication[,1], MY=multiplication[,2], MZ=multiplication[,3]) %>% select(MX,MY,MZ)
open3d()
plot3d(imu.data.body.rotated, col = "blue")
spheres3d(c(0,0,0), radius = 229.1845, col = "red", alpha = 0.4)

scatter3D(imu.data.body.rotated$MX, imu.data.body.rotated$MY, imu.data.body.rotated$MZ, colvar = imu.data.body.rotated$MZ, col = NULL, add = FALSE)
plotrgl()
rglwidget()

# calc yaw pitch roll
pitch=asin(-V[3,1])
roll=asin(V[3,2]/cos(pitch))
yaw=asin(V[2,1]/cos(pitch))
pitch/(2*pi)*360
roll/(2*pi)*360
yaw/(2*pi)*360

# calc axis lengths (from X^2/a^2+Y^2/b^2+Z^2/c^2=1) and scale factors
imu.data.mag.poly.rotated <- imu.data.body.rotated %>% select(MX, MY, MZ) %>% mutate(One = 1, MX2 = MX**2, MY2 = MY**2, MZ2=MZ**2) %>% select(One, MX2, MY2, MZ2)
fit_rotated <- lm(One ~ . - 1, imu.data.mag.poly.rotated)
coeff_rotated.hat <- coef(fit_rotated)
ellipsoid_axis = matrix(0,3,3)
diag(ellipsoid_axis) <- c(sqrt(1/coeff_rotated.hat[1]),sqrt(1/coeff_rotated.hat[2]),sqrt(1/coeff_rotated.hat[3]))
scale_factors_rotated = matrix(0,3,3)
sphere_radius = sqrt(mean(diag(ellipsoid_axis)^2))
diag(scale_factors_rotated) <- sphere_radius/diag(ellipsoid_axis)

# plot sphere
imu.data.mag.sphere <- imu.data.body.rotated %>% mutate(multiplication= as.matrix(imu.data.body.rotated[,]) %*% scale_factors_rotated) %>% select(multiplication) %>% mutate(MX=multiplication[,1], MY=multiplication[,2], MZ=multiplication[,3]) %>% select(MX,MY,MZ)
open3d()
plot3d(imu.data.mag.sphere, col = "blue")
spheres3d(c(0,0,0), radius = sphere_radius, col = "red", alpha = 0.4)

#scatter3D(imu.data.mag.sphere$MX, imu.data.mag.sphere$MY, imu.data.mag.sphere$MZ, colvar = imu.data.mag.sphere$MZ, col = NULL, add = FALSE)
#plotrgl()
#rglwidget()

# return to original space (from rotated)
#scale_factors = V %*% scale_factors_rotated %*% t(V)
#imu.data.mag.sphere <- imu.data.body %>% mutate(multiplication=  as.matrix(imu.data.body[,]) %*% scale_factors) %>% select(multiplication) %>% mutate(MX=multiplication[,1], MY=multiplication[,2], MZ=multiplication[,3]) %>% select(MX,MY,MZ)
imu.data.mag.sphere <- imu.data.mag.sphere %>% mutate(multiplication= as.matrix(imu.data.mag.sphere[,]) %*% t(V)) %>% select(multiplication) %>% mutate(MX=multiplication[,1], MY=multiplication[,2], MZ=multiplication[,3]) %>% select(MX,MY,MZ)
open3d()
plot3d(imu.data.mag.sphere, col = "blue")
spheres3d(c(0,0,0), radius = sphere_radius, col = "yellow", alpha = 0.4)

# plot scaled original data
scatter3D(imu.data.mag.sphere$MX, imu.data.mag.sphere$MY, imu.data.mag.sphere$MZ, colvar = imu.data.mag.sphere$MZ, col = NULL, add = FALSE)
plotrgl()
rglwidget()

# recalc model from sphere
imu.data.mag.poly <- imu.data.mag.sphere %>% select(MX, MY, MZ) %>% mutate(One = 1, MX2 = MX**2, MY2 = MY**2, MZ2=MZ**2) %>% select(One, MX2, MY2, MZ2)
fit_spheric <- lm(One ~ . - 1, imu.data.mag.poly)
coeff_spheric.hat <- coef(fit_spheric)


######################################################################################

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

# Implementazione di Symmetric Kronecker Product per n=3
U=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 1, 0, 1, 0, 0, 0, 0, 0,
           0, 0, 1, 0, 0, 0, 1, 0, 0,
           0, 0, 0, 0, 1, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 1, 0, 1, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=6, ncol=9, byrow = TRUE)
symm_kronecker <- function(A,B,U) {
  U %*% kronecker(A,B)
}
x=matrix(X[1,1:3])
xk=symm_kronecker(x,x, U)
xtAx=t(xk) %*% matrix(vechMat(A))
t(x)%*%A%*%x
xtAx
y=c(xk, X[1,1:3], 1)
psiOls=matrix(y) %*% t(matrix(y))
Uno <- matrix(c(1,1,1,1))
Indx <- matrix(c(1,2,3,0))
M = matrix(c(vechMat(Indx %*% t(Uno)), vechMat(Uno %*% t(Indx))), ncol = 2, byrow = FALSE)

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
  x^4 - 3*x^2*variance + 3*variance^2
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
  Uno <- matrix(rep(1,n+1))
  Indx <- matrix(c(c(1:n),0))
  matrix(c(vechMat(Indx %*% t(Uno)), vechMat(Uno %*% t(Indx))), ncol = 2, byrow = FALSE)
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
  D1 <- seq(1,(n+1)*n/2)
  D2 <- seq(1,n)*(seq(1,n)+1)/2
  setdiff(D1, D2)
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
psi_als = make_psi_als(X, 315.70,make_index_matrix(3))

################################################
## Step 7)8): Find eigenvector of min eigenvalue
################################################
calc_eigenvector_psi_als <- function(psi_als) {
  eigen_als = eigen(psi_als)
  eigen_values = eigen_als$values
  eigen_vectors = eigen_als$vectors 
  idx = which.min(eigen_values)
  beta = matrix(eigen_vectors[, idx])
  beta/norm(beta)
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
  c=-0.5*solve(estimate_A(psi_als,n))%*%estimate_b(psi_als,n)
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
  
  v1=matrix(eigen(Ae)$vectors[,1])
  v2=matrix(eigen(Ae)$vectors[,2])
  v3=matrix(eigen(Ae)$vectors[,3])
  
  l1=as.double(eigen(Ae)$values[1])
  l2=as.double(eigen(Ae)$values[2])
  l3=as.double(eigen(Ae)$values[3])

  print("autovalori pre")
  print(l1)
  print(l2)
  print(l3)
  if(l1 < 0 || l2 < 0 || l3 < 0) {
    Ae2 = matrix(rep(0,n*n), nrow = n, ncol=n)
    if(l1>0) {
      Ae2 = Ae2 + l1*(v1%*%t(v1))
      print("l1")
    }
     if(l2>0) {
      Ae2 = Ae2 + l2*(v2%*%t(v2))
      print("l2")
    }
    if(l3>0) {
      Ae2 = Ae2 + l3*(v3%*%t(v3))
      print("l3")
    }
    als$Ae = Ae2
  }  
  l1=as.double(eigen(als$Ae)$values[1])
  l2=as.double(eigen(als$Ae)$values[2])
  l3=as.double(eigen(als$Ae)$values[3])
  
  print("autovalori post")
  print(l1)
  print(l2)
  print(l3)
  als
}

als=estimate_all(psi_als, 3)

for(i in 1:(dim(X)[1])) {
  x=matrix(X[i,1:3])
  val=as.double(t(x)%*%als$A%*%x+t(als$b)%*%x+als$d)
  if(val > 0) {
    print(val)
  }
}

for(i in 1:(dim(X)[1])) {
  x=matrix(X[i,1:3])-matrix(als$c)
  val=as.double(t(x)%*%als$Ae%*%x)
  if(val < 0) {
    print(val)
  }
}
Ae_eigen_values=matrix(eigen(als$Ae)$values)
a=sqrt(as.double(abs(1/Ae_eigen_values[1,1])))
b=sqrt(as.double(abs(1/Ae_eigen_values[2,1])))
c=sqrt(as.double(abs(1/Ae_eigen_values[3,1])))


# T.B.D
Q=AE
k=t(C) %*% Q %*% C -1
u=-2*t(AE) %*% C
h=matrix(X[1,1:3])
t(h) %*% Q %*% (h) + t(u) %*% h + k


