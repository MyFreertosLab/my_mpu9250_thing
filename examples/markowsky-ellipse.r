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

################################################################################
# from: https://stats.stackexchange.com/questions/268974/ellipse-formula-from-points
center <- c(1,2)
axis.main <- c(3,1)
axis.lengths <- c(4/3, 1/2)
sigma <- 1/20               # Error SD in each coordinate
n <- 60                     # Number of points to generate
set.seed(17)
#
# Compute the axes.
#
axis.main <- axis.main / sqrt(crossprod(axis.main))
axis.aux <- c(-axis.main[2], axis.main[1]) * axis.lengths[2]
axis.main <- axis.main * axis.lengths[1]
axes <- cbind(axis.main, axis.aux)
#
# Generate points along the ellipse.
#
s <- seq(0, 2*pi, length.out=n+1)[-1]
s.c <- cos(s)
s.s <- sin(s)
x <- axis.main[1] * s.c + axis.aux[1] * s.s + center[1] + rnorm(n, sd=sigma)
y <- axis.main[2] * s.c + axis.aux[2] * s.s + center[2] + rnorm(n, sd=sigma)
#
X <- matrix(cbind(x, y), ncol=2)


################################################################################

# Implementazione di Symmetric Kronecker Product per n=3
U=matrix(c(1, 0, 0, 0, 
           0, 1, 1, 0, 
           0, 0, 0, 1 ), nrow=3, ncol=4, byrow = TRUE)
symm_kronecker <- function(A,B,U) {
  U %*% kronecker(A,B, make.dimnames = FALSE)
}
A=matrix(c(2,3,3,2), nrow = 2, ncol=2)
x=matrix(X[1,1:2])
xk=symm_kronecker(x,x, U)
xtAx=t(xk) %*% matrix(vechMat(A))
t(x)%*%A%*%x
xtAx

y=c(xk, X[1,1:2], 1)
psiOls=matrix(y) %*% t(matrix(y))
Uno <- matrix(c(1,1,1))
Indx <- matrix(c(1,2,0))
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
#  matrix(c(vechMat(Indx %*% t(Uno)), vechMat(Uno %*% t(Indx))), ncol = 2, byrow = FALSE)
  matrix(c(1,1,
           2,1,
           2,2,
           0,1,
           0,2,
           0,0), nrow=6, ncol=2, byrow=TRUE)
}
#####################################
## Step 3): Form Tensor R 
#####################################
is_equal <- function(v1, v2) {
  as.integer(v1 == v2)
}
M <- make_index_matrix(2)
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
psi_als = make_psi_als(X, 1/20,make_index_matrix(2))

################################################
## Step 7)8): Find eigenvector of min eigenvalue
################################################
calc_eigenvector_psi_als <- function(psi_als) {
  eigen_als = eigen(psi_als)
  eigen_values = eigen_als$values
  eigen_vectors = eigen_als$vectors 
  idx = which.min(eigen_values)
  beta = matrix(eigen_vectors[, idx])
  beta = beta/norm(beta)
  beta
}
################################################
## Step 9): Estimate A, b, d
################################################
estimate_A <- function(beta, n) {
  A=xpndMat(beta[1:(n*(n+1)/2)])
  # force order
  an=A[1,2]
  ad=A[2,2]
  A[2,2]=an
  A[1,2]=ad
  A[2,1]=ad
  A
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

# T.B.D.
force_positive <- function(A) {
  
  n=3
  v1=matrix(eigen(A)$vectors[,1])
  v2=matrix(eigen(A)$vectors[,2])
  v3=matrix(eigen(A)$vectors[,3])
  
  l1=as.double(eigen(A)$values[1])
  l2=as.double(eigen(A)$values[2])
  l3=as.double(eigen(A)$values[3])
  
  print("autovalori pre")
  print(l1)
  print(l2)
  print(l3)
  if(l1 < 0 || l2 < 0 || l3 < 0) {
    A2 = matrix(rep(0,n*n), nrow = n, ncol=n)
    if(l1>0) {
      A2 = A2 + l1*(v1%*%t(v1))
      print("l1")
    }
    if(l2>0) {
      A2 = A2 + l2*(v2%*%t(v2))
      print("l2")
    }
    if(l3>0) {
      A2 = A2 + l3*(v3%*%t(v3))
      print("l3")
    }
    A2
  }  
  l1=as.double(eigen(als$Ae)$values[1])
  l2=as.double(eigen(als$Ae)$values[2])
  l3=as.double(eigen(als$Ae)$values[3])
  
  print("autovalori post")
  print(l1)
  print(l2)
  print(l3)
  
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

als=estimate_all(psi_als, 2)
evaluate <- function(x, y, beta, d) {
  if (missing(y)) {
    y <- x[, 2]; x <- x[, 1]
  }
  as.vector(cbind(x, y, x*y, x^2, y^2) %*% beta - d)
}

for(i in 1:(dim(X)[1])) {
  x=matrix(X[i,1:2])
  val=as.double(t(x)%*%als$A%*%x+t(als$b)%*%x+als$d)
  if(val > 0) {
    print(val)
  }
}

for(i in 1:(dim(X)[1])) {
  x=matrix(X[i,1:2])-matrix(als$c)
  val=as.double(t(x)%*%als$Ae%*%x)
  if(val < 0) {
    print(val)
  }
}

Ae_eigen_values=matrix(eigen(als$Ae)$values)
a=sqrt(as.double(abs(1/Ae_eigen_values[1,1])))
b=sqrt(as.double(abs(1/Ae_eigen_values[2,1])))

#beta_hat = c(als$b[1], als$b[2],als$Ae[1,1],als$Ae[2,2],2*als$Ae[2,1]) 
beta_hat = c(als$b[1], als$b[2],als$Ae[1,1],2*als$Ae[2,1], als$Ae[2,2]) 
x=X[,1]
y=X[,2]
e.x <- diff(range(x)) / 40
e.y <- diff(range(y)) / 40
n.x <- 1000
n.y <- 1000
u <- seq(min(x)-100*e.x, max(x)+e.x, length.out=n.x)
v <- seq(min(y)-100*e.y, max(y)+e.y, length.out=n.y)

z <- matrix(evaluate(as.matrix(expand.grid(u, v)), beta=beta_hat, d=als$d), n.x)
contour(u, v, z, levels=0, lwd=2, xlab="x", ylab="y", asp=1)
arrows(center[1], center[2], axis.main[1]+center[1], axis.main[2]+center[2],
       length=0.15, angle=15)
arrows(center[1], center[2], axis.aux[1]+center[1], axis.aux[2]+center[2],
       length=0.15, angle=15)
points(center[1], center[2])
points(x,y, pch=19, col="Red")

