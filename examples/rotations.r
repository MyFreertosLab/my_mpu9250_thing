toRad=pi/180
roll=190*toRad
pitch=-90*toRad
yaw=0*toRad
m <- c(cos(yaw)*cos(pitch), cos(yaw)*sin(pitch)*sin(roll) - sin(yaw)*cos(roll), cos(yaw)*sin(pitch)*cos(roll)+sin(yaw)*sin(roll), 
       sin(yaw)*cos(pitch), sin(yaw)*sin(pitch)*sin(roll) + cos(yaw)*cos(roll), sin(yaw)*sin(pitch)*cos(roll)-cos(yaw)*sin(roll),
       -sin(pitch),          cos(pitch)*sin(roll),                                cos(pitch)*cos(roll))
R <- matrix(m, nrow=3, ncol=3, byrow=TRUE)
Rz <- matrix(c(0,0,1,0,1,0,-1,0,0), nrow=3, ncol=3, byrow=TRUE)
Ux=as.matrix(c(1,0,0))
Uy=as.matrix(c(0,1,0))
Uz=as.matrix(c(0,0,1))
Uk=as.matrix(c(1,2,3))

BUx <- R%*%(Ux/norm(Ux, "2"))
BUy <- R%*%(Uy/norm(Uy, "2"))
BUz <- R%*%(Uz/norm(Uz, "2"))
BUk<- R%*%Uk/norm(Uk, "2")

print(c("BUx pitch: ", -asin(BUx[3])/toRad, " must be: ", pitch/toRad), quote = FALSE)
print(c("BUx yaw: ", atan2(BUx[2], BUx[1])/toRad, " must be: ", yaw/toRad), quote = FALSE)

# Questo è corretto in assenza di yaw e se calcolato dai dati accelerometro

if(BUz[3] >= 0) {
  roll_calc = -asin(BUz[2])
  pitch_calc = acos(BUz[3]/cos(roll_calc))
  if(BUz[1] < 0) {
    pitch_calc = -pitch_calc
  }
} else {
  roll_calc = -pi+asin(BUz[2])
  pitch_calc = acos(BUz[3]/cos(roll_calc))
  if(BUz[1] >= 0) {
    pitch_calc = -pitch_calc
  }
}
if(roll_calc > pi) {
  roll_calc = roll_calc - 2*pi
} else if(roll_calc < -pi) {
  roll_calc = roll_calc + 2*pi
}
print(c("BUz roll: ", roll_calc/toRad, " must be: ", roll/toRad), quote = FALSE)
print(c("BUz pitch: ", pitch_calc/toRad, " must be: ", pitch/toRad), quote = FALSE)


