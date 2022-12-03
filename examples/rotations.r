#########################################################################
### Author: Andrea Lambruschini                                       ###
### Date .: 03 December 2022                                          ###
### Calculus of Roll,Pitch (as Axis Rotation) from accelerometer Data ###
###      z                                                            ###
###      | /x                                                         ###
###  y___|/                                                           ###
#########################################################################
#rm(list = ls())
toRad=pi/180
roll=-1.600485e+00*toRad
pitch=7.496511*toRad

# Range Roll [-pi, pi]
# Range Pitch [-pi/2, pi/2]
if(pitch > pi/2) {
  pitch = pitch - pi
  roll = roll + pi
  print(roll)
} else if(pitch < -pi/2) {
  pitch = pitch + pi
  roll = roll + pi
  print(roll)
}
if(roll > pi) {
  roll = roll - 2*pi
} else if(roll < -pi) {
  roll = roll + 2*pi
}
# Roll,Pitch rotation of a vector
m <- c(cos(pitch),  sin(pitch)*sin(roll), sin(pitch)*cos(roll), 
       0,           cos(roll),            -sin(roll),
       -sin(pitch), cos(pitch)*sin(roll), cos(pitch)*cos(roll))
R <- matrix(m, nrow=3, ncol=3, byrow=TRUE)
# Roll, Pitch rotation of axis
R <- t(R)

# Uz = -g (-gravity) in Inertial Frame
Uz=as.matrix(c(0,0,1))

# Uz in Body Frame
BUz <- R%*%(Uz/norm(Uz, "2"))

# Calculate roll, pitch
pitch_calc = -asin(BUz[1])
roll_calc = acos(BUz[3]/cos(pitch_calc))
if(BUz[3] >= 0) {
  if(BUz[2] < 0) {
    roll_calc = -roll_calc
  }
} else {
  if(BUz[2] < 0) {
    roll_calc = -roll_calc
  }
}
if(roll_calc > pi) {
  roll_calc = roll_calc - 2*pi
} else if(roll_calc < -pi) {
  roll_calc = roll_calc + 2*pi
}
print(c("BUz roll: ", roll_calc/toRad, " must be: ", roll/toRad), quote = FALSE)
print(c("BUz pitch: ", pitch_calc/toRad, " must be: ", pitch/toRad), quote = FALSE)
