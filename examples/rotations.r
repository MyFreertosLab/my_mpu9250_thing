toRad=pi/180
roll=190*toRad
pitch=-90*toRad

# Range Roll [-pi, pi]
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
Rz <- matrix(c(0,0,1,0,1,0,-1,0,0), nrow=3, ncol=3, byrow=TRUE)

Uz=as.matrix(c(0,0,1))

BUz <- R%*%(Uz/norm(Uz, "2"))

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



