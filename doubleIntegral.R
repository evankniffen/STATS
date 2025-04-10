### This is OPTIONAL code for exploring how to do a double integral in R
### It is based on the example on slide 15 of joint/conditional notes
f_XY = function(x,y){
  12/7*(x^2 + x*y)
}

doubleIntegralF = function(xLower,xUpper,yLower,yUpper){
  innerF = function(y){sapply(y, function(z) { integrate(f_XY, xLower, xUpper,z)$value })}
  integrate(innerF,yLower,yUpper)$value
}
doubleIntegralF(0,.6,0,.4)

#Note: this won't work for triangular regions... sorry!

fY = function(y){y*((2*.3 + 4*y)/(2*.3 + 2))}

integrate(fY,0,1)
