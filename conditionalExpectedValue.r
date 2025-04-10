### This code goes through a numeric solution to the exercise in class.  Note that nothing
### in this code is something that you "need" to know for the class; it just may be a helpful
### alternative to the analytical solution we discussed in lecture

## conditional pdf (X | Y = 3)
fX = function(x){
  exp(-x/3)/3
}
integrate(f,1,Inf) # P( X > 1 | Y = 3)
##

##
# Expected Value
##
fX = function(x){
  x * exp(-x/3)/3
}
integrate(fX,0,Inf) 
