# The take away from this code: 
#    you need to be able to compute a confidence interval

x = c(35.84, 35.81, 35.62, 35.88, 
      35.95, 35.31, 36.15, 35.62, 
      35.94, 36.42, 35.86, 36.18, 
      36.51, 36.65, 36.15, 36.36)

# Get the necessary components:
n     = length(x)
xBar  = mean(x)
s_X   = sd(x)

alpha = 1-97.2/100
qt(1-alpha/2,n-1)

xBar - qt(1-alpha/2,n-1)*s_X/sqrt(n)
xBar + qt(1-alpha/2,n-1)*s_X/sqrt(n)

###### or alternatively
alpha = 1-97.2/100
t.test(x,conf.level = 1 - alpha)$conf.int
