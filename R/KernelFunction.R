########################
##                    ##
##  Kernel functions  ##
##                    ##
########################

### Epanechnikov kerenl

K2_Epanechnikov <- function(u)
{
  3/4*(1-u^2)*(abs(u)<=1)
}

dK2_Epanechnikov <- function(u)
{
  -3/2*u*(abs(u)<=1)
}

d2K2_Epanechnikov <- function(u)
{
  -3/2*(abs(u)<=1)
}

### 2nd-order biweight kerenl

K2_Biweight <- function(u)
{
  15/16*(1-u^2)^2*(abs(u)<=1)
}

dK2_Biweight <- function(u)
{
  -15/4*u*(1-u^2)*(abs(u)<=1)
}

d2K2_Biweight <- function(u)
{
  -15/4*(1-3*u^2)*(abs(u)<=1)
}

### 4th-order biweight kerenl

K4_Biweight <- function(u)
{
  (105/64)*(1-3*(u^2))*((1-u^2)^2)*(abs(u)<=1)
}

dK4_Biweight <- function(u)
{
  (105/32)*u*(1-u^2)*(9*u^2-5)*(abs(u)<=1)
}

d2K4_Biweight <- function(u)
{
  (105/32)*(-5+42*u^2-45*u^4)*(abs(u)<=1)
}

### Gaussian kernel

K2_Gaussian <- function(u)
{
  dnorm(u,0,1)
}

dK2_Gaussian <- function(u)
{
  -u*dnorm(u,0,1)
}

d2K2_Gaussian <- function(u)
{
  (u^2-1)*dnorm(u,0,1)
}
