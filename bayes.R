help(persp) # help

? persp # help

rnorm(100) # Geneates 100 random numbers

source("c:/work/test.R")

sink("c:/work/log_150605.txt") # start a log
sink() # stop a log

help(sqrt)

# SQRT Example
require(stats) # for spline
require(graphics)
xx <- -9:9
plot(xx, sqrt(abs(xx)),  col = "red")
lines(spline(xx, sqrt(abs(xx)), n=101), col = "pink")

sqrt(-4) # Error
sqrt(-4+0i)

x = seq(from=0, to =1, by=0.01) # make sequence
x = seq(from=0, to=1, length=101)
rm(x)
rep(2, 5) # replication "2" 5 times
x = c(1:5); z=x>3 # Logical Vector
sum(z) # Summerize Logical Vector
# NA : not available
# NaN : Not a number
y=c(-3, -2, -1,0,1,2)
y[y<0] = -y[y<0] # make negative value to positive

z=array(1:20, dim=c(4,5)) # array, a kind of 2nd vector
z[2,3]
i=array(c(1,2,3,3,2,1), dim=c(3,2)) # array
z[i] # z[1,3], z[2,2], z[3,1]

matrix(c(1:8), 2, 4) # matrix(data, nrow, ncole)
matrix(3, 2, 4) # matrix of 3, row is 2, column is 4


# Column Bind(Cbind()) and Row Bind(rbind())
x = c(1, 2, 3); y = c(4, 5, 6)
A = cbind(x, y)
A[1, 2]
B = rbind(x, y)
B[1, 3]

# Outer Product 
x = c(1, 2, 3); y = c(1, 2, 3)
x%o%y # Outer Product
outer(x, y, "*") # Outer Product
func = function(x,y)
  return(cos(y)/(1+x^2)) # func definition
z = outer(x, y, func)

# matrix transpose : t()
A = matrix(0, 3, 2)
t(A)

# operation(*) and matrix operation(%*%)
A = matrix(1, 3, 2)
B = matrix(2, 3, 2)
C = matrix(2, 2, 2)
A*B
A%*%C

# List = list(name=l=object_1, ..., name_m=object_m)
Lst = list(name = "Sam", age = 34, n.child=2)
Lst[2] # age
Lst[[2]] # 34, same as Lst$age
Lst$age

# Read Data File
x = scan()
setwd("c:/work")
data = scan("input.txt", list("", 0, 0))
# format of "input.txt"
# 1
# 2
# 3
# 4
A = matrix(scan("input.txt"), ncol=4, byrow=TRUE)
data = read.table("input.txt", header=T, sep=",")
data = read.csv("input.csv", header=T, sep=",")

#Probabilit Distribution
dnorm(0, mean = 0, sd = 1) # N(0, 1), 0 to density function
pnorm(0, mean = 0, sd = 1) # N(0, 1), 0 to aggregated density function
qnorm(0.5, mean = 0, sd = 1) # N(0, 1), 50% fractile
rnorm(5, mean=0, sd=1)

# if condition
x = 1:4
if(sum(x)>12) {
  y = 1
  z = 0
} else {
  y = 2
  z = 1
}
# for statement
s = 0; p = 1
for(i in 1:10) {
  s = s+i
  p = p*i
}
# while statement
i=1 
while(i <= 10){
  print("Hello, World!")
  i = i+1
}

# funtion definition
fcn = function(x,y)
{
  ans = sqrt(x*x+y*y)
  return(ans)
}
fcn(3,5)
fcn(3,4)

# 2D Graphics
x = seq(from=0, to=2*pi, by=0.01);
y = sin(x); z = cos(x)
plot(x,y)
lines(x,z)
rm(x,y,z)
sin(0); sin(pi/2); sin(pi)
cos(0); cos(pi/2); cos(pi)

# 3D Graphics , Perspective Plots
x = seq(-10, 10, length=30); y = x
f = function(x,y){
  r = sqrt(x^2+y^2);
  r = 10*sin(r)/r
}
z = outer(x,y,f) # Outer Product of Arrarys
persp(x,y,z, theta=30, phi=30, expand=0.5)
?persp

# multi graph
par(mfrow=c(1,3)) # row = 1, column = 3
x = seq(from = -pi, to = pi, by = 0.01)
y = sin(x)
plot(x, y, type = "l", main="sine curve") #type="l" : line
y = cos(x)
plot(x, y, type = "l", main="cosine curve") #type="l" : line
x = seq(from=-1.5, to=1.5, by=0.01)
y = tan(x)
plot(x, y, type = "l", main="tangant curve") #type="l" : line

# save the graph(png)
setwd("c:/work")
png(file="tri.png", width=400, height=350)
par(mfrow=c(1,3)) # row = 1, column = 3
x = seq(from = -pi, to = pi, by = 0.01)
y = sin(x)
plot(x, y, type = "l", main="sine curve") #type="l" : line
y = cos(x)
plot(x, y, type = "l", main="cosine curve") #type="l" : line
x = seq(from=-1.5, to=1.5, by=0.01)
y = tan(x)
plot(x, y, type = "l", main="tangant curve") #type="l" : line
dev.off()

# save the graph(pdf)
setwd("c:/work")
pdf(file="tri.pdf")
par(mfrow=c(1,3)) # row = 1, column = 3
x = seq(from = -pi, to = pi, by = 0.01)
y = sin(x)
plot(x, y, type = "l", main="sine curve") #type="l" : line
y = cos(x)
plot(x, y, type = "l", main="cosine curve") #type="l" : line
x = seq(from=-1.5, to=1.5, by=0.01)
y = tan(x)
plot(x, y, type = "l", main="tangant curve") #type="l" : line
dev.off()

# save the graph(ps : postscript)
setwd("c:/work")
postscript(file="tri.ps")
par(mfrow=c(1,3)) # row = 1, column = 3
x = seq(from = -pi, to = pi, by = 0.01)
y = sin(x)
plot(x, y, type = "l", main="sine curve") #type="l" : line
y = cos(x)
plot(x, y, type = "l", main="cosine curve") #type="l" : line
x = seq(from=-1.5, to=1.5, by=0.01)
y = tan(x)
plot(x, y, type = "l", main="tangant curve") #type="l" : line
dev.off()

# save the graph(jpeg)
setwd("c:/work")
jpeg(file="tri.jpg", width=400, height=350)
par(mfrow=c(1,3)) # row = 1, column = 3
x = seq(from = -pi, to = pi, by = 0.01)
y = sin(x)
plot(x, y, type = "l", main="sine curve") #type="l" : line
y = cos(x)
plot(x, y, type = "l", main="cosine curve") #type="l" : line
x = seq(from=-1.5, to=1.5, by=0.01)
y = tan(x)
plot(x, y, type = "l", main="tangant curve") #type="l" : line
dev.off()
