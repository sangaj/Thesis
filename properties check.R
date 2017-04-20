### Study Properties on Simulated Data

## load packages
library(ggplot2)
library(ggmap)
library(truncnorm)
library(dplyr)
library(plyr)
library(lubridate)
library(MASS)
library(corpcor)
library(gdata)
library(Matrix)
library(reshape2)
library(gridExtra)

#-------------------------------------#

set.seed(2016)
timelength <- 120
myts <- strptime("2016-01-01 00:00:00","%Y-%m-%d %H:%M:%S")+300*seq(0,timelength-1)
loc <-  sprintf("%02d", 1:9)
leng <- length(loc)
numtl <- timelength*leng
mass <- rnorm(numtl,100,10)
myts <- rep(myts,each=leng)
loc <-  rep(loc,timelength)
data <-  data.frame(myts,loc,mass)
data <-  arrange(data,myts)



#---------------------------------
# 1.monotonically increasse comparsion
# initialize A,B,Q,R,pi_1,V_1
# 9 grid test
A_9 <- matrix(c(1,1,0,1,1,1,0,1,1),nrow=3,byrow = T) 
zero <- matrix(0,nrow=3,ncol=3)
diag <- diag(1,nrow=3,ncol=3)
A_9.1 <- cbind(A_9,diag,zero)
A_9.2 <- cbind(diag,A_9,zero)
A_9.3 <- cbind(zero,diag,A_9)
#split and order into list
list <- lapply(split(data, data$myts),function(x) arrange(x,loc))
num <- length(list)
# initialized value through greedy algorithm
newA <- rbind(A_9.1,A_9.2,A_9.3)
newQ <- diag(1,nrow=leng,ncol=leng)
newR <- diag(0.01,nrow=leng,ncol=leng)
pi1 <- matrix(aggregate(data$mass, list(data$loc), mean)[,2])
V1 <-  diag(0.1,nrow=leng,ncol=leng)
jj <- lapply (1 : nrow(newA), function (x) diag(newA[x,]))

# EM settings
max.iter <- 15
tol <- 0.001
likelihood <-  matrix(0, max.iter, 1)
A <- Q <- R <- Vv1 <- array(NA, dim = c(leng, leng, max.iter))
cvg <-  1 + tol


# 1.1 greedy EM algorithm 
for (iter in 1:max.iter) {
  y.updated <- pi1
  V.updated <- V1
  K <- array(NA, dim = c(leng, leng, num))
  y <- ynew <- array(NA, dim = c(leng, 1, num))
  V <- Vnew <- array(NA, dim = c(leng, leng, num))
  y[,,1] <- ynew[,,1] <- y.updated
  V[,,1] <- Vnew[,,1] <- V.updated
  
  for (i in 2:num){
    y.new <- newA %*% y.updated 
    V.new <- newA %*% V.updated %*% t(newA) + newQ
    upperTriangle(V.new) <- lowerTriangle(V.new, diag=FALSE,byrow=TRUE)
    Kn <- V.new %*% solve(newR + V.new)
    eye <- diag(1,leng)
    z <- matrix(list[[i]]$mass)
    y.updated <- y.new + Kn %*% ( z - y.new)
    V.updated <- (eye - Kn) %*% V.new %*% t(eye - Kn) + Kn %*% newR %*% t(Kn)
    K[, ,i] <- Kn
    y[, ,i] <- y.updated
    V[, ,i] <- V.updated
    Vnew[, ,i] <- V.new 
    ynew[,,i] <- y.new
  }
  
  xs <-  array(NA, dim = c(leng, 1, num))
  Ps <-  array(NA, dim = c(leng, leng, num))
  Pcs <- array(NA, dim = c(leng, leng, num))
  J <-  array(NA, dim = c(leng, leng, num))
  xs[, , num] = y[, , num]
  Ps[, , num] = V[, , num]
  eye <- diag(1, leng)
  for (k in num:2) {
    J[, , k - 1] <-  (V[, , k - 1] %*% t(newA)) %*% solve(Vnew[, , k]) 
    xs[, , k - 1] <-  y[, , k - 1] + J[, , k - 1] %*% (xs[, , k] - ynew[,,k])
    Ps[, , k - 1] <-  V[, , k - 1] + J[, , k - 1] %*% (Ps[, , k] - Vnew[, , k]) %*% t(J[, , k - 1])
  }
  
  Pcs[, , num] <-  (eye - K[,,num]) %*% newA %*% V[, , num - 1]
  
  for (k in num:3) {
    Pcs[, , k - 1] <-  V[,,k-1] %*% t(J[,,k-2]) + J[,,k-1] %*% (Pcs[,,k]-newA %*% V[,,k-1]) %*% t(J[,,k-2])
  }
  
  # function for A
  Aa <- array(NA, dim = c(1,ncol(newA), nrow(newA))) 
  sum3 <- sum4  <- 0
  for (i in 1:nrow(newA)){
    for (j in 2:num){
      sum3 <- sum3 + (Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1]))[i,]
      sum4 <- sum4 + Ps[, , j-1] + xs[, , j - 1] %*% t(xs[, , j - 1])
      j <- j + 1
    }
    Aa[,,i] <- sum3 %*% ginv(jj[[i]]%*%sum4)
  }
  newA <- apply(Aa, 2, I)
  
  sum5 <- 0
  for (j in 2:num){
    sum5 <- sum5 + Ps[, , j] + xs[, , j] %*% t(xs[, , j]) + (Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1])) %*% t(newA) 
    + newA %*% t((Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1]))) - newA %*%  (Ps[, , j-1] + xs[, , j-1] %*% t(xs[, , j-1])) %*% t(newA) 
    j <- j + 1
  }
  newQ <- sum5 * (1 - num)^{-1}
  upperTriangle(newQ) <- lowerTriangle(newQ, diag=FALSE,byrow=TRUE)
  newQ <- make.positive.definite(newQ)
  
  sum6 <- 0
  for (j in 1:num){
    zt <-  matrix(list[[j]][,3])
    sum6 <- sum6 + zt %*% t(zt) - zt %*% t(xs[, , j]) - xs[, , j] %*% t(zt) + Ps[, , j] + xs[, , j] %*% t(xs[, , j])
    j <- j + 1
  }
  newR <- sum6 * (num)^{-1}
  upperTriangle(newR) <- lowerTriangle(newR, diag=FALSE,byrow=TRUE)
  newR <- make.positive.definite(newR)
  
  
  pi1 <- pi1 + J[, , 1] %*% (xs[, , 2] - newA %*% y[, , 1])
  
  
  
  # calculate likelihood
  likesum1 <- 0
  
  for (j in 2:num){
    likesum1 <- likesum1 + 0.5 * t((xs[, , j] - newA %*% xs[, , j-1])) %*% solve(newQ) %*% (xs[, , j] - newA %*% xs[, , j-1])
    j <- j + 1
  }
  
  likesum2 <- 0
  for (j in 1:num){
    zt <- matrix(list[[j]][,3])
    likesum2 <- likesum2 + 0.5 * t((zt - xs[, , j])) %*% solve(newR) %*% (zt - xs[, , j])
    j <- j + 1
  }
  
  A[, ,iter] <- newA
  Q[, ,iter] <- newQ
  R[, ,iter] <- newR
  likelihood[iter] <- -0.5 * t((xs[, , 1] - pi1)) %*% solve(V1) %*% (xs[, , 1] - pi1) 
  - likesum1 - likesum2 - 0.5 * log(det(V1)) - 0.5 * num * log(det(newR)) - 0.5 * (num -1) * log(det(newQ)) - num * log(2*pi)
  
  
  if (iter > 1) 
    cvg <-  (likelihood[iter - 1] - likelihood[iter])/abs(likelihood[iter - 1])
  # if (cvg < 0)
  #  stop ("not increase")  
  if (abs(cvg) < tol) 
    break ("converge")
}
likegreedy <- data.frame(likelihood[1:iter])
colnames(likegreedy) <- "likegreedy"

# 1.2 normal EM algorithm
newA <- rbind(A_9.1,A_9.2,A_9.3)
newQ <- diag(1,nrow=leng,ncol=leng)
newR <- diag(0.01,nrow=leng,ncol=leng)
pi1 <- matrix(aggregate(data$mass, list(data$loc), mean)[,2])
V1 <-  diag(0.1,nrow=leng,ncol=leng)

# EM settings
max.iter <- 15
tol <- 0.001
likelihood <-  matrix(0, max.iter, 1)
A <- Q <- R <- Vv1 <- array(NA, dim = c(leng, leng, max.iter))
cvg <-  1 + tol

# kalman filter
for (iter in 1:max.iter) {
  y.updated <- pi1
  V.updated <- V1
  K <- array(NA, dim = c(leng, leng, num))
  y <- ynew <- array(NA, dim = c(leng, 1, num))
  V <- Vnew <- array(NA, dim = c(leng, leng, num))
  y[,,1] <- ynew[,,1] <- y.updated
  V[,,1] <- Vnew[,,1] <- V.updated
  
  for (i in 2:num){
    y.new <- newA %*% y.updated 
    V.new <- newA %*% V.updated %*% t(newA) + newQ
    upperTriangle(V.new) <- lowerTriangle(V.new, diag=FALSE,byrow=TRUE)
    Kn <- V.new %*% solve(newR + V.new)
    eye <- diag(1,leng)
    z <- matrix(list[[i]]$mass)
    y.updated <- y.new + Kn %*% ( z - y.new)
    V.updated <- (eye - Kn) %*% V.new %*% t(eye - Kn) + Kn %*% newR %*% t(Kn)
    K[, ,i] <- Kn
    y[, ,i] <- y.updated
    V[, ,i] <- V.updated
    Vnew[, ,i] <- V.new 
    ynew[,,i] <- y.new
  }
  # Kalman smoother
  xs <-  array(NA, dim = c(leng, 1, num))
  Ps <-  array(NA, dim = c(leng, leng, num))
  Pcs <- array(NA, dim = c(leng, leng, num))
  J <-  array(NA, dim = c(leng, leng, num))
  xs[, , num] = y[, , num]
  Ps[, , num] = V[, , num]
  eye <- diag(1, leng)
  for (k in num:2) {
    J[, , k - 1] <-  (V[, , k - 1] %*% t(newA)) %*% solve(Vnew[, , k]) 
    xs[, , k - 1] <-  y[, , k - 1] + J[, , k - 1] %*% (xs[, , k] - ynew[,,k])
    Ps[, , k - 1] <-  V[, , k - 1] + J[, , k - 1] %*% (Ps[, , k] - Vnew[, , k]) %*% t(J[, , k - 1])
  }
  
  Pcs[, , num] <-  (eye - K[,,num]) %*% newA %*% V[, , num - 1]
  
  for (k in num:3) {
    Pcs[, , k - 1] <-  V[,,k-1] %*% t(J[,,k-2]) + J[,,k-1] %*% (Pcs[,,k]-newA %*% V[,,k-1]) %*% t(J[,,k-2])
  }
  
  # function for A
  sum3 <- sum4  <- 0
  for (j in 2:num){
    sum3 <- sum3 + (Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1]))
    sum4 <- sum4 + Ps[, , j-1] + xs[, , j - 1] %*% t(xs[, , j - 1])
    j <- j + 1
  }
  newA <- sum3 %*% solve(sum4)
  
  
  sum5 <- 0
  for (j in 2:num){
    sum5 <- sum5 + Ps[, , j] + xs[, , j] %*% t(xs[, , j]) + (Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1])) %*% t(newA) 
    + newA %*% t((Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1]))) - newA %*%  (Ps[, , j-1] + xs[, , j-1] %*% t(xs[, , j-1])) %*% t(newA) 
    j <- j + 1
  }
  newQ <- sum5 * (1 - num)^{-1}
  upperTriangle(newQ) <- lowerTriangle(newQ, diag=FALSE,byrow=TRUE)
  newQ <- make.positive.definite(newQ)
  
  sum6 <- 0
  for (j in 1:num){
    zt <-  matrix(list[[j]][,3])
    sum6 <- sum6 + zt %*% t(zt) - zt %*% t(xs[, , j]) - xs[, , j] %*% t(zt) + Ps[, , j] + xs[, , j] %*% t(xs[, , j])
    j <- j + 1
  }
  newR <- sum6 * (num)^{-1}
  upperTriangle(newR) <- lowerTriangle(newR, diag=FALSE,byrow=TRUE)
  newR <- make.positive.definite(newR)
  
  
  pi1 <- pi1 + J[, , 1] %*% (xs[, , 2] - newA %*% y[, , 1])
  
  
  # calculate likelihood
  likesum1 <- 0
  
  for (j in 2:num){
    likesum1 <- likesum1 + 0.5 * t((xs[, , j] - newA %*% xs[, , j-1])) %*% solve(newQ) %*% (xs[, , j] - newA %*% xs[, , j-1])
    j <- j + 1
  }
  
  likesum2 <- 0
  for (j in 1:num){
    zt <- matrix(list[[j]][,3])
    likesum2 <- likesum2 + 0.5 * t((zt - xs[, , j])) %*% solve(newR) %*% (zt - xs[, , j])
    j <- j + 1
  }
  
  A[, ,iter] <- newA
  Q[, ,iter] <- newQ
  R[, ,iter] <- newR
  likelihood[iter] <- -0.5 * t((xs[, , 1] - pi1)) %*% solve(V1) %*% (xs[, , 1] - pi1) 
  - likesum1 - likesum2 - 0.5 * log(det(V1)) - 0.5 * num * log(det(newR)) - 0.5 * (num -1) * log(det(newQ)) - num * log(2*pi)
  
  
  if (iter > 1) 
    cvg <-  (likelihood[iter - 1] - likelihood[iter])/abs(likelihood[iter - 1])
  # if (cvg < 0)
  #  stop ("not increase")  
  if (abs(cvg) < tol) 
    break ("converge")
}

# plot it

likenormal <- data.frame(likelihood[1:iter])
colnames(likenormal) <- "likenormal"
like <- cbind(likegreedy,likenormal)
like$ID<-seq.int(nrow(like))

p1 <- ggplot(data=like, aes(x=ID, y=likenormal)) + geom_line(color="blue") + scale_x_continuous(breaks = 1:15) + labs(x="Iteration",y="Normal LogLikelihood",title="Comparison of Two Loglikelihood in Simulated data")  + theme(axis.line = element_line(size=1, colour = "black"),
                                                                                                                                                                                                                                 panel.grid.major = element_line(colour = "#d3d3d3"),
                                                                                                                                                                                                                                 panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                 panel.border = element_blank(), panel.background = element_blank(),
                                                                                                                                                                                                                                 axis.text.x=element_text(colour="black", size = 9),
                                                                                                                                                                                                                                 axis.text.y=element_text(colour="black", size = 9),
                                                                                                                                                                                                                                 legend.title = element_blank())
p2 <- ggplot(data=like, aes(x=ID, y=likegreedy)) +coord_cartesian(ylim=c(-700000, 0)) + geom_line(color="green") + scale_x_continuous(breaks = 1:15)+ labs(x="Iteration",y="Greedy LogLikelihood")  + theme(axis.line = element_line(size=1, colour = "black"),
                                                                                                                                                                                                            panel.grid.major = element_line(colour = "#d3d3d3"),
                                                                                                                                                                                                            panel.grid.minor = element_blank(),
                                                                                                                                                                                                            panel.border = element_blank(), panel.background = element_blank(),
                                                                                                                                                                                                            axis.text.x=element_text(colour="black", size = 9),
                                                                                                                                                                                                            axis.text.y=element_text(colour="black", size = 9),
                                                                                                                                                                                                            legend.title = element_blank())

# Multiple plot function
# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
multiplot(p1, p2, layout=matrix(c(1,1,2,2), nrow=2, byrow=TRUE))


##-----------------------------------------------

# initialization comparsion

A_9 <- matrix(c(1,1,0,1,1,1,0,1,1),nrow=3,byrow = T) 
zero <- matrix(0,nrow=3,ncol=3)
diag <- diag(1,nrow=3,ncol=3)
A_9.1 <- cbind(A_9,diag,zero)
A_9.2 <- cbind(diag,A_9,zero)
A_9.3 <- cbind(zero,diag,A_9)
#split and order into list
list <- lapply(split(data, data$myts),function(x) arrange(x,loc))
num <- length(list)
# initialized value through greedy algorithm

invalues <- c(0.01,0.1,0.2,0.5,1)
maxlike <- maxsteps <- numeric(0)
for (kk in 1:length(invalues)){
  newA <- rbind(A_9.1,A_9.2,A_9.3)
  newQ <- diag(aggregate(data$mass, list(data$loc), var)[,2])
  newR <- diag(invalues[kk],nrow=leng,ncol=leng)
  pi1 <- matrix(aggregate(data$mass, list(data$loc), mean)[,2])
  #pi1 <- matrix(10,nrow=9,ncol=1)
  V1 <-  diag(1, nrow=leng, ncol=leng)
  jj <- lapply (1 : nrow(newA), function (x) diag(newA[x,]))
  
  # EM settings
  max.iter <- 50
  tol <- 0.001
  likelihood <-  matrix(0, max.iter, 1)
  A <- Q <- R <- Vv1 <- array(NA, dim = c(leng, leng, max.iter))
  cvg <-  1 + tol
  # kalman filter
  for (iter in 1:max.iter) {
    y.updated <- pi1
    V.updated <- V1
    K <- array(NA, dim = c(leng, leng, num))
    y <- ynew <- array(NA, dim = c(leng, 1, num))
    V <- Vnew <- array(NA, dim = c(leng, leng, num))
    y[,,1] <- ynew[,,1] <- y.updated
    V[,,1] <- Vnew[,,1] <- V.updated
    
    for (i in 2:num){
      y.new <- newA %*% y.updated 
      V.new <- newA %*% V.updated %*% t(newA) + newQ
      upperTriangle(V.new) <- lowerTriangle(V.new, diag=FALSE,byrow=TRUE)
      Kn <- V.new %*% solve(newR + V.new)
      eye <- diag(1,leng)
      z <- matrix(list[[i]]$mass)
      y.updated <- y.new + Kn %*% ( z - y.new)
      V.updated <- (eye - Kn) %*% V.new %*% t(eye - Kn) + Kn %*% newR %*% t(Kn)
      K[, ,i] <- Kn
      y[, ,i] <- y.updated
      V[, ,i] <- V.updated
      Vnew[, ,i] <- make.positive.definite(V.new) 
      ynew[,,i] <- y.new
    }
    # kalman smoother
    xs <-  array(NA, dim = c(leng, 1, num))
    Ps <-  array(NA, dim = c(leng, leng, num))
    Pcs <- array(NA, dim = c(leng, leng, num))
    J <-  array(NA, dim = c(leng, leng, num))
    xs[, , num] = y[, , num]
    Ps[, , num] = V[, , num]
    eye <- diag(1, leng)
    for (k in num:2) {
      J[, , k - 1] <-  (V[, , k - 1] %*% t(newA)) %*% solve(Vnew[, , k]) 
      xs[, , k - 1] <-  y[, , k - 1] + J[, , k - 1] %*% (xs[, , k] - ynew[,,k])
      Ps[, , k - 1] <-  V[, , k - 1] + J[, , k - 1] %*% (Ps[, , k] - Vnew[, , k]) %*% t(J[, , k - 1])
    }
    
    Pcs[, , num] <-  (eye - K[,,num]) %*% newA %*% V[, , num - 1]
    
    for (k in num:3) {
      Pcs[, , k - 1] <-  V[,,k-1] %*% t(J[,,k-2]) + J[,,k-1] %*% (Pcs[,,k]-newA %*% V[,,k-1]) %*% t(J[,,k-2])
    }
    
    # function for A
    Aa <- array(NA, dim = c(1,ncol(newA), nrow(newA))) 
    sum3 <- sum4  <- 0
    for (i in 1:nrow(newA)){
      for (j in 2:num){
        sum3 <- sum3 + (Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1]))[i,]
        sum4 <- sum4 + Ps[, , j-1] + xs[, , j - 1] %*% t(xs[, , j - 1])
        j <- j + 1
      }
      Aa[,,i] <- sum3 %*% tryCatch(solve(jj[[i]]%*%sum4), error=function(e) ginv(jj[[i]]%*%sum4))
    }
    newA <- apply(Aa, 2, I)
    
    sum5 <- 0
    for (j in 2:num){
      sum5 <- sum5 + Ps[, , j] + xs[, , j] %*% t(xs[, , j]) + (Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1])) %*% t(newA) 
      + newA %*% t((Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1]))) - newA %*%  (Ps[, , j-1] + xs[, , j-1] %*% t(xs[, , j-1])) %*% t(newA) 
      j <- j + 1
    }
    newQ <- sum5 * (1 - num)^{-1}
    upperTriangle(newQ) <- lowerTriangle(newQ, diag=FALSE,byrow=TRUE)
    newQ <- make.positive.definite(newQ)
    
    sum6 <- 0
    for (j in 1:num){
      zt <-  matrix(list[[j]][,3])
      sum6 <- sum6 + zt %*% t(zt) - zt %*% t(xs[, , j]) - xs[, , j] %*% t(zt) + Ps[, , j] + xs[, , j] %*% t(xs[, , j])
      j <- j + 1
    }
    newR <- sum6 * (num)^{-1}
    upperTriangle(newR) <- lowerTriangle(newR, diag=FALSE,byrow=TRUE)
    newR <- make.positive.definite(newR)
    
    
    pi1 <- pi1 + J[, , 1] %*% (xs[, , 2] - newA %*% y[, , 1])
    
    
    
    # calculate likelihood
    likesum1 <- 0
    
    for (j in 2:num){
      likesum1 <- likesum1 + 0.5 * t((xs[, , j] - newA %*% xs[, , j-1])) %*% solve(make.positive.definite(newQ)) %*% (xs[, , j] - newA %*% xs[, , j-1])
      j <- j + 1
    }
    
    likesum2 <- 0
    for (j in 1:num){
      zt <- matrix(list[[j]][,3])
      likesum2 <- likesum2 + 0.5 * t((zt - xs[, , j])) %*% solve(newR) %*% (zt - xs[, , j])
      j <- j + 1
    }
    
    A[, ,iter] <- newA
    Q[, ,iter] <- newQ
    R[, ,iter] <- newR
    likelihood[iter] <- -0.5 * t((xs[, , 1] - pi1)) %*% solve(V1) %*% (xs[, , 1] - pi1) 
    - likesum1 - likesum2 - 0.5 * log(det(V1)) - 0.5 * num * log(det(newR)) - 0.5 * (num -1) * log(det(newQ)) - num * log(2*pi)
    
    
    if (iter > 1) 
      cvg <-  likelihood[iter] - likelihood[iter - 1]
    # if (cvg < 0)
    #  stop ("not increase")  
    if (abs(cvg) < tol) 
      break ("converge")
  }
  maxlike[kk] <- max(likelihood)
  maxsteps[kk] <- which.max(likelihood)
}