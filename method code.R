### Methods application on Schiphol sample data

# loading Package
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
library(stringr)
library(reshape2)
library(caret)

## Data Preprocessing

# join two tables and extract the information we need
track <- read.table("track.csv",header=T,sep=",")
trackestimate <- read.table("trackestimate.csv",header=T,sep=",")
trackestimate$timestamp <- ymd_hms(substr(trackestimate$timestamp,1,19))
track <- subset(track,classification_id != 1,select=c(id,classification_id,distance_travelled,score))
colnames(track)[1] <- "track_id"
trackestimate$st_astext <- gsub("POINT Z ","",trackestimate$st_astext) 
trackestimate$st_astext <- gsub("[()]","",trackestimate$st_astext) 
location <- as.data.frame(str_split_fixed(trackestimate$st_astext, " ", 3))
names(location) <- c("longitude", "latitude", "altitude")
trackestimate <- cbind(subset(trackestimate, select=-st_astext),subset(location, select=c(latitude,longitude)))
bird <- join(trackestimate,track,by="track_id")
bird <- subset(bird,select=-c(id,track_id,classification_id))
bird[is.na(bird)] <- 0

# define grid
## runway coordinates
l_lon <-  4.706
r_lat <-  52.368
r_lon <-  4.717
l_lat <-  52.325
## grid numbers and assign each location into cells
rowN <- 6
colN <- 4
cellN <- rowN*colN
lat_det <-  2/111
lon_det <-  2/(cos(r_lat*pi/180)*111.321)
bound_x <-  c(l_lon-lon_det,l_lat-lat_det)
bound_y <-  c(r_lon+lon_det,r_lat+lat_det)
gridx <- seq(l_lon-lon_det,r_lon+lon_det,length.out=colN+1)
gridy <- seq(l_lat-lat_det,r_lat+lat_det,length.out=rowN+1)
bird$latitude <- as.numeric(as.character(bird$latitude))
bird$longitude <- as.numeric(as.character(bird$longitude))
bird$yy <- findInterval(bird$latitude, gridy, rightmost.closed = TRUE)
bird$xx <- findInterval(bird$longitude, gridx, rightmost.closed = TRUE)
bird <- subset(bird, !(xx %in% c(0,colN+1)))
bird <- subset(bird, !(yy %in% c(0,rowN+1)))
for( i in 1:rowN-1){
  bird$yy <-  replace(bird$yy,bird$yy== i, 2*rowN-i)
}
bird$xy <- with(bird, paste(yy,xx, sep=""))
bird$xy <- as.numeric(as.character(bird$xy))
cell <- data.frame(xy=unique(bird$xy))
cell$xy <- as.numeric(as.character(cell$xy))
cell <- arrange(cell,xy)
cell$cell<-seq.int(nrow(cell))
birds <- left_join(bird,cell,by="xy")
birds <- subset(birds,select=-c(yy,xx,xy))
# aggregate each variable into minutes
birds$timestamp <- trunc(birds$timestamp, "min")
birds$timestamp <- ymd_hms(birds$timestamp)
#birds_melt <- melt(birds, id = c("timestamp","cell"))
#birds <- dcast(birds_melt, timestamp + cell ~ variable, mean) 
birds <- recast(birds, timestamp + cell ~ variable, mean, id.var = c("timestamp", "cell"))

# clean weather data
weather <- read.table("weather_2016.csv",sep=",",head=T)
weatherdate <- data.frame(do.call(rbind, strsplit(as.vector(weather$DATE.LT), split = "-"))) %>% mutate(date=ymd(paste(X3,X2,X1)))
weather <- cbind.data.frame(weather,ten=weatherdate$date) %>% mutate(ten= paste(ten,TIME.LT)) 

start <- ymd_hms("2016-04-07 00:00:00")
end <- ymd_hms("2016-04-08 23:59:59")
diffday <- as.numeric(round(difftime(end,start,units="days")))
diffmin <- as.numeric(round(difftime(end,start,units="mins")))
timestamp <- rep(seq(from=start, by=60, to=end),each=cellN)
ten <-cut(timestamp,breaks=6*24*diffday) 
times <- rep(c(1:diffmin),each=cellN)
cell <- rep(c(1:cellN),times=diffmin)
Hours <- cbind.data.frame(timestamp,cell,ten)

weather.join <- join(Hours,weather,by=c("ten"),type="left") %>% subset(select=-c(TIME.LT,DATE.LT))
Bird <- join(weather.join,birds,by=c("timestamp","cell"),type="left")

# only extract the active time period 
Birds <- subset(Bird,lubridate::hour(Bird$timestamp)==10 | lubridate::hour(Bird$timestamp)==11,select=-ten)
Birds[is.na(Birds)] <- 0
save.image(file="data.RData")
load("data.RData")



# features reduction 
library(factoextra)
library(FactoMineR)
# change column order and drop same value column 
Birds <- subset(Birds,select=-c(X.VIS,X.GML,X.SHWR,heading_vertical)) %>% dplyr::select(timestamp, cell, mass, X.CEIL:peak_mass, mass_correction:score)
# divided into train and test set 
train <- subset(Birds,lubridate::day(Birds$timestamp)==7)
test <- subset(Birds,lubridate::day(Birds$timestamp)==8)
# drop same value column in train set and then drop same features in test set. 
train <- subset(train, select=-c(X.CEIL,X.WDIR,X.WSPD,X.WGUS))
test <- subset(test, select=-c(X.CEIL,X.WDIR,X.WSPD,X.WGUS))
# PCA
res.pca <- PCA(train[,-c(1:2)], graph =FALSE,ncp=10,scale.unit=TRUE,quanti.sup=1)
fviz_screeplot(res.pca, ncp=10)
fviz_contrib(res.pca, choice = "var", axes = 1:5,sort.val="desc") # variance>90
fviz_pca_var(res.pca, col.var="cos2") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=0.5) + theme_minimal()


# keep lon,lat,heading,mass_correction,velocity interms of 90% explained variation

train <- subset(train,select=c(timestamp,cell,mass,latitude,longitude,heading,velocity,mass_correction))
test <-  subset(test,select=c(timestamp,cell,mass,latitude,longitude,heading,velocity,mass_correction))

# initialize A,B,Q,R,pi_1,V_1
subA <- matrix(c(1,1,0,0,1,1,1,0,0,1,1,1,0,0,1,1),nrow=4,byrow = T) 
zeroM <- matrix(0,nrow=4,ncol=4,byrow=T)


# 24 grid setting
A.1 <-  cbind(matrix(rep(subA,2), nrow = nrow(subA) , byrow = F),matrix(rep(zeroM,4), nrow = nrow(zeroM) , byrow = F))
A.2 <-  cbind(matrix(rep(subA,3), nrow = nrow(subA) , byrow = F),matrix(rep(zeroM,3), nrow = nrow(zeroM) , byrow = F))
A.3 <-  cbind(zeroM,matrix(rep(subA,3), nrow = nrow(subA) , byrow = F),matrix(rep(zeroM,2), nrow = nrow(zeroM) , byrow = F))
A.4 <-  cbind(matrix(rep(zeroM,2), nrow = nrow(zeroM) , byrow=F),matrix(rep(subA,3), nrow = nrow(subA) , byrow = F) , zeroM)
A.5 <-  cbind(matrix(rep(zeroM,3), nrow = nrow(zeroM) , byrow=F),matrix(rep(subA,3), nrow = nrow(subA) , byrow = F))
A.6 <-  cbind(matrix(rep(zeroM,4), nrow = nrow(zeroM) , byrow = F),matrix(rep(subA,2), nrow = nrow(subA) , byrow = F))


# assign different initial value
value <- c(0.0001,0.0002,0.0005,0.001,0.002)
RMSE1 <- RMSE2 <-  numeric(0)
MAE1 <- MAE2 <- numeric(0)
likely <- numeric(0)

### Kalman Model
for (kkk in 1:1){
    subtrain <- train
    procValues <- preProcess(subtrain[,-c(1:2)], method = c("center", "scale"))
    scaledTrain <-  predict(procValues, subtrain)
    scaledTest <-  predict(procValues, test)
    list <- lapply(split(scaledTrain, scaledTrain$timestamp),function(x) arrange(x,cell))
    num <- length(list)
    featureN <- ncol(scaledTrain)-3
    
    # initialized value
    
    leng <- length(unique(scaledTrain$cell))
    newA <- rbind(A.1,A.2,A.3,A.4,A.5,A.6)
    newB <- do.call(cbind, replicate(featureN, newA, simplify=FALSE))
    #newQ <- diag(0.05,nrow=leng,ncol=leng)
    newQ  <-  diag(aggregate(scaledTrain$mass, list(scaledTrain$cell), var)[,2])
    
    newR <- diag(value[kkk],nrow=leng,ncol=leng)
    pi1 <- matrix(aggregate(scaledTrain$mass, list(scaledTrain$cell), mean)[,2])
    #pi1 <- matrix(list[[1]]$mass)
    #V1 <-  diag(aggregate(scaledTrain$mass, list(scaledTrain$cell), var)[,2])
    V1 <- diag(0.1,nrow=leng,ncol=leng)
    
    # indicator matrix
    
    jj <- lapply (1 : nrow(newA), function (x) diag(newA[x,]))
    
    kk <- lapply (1 : nrow(newB), function (x) diag(newB[x,]))
    
    # EM settings
    max.iter <- 50
    tol <- 0.001
    likelihood <-  matrix(0, max.iter, 1)
    A <- Q <- R <- Vv1 <- array(NA, dim = c(leng, leng, max.iter))
    Pi1<- array(NA, dim = c(leng, 1, max.iter))
    B <- array(NA, dim = c(leng, leng * featureN, max.iter))
    cvg <-  1 + tol
    
    # kalman filter and smoother
    for (iter in 1:max.iter) {
      y.updated <- pi1
      V.updated <- V1
      K <- array(0, dim = c(leng, leng, num))
      y <- ynew <- array(0, dim = c(leng, 1, num))
      V <- Vnew <- array(0, dim = c(leng, leng, num))
      y[,,1] <- ynew[,,1] <- y.updated
      V[,,1] <- Vnew[,,1] <- V.updated
      
      for (i in 2:num){
        u <-  matrix(stack(list[[i]][,-(1:3)])[,1])
        y.new <- newA %*% y.updated + newB %*% u
        V.new <- newA %*% V.updated %*% t(newA) + newQ
        upperTriangle(V.new) <- lowerTriangle(V.new, diag=FALSE,byrow=TRUE)
        Kn <- V.new %*% solve(newR + V.new)
        eye <- diag(1,leng)
        z <- matrix(list[[i]]$mass)
        y.updated <- y.new + Kn %*% ( z - y.new)
        V.updated <- (eye - Kn) %*% V.new %*% t(eye - Kn) + Kn %*% newR %*% t(Kn)
        y[, ,i] <- y.updated
        V[, ,i] <- V.updated
        ynew[,,i] <- y.new
        Vnew[,,i] <- V.new
        K[, ,i] <- Kn
      }
      
      
      xs <-  array(NA, dim = c(leng, 1, num))
      Ps <-  array(NA, dim = c(leng, leng, num))
      Pcs <- array(NA, dim = c(leng, leng, num))
      J <-  array(NA, dim = c(leng, leng, num))
      xs[, , num] = y[, , num]
      Ps[, , num] = V[, , num]
      eye <- diag(1, leng)
      for (k in num:2) {
        u <-  matrix(stack(list[[k-1]][,-(1:3)])[,1])
        J[, , k - 1] <-   V[, , k - 1] %*% t(newA) %*% solve(Vnew[, , k]) 
        xs[, , k - 1] <-  y[, , k - 1] + J[, , k - 1] %*% (xs[, , k] - ynew[,,k])
        Ps[, , k - 1] <-  V[, , k - 1] + J[, , k - 1] %*% (Ps[, , k] - Vnew[, , k]) %*% t(J[, , k - 1])
      }
      Pcs[, , num] <-  (eye - K[,,num]) %*% newA %*% V[, , num - 1]
      
      for (k in num:3) {
        Pcs[, , k - 1] <-  V[,,k-1] %*% t(J[,,k-2]) + J[,,k-1] %*% (Pcs[,,k]-newA %*% V[,,k-1]) %*% t(J[,,k-2])
      }
      
      # equation for B
      Bb <- array(NA, dim = c(1,ncol(newB), nrow(newB))) 
      sum1 <- sum2  <- 0
      for(i in 1:nrow(newB)){
        for (j in 2:num){
          u <- matrix(stack(list[[j-1]][,-(1:3)])[,1])
          #sum1 <- sum1 + t(u) %*% u
          sum1 <- sum1 + u %*% t(u)
          sum2 <- sum2 + xs[, , j][i]%*%t(kk[[i]] %*% u) - newA[i,] %*% jj[[i]] %*% xs[, , j-1]%*% t(kk[[i]] %*% u)
          j <- j + 1
        }
        Bb[,,i] <- sum2 %*% ginv( kk[[i]] %*% sum1 %*% t(kk[[i]])) #psudeo inverse
      }
      newB <- apply(Bb, 2, I) 
      
      # equation for A
      Aa <- array(NA, dim = c(1,ncol(newA), nrow(newA))) 
      sum3 <- sum4  <- 0
      for (i in 1:nrow(newA)){
        for (j in 2:num){
          u <-  matrix(stack(list[[j-1]][,-(1:3)])[,1])
          sum3 <- sum3 + (Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1]))[i,]  - newB[i,] %*% kk[[i]] u %*% t(xs[, , j-1])
          sum4 <- sum4 + Ps[, , j-1] + xs[, , j - 1] %*% t(xs[, , j - 1])
          j <- j + 1
        }
        Aa[,,i] <- sum3 %*% ginv(jj[[i]]%*%sum4)  #psudeo inverse
      }
      newA <- apply(Aa, 2, I)
      
      
      # equation for Q
      sum5 <- 0
      for (j in 2:num){
        u <-  matrix(stack(list[[j-1]][,-(1:3)])[,1])
        sum5 <- sum5 + ( Ps[, , j] + xs[, , j] %*% t(xs[, , j])) - (Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1])) %*% t(newA) -
          xs[, , j] %*% t(u) %*% t(newB) - newA %*% t((Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1]))) + newA %*%  (Ps[, , j-1] + xs[, , j-1] %*% t(xs[, , j-1])) %*% t(newA) +
          newA %*% xs[, , j-1] %*% t(u) %*% t(newB) - newB %*% u %*% t(xs[, , j]) + newB %*% u %*% t(xs[, , j-1]) %*% t(newA) + newB %*% u %*% t(u) %*% t(newB)
        j <- j + 1
      }
      newQ <- sum5 * (num - 1)^{-1}
      upperTriangle(newQ) <- lowerTriangle(newQ, diag=FALSE,byrow=TRUE)
      newQ <- make.positive.definite(newQ)
      
      # equation for R
      sum6 <- 0
      for (j in 1:num){
        zt <-  matrix(list[[j]][,3])
        sum6 <- sum6 + zt %*% t(zt) - zt %*% t(xs[, , j]) - xs[, , j] %*% t(zt) + Ps[, , j] + xs[, , j] %*% t(xs[, , j])
        j <- j + 1
      }
      newR <- sum6 * (num)^{-1}
      upperTriangle(newR) <- lowerTriangle(newR, diag=FALSE,byrow=TRUE)
      newR <- make.positive.definite(newR)
      
      # equation for pi1 and V1
      u1 <-  matrix(stack(list[[1]][,-(1:3)])[,1])
      pi1 <- pi1 + J[, , 1] %*% (xs[, , 2] - newA %*% y[, , 1] - newB %*% u1)
      # if we already know it, don't need to estimate it
      V1 <- V1 + J[, , 1] %*%(Ps[, , 2] - Vnew[, , 2])%*%t(J[, , 1])
      upperTriangle(V1) <- lowerTriangle(V1, diag=FALSE,byrow=TRUE)
      V1 <- make.positive.definite(V1)
      # check
      if (suppressWarnings(is.nan(log(det(V1)))==TRUE)){
        V1 <- make.positive.definite(V1)
      } else{
        V1 <- V1
      }
      
      
      
      # calculate likelihood
      likesum1 <- 0
      
      for (j in 2:num){
        u <- matrix(stack(list[[j-1]][,-(1:3)])[,1])
        likesum1 <- likesum1 + 0.5 * t((xs[, , j] - newA %*% xs[, , j-1] - newB %*% u)) %*% solve(newQ) %*% (xs[, , j] - newA %*% xs[, , j-1] - newB %*% u)
        j <- j + 1
      }
      
      likesum2 <- 0
      for (j in 1:num){
        zt <- matrix(list[[j]][,3])
        likesum2 <- likesum2 + 0.5 * t((zt - xs[, , j])) %*% solve(newR) %*% (zt - xs[, , j])
        j <- j + 1
      }
      
      A[, ,iter] <- newA
      B[, ,iter] <- newB
      Q[, ,iter] <- newQ
      R[, ,iter] <- newR
      Vv1[, ,iter] <- V1
      Pi1[, ,iter] <- pi1
      likelihood[iter] <- -0.5 * t((xs[, , 1] - pi1)) %*% solve(V1) %*% (xs[, , 1] - pi1) 
      - likesum1 - likesum2 - 0.5 * log(det(V1)) - 0.5 * num * log(det(newR)) - 0.5 * (num -1) * log(det(newQ)) - num * log(2*pi)
      
      
      if (iter > 1) 
        cvg <-  likelihood[iter] - likelihood[iter - 1]
      # if (cvg < 0)
      #  stop ("not increase")  
      if (abs(cvg) < tol) 
        break ("converge")
    }
    
    
    # train error
    stepmax <- which.max(likelihood[1:iter])
    likely[kkk] <- max(likelihood[1:iter])
    As <- A[,,stepmax]
    Bs <- B[,,stepmax]
    Qs <- Q[,,stepmax]
    Rs <- R[,,stepmax]
    V1s <- Vv1[,,stepmax]
    pi1 <- Pi1[, ,stepmax]
    ytest <- pi1
    Vtest <- V1s 
    yts <- zts <- array(0, dim = c(leng, 1, num))
    
    for (i in 2:num){
      u <-  matrix(stack(list[[i]][,-(1:3)])[,1])
      ytsnew <- As %*% ytest + Bs %*% u
      Vtsnew <- As %*% Vtest %*% t(As) + Qs
      upperTriangle(Vtsnew) <- lowerTriangle(Vtsnew, diag=FALSE,byrow=TRUE)
      Ktest <- Vtsnew %*% solve(Rs + Vtsnew)
      z <- matrix(list[[i]]$mass)
      ytest <- ytsnew + Ktest %*% ( z - ytsnew)
      Vtest <- (eye - Ktest) %*% Vtsnew %*% t(eye - Ktest) + Ktest %*% Rs %*% t(Ktest)
      upperTriangle(Vtest) <- lowerTriangle(Vtest, diag=FALSE,byrow=TRUE)
      yts[, ,i] <- ytsnew
      zts[,,i] <- z
    }
    
    error <- yts-zts
    dim(error) <- c(leng,num)
    error <- error[,-1]
    RMSE1[kkk] <- sqrt(sum(error^2)/(nrow(error)*ncol(error)))
    MAE1[kkk] <-  sum(abs(error))/(nrow(error)*ncol(error))
  
}













# plot and prediction on test set

RMSE <- MAE <- numeric(0)
subtrain <- train
procValues <- preProcess(subtrain[,-c(1:2)], method = c("center", "scale"))
scaledTrain <-  predict(procValues, subtrain)
scaledTest <-  predict(procValues, test)
list <- lapply(split(scaledTrain, scaledTrain$timestamp),function(x) arrange(x,cell))
num <- length(list)
featureN <- ncol(scaledTrain)-3

# initialized value

leng <- length(unique(scaledTrain$cell))
newA <- rbind(A.1,A.2,A.3,A.4,A.5,A.6)
newB <- do.call(cbind, replicate(featureN, newA, simplify=FALSE))
newQ <- diag(0.05,nrow=leng,ncol=leng)
newR <- diag(0.0001,nrow=leng,ncol=leng)
pi1 <- matrix(aggregate(scaledTrain$mass, list(scaledTrain$cell), mean)[,2])
#pi1 <- matrix(list[[1]]$mass)
V1 <-  diag(aggregate(scaledTrain$mass, list(scaledTrain$cell), var)[,2])
#V1 <- diag(0.5,nrow=leng,ncol=leng)

# indicator matrix

jj <- lapply (1 : nrow(newA), function (x) diag(newA[x,]))

kk <- lapply (1 : nrow(newB), function (x) diag(newB[x,]))

# EM settings
max.iter <- 50
tol <- 0.001
likelihood <-  matrix(0, max.iter, 1)
A <- Q <- R <- Vv1 <- array(NA, dim = c(leng, leng, max.iter))
Pi1<- array(NA, dim = c(leng, 1, max.iter))
B <- array(NA, dim = c(leng, leng * featureN, max.iter))
cvg <-  1 + tol
# kalman filter
# give the u value

for (iter in 1:max.iter) {
  y.updated <- pi1
  V.updated <- V1
  K <- array(0, dim = c(leng, leng, num))
  y <- ynew <- array(0, dim = c(leng, 1, num))
  V <- Vnew <- array(0, dim = c(leng, leng, num))
  y[,,1] <- ynew[,,1] <- y.updated
  V[,,1] <- Vnew[,,1] <- V.updated
  
  for (i in 2:num){
    u <-  matrix(stack(list[[i]][,-(1:3)])[,1])
    y.new <- newA %*% y.updated + newB %*% u
    V.new <- newA %*% V.updated %*% t(newA) + newQ
    upperTriangle(V.new) <- lowerTriangle(V.new, diag=FALSE,byrow=TRUE)
    Kn <- V.new %*% solve(newR + V.new)
    eye <- diag(1,leng)
    z <- matrix(list[[i]]$mass)
    y.updated <- y.new + Kn %*% ( z - y.new)
    V.updated <- (eye - Kn) %*% V.new %*% t(eye - Kn) + Kn %*% newR %*% t(Kn)
    y[, ,i] <- y.updated
    V[, ,i] <- V.updated
    ynew[,,i] <- y.new
    Vnew[,,i] <- V.new
    K[, ,i] <- Kn
  }
  
  
  xs <-  array(NA, dim = c(leng, 1, num))
  Ps <-  array(NA, dim = c(leng, leng, num))
  Pcs <- array(NA, dim = c(leng, leng, num))
  J <-  array(NA, dim = c(leng, leng, num))
  xs[, , num] = y[, , num]
  Ps[, , num] = V[, , num]
  eye <- diag(1, leng)
  for (k in num:2) {
    u <-  matrix(stack(list[[k-1]][,-(1:3)])[,1])
    J[, , k - 1] <-   V[, , k - 1] %*% t(newA) %*% solve(Vnew[, , k]) 
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
      u <-  matrix(stack(list[[j-1]][,-(1:3)])[,1])
      sum3 <- sum3 + (Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1]))[i,]  - newB[i,] %*% kk[[i]] u %*% t(xs[, , j-1])
      sum4 <- sum4 + Ps[, , j-1] + xs[, , j - 1] %*% t(xs[, , j - 1])
      j <- j + 1
    }
    Aa[,,i] <- sum3 %*% ginv(jj[[i]]%*%sum4)
  }
  newA <- apply(Aa, 2, I)
  
  
  
  Bb <- array(NA, dim = c(1,ncol(newB), nrow(newB))) 
  sum1 <- sum2  <- 0
  for(i in 1:nrow(newB)){
    for (j in 2:num){
      u <- matrix(stack(list[[j-1]][,-(1:3)])[,1])
      sum1 <- sum1 + u %*% t(u)
      sum2 <- sum2 + xs[, , j][i]%*%t(kk[[i]] %*% u) - newA[i,] %*% jj[[i]] %*% xs[, , j-1]%*% t(kk[[i]] %*% u)
      j <- j + 1
    }
    Bb[,,i] <- sum2 %*% ginv(kk[[i]] %*% sum1 %*% t(kk[[i]]))  
  }
  newB <- apply(Bb, 2, I) 
  
  
  sum5 <- 0
  for (j in 2:num){
    u <-  matrix(stack(list[[j-1]][,-(1:3)])[,1])
    sum5 <- sum5 + ( Ps[, , j] + xs[, , j] %*% t(xs[, , j])) - (Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1])) %*% t(newA) -
      xs[, , j] %*% t(u) %*% t(newB) - newA %*% t((Pcs[, , j] + xs[, , j] %*% t(xs[, , j - 1]))) + newA %*%  (Ps[, , j-1] + xs[, , j-1] %*% t(xs[, , j-1])) %*% t(newA) +
      newA %*% xs[, , j-1] %*% t(u) %*% t(newB) - newB %*% u %*% t(xs[, , j]) + newB %*% u %*% t(xs[, , j-1]) %*% t(newA) + newB %*% u %*% t(u) %*% t(newB)
    j <- j + 1
  }
  newQ <- sum5 * (num - 1)^{-1}
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
  
  u1 <-  matrix(stack(list[[1]][,-(1:3)])[,1])
  pi1 <- pi1 + J[, , 1] %*% (xs[, , 2] - newA %*% y[, , 1] - newB %*% u1)
  # if we already know it, don't need to estimate it
  V1 <- V1 + J[, , 1] %*%(Ps[, , 2] - Vnew[, , 2])%*%t(J[, , 1])
  upperTriangle(V1) <- lowerTriangle(V1, diag=FALSE,byrow=TRUE)
  V1 <- make.positive.definite(V1)
  # check
  if (suppressWarnings(is.nan(log(det(V1)))==TRUE)){
    V1 <- make.positive.definite(V1)
  } else{
    V1 <- V1
  }
  
  
  
  # calculate likelihood
  likesum1 <- 0
  
  for (j in 2:num){
    u <- matrix(stack(list[[j-1]][,-(1:3)])[,1])
    likesum1 <- likesum1 + 0.5 * t((xs[, , j] - newA %*% xs[, , j-1] - newB %*% u)) %*% solve(newQ) %*% (xs[, , j] - newA %*% xs[, , j-1] - newB %*% u)
    j <- j + 1
  }
  
  likesum2 <- 0
  for (j in 1:num){
    zt <- matrix(list[[j]][,3])
    likesum2 <- likesum2 + 0.5 * t((zt - xs[, , j])) %*% solve(newR) %*% (zt - xs[, , j])
    j <- j + 1
  }
  
  A[, ,iter] <- newA
  B[, ,iter] <- newB
  Q[, ,iter] <- newQ
  R[, ,iter] <- newR
  Vv1[, ,iter] <- V1
  Pi1[, ,iter] <- pi1
  likelihood[iter] <- -0.5 * t((xs[, , 1] - pi1)) %*% solve(V1) %*% (xs[, , 1] - pi1) 
  - likesum1 - likesum2 - 0.5 * log(det(V1)) - 0.5 * num * log(det(newR)) - 0.5 * (num -1) * log(det(newQ)) - num * log(2*pi)
  
  
  # prediction in train
  ytest <- pi1
  Vtest <- V1
  yts <- zts <- array(0, dim = c(leng, 1, num))
  
  for (i in 2:num){
    u <-  matrix(stack(list[[i]][,-(1:3)])[,1])
    ytsnew <- newA %*% ytest + newB %*% u
    Vtsnew <- newA %*% Vtest %*% t(newA) + newQ
    upperTriangle(Vtsnew) <- lowerTriangle(Vtsnew, diag=FALSE,byrow=TRUE)
    Ktest <- Vtsnew %*% solve(newR + Vtsnew)
    z <- matrix(list[[i]]$mass)
    ytest <- ytsnew + Ktest %*% ( z - ytsnew)
    Vtest <- (eye - Ktest) %*% Vtsnew %*% t(eye - Ktest) + Ktest %*% newR %*% t(Ktest)
    upperTriangle(Vtest) <- lowerTriangle(Vtest, diag=FALSE,byrow=TRUE)
    yts[, ,i] <- ytsnew
    zts[,,i] <- z
  }
  
  error <- yts-zts
  dim(error) <- c(leng,num)
  error <- error[,-1]
  RMSE[iter] <- sqrt(sum(error^2)/(nrow(error)*ncol(error)))
  MAE[iter] <-  sum(abs(error))/(nrow(error)*ncol(error))
  
  if (iter > 1) 
    cvg <-  likelihood[iter] - likelihood[iter - 1]
  # if (cvg < 0)
  #  stop ("not increase")  
  if (abs(cvg) < tol) 
    break ("converge")
}


# prediction in test
stepmax <- which.max(likelihood[1:iter])
As <- A[,,stepmax]
Bs <- B[,,stepmax]
Qs <- Q[,,stepmax]
Rs <- R[,,stepmax]  
V1s <- Vv1[,,stepmax]
testlist <- lapply(split(scaledTest, scaledTest$timestamp),function(x) arrange(x,cell))
subtestlist <- testlist[1:120]
testnum <- length(subtestlist)
pi1 <- matrix(subtestlist [[1]]$mass)
ytest <-  pi1
Vtest <- V1s 
yts <- zts <- array(0, dim = c(leng, 1, testnum))

for (i in 2:testnum){
  u <-  matrix(stack(subtestlist[[i]][,-(1:3)])[,1])
  ytsnew <- As %*% ytest + Bs %*% u
  Vtsnew <- As %*% Vtest %*% t(As) + Qs
  upperTriangle(Vtsnew) <- lowerTriangle(Vtsnew, diag=FALSE,byrow=TRUE)
  Ktest <- Vtsnew %*% solve(Rs + Vtsnew)
  z <- matrix(subtestlist[[i]]$mass)
  ytest <- ytsnew + Ktest %*% ( z - ytsnew)
  Vtest <- (eye - Ktest) %*% Vtsnew %*% t(eye - Ktest) + Ktest %*% Rs %*% t(Ktest)
  upperTriangle(Vtest) <- lowerTriangle(Vtest, diag=FALSE,byrow=TRUE)
  yts[, ,i] <- ytsnew
  zts[,,i] <- z
}

error <- yts-zts
dim(error) <- c(leng,testnum)
error <- error[,-1]
RMSEtest <- sqrt(sum(error^2)/(nrow(error)*ncol(error)))
MAEtest <-  sum(abs(error))/(nrow(error)*ncol(error))
# only in core region
coreerror <- error[c(6,7,10,11,14,15,18,19),]
RMSEcore <- sqrt(sum(coreerror^2)/(nrow(coreerror)*ncol(coreerror)))
MAEcore <-  sum(abs(coreerror))/(nrow(coreerror)*ncol(coreerror))






library(ggplot2)
library(gtable)
library(grid)  
like <- data.frame(exp(likelihood[1:iter]))
colnames(like) <- "like"
like$ID<-seq.int(nrow(like))
like$RMSE <- RMSE
like <- mutate(like,color = (max(like) == like))

plot1 <- ggplot(like,aes(ID,like))+geom_line(color = "blue") 
plot1 <- plot1  + theme_bw() + geom_point(aes(color = color)) + scale_color_manual(values = c(NA, "red"),guide = FALSE) + labs(x="Iteration",y="Likelihood")  + theme(axis.line = element_line(size=1, colour = "black"),
                                                                                                                                                                      panel.grid.major = element_line(colour = "#d3d3d3"),
                                                                                                                                                                      panel.grid.minor = element_blank(),
                                                                                                                                                                      panel.border = element_blank(), panel.background = element_blank(),
                                                                                                                                                                      axis.text.x=element_text(colour="black", size = 9),
                                                                                                                                                                      axis.text.y=element_text(colour="black", size = 9),
                                                                                                                                                                      legend.title = element_blank())
plot2 <- ggplot(like,aes(ID,RMSE))+geom_line(color = "green")
plot2 <- plot2  + theme_bw() + labs(x="Iteration",y="RMSE")  + theme(axis.line = element_line(size=1, colour = "black"),
                                                                     panel.grid.major = element_blank(),
                                                                     panel.grid.minor = element_blank(),
                                                                     panel.border = element_blank(), panel.background = element_blank(),
                                                                     axis.text.x=element_text(colour="black", size = 9),
                                                                     axis.text.y=element_text(colour="black", size = 9),
                                                                     legend.title = element_blank())
# Get the ggplot grobs
g1 <- ggplotGrob(plot1)
g2 <- ggplotGrob(plot2)


# Get the location of the plot panel in g1.
# These are used later when transformed elements of g2 are put back into g1
pp <- c(subset(g1$layout, name == "panel", se = t:r))

# Overlap panel for second plot on that of the first plot
g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)



hinvert_title_grob <- function(grob){
  
  # Swap the widths
  widths <- grob$widths
  grob$widths[1] <- widths[3]
  grob$widths[3] <- widths[1]
  grob$vp[[1]]$layout$widths[1] <- widths[3]
  grob$vp[[1]]$layout$widths[3] <- widths[1]
  
  # Fix the justification
  grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
  grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
  grob$children[[1]]$x <- unit(1, "npc") - grob$children[[1]]$x
  grob
}

# Get the y axis title from g2
index <- which(g2$layout$name == "ylab-l") # Which grob contains the y axis title?
ylab <- g2$grobs[[index]]                # Extract that grob
ylab <- hinvert_title_grob(ylab)         # Swap margins and fix justifications

# Put the transformed label on the right side of g1
g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = "off", name = "ylab-r")

# Get the y axis from g2 (axis line, tick marks, and tick mark labels)
index <- which(g2$layout$name == "axis-l")  # Which grob
yaxis <- g2$grobs[[index]]                  # Extract the grob

# yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
# The relevant grobs are contained in axis$children:
#   axis$children[[1]] contains the axis line;
#   axis$children[[2]] contains the tick marks and tick mark labels.

# First, move the axis line to the left
yaxis$children[[1]]$x <- unit.c(unit(0, "npc"), unit(0, "npc"))

# Second, swap tick marks and tick mark labels
ticks <- yaxis$children[[2]]
ticks$widths <- rev(ticks$widths)
ticks$grobs <- rev(ticks$grobs)

# Third, move the tick marks
ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, "npc") + unit(3, "pt")

# Fourth, swap margins and fix justifications for the tick mark labels
ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])

# Fifth, put ticks back into yaxis
yaxis$children[[2]] <- ticks

# Put the transformed yaxis on the right side of g1
g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = "off", name = "axis-r")



# Construct the four grobs - two symbols and two labels
L1 = linesGrob(x = unit(c(.7, .25), "npc"), y = unit(c(.5, .5), "npc"), gp = gpar(col = "blue"))
L2 = linesGrob(x = unit(c(.7, .25), "npc"), y = unit(c(.5, .5), "npc"), gp = gpar(col = "green"))
T1 = textGrob("Likelihood", x = .02, just = "left",gp = gpar(cex=0.55))
T2 = textGrob("RMSE", x = .02, just = "left",gp = gpar(cex=0.55))


# Construct a gtable - 2 columns X 4 rows
leg = gtable(width = unit(c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5), "cm"), height = unit(c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8), "cm"))
leg = gtable_add_grob(leg, rectGrob(gp = gpar(fill = NA, col = "white")), t=6,l=7,b=7,r=8)

# Place the 4 grob into the table
leg = gtable_add_grob(leg, L1, t=6, l=7)
leg = gtable_add_grob(leg, L2, t=7, l=7)
leg = gtable_add_grob(leg, T1, t=6, l=8)
leg = gtable_add_grob(leg, T2, t=7, l=8)
pos = g1$layout[grepl("panel", g1$layout$name), c('t', 'l')]
gg = gtable_add_grob(g1, leg, t = pos$t+0.5, l = pos$l+0.5)



grid.newpage()
grid.draw(gg)


##-----------------------##
## penalized  regression ##
### data preparation #
load("data.RData")
Birds.r <- subset(Birds,select=c(timestamp,cell,mass,latitude,longitude,heading,velocity,mass_correction))

train.r <- subset(Birds.r,lubridate::day(Birds.r$times)==7)
test.r <- subset(Birds.r,lubridate::day(Birds.r$times)==8)
procValues <- preProcess(train.r[,-c(1:2)], method = c("center", "scale"))
scaledTrain <-  predict(procValues, train.r)
scaledtest <-  predict(procValues, test.r)


trans <- function(mat){
  grid <- alply(mat,1,.dims=TRUE)
  gridmatrix <- lapply(grid, function(x) matrix(x,nrow=6,byrow=TRUE))
  gridmatrix <- lapply(gridmatrix,function(x) replace(x,is.na(x),0))
  
#  names <- c("loc22","loc32","loc42","loc52",
#             "loc23","loc33","loc43","loc53")
   names <- c("loc11","loc21","loc31","loc41","loc51","loc61","loc12","loc22","loc32","loc42","loc52","loc62",
              "loc13","loc23","loc33","loc43","loc53","loc63","loc14","loc24","loc34","loc44","loc54","loc64")
  func <- function(x) {
    mat <- matrix(unlist(x),6,4)
    w <-  which(mat==mat, arr.ind=TRUE)
    d <- as.matrix(dist(w, "maximum", diag=TRUE, upper=TRUE))
    d[d==0]=1
    a <- apply(d, 1, function(i) mat[i == 1] )
#    length <- lapply(a,length) == 9
    a[[1]] <- c(0,0,0,0,a[[1]][1:2],0,a[[1]][3:4])
    a[[2]] <- c(0,a[[2]][1:2],0,a[[2]][3:4],0,a[[2]][5:6])
    a[[3]] <- c(0,a[[3]][1:2],0,a[[3]][3:4],0,a[[3]][5:6])
    a[[4]] <- c(0,a[[4]][1:2],0,a[[4]][3:4],0,a[[4]][5:6])
    a[[5]] <- c(0,a[[5]][1:2],0,a[[5]][3:4],0,a[[5]][5:6])
    a[[6]] <- c(0,a[[6]][1:2],0,a[[6]][3:4],0,0,0)
    a[[7]] <- c(0,0,0,a[[7]])
    a[[12]] <- c(a[[12]],0,0,0)
    a[[13]] <- c(0,0,0,a[[13]])
    a[[18]] <- c(a[[18]],0,0,0)
    a[[19]] <- c(0,0,0,a[[19]][1:2],0,a[[19]][3:4],0)
    a[[20]] <- c(a[[20]][1:2],0,a[[20]][3:4],0,a[[20]][5:6],0)
    a[[21]] <- c(a[[21]][1:2],0,a[[21]][3:4],0,a[[21]][5:6],0)
    a[[22]] <- c(a[[22]][1:2],0,a[[22]][3:4],0,a[[22]][5:6],0)
    a[[23]] <- c(a[[23]][1:2],0,a[[23]][3:4],0,a[[23]][5:6],0)
    a[[24]] <- c(a[[24]][1:2],0,a[[24]][3:4],0,0,0,0)
#    b <- a[length]
#    names(b) <- names
#    b
    names(a) <- names
    a
  }
  listt <-lapply(gridmatrix, func)
#  .id <- rep(c(1:120),each=8)
   .id <- rep(c(1:120),each=24)
#  location <- rep(c("loc22","loc32","loc42","loc52",
#                    "loc23","loc33","loc43","loc53"),times=120)
   location <- rep(c("loc11","loc12","loc13","loc14","loc21","loc22","loc23","loc24","loc31","loc32","loc33","loc34",
                       "loc41","loc42","loc43","loc44","loc51","loc52","loc53","loc54","loc61","loc62","loc63","loc64"),times=120)
  
  prehead <-  as.data.frame(cbind(.id,location))
  
  unlistdf <- ldply(listt, function(x){
    df <- ldply(x, function(z) as.data.frame(matrix(z,ncol=9,byrow=T)))
    names(df)[1] <- "location"; 
    df
  })
  unlistdf
}
aqm <- melt(scaledtest, id=c("timestamp", "cell"), na.rm=TRUE)
acast <- acast(aqm, timestamp ~ cell ~ variable)
attply <- alply(acast,3,.dims = TRUE)
ans <- lapply(attply,trans)
allattribute <- do.call(cbind, ans)
colnames(allattribute)[which(names(allattribute) == "mass..id")] <- "times"
colnames(allattribute)[which(names(allattribute) == "mass.location")] <- "cell"
pattern <- "id|location"
alldata <- allattribute[,-grep(pattern,names(allattribute))]
n.features <- length(names(alldata)) -2 
addframe <- as.data.frame(matrix(0,ncol=n.features,nrow=8))
names(addframe) <- names(alldata )[-c(1:2)]
pretime.fea <- rbind(do.call("rbind", replicate(1, addframe, simplify = FALSE)),alldata[c(1:(nrow(alldata)-8*1)),-c(1:2)])
names(pretime.fea) <- paste(names(addframe),'.p1',sep="")
scaledtest <- cbind(alldata[,c(1,2,7)],pretime.fea)
scaledtest  <- scaledtest[-c(1:24),]
colnames(scaledtest)[which(names(scaledtest ) == "mass.V5")] <- "mass"
scaledtest$cell <- as.factor(scaledtest$cell) 
#train set
aqm <- melt(scaledTrain, id=c("timestamp", "cell"), na.rm=TRUE)
acast <- acast(aqm, timestamp ~ cell ~ variable)
attply <- alply(acast,3,.dims = TRUE)
ans <- lapply(attply,trans)
allattribute <- do.call(cbind, ans)
colnames(allattribute)[which(names(allattribute) == "mass..id")] <- "times"
colnames(allattribute)[which(names(allattribute) == "mass.location")] <- "cell"
pattern <- "id|location"
alldata <- allattribute[,-grep(pattern,names(allattribute))]
n.features <- length(names(alldata)) -2 
addframe <- as.data.frame(matrix(0,ncol=n.features,nrow=8))
names(addframe) <- names(alldata )[-c(1:2)]
pretime.fea <- rbind(do.call("rbind", replicate(1, addframe, simplify = FALSE)),alldata[c(1:(nrow(alldata)-8*1)),-c(1:2)])
names(pretime.fea) <- paste(names(addframe),'.p1',sep="")
scaledTrain <- cbind(alldata[,c(1,2,7)],pretime.fea)
scaledTrain  <- scaledTrain[-c(1:24),]
colnames(scaledTrain)[which(names(scaledTrain ) == "mass.V5")] <- "mass"
scaledTrain$cell <- as.factor(scaledTrain$cell) 

### lasso regression #


# assign fold
library(h2o)
h2o.init(max_mem_size = "8g",nthreads = -1)
rmse.r <- mae.r <- RMSE.r <-  numeric(0)
fold <- c(44,59,74,89,104,119)
lamb <- seq(0.005,0.25,by=0.005)
for (j in 1:length(lamb)){
for (ii in 1:5){
subtrain.r <- scaledTrain[c(1:(fold[ii]*8)),]
validation.r <- scaledTrain[c((fold[ii]*8+1):(fold[ii+1]*8)),]
birds.models <- h2o.glm(x = c(2,4:57), y = 3, training_frame = as.h2o(subtrain.r),
                        family = "gaussian",lambda = lamb[j],
                        validation_frame=as.h2o(validation.r),
                        alpha = 1)
rmse.r[ii] <- sqrt(h2o.mse(birds.models,valid = TRUE))
}
RMSE.r[j] <- mean(rmse.r)
}
bestlamb <- lamb[which.min(RMSE.r)]
model.r <- h2o.glm(x = c(2,4:57), y = 3, training_frame = as.h2o(scaledTrain),
                   family = "gaussian",lambda = bestlamb,
                   alpha = 1)
predictions <- h2o.predict(model.r, as.h2o(scaledTrain))
RMSE.train <- sqrt(mean((scaledTrain[,3] - as.data.frame(predictions$predict)[,1])^2))
predictions.test <- h2o.predict(model.r, as.h2o(scaledtest[1:1416,]))
RMSE.test <- sqrt(mean((scaledtest[1:1416,][,3] - as.data.frame(predictions.test$predict)[,1])^2))
MAE.test <- mean(abs(scaledtest[1:1416,][,3] - as.data.frame(predictions.test$predict)[,1]))
## core region
coreregion <- c("loc22","loc32","loc42","loc52","loc23","loc33","loc43","loc53")
coretest <- subset(scaledtest, cell %in% coreregion)
predictions.test <- h2o.predict(model.r, as.h2o(coretest[1:472,]))
RMSE.test <- sqrt(mean((coretest[1:472,][,3] - as.data.frame(predictions.test$predict)[,1])^2))
MAE.test <- mean(abs(coretest[1:472,][,3] - as.data.frame(predictions.test$predict)[,1]))


