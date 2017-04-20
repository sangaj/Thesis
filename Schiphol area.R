### Schiphol Area Figure

library(ggplot2)
library(ggmap)

map <- qmap("2141CN",zoom = 12, maptype="toner-lite")
l_lon <-  4.706
r_lat <-  52.368
r_lon <-  4.717
l_lat <-  52.325
#map <- map+geom_rect(aes(xmin=l_lon, xmax=r_lon, ymin=l_lat ,ymax=r_lat), fill="lightblue", alpha=0.2)
lat_det <-  2/111
lon_det <-  2/(cos(r_lat*pi/180)*111.321)
bound_x <-  c(l_lon-lon_det,l_lat-lat_det)
bound_y <-  c(r_lon+lon_det,r_lat+lat_det)
#map <- map+geom_rect(aes(xmin=l_lon-lon_det, xmax=r_lon+lon_det, ymin=l_lat-lat_det ,ymax=r_lat+lat_det), fill="lightblue", alpha=0.2)

gridx <- seq(l_lon-lon_det,r_lon+lon_det,length.out=4)
gridy <- seq(l_lat-lat_det,r_lat+lat_det,length.out=6)
grid <- expand.grid(x = gridx, y = gridy)
grid$XX <- rep(1:4, times=6)
grid$YY <- rep(1:6, each=4)
subgrid <- subset(grid, YY > 1 & YY < 6 & XX < 4 & XX > 1)
map + geom_tile(data=grid,aes(x, y),color="lightblue",alpha=0.1) + geom_tile(data=subgrid,aes(x, y),color="orangered",alpha=0.1)
