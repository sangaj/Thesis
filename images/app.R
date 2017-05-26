#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(reshape2)
library(plyr)
library(ggmap)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
   titlePanel("Prediction Error in Polderbaan"),
  
   sliderInput("Times",
               "Choose times on 8th April 2016:",
                min = strptime("2016-04-08 10:01:00","%Y-%m-%d %H:%M"),
                max = strptime("2016-04-08 10:59:00","%Y-%m-%d %H:%M"),
                value=c(strptime("2016-04-08 10:01:00","%Y-%m-%d %H:%M")),timeFormat="%a %H:%M",ticks = F, animate = T,step = 60),

  plotOutput("Plot1")
)


# Define server logic required to draw a histogram
server <- function(input, output) {
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
  data <- read.csv("yts.csv",header=TRUE)
  data2 <- read.csv("zts.csv",header=TRUE)
  error <- abs(data - data2)
  subdata <- error[c(6,7,10,11,14,15,18,19),]
  subdata$X <- c("22","23","32","33","42","43","52","53")
  meltsubdata <- melt(subdata,id="X")
  meltsubdata$time <- rep(seq(from = as.POSIXct("2016-04-08 10:00"),to = as.POSIXct("2016-04-08 10:59"), by =60),each=8)
  subgrid$X <- paste(subgrid$YY,subgrid$XX,sep="")
  joindata <- join(meltsubdata,subgrid,by="X")
  plotdata <- subset(joindata,select=c("value","time","x","y"))
  map <- qmap("2141CN",zoom = 12, maptype="toner-lite")
  output$Plot1 <- renderPlot({
     subplotdata <- subset(plotdata, time==input$Times)
     map + geom_tile(data = subplotdata, aes(x = x, y = y, fill = value)) + scale_fill_gradient(low="white",high="blue",limit=c(0,5))
      
   })
 
}

# Run the application 
shinyApp(ui = ui, server = server)

