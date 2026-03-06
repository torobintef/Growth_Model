  # gam ---------------------------------------------------------------------
  
  # renv::install("mgcv)
  # renv::install("gratia")
  # renv::install("gratia")
  
  # devtools::install_github("BlakeRMills/MetBrewer") 
  library(mgcv) #gam
  library(gratia)
  library(dplyr)
  
  LT <- "DB"
  
  data <- read.csv("result/glmm/distance.csv",header = T)
  d1 <- subset(data,data$group == LT & data$SP == "canopy" & data$layer == "understory")
  d2 <- subset(data,data$group == LT & data$SP == "canopy" & data$layer == "overstory")
  d3 <- subset(data,data$group == LT & data$SP == "non_canopy" & data$layer == "understory")
  
  col1 <- MetBrewer::met.brewer("Egypt", 3)[1]
  col2 <- MetBrewer::met.brewer("Egypt", 3)[2]
  col3 <- MetBrewer::met.brewer("Egypt", 3)[3]
  
  
  fit1 <- gam(deltaAIC ~ s(distance),data = d1)
  fit2 <- gam(deltaAIC ~ s(distance),data = d2)
  fit3 <- gam(deltaAIC ~ s(distance),data = d3)
  deriv1 <- derivatives(fit1, select = "s(distance)")
  deriv2 <- derivatives(fit2, select = "s(distance)")
  deriv3 <- derivatives(fit3, select = "s(distance)")
  plot(deriv1$distance,deriv1$.derivative,type = "l",col = col1)
  lines(deriv2$distance,deriv2$.derivative,type = "l",col = col2)
  lines(deriv2$distance,deriv3$.derivative,type = "l",col = col3)
  deriv1$zero_in_ci <- deriv1$.lower_ci <= 0 & deriv1$.upper_ci >= 0
  deriv2$zero_in_ci <- deriv2$.lower_ci <= 0 & deriv2$.upper_ci >= 0
  deriv3$zero_in_ci <- deriv3$.lower_ci <= 0 & deriv3$.upper_ci >= 0
  saturation_point1 <- min(deriv1$distance[deriv1$zero_in_ci])
  saturation_point1
  saturation_point2 <- min(deriv2$distance[deriv2$zero_in_ci])
  saturation_point2
  saturation_point3 <- min(deriv3$distance[deriv3$zero_in_ci])
  saturation_point3
  
  
  pred1 <- predict(fit1)
  pred2 <- predict(fit2)
  pred3 <- predict(fit3)
  y1 <- approx(d1$distance,pred1,xout = saturation_point1)$y
  y2 <- approx(d2$distance,pred2,xout = saturation_point2)$y
  y3 <- approx(d3$distance,pred3,xout = saturation_point3)$y
  
  svg(paste0("figure/distance_",LT,".svg"),width = 5.7,height = 5,family = "Aerial")
  par(oma = c(4,6,2,2),ps = 20,pin = c(4,3.5),bg = "#FEF7DB")
  plot(pred1,type = "n",xlab = "",axes = F,
       ylab = "",lwd = 3,ylim = c(0,300))
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "white")
  lines(d1$distance,pred1,col = col1,lwd = 3)
  lines(d1$distance,pred2,col = col2,lwd = 3)
  lines(d1$distance,pred3,col = col3,lwd = 3)
  lines(d1$distance,d1$deltaAIC,col = adjustcolor("black",alpha.f = 0.5),lwd = 1.5)
  lines(d1$distance,d2$deltaAIC,col = adjustcolor("black",alpha.f = 0.5),lwd = 1.5)
  lines(d1$distance,d3$deltaAIC,col = adjustcolor("black",alpha.f = 0.5),lwd = 1.5)
  segments(
    x0 = saturation_point1, x1 = saturation_point1,
    y0 = par("usr")[3],y1 = y1,
    lty = 2,lwd = 3,col = adjustcolor(col1, alpha.f = 1)
  )
  segments(
    x0 = saturation_point2, x1 = saturation_point2,
    y0 = par("usr")[3],y1 = y2,
    lty = 2,lwd = 3,col = adjustcolor(col2, alpha.f = 1)
  )
  segments(
    x0 = saturation_point3, x1 = saturation_point3,
    y0 = par("usr")[3],y1 = y3,
    lty = 2,lwd = 3,col = adjustcolor(col3, alpha.f = 1)
  )
  legend("topleft",legend = c("canopy species(DBH<30)","canopy species(DBH>=30)","non-canopy species"),col = c(col1,col2,col3),
         lty = 1,lwd = 3,bty = "n",y.intersp = 1.5)
  points(saturation_point1,y1,col = col1,cex = 1.5)
  points(saturation_point2,y2,col = col2,cex = 1.5)
  points(saturation_point3,y3,col = col3,cex = 1.5)
  box(lwd = 2)
  axis(side = 1,at = seq(0,50,by = 10),lwd = 1.5,tck = -0.03,padj = 0)
  axis(side = 1,at = seq(5,45,by = 10),lwd = 1.5,tck = -0.015,padj = 1,labels = F)
  # axis(side = 2,at = seq(0,300,by = 100),lwd = 1.5,las = 1,tck = -0.05,hadj = 1.2,labels = c("   0","100","200","300"))
  axis(side = 2,at = seq(0,300,by = 100),lwd = 1.5,las = 1,tck = -0.05,hadj = 1.2,labels = F)
  mtext("（B）落葉広葉樹",side = 3,line = 1,cex = 1,family = "BIZ UDPGothic")
  mtext("近接距離",side = 1,line = 3.3,cex = 1,family = "BIZ UDPGothic")
  # mtext("(m)",side = 1,line = 1,at = 55,cex = 0.8)
  # mtext(expression(paste(Delta,AIC,sep = "")),side = 2,line = 5,cex = 1)
  dev.off()
  
