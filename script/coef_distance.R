library(tidyverse)

functional <- "understory"
leaf <- "DB"
df <- read.csv("result/summary/parameters.csv",header = T) %>%
  subset(distance == 5 | distance == 10 | distance == 15 | distance == 20 | distance == 25 |
           distance == 30 | distance == 35 | distance == 40 | distance == 45 | distance == 50 ) %>%
  subset(LT == leaf & FT == functional)
col1 <- MetBrewer::met.brewer("Egypt", 2)[1]
col2 <- MetBrewer::met.brewer("Egypt", 2)[2]

# 例: dfの列
# df$distance   (5,10,...,50)
# df$Estimate       回帰係数
# df$p          P値
# df$CI_low     95%CI下限
# df$CI_high    95%CI上限

    # 有意かどうか
    sig <- df$Pr...t.. < 0.05
    # 色
    cols <- ifelse(df$Row.names == "two_DB", col1, col2)
    shape <- ifelse(sig,19,1)
    # 空のプロット
    
    svg(paste0("figure/coef_distance",leaf,"_",functional,".svg"),width = 8,height = 3.5,family = "Aerial")
    par(oma = c(4,6,0,2),ps = 20,pin = c(6,2),bg = "#FEF7DB")
    plot(df$distance, df$Estimate,
         type = "n",
         xlab = "",
         ylab = "",
         xlim = c(5, 50),
         ylim = c(-1,0.5),
         axes = F)
    rect(par("usr")[1], par("usr")[3],
         par("usr")[2], par("usr")[4],
         col = "white")
    abline(h = 0,lty = 2)
    # エラーバー
    arrows(df$distance, df$CI_low,
           df$distance, df$CI_high,
           angle = 90, code = 3, length = 0.05,
           col = cols)
    
    # 回帰係数の点
    points(df$distance, df$Estimate,
           pch = shape,
           col = cols)
    
    axis(side = 1,at = seq(5,45,by = 10),lwd = 1.5,tck = -0.03,padj = 1,labels = F)
    # axis(side = 1,at = seq(10,50,by = 10),lwd = 1.5,tck = -0.06,padj = 0)
    axis(side = 1,at = seq(0,50,by = 10),lwd = 1.5,tck = -0.06,padj = 0,labels = F)
    # axis(side = 2,at = seq(-1,0.5,by = 0.5),lwd = 1.5,las = 1,tck = -0.05,hadj = 1.2)
    axis(side = 2,at = seq(-1,0.5,by = 0.5),lwd = 1.5,las = 1,tck = -0.05,hadj = 1.2,labels = F)
    text(4,0.4,adj = 0,labels = "（B）落葉広葉樹－林冠構成種",family = "BIZ UDPGothic")
    # mtext("近接距離",side = 1,line = 3.3,cex = 1,family = "BIZ UDPGothic")
    box(lwd = 1.5)
    # mtext("(m)",side = 1,line = 1,at = 53,cex = 0.8)
    # mtext("Coefficient",side = 2,line = 5,cex = 1)
    dev.off()

