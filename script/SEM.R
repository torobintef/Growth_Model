
# packages ----------------------------------------------------------------

# install.packages("glmmML")
# install.packages("lme4")
# install.packages("performance")
# install.packages("broom")
# install.packages("MuMIn")
# install.packages("ggeffects")

library(tidyverse)
library(performance)
library(lme4)
library(broom)
library(MuMIn)
options(na.action="na.fail")
library(ggeffects)
library(parameters)
library(spdep)
library(spatialreg)
library(car)
library(mgcv) #gam


# repeat ------------------------------------------------------------------


mname <- "understory"

# 空のベクトルまたはリストを作成
LeafType <- c("EB","DB")
model_results <- list()
distance <- vector()
best.R2 <- vector()
converge <- vector()
best.converge <- vector()
varcor.SP <- vector()
varcor.MESH <- vector()
best.varcor.SP <- vector()
best.varcor.MESH <- vector()
R2 <- vector()
d600 <- read.csv("KN600.csv")[c("STEM_ID","ALT","slope","TPI")]
d700 <- read.csv("KN700.csv")[c("STEM_ID","ALT","slope","TPI")]
d800 <- read.csv("KN800.csv")[c("STEM_ID","ALT","slope","TPI")]

for(j in 1:2){
  for(r in 1:50){

    distance[r] <- r
    d1 <- read.csv(paste0("Dd_20_edge/Dd600_r",r,"_1.csv"),header = T) %>%
      left_join(d600,by = "STEM_ID") %>%
      group_by(IND_ID) %>%
      mutate(DBH = max(DBH)) %>%
      distinct(IND_ID,.keep_all = TRUE) %>%
      ungroup()
    
    d2 <- read.csv(paste0("Dd_20_edge/Dd700_r",r,"_1.csv"),header = T) %>%
      left_join(d700,by = "STEM_ID") %>%
      group_by(IND_ID) %>%
      mutate(DBH = max(DBH)) %>%
      distinct(IND_ID,.keep_all = TRUE) %>%
      ungroup()
    
    d3 <- read.csv(paste0("Dd_20_edge/Dd800_r",r,"_1.csv"),header = T) %>%
      left_join(d800,by = "STEM_ID") %>%
      group_by(IND_ID) %>%
      mutate(DBH = max(DBH)) %>%
      distinct(IND_ID,.keep_all = TRUE) %>%
      ungroup()
    
    d1$PLOT <- factor("ALT600")
    d2$PLOT <- factor("ALT700")
    d3$PLOT <- factor("ALT800")
    d1$X <- d1$X + 0
    d1$Y <- d1$Y + 0
    d2$X <- d2$X + 1000
    d2$Y <- d2$Y + 1000
    d3$X <- d3$X + 2000
    d3$Y <- d3$Y + 2000
    
    td1 <- rbind(d1,d2,d3) %>%
      rename(one_all = ci1, one_EB = ci2, one_DB = ci3, one_EC = ci4)
    
    d4 <- read.csv(paste0("Dd_20_edge/Dd600_r",r,"_2.csv"),header = T) %>%
      group_by(IND_ID) %>%
      mutate(DBH = max(DBH)) %>%
      distinct(IND_ID,.keep_all = TRUE) %>%
      ungroup()
    d5 <- read.csv(paste0("Dd_20_edge/Dd700_r",r,"_2.csv"),header = T) %>%
      group_by(IND_ID) %>%
      mutate(DBH = max(DBH)) %>%
      distinct(IND_ID,.keep_all = TRUE) %>%
      ungroup()
    d6 <- read.csv(paste0("Dd_20_edge/Dd800_r",r,"_2.csv"),header = T) %>%
      group_by(IND_ID) %>%
      mutate(DBH = max(DBH)) %>%
      distinct(IND_ID,.keep_all = TRUE) %>%
      ungroup()
    
    td2 <- rbind(d4,d5,d6) %>%
      rename(two_all = ci1, two_EB = ci2, two_DB = ci3, two_EC = ci4)
    
    td2_unique <- td2 %>%
      dplyr::select(-intersect(names(td1), names(td2)))
    
    d <- bind_cols(td1, td2_unique) %>%
      subset(!is.na(RGR20) & !is.na(DBH05) &
               !is.na(one_EB) & !is.na(one_DB) & !is.na(one_EC) &
               !is.na(two_EB) & !is.na(two_DB) & !is.na(two_EC)) %>%
      subset(SP != "アカガシ" & SP != "イヌガシ" & SP != "ウラジロガシ" &
               SP != "ブナ" & SP != "ケヤキ" & SP != "ヒメシャラ" &
               SP != "イヌシデ" & SP != "オオモミジ" & SP != "イタヤカエデ" & SP != "カジカエデ") %>%
      subset(LT == LeafType[j] & LD == 0 & DBH05 < 30 & RGR20 >= 0)
    
    # データの対数変換と正規化(z-score)
    d$RGR20 <- as.numeric(scale(d$RGR20))
    d$DBH05 <- as.numeric(scale(d$DBH05))
    d$one_EB <- as.numeric(scale(d$one_EB + 0.01))
    d$one_DB <- as.numeric(scale(d$one_DB + 0.01))
    d$two_EB <- as.numeric(scale(d$two_EB + 0.01))
    d$two_DB <- as.numeric(scale(d$two_DB + 0.01))
    d$ALT <- as.numeric(scale(d$ALT))
    d$slope <- as.numeric(scale(d$slope))
    d$TPI <- as.numeric(scale(d$TPI))
    
    # 近接行列
    coords <- as.matrix(d[, c("X", "Y")])
    # sf_d <- st_as_sf(d, coords = c("X", "Y"), crs = NA)
    nb <- dnearneigh(coords, d1 = 0, d2 = 30)
    lw <- nb2listw(nb, style = "W",zero.policy = T)
    # lw <- nb2listwdist(nb, x = sf_d, type="idw", style="W",
    #                    alpha = 1, dmax = NULL, longlat = NULL, zero.policy=TRUE)
 
    # model_summary
    print(paste0("---Now caliculating ",mname," model of ",LeafType[j]," for r = ",r,"----------"))
    
    formula <- RGR20 ~ DBH05 + two_EB + two_DB + ALT + slope + TPI
    base_formula <- RGR20 ~ DBH05 + ALT + slope + TPI
    lmod <- lm(formula,data = d)
    model <- spatialreg::errorsarlm(formula,data = d,listw = lw,zero.policy = T)
    base_model <- spatialreg::errorsarlm(base_formula,data = d,listw = lw,zero.policy = T)
    model_summary <- summary(model)
    res <- residuals(model)
    mtest <- moran.test(res,lw)
    
    # FE(各固定効果のP値とVIF値)
    coef <- as.data.frame(model_summary$Coef)
    vif <- as.data.frame(vif(lmod))
    FE <- merge(coef,vif, by = "row.names", all = TRUE)
    
    # result(モデルの回帰係数とAIC,logLik)
    AIC <- AIC(model)
    deltaAIC <- AIC(base_model) - AIC(model)
    logLik <- logLik(model)
    lambda <- model_summary$lambda
    M_stat <- mtest[["statistic"]]
    M_p.value <- mtest[["p.value"]]
    
    
    # 回帰係数をデータフレームに格納し、行と列を反転
    result <- as.data.frame(cbind(AIC,deltaAIC,logLik,lambda,M_stat,M_p.value))
    
    
    # result（best.result），result（best.result）の結合
    if(r == 1){
      total_FE <- FE
      total_results <- result
      sd_ci <- sd_value
    }else{
      total_FE <- bind_rows(total_FE,FE)
      sd_ci <- bind_rows(sd_ci,sd_value)
      
      common_col <- intersect(names(total_results), names(result))
      for (col in common_col){
        # 列の型が異なる場合
        if (class(total_results[[col]])[1] != class(result[[col]])[1]) {
          total_results[[col]] <- as.double(as.character(total_results[[col]]))  # total_best.results の factor を character に変換
          result[[col]] <- as.double(as.character(result[[col]]))  # best.result の factor を character に変換
        }
      }
      total_results <- bind_rows(total_results,result)
    }
    
    # モデルが正常に作成された場合
    if (!is.null(model)) {
      # 収束状況の確認とconvergeの格納
      if (length(model_summary[["optinfo"]][["conv"]][["lme4"]]) == 0){
        print("full-Model has converged.")
        converge[r] <- TRUE
      } else {
        print("full-Model did not converge.")
        converge[r] <- FALSEsc
      }
    }
  }
  
  fit <- cbind(distance,total_results,sd_ci,converge)
  
  write.csv(fit,paste0("Model_Results/SEM/All_Plots/Summary_",mname,"_",LeafType[j],".csv"),fileEncoding = "Shift-jis",row.names = F)
  write.csv(total_FE,paste0("Model_Results/SEM/All_Plots/FE_",mname,"_",LeafType[j],".csv"),fileEncoding = "Shift-jis",row.names = F)
}

# predict
# g <- ggpredict(model,terms = "two_EB")
# ggplot() + 
#   geom_point(data = d, aes(x = two_DB, y = RGR20)) + 
#   geom_line(data = g, aes(x = x, y = predicted)) + 
#   geom_ribbon(data = g, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.1)

