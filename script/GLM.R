
# packages ----------------------------------------------------------------

# renv::install("tidyverse")
# renv::install("lme4")
# renv::install("car")
# renv::install("spdep")
# renv::install("performance")
# renv::install("parameters")

library(tidyverse)
library(lme4)
library(car)
library(spdep)
library(performance)
library(parameters)
# library(broom)
# library(MuMIn)
# options(na.action="na.fail")
# library(ggeffects)
# library(mgcv) #gam


# repeat ------------------------------------------------------------------


mname <- "canopy"

# 空のベクトルまたはリストを作成
LeafType <- c("EB","DB")
distance <- vector()
best.R2 <- vector()
converge <- vector()
best.converge <- vector()
varcor.SP <- vector()
varcor.MESH <- vector()
best.varcor.SP <- vector()
best.varcor.MESH <- vector()
R2 <- vector()
d600 <- read.csv("data/KN600.csv")[c("STEM_ID","ALT","slope","TPI")]
d700 <- read.csv("data/KN700.csv")[c("STEM_ID","ALT","slope","TPI")]
d800 <- read.csv("data/KN800.csv")[c("STEM_ID","ALT","slope","TPI")]

for(j in 1:2){
  for(r in 1:50){
    
    distance[r] <- r
    d1 <- read.csv(paste0("Dd_20_edge/Dd600_r",r,"_1.csv"),header = T) %>%
      mutate(RGR20 = DBH/DBH05) %>%
      mutate(RGRa = DBH15/DBH05) %>%
      mutate(RGRb = DBH/DBH15) %>%
      mutate(t = 19) %>%
      mutate(IDBH = DBH-DBH05) %>%
      left_join(d600,by = "STEM_ID") %>%
      group_by(IND_ID) %>%
      mutate(DBH = max(DBH)) %>%
      distinct(IND_ID,.keep_all = TRUE) %>%
      ungroup()
    d1$PLOT <- factor("ALT600")
    d2 <- read.csv(paste0("Dd_20_edge/Dd700_r",r,"_1.csv"),header = T) %>%
      mutate(RGR20 = DBH/DBH05) %>%
      mutate(RGRa = DBH15/DBH05) %>%
      mutate(RGRb = DBH/DBH15) %>%
      mutate(t = 17) %>%
      mutate(IDBH = DBH-DBH05) %>%
      left_join(d700,by = "STEM_ID") %>%
      group_by(IND_ID) %>%
      mutate(DBH = max(DBH)) %>%
      distinct(IND_ID,.keep_all = TRUE) %>%
      ungroup()
    d2$PLOT <- factor("ALT700")
    d3 <- read.csv(paste0("Dd_20_edge/Dd800_r",r,"_1.csv"),header = T) %>%
      mutate(RGR20 = DBH/DBH05) %>%
      mutate(RGRa = DBH15/DBH05) %>%
      mutate(RGRb = DBH/DBH15) %>%
      mutate(t = 19) %>%
      mutate(IDBH = DBH-DBH05) %>%
      left_join(d800,by = "STEM_ID") %>%
      group_by(IND_ID) %>%
      mutate(DBH = max(DBH)) %>%
      distinct(IND_ID,.keep_all = TRUE) %>%
      ungroup()
    d3$PLOT <- factor("ALT800")
    
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
      subset(!is.na(RGR20) & !is.na(DBH05)& !is.na(ALT) & !is.na(IDBH) &
               !is.na(one_EB) & !is.na(one_DB) & !is.na(one_EC) &
               !is.na(two_EB) & !is.na(two_DB) & !is.na(two_EC)) %>%
      subset(SP == "アカガシ" | SP == "イヌガシ" | SP == "ウラジロガシ" |
               SP == "ブナ" | SP == "ケヤキ" | SP == "ヒメシャラ" |
               SP == "イヌシデ" | SP == "オオモミジ" | SP == "イタヤカエデ" | 
               SP == "カジカエデ") %>%
      subset(LT == LeafType[j] & LD == 0 & DBH05 >= 30 & IDBH >= 0)
    sd_value <- data.frame(sd.one_DB = sd(d$one_EB),sd.one_EB = sd(d$one_DB),
                           sd.two_DB = sd(d$two_EB),sd.two_EB = sd(d$two_DB))
    
    # データの正規化(z-score)
    d$IDBH <- as.numeric(d$IDBH + 0.01)
    d$RGR20 <- as.numeric(d$RGR20)
    d$DBH05 <- as.numeric(scale(log(d$DBH05)))
    d$one_EB <- as.numeric(scale(log(d$one_EB + 0.01)))
    d$one_DB <- as.numeric(scale(log(d$one_DB + 0.01)))
    d$one_all <- as.numeric(scale(log(d$one_all + 0.01)))
    d$two_EB <- as.numeric(scale(log(d$two_EB + 0.01)))
    d$two_DB <- as.numeric(scale(log(d$two_DB + 0.01)))
    d$two_all <- as.numeric(scale(log(d$two_all + 0.01)))
    d$ALT <- as.numeric(scale(log(d$ALT)))
    d$slope <- as.numeric(scale(log(d$slope)))
    d$TPI <- as.numeric(scale(log(d$TPI+1)))
    
    # 近接行列
    coords <- as.matrix(d[, c("X", "Y")])
    nb <- dnearneigh(coords, d1 = 0, d2 = 20)
    lw <- nb2listw(nb, style = "W",zero.policy = T)
    
    
    # model_summary
    print(paste0("---Now caliculating ",mname," model of ",LeafType[j]," for r = ",r,"----------"))
    
    base_formula <- IDBH ~ DBH05 + slope + TPI + PLOT
    formula <- IDBH ~ DBH05 + TPI + slope + PLOT + two_EB + two_DB
    
    base_model <- glm(base_formula, offset = log(t), family = Gamma(link = "log"), data = d)
    model <- glm(formula, offset = log(t), family = Gamma(link = "log"), data = d)
    
    model_summary <- summary(model)
    res <- residuals(model)
    mtest <- moran.test(res,lw)
    
    # FE(各固定効果のP値とVIF値)
    fixed <- as.data.frame(model_summary$coefficients)
    colltest <- check_collinearity(model)
    rownames(colltest) <- colltest[[1]]
    colltest <- colltest[,-1]
    param <- as.data.frame(model_parameters(model))
    rownames(param) <- param[[1]]
    param <- param[,-1]
    FE <- merge(fixed,colltest, by = "row.names", all = TRUE)
    rownames(FE) <- FE[[1]]
    FE <- FE[,-1]
    FE <- merge(FE,param, by = "row.names", all = TRUE)
    
    # result(モデルの回帰係数とAIC,logLik)
    AIC <- AIC(model)
    deltaAIC <- AIC(base_model) - AIC(model)
    logLik <- logLik(model)
    M_stat <- mtest[["statistic"]]
    M_p.value <- mtest[["p.value"]]
    
    
    # 回帰係数をデータフレームに格納し、行と列を反転
    result <- as.data.frame(cbind(AIC,deltaAIC,logLik,M_stat,M_p.value))
    
    # result（best.result），result（best.result）の結合
    if(r == 1){
      total_FE <- FE
      total_results <- result
    }else{
      total_FE <- bind_rows(total_FE,FE)
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
    
    # 決定係数R2の格納(conR = 固定効果とランダム効果の分散を考慮，marR = 固定効果の分散のみ)
    
    # R2[r] <- model_summary[["r.sq"]]
    # best.R2[r] <- best.summary[["r.sq"]]
    R2[r] <- r2(model)[["R2_Nagelkerke"]][["Nagelkerke's R2"]]
    
    # モデルが正常に作成された場合
    if (!is.null(model)) {
      # 収束状況の確認とconvergeの格納
      if (length(model_summary[["optinfo"]][["conv"]][["lme4"]]) == 0){
        print("full-Model has converged.")
        converge[r] <- TRUE
      } else {
        print("full-Model did not converge.")
        converge[r] <- FALSE
      }
    }
  }
  
  fit <- cbind(distance,total_results,R2,converge)
  
  write.csv(fit,paste0("result/glmm/Summary_",mname,"_",LeafType[j],".csv"),fileEncoding = "Shift-jis",row.names = F)
  write.csv(total_FE,paste0("result/glmm/FE_",mname,"_",LeafType[j],".csv"),fileEncoding = "Shift-jis",row.names = F)
}
