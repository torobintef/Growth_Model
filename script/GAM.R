
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

for(p in c("ALT600","ALT700","ALT800")){
  for(j in 1:2){
    for(r in 1:50){
      distance[r] <- r
      d1 <- read.csv(paste0("Dd_20_edge/Dd600_r",r,"_1.csv"),header = T) %>%
        mutate(RGR20 = (log(DBH/DBH05))/19*100) %>%
        mutate(RGRa = (log(DBH15/DBH05))/10*100) %>%
        mutate(RGRb = (log(DBH/DBH15))/9*100) %>%
        left_join(d600,by = "STEM_ID") %>%
        group_by(IND_ID) %>%
        mutate(DBH = max(DBH)) %>%
        distinct(IND_ID,.keep_all = TRUE) %>%
        ungroup()
      d1$PLOT <- factor("ALT600")
      d2 <- read.csv(paste0("Dd_20_edge/Dd700_r",r,"_1.csv"),header = T) %>%
        mutate(RGR20 = (log(DBH/DBH05))/17*100) %>%
        mutate(RGRa = (log(DBH15/DBH05))/10*100) %>%
        mutate(RGRb = (log(DBH/DBH15))/7*100) %>%left_join(d700,by = "STEM_ID") %>%
        group_by(IND_ID) %>%
        mutate(DBH = max(DBH)) %>%
        distinct(IND_ID,.keep_all = TRUE) %>%
        ungroup()
      d2$PLOT <- factor("ALT700")
      d3 <- read.csv(paste0("Dd_20_edge/Dd800_r",r,"_1.csv"),header = T) %>%
        mutate(RGR20 = (log(DBH/DBH05))/19*100) %>%
        mutate(RGRa = (log(DBH15/DBH05))/10*100) %>%
        mutate(RGRb = (log(DBH/DBH15))/9*100) %>%
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
        select(-intersect(names(td1), names(td2)))
      
      d <- bind_cols(td1, td2_unique) %>%
        subset(!is.na(RGR20) & !is.na(DBH05)& !is.na(ALT) &
                 !is.na(one_EB) & !is.na(one_DB) & !is.na(one_EC) &
                 !is.na(two_EB) & !is.na(two_DB) & !is.na(two_EC)) %>%
        subset(SP == "アカガシ" | SP == "イヌガシ" | SP == "ウラジロガシ" |
                 SP == "ブナ" | SP == "ケヤキ" | SP == "ヒメシャラ" |
                 SP == "イヌシデ" | SP == "オオモミジ" | SP == "イタヤカエデ" | SP == "カジカエデ") %>%
        subset(LT == LeafType[j] & LD == 0 & DBH05 < 30 & RGR20 >= 0 & PLOT == p)
      sd_value <- data.frame(sd.one_DB = sd(log(d$one_EB + 0.01)),sd.one_EB = sd(log(d$one_DB + 0.01)),
                             sd.two_DB = sd(log(d$two_EB + 0.01)),sd.two_EB = sd(log(d$two_DB + 0.01)))
      
      # データの対数変換と正規化(z-score)
      d$RGR20 <- as.numeric(scale(d$RGR20))
      d$DBH05 <- as.numeric(scale(log(d$DBH05)))
      d$one_EB <- as.numeric(scale(log(d$one_EB + 0.01)))
      d$one_DB <- as.numeric(scale(log(d$one_DB + 0.01)))
      d$two_EB <- as.numeric(scale(log(d$two_EB + 0.01)))
      d$two_DB <- as.numeric(scale(log(d$two_DB + 0.01)))
      d$ALT <- as.numeric(scale(d$ALT))
      d$slope <- as.numeric(scale(d$slope))
      d$TPI <- as.numeric(scale(d$TPI))
      
      
      
      # model_summary
      print(paste0("---Now caliculating ",mname," model of ",LeafType[j]," at ",p," for r = ",r,"----------"))
      
      model <- gam(RGR20 ~ s(DBH05,k = 10) + s(two_EB,k = 10) + s(two_DB,k = 10) +
                     s(ALT,k = 10) + s(slope,k = 10) + s(TPI,k = 10) + s(X,Y,by = PLOT,k = 10), data = d)
      base_model <- gam(RGR20 ~ s(DBH05,k = 10) + s(ALT,k = 10) + 
                          s(slope,k = 10) + s(TPI,k = 10),data = d)
      model_summary <- summary(model)
      
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
      # AIC <- model_summary[["aic"]]
      AIC <- AIC(model)
      deltaAIC <- AIC(base_model) - AIC(model)
      logLik <- logLik(model)
      # 回帰係数をデータフレームに格納し、行と列を反転
      result <- as.data.frame(cbind(AIC,deltaAIC,logLik))
      
      # best.model（収束した中で最もAICが低い最良モデル）
      # best.result（最良モデルの回帰係数とAIC，logLik）
      # best.FE（最良モデルにおける各説明変数のP値とVIF値）
      # dredge_test <- dredge(model,rank = "AIC")
      # t <- 0
      # s <- 1
      # while(s != 0){
      #   t <- t + 1
      #   model.t <- get.models(dredge_test, t)[[1]]
      #   model.t_summary <- summary(model.t)
      #   s <- length(model.t_summary[["optinfo"]][["conv"]][["lme4"]])
      # }
      # best.model <- get.models(dredge_test, t)[[1]]
      # best.result <- as.data.frame(dredge_test[t,])
      # best.summary <- model.t_summary
      # best.fixed <- as.data.frame(best.summary$coefficients)
      # best.colltest <- check_collinearity(best.model)
      # rownames(best.colltest) <- best.colltest[[1]]
      # best.colltest <- best.colltest[,-1]
      # param <- as.data.frame(model_parameters(best.model))
      # param <- head(param, n = nrow(param) - 3)
      # rownames(param) <- param[[1]]
      # param <- param[,-1]
      # best.FE <- merge(best.fixed,best.colltest, by = "row.names", all = T)
      # rownames(best.FE) <- best.FE[[1]]
      # best.FE <- best.FE[,-1]
      # best.FE <- merge(best.FE,param,by = "row.names",all = T)
      # 
      # result（best.result），result（best.result）の結合
      if(r == 1){
        total_FE <- FE
        # total_best.FE <- best.FE
        total_results <- result
        # total_best.results <- best.result
        sd_ci <- sd_value
      }else{
        total_FE <- bind_rows(total_FE,FE)
        # total_best.FE <- bind_rows(total_best.FE,best.FE)
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
        
        # best.common_col <- intersect(names(total_best.results), names(best.result))
        # for (col in best.common_col){
        #   # 列の型が異なる場合
        #   if (class(total_best.results[[col]])[1] != class(best.result[[col]])[1]) {
        #     total_best.results[[col]] <- as.double(as.character(total_best.results[[col]]))  # total_best.results の factor を character に変換
        #     best.result[[col]] <- as.double(as.character(best.result[[col]]))  # best.result の factor を character に変換
        #   }
        # }
        # total_best.results <- bind_rows(total_best.results,best.result)
      }
      
      # 決定係数R2の格納(conR = 固定効果とランダム効果の分散を考慮，marR = 固定効果の分散のみ)
      
      R2[r] <- model_summary[["r.sq"]]
      # best.R2[r] <- best.summary[["r.sq"]]
      # R2[r] <- r2(model)[["R2"]][1]
      # best.R2[r] <- r2(best.model)[["R2"]][1]
      
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
      # if (!is.null(best.model)) {
      #   # 収束状況の確認とconvergeの格納
      #   if (length(best.summary[["optinfo"]][["conv"]][["lme4"]]) == 0){
      #     print("Best Model has converged.")
      #     best.converge[r] <- TRUE
      #   } else {
      #     print("Best Model did not converge.")
      #     best.converge[r] <- FALSE
      #   }
      # }
      # # ランダム切片の分散varcorの格納
      # varcor <- model_summary$varcor
      # best.varcor <- best.summary$varcor
      # varcor.SP[r] <- varcor[["SP"]][1][1]
      # # varcor.MESH[r] <- varcor[["MESH"]][1][1]
      # best.varcor.SP[r] <- best.varcor[["SP"]][1][1]
      # # best.varcor.MESH[r] <- best.varcor[["MESH"]][1][1]
    }
    
    fit <- cbind(distance,total_results,sd_ci,R2,converge)
    # best.fit <- cbind(distance,total_best.results,sd_ci,best.R2,best.converge)
    
    write.csv(fit,paste0("Model_Results/GAM/",p,"/Summary_",mname,"_",LeafType[j],".csv"),fileEncoding = "Shift-jis",row.names = F)
    # write.csv(best.fit,paste0("260201result/gam2/Best_Summary_",mname,"_",LeafType[j],".csv"),fileEncoding = "Shift-jis",row.names = F)
    write.csv(total_FE,paste0("Model_Results/GAM/",p,"/FE_",mname,"_",LeafType[j],".csv"),fileEncoding = "Shift-jis",row.names = F)
    # write.csv(total_best.FE,paste0("260201result/gam2/Best_FE_",mname,"_",LeafType[j],".csv"),fileEncoding = "Shift-jis",row.names = F)
  }
}


