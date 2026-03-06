functional <- c("understory","canopy","non-canopy")
SP <- c("canopy","canopy","non-canopy")
layer <- c("understory","overstory","understory")
leaf <- c("EB","DB")

for(i in 1:3){
  for(j in 1:2){
    data <- read.csv(paste0("result/glmm/Summary_",functional[i],"_",leaf[j],".csv"),header = T) %>%
      mutate(distance = rep(1:50, each = 1),LT = leaf[j],FT = functional[i],SP = SP[i],layer = layer[i])
    if(i == 1 & j == 1){
      result1 <- data
    }else{
      result1 <- bind_rows(result1,data)
    }
  }
}

SP <- c("Aa","Qa")
species <- c("Aamoenum","Qacuta")
layer <- c("understory","overstory")
functional <- c("understory","canopy")
for(i in 1:2){
  for(j in 1:2){
    data <- read.csv(paste0("result/glmm/Summary_",species[i],"_",functional[j],".csv"),header = T) %>%
      subset(Row.names == "two_DB" | Row.names == "two_EB") %>%
      mutate(distance = rep(1:50, each = 2),LT = leaf[i],FT = functional[j],SP = SP[i],layer = layer[j])
    if(i == 1 & j == 1){
      result2 <- data
    }else{
      result2 <- bind_rows(result2,data)
    }
  }
}
result <- bind_rows(result1,result2)
write.csv(result,"result/summary/distance.csv",row.names = F,fileEncoding = "utf-8")
