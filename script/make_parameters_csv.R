
functional <- c("understory","non-canopy")
leaf <- c("EB","DB")

for(i in 1:2){
  for(j in 1:2){
    data <- read.csv(paste0("result/glmm/FE_",functional[i],"_",leaf[j],".csv"),header = T) %>%
      subset(Row.names == "two_DB" | Row.names == "two_EB") %>%
      mutate(distance = rep(1:50, each = 2),LT = leaf[i],FT = functional[j])
    if(i == 1 & j == 1){
      result <- data
    }else{
      result <- bind_rows(result,data)
    }
  }
}
write.csv(result,"result/summary/parameters.csv",row.names = F,fileEncoding = "utf-8")
