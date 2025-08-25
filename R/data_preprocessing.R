#set directory
getwd()
setwd("C:/Users/USER/Documents/pancreas")

#load libraries
library(data.table)
library(dplyr)

#load data
list.files(paste0(getwd(), "/input"))
days2= read.csv("input/FOR_METABO_FINAL_2days.csv")[-1,]
uniques= read_xlsx("input/Stat Panel MO last.xlsx" ,sheet= "Uniques Passing 25% CUTOFF") %>% as.data.frame()
names(uniques)[1]= 'Metabolite.name'
dim(uniques)

sum(is.na(uniques))
sum(uniques == "NA")
# Replace "NA" strings with NA values
uniques <- uniques %>% mutate_all(~na_if(., "NA"))

#impute uniques
impute_me= function(preImputed,
                    groups=data.frame(group1=c(1,10), group2=c(11,20))){
  names(preImputed)[1]= 'ID'
  imputed= as.data.frame(preImputed[1])
  for (i in groups) {
    sampleColumns <- seq(as.numeric(i[1]+1), as.numeric(i[2]+1))
    singleGroup <- preImputed[, c(1, sampleColumns)] %>% as.data.frame()
    for(i in 1:nrow(singleGroup[, -1])) {
      for(j in 2:ncol(singleGroup)){
        if (sum(!is.na(singleGroup[i, -1]))==2){
          imp= mean(as.numeric(unlist(singleGroup[i, -1])), na.rm = TRUE)
          if(is.na(singleGroup[i,j])){
            singleGroup[i,j] <- imp
          }
        }
        else if (sum(is.na(singleGroup[i, -1]))/length(singleGroup[i, -1]) != 1) { 
          medianVal <- median(as.numeric(unlist(singleGroup[i, -1])), na.rm = T)
          imp <- runif(1, medianVal - (medianVal*0.01), medianVal + (medianVal*0.01))
          if(is.na(singleGroup[i,j])){
            singleGroup[i,j] <- imp}
        }}}
    
    imputed <- merge(imputed, singleGroup, by = 'ID', all.x = TRUE)
  }
  return(imputed)
}
#run
imputed_uniques= impute_me(uniques)
imputed_uniques[is.na(imputed_uniques)]= 0
imputed_uniques %>% names(.)


days2= read.csv("input/FOR_METABO_FINAL_2days.csv")[-1,]
days2_ordered= cbind(ID= days2$Sample,
                     days2[,grep("B", names(days2))],
                     days2[,grep("A", names(days2))])
days2_ordered %>% names(.)
days2_all= rbind(days2_ordered,imputed_uniques )

days7= read.csv("input/FOR_METABO_FINAL_7days.csv" )[-1,]
names(days7)[1]= "ID"
data= Reduce(function(x,y) inner_join(x, y, by= "ID"), list(days2_all, days7))
data[is.na(data)]=0

data %>% names()

write.csv(data, "input/pancreas_alldata.csv", row.names = F)


