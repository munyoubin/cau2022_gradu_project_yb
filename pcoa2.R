library(ggplot2)
library(dplyr)
library(gridExtra)
library(pheatmap)
library(lmtest)
library(corrplot)
library(heatmaply)
library(vegan)

m <- read.csv("bray_m.csv",check.names=FALSE, row.names=1)
m <- m[c("Sample", "Env_Site","Organism", "Env", "Site")]
m<-m[-92,]

total = read.csv("total_re.csv", check.names=FALSE, row.names=1)
total_name = read.csv("total_gene.csv", check.names = F, row.names = 1)
names(total) <- names(total_name)
re_name = names(sample)

total_name_sort <- total[re_name]
ani_env_m <- m[m$Env == "Swine" | m$Env == "Env", c("Sample", "Env")]
ani_hu_m <- m[m$Env == "Swine" | m$Env == "Human", c("Sample", "Env")]
env_hu_m <- m[m$Env == "Env" | m$Env == "Human", c("Sample", "Env")]

row.names(df2)


ani_efdi2 <- vegdist(ani_env_df, method="bray")

total_del <- total_name_sort[,-c(93, 94)] #char 삭제해야지 제대로 num로 traspose 됨
df2 <- t(total_del)

ani_endf<- df[1:58,]
ani_endi <- vegdist(ani_endf, method="bray")
anienvdist.pcoa <- cmdscale(ani_endi, eig=TRUE)
ani_env <- anienvdist.pcoa$points

ani_hudf<- df[c(1:45, 59:92),]
ani_hudi <- vegdist(ani_hudf, method="bray")

ani_endi <- vegdist(ani_endf, method="bray")
anienvdist.pcoa <- cmdscale(ani_endi, eig=TRUE)
ani_env <- anienvdist.pcoa$points
ani_hu_m <- ani_hu_m[-91,]

env_hu <- df[47: 92,]
env_hudi <- vegdist(env_hu, method="bray")
anosim(env_hudi, grouping = env_hu_m$Env)
anosim(ani_hudi, grouping = ani_hu_m$Env)
anosim(ani_endi, grouping = ani_env_m$Env)

ani_hu_m <- ani_hu_m[-80, ]

