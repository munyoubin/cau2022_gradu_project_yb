library(ggplot2)
library(dplyr)
library(vegan)

# 언니가 준 걸로 그리기
bray_pc_m <- read.csv('bray_m.csv', row.names = 1)

## 묶음 별로 뽑기
swine_nasal <- bray_pc_m %>% filter(Env=="Swine_Nasal")
swine_fecal <- bray_pc_m %>% filter(Env=="Swine_Fecal")
human_na_fe <- bray_pc_m %>% filter(Env=="Human_Nasal"|Env=="Human_Fecal")
nasal <- bray_pc_m %>% filter(Env=="Human_Nasal"|Env=="Swine_Nasal")
fecal <- bray_pc_m %>% filter(Env=="Human_Fecal"|Env=="Swine_Fecal")
finishing <- bray_pc_m %>%  filter(Env=="Floor_Finishing"|Site=="Finishing_Fecal"|Site=="Finishing_Nasal")
lactating <- bray_pc_m %>%  filter(Env=="Floor_Lactating"|Site=="Lactating_Fecal"|Site=="Lactating_Nasal")
swine_na_fe <- bray_pc_m %>%  filter(Env=="Swine_Nasal"|Env=="Swine_Fecal")

##anosim
anosim(distance, grouping=bray_pc_m$Env_Site, permutations = 999, distance="bray") # 0.001
anosim(distance, grouping=bray_pc_m$Site, permutations = 999, distance="bray") # 0.001
anosim(distance, grouping=bray_pc_m$Organism, permutations = 999, distance="bray") # 0.001
anosim(distance, grouping=bray_pc_m$Env, permutations = 999, distance="bray") # 0.001

#env_site
ggplot(data=bray_pc_m, aes(x=PC1, y=PC2, color=Env_Site))+
  geom_point(size=2)+
  stat_ellipse(geom = "polygon",
               aes(fill = Env_Site), 
               alpha = 0.25)+
  theme_classic()+
  theme(axis.text.x = element_text(size=14, face="bold"), 
        axis.text.y = element_text(size=14, face="bold"))+
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))


#site
ggplot(data=bray_pc_m, aes(x=PC1, y=PC2, color=Site))+
  geom_point(size=2)+
  stat_ellipse(geom = "polygon",
               aes(fill = Site), 
               alpha = 0.25)+
  theme_classic()+
  theme(axis.text.x = element_text(size=14, face="bold"), 
        axis.text.y = element_text(size=14, face="bold"))+
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))

#Organism
ggplot(data=bray_pc_m, aes(x=PC1, y=PC2, color=Organism))+
  geom_point(size=2)+
  stat_ellipse(geom = "polygon",
               aes(fill = Organism), 
               alpha = 0.25)+
  theme_classic()+
  theme(axis.text.x = element_text(size=14, face="bold"), 
        axis.text.y = element_text(size=14, face="bold"))+
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))

#Env
ggplot(data=bray_pc_m, aes(x=PC1, y=PC2, color=Env))+
  geom_point(size=2)+
  stat_ellipse(geom = "polygon",
               aes(fill = Env), 
               alpha = 0.25)+
  theme_classic()+
  theme(axis.text.x = element_text(size=14, face="bold"), 
        axis.text.y = element_text(size=14, face="bold"))+
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))

