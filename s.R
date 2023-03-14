phy = read.csv("sort_phylum.csv", check.names = F, row.names = 1)
phy = t(phy)
swine_f = phy[1:24,]
swine_n = phy[25:46,]
env = phy[47:58,]
hum_f = phy[59:66,]
hum_n = phy[67:74,]
f_bind = rbind(swine_f, hum_f)
install.packages("vegan")
library(ggplot2)

library(vegan)
help(vegan)
distance <- vegdist(f_bind, method="bray")
head(distance)
dist.pcoa <- cmdscale(distance, eig=TRUE)
head(dist.pcoa)
dist.wa <- wascores(dist.pcoa$points[,1:2], df.t)
summary(dist.pcoa)

bray_pc <- dist.pcoa$points
ggplot(data=bray_pc, aes(x=V1, y=V2))+
  geom_point(size=2)+
  theme(axis.text.x = element_text(size=14, face="bold"), 
        axis.text.y = element_text(size=14, face="bold"))+
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))+
  ggtitle("Jaccard PCoA")
