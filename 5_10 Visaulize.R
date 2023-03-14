library(ggplot2)
library(dplyr)
library(gridExtra)
library(pheatmap)
library(lmtest)
library(corrplot)
library(heatmaply)
install.packages("devtools")
library(devtools)
install.packages("ggbiplot")
library(ggbiplot)
install.packages('pkgload')
install.packages('ggfortify')
library(ggfortify)

sample = read.csv("finish2.csv", check.names=FALSE, row.names = 1)
class_sample = read.csv("all.csv", row.names=1, stringsAsFactors = TRUE)
total_name = read.csv("total_gene.csv", check.names = F, row.names = 1)
total = read.csv("total_re.csv", check.names=FALSE, row.names=1)
full = read.csv("full.csv", row.names=1)

names(total) = names(total_name)


# animal, env, human 3개로 그룹화
sample.a = subset(class_sample, class_sample$class=="animal")
sample.e = subset(class_sample, class_sample$class=="environment")
sample.h = subset(class_sample, class_sample$class=="human")

table(class_sample$class)

qqnorm(sample.a$MDR)
qqline(sample.a$MDR)

dwtest(data=class_sample, MDR~class)
dwtest(data=class_sample, MLS~class)
dwtest(data=class_sample, aminoglycoside.antibiotic~class)


sha.a = lapply(sample.a[,2:20], shapiro.test)
unlist(lapply(sha.a, function(x) x$p.value)) # 정규분포x

sample_cor = round(cor(class_sample[,2:20], method="spearman"),2)
sample.a_cor = round(cor(sample.a[,2:20], method="spearman"),2) 
total_cor = round(cor(total[1:92], method="spearman"),2)

corrplot(total_cor, order="hclust", 
         tl.col="black", cl.pos="b",
         number.cex = 0.2, tl.cex=0.35, addrect=5,rect.lwd = 4,
         method="color", title="All samples corrplot", mar=c(0,0,1,1))

sample.a_cor
corrplot(sample_cor, type="lower", order="hclust", 
         addCoef.col="black", tl.col="black",
         tl.srt=45, number.cex = 0.3, tl.cex=0.8, 
         method="color", title="All samples corrplot", mar=c(0,0,0.5,0))
corrplot(sample.a_cor, type="lower", order="hclust", 
         addCoef.col="black", tl.col="black",
         tl.srt=45, number.cex = 0.3, tl.cex=0.8, 
         method="color", title="Animal samples corrplot", mar=c(0,0,0.5,0))

library(factoextra)
library(FactoMineR)


Principal_Component = PCA(total[1:92],graph = FALSE)

df <- total[1:92]
pca = prcomp(df, center = T, scale=T, retx = T)
summary(pca)
autoplot(pca, data=total)
options(ggrepel.max.overlaps = Inf)
pca_numb = pca$rotation[,c(1,2)]
pca_12 = as.data.frame(pca_numb)

df <- full[,c(2:526)]
new_pca <- prcomp(df, center =T, scale. =T)
summary(new_pca)
autoplot(new_pca, data=full, colour='Site')

#pca csv 만들기
rownames(pca_12)
pca_12$samples = rownames(pca_12)
rownames(pca_12) <- NULL
pca_12 = pca_12[,c(3,1,2)]
pca_meta = read.csv('pca_meta.csv', row.names = 1)

write.csv(pca_12, "pca.csv", row.names = T)
pca_numb[,1]
fviz_pca_var(Principal_Component, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T, tl.cex=0.05
)

plot(pca$x[,c(1,2)])

# 언니가 주신거
head(total)
str(total)
head(total)
head(df)
total_del <- total[,-c(93, 94)] #char 삭제해야지 제대로 num로 traspose 됨
df <- t(total_del)
str(df)
head(total)
library(vegan)
distance <- vegdist(df, method="bray")
head(distance)
dist.pcoa <- cmdscale(distance, eig=TRUE)
head(dist.pcoa)
dist.wa <- wascores(dist.pcoa$points[,1:2], df.t)
summary(dist.pcoa)
bray_pc <- dist.pcoa$points

write.csv(bray_pc, "bray_pc.csv", row.names = T)

distance2 <- vegdist(df, method="jaccard")
dist2.pcoa <- cmdscale(distance2, eig=T)
jaccard_pc <- dist2.pcoa$points

write.csv(jaccard_pc, "jaccard_pc.csv", row.names =T)

# 언니가 준 걸로 그리기
bray_pc_m <- read.csv('bray_m.csv', row.names = 1)
ggplot(data=bray_pc_m, aes(x=PC1, y=PC2, color=Site))+
  geom_point(size=2)+
  theme(axis.text.x = element_text(size=14, face="bold"), 
        axis.text.y = element_text(size=14, face="bold"))+
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))+
  ggtitle("Bray PCoA", size=16)

jaccard_pc_m <- read.csv('jaccard_m.csv', row.names=1)
ggplot(data=jaccard_pc_m, aes(x=PC1, y=PC2, color=Site))+
  geom_point(size=2)+
  theme(axis.text.x = element_text(size=14, face="bold"), 
        axis.text.y = element_text(size=14, face="bold"))+
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))+
  ggtitle("Jaccard PCoA")

# pca plot 그리기 + metadata
ggplot(data=pca_meta, aes(x=PC1, y=PC2, color=Env), group=Site)+
  geom_point(size=2)+
  theme(axis.text.x = element_text(size=14, face="bold"), 
        axis.text.y = element_text(size=14, face="bold"))+
  theme(axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))

autoplot(pca)
ggbiplot(pca_meta[,c(2,3)], ellips=T, labels = pca_meta$Sample, groups=pca_meta$Site)

? pheatmap
my_palette <- colorRampPalette(c("#f1eef6","#d4b9da",
                                 "#c994c7","#df65b0", 
                                 "#e7298a", "#ce1256", 
                                 "#91003f"))(n = 699)

col_breaks = c(seq(0.0,0.5, length=100),  # for red
               seq(0.5000000001,5, length=100),
               seq(5.0000000001,20, length=100),
               seq(20.0000000001, 80, length=100),
               seq(80.00000000001, 160, length=100),
               seq(160.0000000001, 300, length=100),# for yellow
               seq(300.00000001, 1630, length=100))          # for green

pheatmap(sample, dendrogram ="none", cluster_rows = FALSE, cluster_cols = FALSE, fontsize_row = 8, fontsize_col = 7,scale="none", na.color="grey", border_color = "black",
         trace = "none", col = my_palette, breaks = col_breaks, gaps_col = c(46, 58))

t_s = as.data.frame(t(sample)) # transpose

boxplot(t_s$MDR)$stats
names(t_s)
