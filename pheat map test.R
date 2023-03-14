library(ggplot2)
library(dplyr)
library(gridExtra)
library(pheatmap)
library(lmtest)
library(corrplot)
install.packages("heatmaply")
library(heatmaply)
sample = read.csv("finish2.csv", check.names=FALSE, row.names = 1)
class_sample = read.csv("all.csv", row.names=1, stringsAsFactors = TRUE)
total_name = read.csv("total_gene.csv", check.names = F, row.names = 1)
total = read.csv("total_re.csv", check.names=FALSE, row.names=1)

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

