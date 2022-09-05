df = read.delim("@@File@@", header=FALSE)
head(df)
l= df$V1
k= df$V2
df = data.frame(l,k)
names(df) = c("length","count")
library("ggpubr")
library(svglite)
res <- 144
svglite("@@Output@@", width = 1080/res, height = 720/res)
ggplot(data = df, aes(x = length, y = count)) + geom_point(color = "dark blue") + xlab("length") + ylab("count") + theme_classic()
dev.off()