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
df_ <- aggregate(df$length, by=list(count=df$count), FUN="sum")
x <- df_[df_$count == max(df_$count),]$x
y <- df_[df_$count == max(df_$count),]$count
ggplot(data=df, aes(x=length, y=count)) +
       geom_vline(aes(xintercept=x)) +
       annotate("text", x=x-250, y=y, label=x) +
       scale_x_continuous(limits=c(0, x+1000)) +
       geom_line(color="blue") + xlab("Sequenced Read Length (nt)") + ylab("Read Count") + theme_classic()
dev.off()