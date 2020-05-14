
#Plot dn:ds p-values across gene data from hyphy is read from csv files in 4_ file
#For FEL plot, should remove sites with low dn:ds, since p-value could represent high or low dn:ds

plot(data=all_fel_table, (1-p.value)~Site, ylab="(1-p value)", xlab='Aligned Codon Site', ylim=c(0.85,1), cex=0.3, col="dark grey", type="h", pch=19)
par(new=TRUE)
plot(data=all_meme_table, (1-p.value)~Site, ylab="(1-p value)",xlab='Aligned Codon Site', ylim=c(0,1), cex=0.3, col="blue", type="h", pch=19)

all_meme_table <- read.csv("~/Documents/GitHub/Cypridinidae_Emission/LuciferaseTree_dNds/results/hyphy/lucclade.meme.allsites.csv")
colnames(all_meme_table) <- c("site","partition","alpha_0","beta_0","p_val_minus","beta_1","p_val_plus","LRT","p_value_MEME","branches","b_length")

#read in fel analyses
all_fel_table <- read.csv("~/Documents/GitHub/Cypridinidae_Emission/LuciferaseTree_dNds/results/hyphy/lucclade.fel.allsites.csv")
colnames(all_fel_table)[1] <- "site"
colnames(all_fel_table)[8] <- "p_value_FEL"

#read in Multihit table
multihit <- read.csv("~/Desktop/hyphy_more/datamonkey-table-nultihit.csv",col.names = c("site","three","three_island_two","three_v_three","two"))

#combine hyphy MEME and FEL results into one table
sel_dat <- merge(all_meme_table,all_fel_table,by= "site")
sel_dat <- melt(sel_dat,id.vars = "site",measure.vars = c("p_value_MEME","p_value_FEL"),value.name = "p_value")

#merge analyses by site for plotting
hyphy_dat <- merge(sel_dat,multihit,by= "site")

#plot hyphy data together
library(ggplot2); library(reshape2); library(gridExtra)
coeff <- 1
ggplot(data=hyphy_dat,aes(x=site)) + geom_line(aes(y=(1-p_value)),size=0.5) +
  geom_hline(yintercept = 0.95,col="red",lty=2) +
  geom_point(aes(y=two/coeff),size=0.5,alpha=0.5,color="blue") +
  scale_y_continuous(name="1 - P-value",sec.axis = sec_axis(~.*coeff, name="Evidence ratio of Two MNS"),limits = c(0,400))

df <- melt(hyphy_dat,id.vars="site",measure.vars=c("three","two"))
ggplot(data=df,aes(x=site,y=value)) + geom_line() + facet_grid(.~variable)

hyphy_cols <- c("blue","grey")
p1 <- ggplot(data=hyphy_dat,aes(x=site)) + #geom_point(aes(y=(1-p_value)),size=0.5) +
  geom_linerange(aes(x=site,ymax=(1-p_value),ymin=0,color=variable)) +
  geom_hline(yintercept = 0.95,col="red",lty=2) + scale_color_manual(values=hyphy_cols) +
  theme(legend.position = "none") + ylab("1 - P-value")
p2 <- ggplot(data=df,aes(x=site)) + geom_line(aes(y=log(value),color=as.factor(variable)),size=0.5,alpha=0.75) +
  geom_hline(yintercept = log(5),col="grey",lty=2) +
  geom_hline(yintercept = log(2),col="grey",lty=2) +
  theme(legend.position = "none") + ylab("Natural log (Evidence Ratio for MNS model)")

p1 <- ggplot(data=hyphy_dat,aes(x=site)) +
  geom_hline(yintercept = 0.95,col="red",lty=2) + scale_color_manual(values=hyphy_cols) +
  geom_point(aes(y=(1-p_value),color=variable),alpha=0.5) +
  geom_point(data=subset.data.frame(hyphy_dat[hyphy_dat$p_value <= 0.05,]),aes(y=(1-p_value),color=variable)) +
  ylab("1 - P-value") + scale_y_continuous(limits=c(0.5,1)) + theme_bw() + theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p2 <- ggplot(data=df,aes(x=site)) +
  geom_vline(xintercept = c(207,211,406,435,102,142,189,41,261),lty=1,col="grey") +
  geom_hline(yintercept = log(5),col="black",lty=2) +
  geom_hline(yintercept = log(2),col="black",lty=2) +
  geom_point(aes(y=log(value),color=as.factor(variable)),alpha=0.25) +
  geom_point(data=subset.data.frame(df[df$value > 1,]),aes(y=log(value),color=as.factor(variable))) + ylab("Natural log (Evidence Ratio for MNS model)") +
  theme_bw() + theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))

grid.arrange(p1,p2,ncol=1)

###VWD domains according to Pfam are from sites 78 - 231 & 317 - 466
cyp2aligned(c(78,231,317,466))
#106-260 & 348-498



