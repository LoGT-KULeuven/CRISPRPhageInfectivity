library(ggplot2)
library(UpSetR)
library(reshape2)
library(scales)
library(ggpubr)
library(rstatix)
library(mltools)

# Adjust to your settings
setwd(dir = "/mnt/sequencing/projects/ISMEJCrispr/R/R_latest/")


##########################################################
## Part I. Dataset with phage infectivity (n=173)       ##
## Link between defense systems and phage infectivity   ##
##########################################################

crispr.pi <- read.csv("crispr_phage_infectivity.csv")
crispr.pi.nooutlier <- subset(crispr.pi, TotalDefenseGenes>=130) # 1 outlier strain has very few defense systems
ggplot(crispr.pi.nooutlier, aes(y=PhageProd, x=TotalDefenseGenes, group=Presence)) +
  geom_point(aes(colour=Presence)) +
  geom_smooth(method="lm", aes(colour=Presence)) +  
  theme_classic() +
  theme(strip.text.x=element_blank(), axis.ticks.x = element_blank(),
        axis.title=element_blank()) +
  scale_fill_manual(values=c("#7570B3","#1B9E77")) +
  scale_colour_manual(values=c("#7570B3","#1B9E77"))

# calculation of linear models
crispr.pi.nooutlier.nocrispr <- subset(crispr.pi.nooutlier, Presence == "Not Present")
crispr.pi.nooutlier.crispr <- subset(crispr.pi.nooutlier, Presence == "Present")

lm.nocrispr <- lm(PhageProd ~ TotalDefenseGenes, data=crispr.pi.nooutlier.nocrispr)
lm.crispr <- lm(PhageProd ~ TotalDefenseGenes, data=crispr.pi.nooutlier.crispr)

print(lm.nocrispr)
print(lm.crispr)

# Link presence CRISPR and Defense systems
ds <- c("Abi", "BREX", "DISARM",
        "DND", "Druantia", "Gabija", "Hachiman",
        "Lamassu", "Pagos", "RM", "Septu",
        "Shedu", "TA", "Thoeris", "Wadjet", "Zorya")

crispr.melt <- melt(crispr.pi, id.var="Presence", measure.vars=ds)

ggplot(data=crispr.melt, aes(x=variable, y=value, fill=Presence)) +
  facet_wrap(~variable, scales="free") +
  geom_boxplot(aes(fill=Presence), 
               outlier.color = NA, 
               position=position_dodge(width=1), 
               alpha=0.7,
               width=0.5) +
  geom_point(aes(colour=Presence), 
             position=position_jitterdodge(dodge.width = 1, jitter.width = 0.05, jitter.height=0.15), 
             shape=21, size=1, alpha=0.5) +
  scale_fill_manual(values=c("#7570B3","#1B9E77")) +
  scale_colour_manual(values=c("#7570B3","#1B9E77")) +
  theme_classic() +
  theme(strip.text.x=element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.title=element_blank()) +
  scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
  stat_compare_means(aes(label=..p.signif..)) 

crispr.yes <- subset(crispr.pi, Presence == "Present")
crispr.no <- subset(crispr.pi, Presence == "Not Present")

#####################################################
### Link between CRISPR+ and Productive infection ###
#####################################################

ggplot(data = crispr.pi, aes(x=Presence, y=PhageProd)) +
  geom_boxplot(aes(fill=Presence), 
               outlier.color = NA, 
               alpha=0.7,
               width=0.4) +
  geom_jitter(aes(colour=Presence), position=position_jitter(0.1), alpha=.4) +
  theme_classic() +
  theme(strip.text.x=element_blank(), axis.ticks.x = element_blank(),
        axis.title=element_blank()) +
  scale_fill_manual(values=c("#7570B3","#1B9E77")) +
  scale_colour_manual(values=c("#7570B3","#1B9E77")) +
  stat_compare_means(aes(label=..p.signif..)) 


#####################################################
## Part II: extended analysis of CRISPR/Population ##
##  Including resolved genomes (n=468)             ##
#####################################################
crispr.nopi <- read.csv("crispr_no_phage_infectivity.csv")

# Entire dataset
crispr.upset <- read.csv("crispr_upset.csv")
names(crispr.upset) <- c("Type I-E", "Type I-F", "Type I-C", "No CRISPR")
upset(crispr.upset,
      sets.bar.color = "#FF8400", 
      mainbar.y.label = "CRISPR presence per type", 
      sets.x.label = "CRISPR type",
      main.bar.color = "#0069FF",
      text.scale = c(1.3, 1.3, 1.3, 1.3, 1.3, 1.3))

# Upset population group 1
crispr.upset.C1 <- read.csv("crispr_upset_C1.csv")
names(crispr.upset.C1) <- c("Type I-E", "Type I-F", "Type I-C", "No CRISPR")
upset(crispr.upset.C1,
      sets.bar.color = "#69696979", 
      mainbar.y.label = "CRISPR presence per type", 
      sets.x.label = "CRISPR type",
      main.bar.color = "#e41a1c79",
      text.scale = c(1.3, 1.3, 1.3, 1.3, 1.3, 1.3))

# Upset population group 2
crispr.upset.C2 <- read.csv("crispr_upset_C2.csv")
names(crispr.upset.C2) <- c("Type I-E", "Type I-F", "Type I-C", "No CRISPR")
upset(crispr.upset.C2,
      sets.bar.color = "#69696979", 
      mainbar.y.label = "CRISPR presence per type", 
      sets.x.label = "CRISPR type",
      main.bar.color = "#377eb879",
      text.scale = c(1.3, 1.3, 1.3, 1.3, 1.3, 1.3))

########################################
### Genome size vs. population group ###
########################################
groups <- subset(crispr.nopi, Clade == 1 | Clade ==2)
groups[,"Clade"] <- factor(groups[,"Clade"])

groupColors <- c("1"="#e41a1c79", "2"="#377eb879")
ggplot(data = groups, aes(x=Length, fill=Clade)) +
  geom_density(alpha=.4)+
  scale_fill_manual(values=groupColors) +
  theme_classic() +
  theme(strip.text.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.title=element_blank(),  axis.line.y = element_blank())

# Is there still an association between genome size and CRISPR when looking within groups?
# group 1
group1 <- subset(crispr.nopi, Clade == 1)
group2 <- subset(crispr.nopi, Clade == 2)

# normality test (shapiro)
with(group1, shapiro.test(Length[Presence=="Present"])) # not normally distributed
with(group1, shapiro.test(Length[Presence=="Not Present"])) # not normally distributed

wilcox.test(Length~Presence, data=group1)
wilcox.test(Length~Presence, data=group2)

with(group2, shapiro.test(Length[Presence=="Present"]))
with(group2, shapiro.test(Length[Presence=="Not Present"]))

ggplot(data = subset(crispr.nopi, Clade == 1), aes(x=Presence, y=Length)) +
  geom_boxplot(aes(fill=Presence), 
               outlier.color = NA, 
               alpha=0.7,
               width=0.4) +
  geom_jitter(aes(colour=Presence), position=position_jitter(0.1), alpha=.4) +
  theme_classic() +
  theme(strip.text.x=element_blank(), axis.ticks.x = element_blank(),
        axis.title=element_blank()) +
  scale_fill_manual(values=c("#7570B3","#1B9E77")) +
  scale_colour_manual(values=c("#7570B3","#1B9E77")) +
  stat_compare_means(aes(label=..p.signif..)) 

# group2
ggplot(data = subset(crispr.nopi, Clade == 2), aes(x=Presence, y=Length)) +
  geom_boxplot(aes(fill=Presence), 
               outlier.color = NA, 
               alpha=0.7,
               width=0.4) +
  geom_jitter(aes(colour=Presence), position=position_jitter(0.1), alpha=.4) +
  theme_classic() +
  theme(strip.text.x=element_blank(), axis.ticks.x = element_blank(),
        axis.title=element_blank()) +
  scale_fill_manual(values=c("#7570B3","#1B9E77")) +
  scale_colour_manual(values=c("#7570B3","#1B9E77")) +
  stat_compare_means(aes(label=..p.signif..))

########################################
## Part III. RefSeq dataset (n=4812)  ##   
## CRISPR-Acr effect on genome size   ##
########################################
ncbi.full <- read.csv("ncbi_full.csv")

# Dataset exploration (how many strains with CRISPR, type)
n <- 4812
n.IC <- sum(ncbi.full$IC > 0 & ncbi.full$IE == 0 & ncbi.full$IF == 0)
n.IE <- sum(ncbi.full$IE > 0 & ncbi.full$IC == 0 & ncbi.full$IF == 0)
n.IF <- sum(ncbi.full$IF > 0 & ncbi.full$IE == 0 & ncbi.full$IC == 0)
n.ICIE <- sum(ncbi.full$IC > 0 & ncbi.full$IE > 0 & ncbi.full$IF == 0)
n.ICIF <- sum(ncbi.full$IC > 0 & ncbi.full$IE == 0 & ncbi.full$IF > 0)
n.IEIF <- sum(ncbi.full$IC == 0 & ncbi.full$IE > 0 & ncbi.full$IF > 0)
n.ICIEIF <- sum(ncbi.full$IC > 0 & ncbi.full$IE > 0 & ncbi.full$IF > 0)
combinations <- n.ICIE + n.ICIF + n.IEIF + n.ICIEIF
none <- n - n.IC - n.IE - n.IF - n.ICIE - n.ICIF - n.IEIF - n.ICIEIF
proportions <- data.frame(
  group = c("I-C", "I-E", "I-F", "combinations", "No CRISPR"),
  value = c(n.IC, n.IE, n.IF, combinations, none)
)

ggplot(proportions, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="white", alpha=0.8) +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_brewer(palette = "Greens")

ncbi.full.IC <- subset(ncbi.full, IE==0 & IF==0 & IC > 0)
ncbi.full.IE <- subset(ncbi.full, IE >0 & IF==0 & IC == 0)
ncbi.full.IF <- subset(ncbi.full, IE==0 & IF > 0 & IC == 0)

# deactivated vs not
attach(ncbi.full.IC)
ncbi.full.IC$deactivated <- ifelse(grepl("AcrIC", Acr_genes_hom), 1, ifelse(gbaMC_IC > 0 | gbaHC_IC > 0, 1, 0))
IC.deactivated <- sum(ncbi.full.IC$deactivated > 0) #38/96 = 39%
# Other Acr found?                     
sum(ncbi.full.IC$acrIF > 0) # 14/96 = 14%
sum(ncbi.full.IC$acrIE > 0) # 15/96 = 15%
detach(ncbi.full.IC)

attach(ncbi.full.IE)
ncbi.full.IE$deactivated <- ifelse(grepl("AcrIE", Acr_genes_hom), 1, ifelse(gbaMC_IE > 0 | gbaHC_IE > 0, 1, 0))
IE.deactivated <- sum(ncbi.full.IE$deactivated > 0) # 208/352 = 59%
# Other Acr found?  
sum(ncbi.full.IE$acrIF > 0) # 190/352 = 53%
detach(ncbi.full.IE)

attach(ncbi.full.IF)
ncbi.full.IF$deactivated <- ifelse(grepl("AcrIF", Acr_genes_hom), 1, ifelse(gbaMC_IF > 0 | gbaHC_IF > 0, 1, 0))
IF.deactivated <- sum(ncbi.full.IF$deactivated > 0) # 700/1657 = 42%
# Other Acr found? 
sum(ncbi.full.IF$acrIE > 0) # 522/1657 = 31%
detach(ncbi.full.IF)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

p1 <- ggplot(ncbi.full.IC, aes(x=Size, color=as.factor(deactivated), fill=as.factor(deactivated)))+
  geom_density(alpha=0.4)+
  ylab("Density") + xlab("Genome size (Mbp)") +
  scale_fill_manual(values=c("#1b9e77","#d95f02")) +
  scale_colour_manual(values=c("#1b9e77","#d95f02")) +
  theme_classic() +
  xlim(6,7.5)

p2 <- ggplot(ncbi.full.IE, aes(x=Size, color=as.factor(deactivated), fill=as.factor(deactivated)))+
  geom_density(alpha=0.4)+
  ylab("Density") + xlab("Genome size (Mbp)") +
  scale_fill_manual(values=c("#1b9e77","#d95f02")) +
  scale_colour_manual(values=c("#1b9e77","#d95f02")) +
  theme_classic() +
  xlim(6,7.5)

p3 <- ggplot(ncbi.full.IF, aes(x=Size, color=as.factor(deactivated), fill=as.factor(deactivated)))+
  geom_density(alpha=0.4)+
  ylab("Density") + xlab("Genome size (Mbp)") +
  scale_fill_manual(values=c("#1b9e77","#d95f02")) +
  scale_colour_manual(values=c("#1b9e77","#d95f02")) +
  theme_classic()+
  xlim(6,7.5)

multiplot(p1, p2, p3, cols=1)

