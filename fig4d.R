####R code for analysis and visulize the simulation and validation result

library(forcats)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyverse)
library(hrbrthemes)
library(psych)
library(plyr)
# library(extrafont)
# font_import()
# fonts()
# loadfonts()


library(forcats)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyverse)
library(hrbrthemes)
library(psych)
library(pheatmap)
library(psych)
library(igraph)
library(ggm)
library(ggpubr)

library(ggsci)

summary_func <- function(x, col){
         c(mean = mean(x[[col]], na.rm=TRUE), sd = sd(x[[col]], na.rm=TRUE)) 
    } 
data_summary <- function(data, varname, groupnames){ 
    data_sum<-ddply(data, groupnames, .fun=summary_func, varname)
    data_sum <- rename(data_sum, c("mean" = varname)) 
    return(data_sum) 
    }


##########################
##########################
##########################
##########################


data<- read.table("kegg_r.csv", header = TRUE, sep = ",",colClasses=c("character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

mae_group<-c(
rep('chondroitin sulfate degradation I (bacterial)',dim(data)[1]),
rep('sucrose degradation III (sucrose invertase)',dim(data)[1]),
rep('superpathway of glucose and xylose degradation',dim(data)[1]),
rep('L-rhamnose degradation I',dim(data)[1]),
rep('lipid IVA biosynthesis',dim(data)[1]),
rep('queuosine biosynthesis',dim(data)[1]),
rep('preQ0 biosynthesis',dim(data)[1]),
rep('TCA cycle VI (obligate autotrophs)',dim(data)[1]),
rep('D-fructuronate degradation',dim(data)[1]),
rep('4-deoxy-L-threo-hex-4-enopyranuronate degradation',dim(data)[1])
)
mae_value<-c(data$X1,data$X2,data$X3,data$X4,data$X5,data$X6,data$X7,data$X8,data$X9,data$X10)

mae_con<-c(rep(data$g,10))
mae_data<-data.frame(mae_group,mae_value,mae_con)
mae_data$mae_con <-factor(mae_data$mae_con,,levels=c('WT','KO'))
mae_data$mae_group<-factor(mae_data$mae_group,, levels = rev(c('chondroitin sulfate degradation I (bacterial)', 'sucrose degradation III (sucrose invertase)', 'superpathway of glucose and xylose degradation', 'L-rhamnose degradation I', 'lipid IVA biosynthesis', 'queuosine biosynthesis', 'preQ0 biosynthesis', 'TCA cycle VI (obligate autotrophs)', 'D-fructuronate degradation', '4-deoxy-L-threo-hex-4-enopyranuronate degradation')))

p3<-ggplot(mae_data, aes(x=mae_group, y=mae_value, fill=mae_con)) + 
geom_boxplot(alpha=0.5)+labs(x="Pathway", y ="Proportion") +scale_color_nejm()+theme_ipsum()+coord_flip()


ggsave(p3,device=cairo_pdf, filename='top10pathwaykegg.pdf', width=25, height=25, units=c("cm"))

