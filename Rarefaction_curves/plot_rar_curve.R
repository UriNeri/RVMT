library(ggplot2)
library(dplyr)
library(extrafont)
nice_theme<-theme(text=element_text(color="black",size=8,family="Source Sans Pro"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans Pro"))

## Overall curve
data_rc_2 <- read.delim("Rar_by_sample-from-mapping_trace.tsv", stringsAsFactors = T)
summary(data_rc_2)
ggplot(data_rc_2) + geom_ribbon(aes(x=N_sample,ymin=Min,ymax=Max),alpha=0.3) + geom_line(aes(x=N_sample,y=N_ANI),col="red") + xlab("N Samples") + ylab("N RvANI90") + nice_theme + scale_x_continuous(labels=scales::comma) + scale_y_continuous(labels=scales::comma) + nice_theme

## By ecosystem
data_rc_3 <- read.delim("Rar_by_sample-from-mapping-spliteco_trace.tsv", stringsAsFactors = T)
summary(data_rc_3)
data_rc_3$Eco<-factor(data_rc_3$Eco,ordered=T,levels=rev(c("Engineered","Terrestrial;Soil","Terrestrial;Plant litter and peat","Host-associated;Plant Phyllosphere","Host-associated;Plant Rhizosphere","Host-associated;Animal microbiome","Aquatic;Freshwater","Aquatic;Marine","Aquatic;Other")))
summary(data_rc_3)
custom_scale <- rev(c("#9e9e9eff","#b4662eff","#db9763ff","#00c997ff","#006149ff","#6cf0cfff","#4f97dcff","#073969ff","#afd4f8ff"))
#                     Engineered   Soil    Plant litter    Phyllosphere Rhizosphere    Animal   Freshwater     Marine    Other
ggplot(data_rc_3) + geom_ribbon(aes(x=N_sample,ymin=Min,ymax=Max,fill=Eco),alpha=0.3) + geom_line(aes(x=N_sample,y=N_ANI,col=Eco)) + xlab("N Samples") + ylab("N RvANI90") + nice_theme + scale_x_continuous(labels=scales::comma) + scale_y_continuous(labels=scales::comma) + scale_color_manual(values=custom_scale) +  scale_fill_manual(values=custom_scale) + nice_theme
