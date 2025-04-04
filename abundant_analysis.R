library(ggstance)
library(dplyr)
library(phyloseq)
library(castor)
library(ggplot2)
library(ggtree)
library(stringr)
library(ggtreeExtra)
library(reshape2)
library(vegan)
library("relaimpo")
library("caper")
library(DescTools)

overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#import data
ps<-readRDS("~/Desktop/WSNS/phyloseq_clean.RData")
ps<- prune_samples(sample_names(ps) %in% rownames(sample_data(ps)[sample_data(ps)$type == "watercolumn",]), ps)
#otu_table(ps)<-otu_table(decostand(otu_table(ps), MARGIN = 2, "total")*100, taxa_are_rows = TRUE)
comm<-otu_table(ps)
ps<-prune_taxa(taxa =  (rowSums(otu_table(ps))/sum(otu_table(ps)))*100>1, ps)
sum(otu_table(ps))/sum(comm)
otu_table(ps)
new.m<-mapply("/", data.frame(otu_table(ps)), colSums(comm))*100
rownames(new.m)<-rownames(otu_table(ps))
otu_table(ps)<-otu_table(new.m, taxa_are_rows = TRUE)
comm<-decostand(comm,MARGIN = 2, method = "total")

ori<-readRDS("~/Desktop/WSNS/origin_asv.RData")
fapro<-read.table("~/Desktop/WSNS/WS_NS_FAPROTAX_func_annot.list", header = TRUE, row.names = 1)
fapro<-fapro[rownames(comm),]
cwm<-FD::functcomp(fapro, t(as.data.frame(comm, MARGIN = 2, method = 'total')), bin.num = TRUE, CWM.type = "all")

bray.func<-vegdist(cwm)
bray.all<-vegdist(t(as.data.frame(comm)))
bray.ab<-vegdist(t(as.data.frame(otu_table(ps))), "bray")
cor.test(bray.all, bray.ab, method = "spearman")
cor.test(bray.func, bray.ab, method = "spearman")

#View(merge(as.data.frame(refseq(ps)), tax_table(ps), by.x = "row.names", by.y = "asv"))
## wind information and merging with data
# import and format
wind<-read.csv('~/Desktop/WSNS/wind_wsns.csv', sep = ';')
wind$time<-as.POSIXct(paste(wind$Date, wind$Time), format="%d/%m/%Y %H:%M:%S",tz = "GMT")
wind$time.t<-as.numeric(wind$time)
# subset watercolumn data and format for merging with wind data
info2<-as.data.frame(as.matrix(sample_data(ps)))
info2$time.t<-as.numeric(as.POSIXct(info2$time))
for (i in 1:nrow(info2)){info2$time.wind[i]<-Closest(a = info2$time.t[i], x = wind$time.t)}
info2$time.wind<-as.POSIXct(info2$time.wind, format="%d/%m/%Y %H:%M:%S",tz = "GMT", origin = "01/01/1970 00:00:00")
# new info file with wind data
info2<-merge(info2, wind, by.x = "time.wind", by.y = "time")

#color palettes
pal<-readRDS( "~/Desktop/WSNS/pal_ori.RData")
names(pal)<-gsub("/Watercolumn","",names(pal))
pal.t<-readRDS( "~/Desktop/WSNS/pal_t.RData")

## 1/ tree + taxonomy
tree<-ps@phy_tree
df=merge(ps@tax_table, ori[, 14:15], by = "asv")
labs<-as.data.frame(t(data.frame(asv.1 = "Planktomarina temperata",asv.3 = "Synechococcus CC9902",asv.4 = "Flavobacteriaceae NS3 MG",asv.5 = "SAR11 Clade Ia",asv.6 = "SAR11 Clade Ia",asv.7 = "Actinomarina",asv.8 = "Roseibacillus",asv.10 = "Cryomorphaceae",asv.11 = "Aurantivirga",asv.12 = "Synechococcus CC9902",asv.14 = "Fluviicola",asv.15 = "Amylibacter",asv.16 = "Flavobacteriaceae NS2b MG",asv.19 = "Flavobacteriaceae NS5 MG",asv.20 = "Lentibacter algarum",asv.22 = "Roseibacillus",asv.24 = "Ascidiaceihabitans")))
df<-merge(df, labs,by.x ="asv" ,by.y = "row.names" )
df$asv<-factor(df$asv, levels = tree$tip.label, ordered = TRUE)
df$Origin<-gsub("/Watercolumn","",str_to_title(df$origin))
df$Origin<-factor(df$Origin, levels = names(pal), ordered = TRUE)
df$asv_names<-as.character(df$asv)

# plot
gt<-ggtree(tree, ladderize=F)+
  #xlim(0,4)+
  geom_tiplab(align = TRUE, aes(label = ""))+
  geom_fruit(data=df,geom=geom_point,mapping=aes(y=asv, color=Phylum),position="identity",size = 4) +
  geom_fruit(data=df,geom=geom_text,mapping=aes(y=asv, label= V1),offset = 0, hjust = 0,size = 4) +
  #geom_fruit(data=df,geom=geom_text,mapping=aes(y=asv, label= asv_names),offset = 0.5, hjust = 0,size = 4) +
  scale_color_manual(values=pal.t)+
  theme_void()+
  theme(legend.text = element_text(size = 12), 
        legend.position = "bottom", 
        legend.direction = "vertical", 
        legend.title = element_text(size = 12, face="bold"));gt

## 2/ Origin
gt<-gt+geom_fruit(data=df,geom=geom_tile,mapping=aes(y=asv, fill = Origin), offset = 0.55, width = 0.25) +
  scale_fill_manual(values=pal, name = "Origin")+
  theme(legend.text = element_text(size = 12), 
        legend.position = "bottom", 
        legend.direction = "vertical", 
        legend.title = element_text(size = 12, face="bold"));gt
  
## 3/ Abundance distribution
m.ps<-melt(otu_table(ps))
m.ps<-merge(m.ps, sample_data(ps), by.x = "Var2", by.y = "row.names")
m.ps$time<-as.POSIXct(m.ps$time, format="%d/%m/%Y %H:%M:%S",tz = "GMT", origin = "01/01/1970 00:00:00")
m.ps$location<-factor(m.ps$location,levels = c("north sea", "wadden sea"), labels =c("NS", "WS") , ordered = TRUE )
m.ps$season<-factor(m.ps$season,levels = c("winter", "spring", "summer","autumn"), labels =c("Winter", "Spring", "Summer", "Autumn") , ordered = TRUE )
m.ps$asv<-factor(m.ps$Var1, levels = rev(tree$tip.label), ordered = TRUE)
m.ps$Depth<-ifelse(m.ps$season == "Summer" & m.ps$location == "NS",str_to_title(m.ps$DEPTH), "" )
m.ps$Depth<-factor(m.ps$Depth , levels = c("Surface", 'Depth', ""), labels = c("Surface", "Bottom", "") , ordered = TRUE)
m.ps<-m.ps[!m.ps$DEPTH %in% c("depth") | !m.ps$Depth %in% c("") ,]
m.ps$id<-m.ps$asv
m.ps<-m.ps[order(m.ps$asv),]
m.ps<-merge(m.ps, labs, by.x = "asv", by.y = "row.names")

top<-m.ps %>% group_by(asv, location) %>% top_n(3,value)
top$asv<-factor(top$asv,levels = rev(tree$tip.label), ordered = TRUE)
top<-top[order(top$asv, decreasing = TRUE),c(1,4,33,7)]

gd<-ggplot(data = m.ps[m.ps$value>0,], 
       aes(x = time, y = asv, size = value, color = value))+
  facet_grid(.~location+season+Depth, scales = "free", space = "free")+
  scale_color_viridis_c(guide = "none")+labs(size='%')+
  geom_point()+
  theme_void()+
  theme(strip.text =  element_text(size = 12))+
  theme(legend.text = element_text(size = 12), 
        legend.position = "bottom", 
        legend.title = element_text(size = 12, face="bold"));gd

##
# Figure 2
##
library(patchwork)
gt+gd # + plot_layout(widths = c(2,1 ))




## 4/ ANOVAS
ano<-merge(t(otu_table(ps)), info2, by.x = "row.names", by.y = "sampleID")
ano$location<-factor(ano$location,levels = c("north sea", "wadden sea"), labels =c("NS", "WS") , ordered = TRUE )
ano$season<-factor(ano$season,levels = c("winter", "spring", "summer","autumn"), labels =c("Winter", "Spring", "Summer", "Autumn") , ordered = TRUE )
ano$Depth<-ifelse(ano$season == "Summer" & ano$location == "NS",str_to_title(ano$DEPTH), "" )
ano$Depth<-factor(ano$Depth , levels = c("Surface", 'Depth', ""), labels = c("Surface", "Bottom", "") , ordered = TRUE)
ano<-ano[!ano$DEPTH %in% c("depth") | !ano$Depth %in% c("") ,]
for (i in grep("time.t.x|light.ac|temp|sal|tide|Winddirection|Windspeed.avgph",colnames(ano)) ){ano[,i]<-as.numeric(ano[,i])}

citation("relaimpo")
list_anova<-list()
for (i in tree$tip.label){
      # linear model
      lm_fit<-lm(ano[[i]] ~ time.t.x+light.ac+temp+sal+tide+Winddirection+Windspeed.avgph , data=ano)
      # relative importance of each variable
      car_relimpo_metrics <- calc.relimp(lm_fit, type = c("car"))
      car_relimpo_metricsDF<-data.frame(cbind(car_relimpo_metrics$car,rownames(car_relimpo_metrics)))
      car_relimpo_metricsDF$R_name<-rownames(car_relimpo_metricsDF)
      colnames(car_relimpo_metricsDF)<-c("car","R_name")
      car_relimpo_metricsDF$asv<-i
      car_relimpo_metricsDF$pv.var<-coefficients(summary(lm_fit))[-1,4]
      car_relimpo_metricsDF$pv<-overall_p(lm_fit)
      car_relimpo_metricsDF$R2<-summary(lm_fit)$r.squared
      car_relimpo_metricsDF$adj.R2<-summary(lm_fit)$adj.r.squared
      list_anova[[i]]<-car_relimpo_metricsDF
    }
list_anova<-bind_rows(list_anova, .id = "asv")
list_anova$asv<-factor(list_anova$asv,levels = rev(tree$tip.label), ordered = TRUE) 
list_anova$R_name<-factor(list_anova$R_name, 
                          levels = (c("time.t.x","light.ac","temp", "sal", "tide", "Winddirection", "Windspeed.avgph")),
                          labels = (c("Time","Accumulated Daylight","Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed")), ordered = TRUE)
list_anova$NS<-ifelse(list_anova$pv.var>0.05, "Non-significant", as.character(list_anova$R_name))
list_anova$NS<-factor(list_anova$NS, levels = c("Non-significant","Time","Accumulated Daylight","Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed"), ordered = TRUE)
list_anova$lab<-"All\nt.s."

list_anova_ori<-list()
list_lag<-list()
for (i in tree$tip.label){
  # subset of origin
  orig=str_split_fixed(df[df$asv == i,"Origin"], "/", n = 3)
  ano.sub<-ano[ano$Sea %in% orig[1] & ano$Season %in% orig[2] & ano$Depth %in% orig[3],] 
  
  # linear model
  lm_fit<-lm(ano.sub[[i]] ~ time.t.x+light.ac+temp+sal+tide+Winddirection+Windspeed.avgph , data=ano.sub)
  # relative importance of each variable
  car_relimpo_metrics <- calc.relimp(lm_fit, type = c("car"))
  car_relimpo_metricsDF<-data.frame(cbind(car_relimpo_metrics$car,rownames(car_relimpo_metrics)))
  car_relimpo_metricsDF$R_name<-rownames(car_relimpo_metricsDF)
  colnames(car_relimpo_metricsDF)<-c("car","R_name")
  car_relimpo_metricsDF$pv.var<-coefficients(summary(lm_fit))[-1,4]
  car_relimpo_metricsDF$pv<-overall_p(lm_fit)
  car_relimpo_metricsDF$R2<-summary(lm_fit)$r.squared
  car_relimpo_metricsDF$adj.R2<-summary(lm_fit)$adj.r.squared
  list_anova_ori[[i]]<-car_relimpo_metricsDF
  
  res<-as.data.frame(matrix(ncol= 6, nrow = 7));colnames(res)<-c("var", "r2", "pv", "lag.r2", "lag", "lag.pv")
  res$var<-c("time.t.x","light.ac","temp","sal","tide","Winddirection", "Windspeed.avgph")
  for(y in res$var){
    res[res$var == y, "r2"]<-cor.test(ano.sub[[i]], ano.sub[,y])$estimate
    res[res$var == y, "pv"]<-cor.test(ano.sub[[i]], ano.sub[,y])$p.value
    corLag<-ccf(ano.sub[[i]],ano.sub[,y])
    res[res$var == y, "lag.r2"]<-max(corLag$acf)
    res[res$var == y, "lag"]<-corLag$lag[which.max(corLag$acf)]
    res[res$var == y, "lag.pv"]<- 2*(1 - pnorm(max(abs(corLag$acf)), mean = 0, sd = 1/sqrt(corLag$n.used)))
  }
  list_lag[[i]]<-res
}

list_anova_ori<-bind_rows(list_anova_ori, .id = "asv")
list_anova_ori$asv<-factor(list_anova_ori$asv,levels = rev(tree$tip.label), ordered = TRUE) 
list_anova_ori$R_name<-factor(list_anova_ori$R_name, 
                          levels = (c("time.t.x","light.ac","temp", "sal", "tide", "Winddirection", "Windspeed.avgph")),
                          labels = (c("Time","Accumulated Daylight","Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed")), ordered = TRUE)
list_anova_ori$NS<-ifelse(list_anova_ori$pv.var>0.1, "Non-significant", as.character(list_anova$R_name))
list_anova_ori$NS<-factor(list_anova_ori$NS, levels = c("Non-significant","Time","Accumulated Daylight","Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed"), ordered = TRUE)
list_anova_ori$lab<-"Origin\nt.s."

pal.env<-c("gray75","black","gold3","#f0624f", "#913e33", "#224282","#71bf9b","#46705d")
names(pal.env)<-c("Non-significant","Time","Accumulated Daylight","Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed")

gano<-ggplot(data=list_anova,aes(y=asv, x = car,  fill = NS ))+
  facet_grid(.~lab)+
  geom_barh(stat = "identity", position = "stack")+
  scale_x_continuous(limits = c(0,1.1), breaks = c(0,0.5,1),labels = c("0", "0.5", "1"))+
  scale_fill_manual(values=pal.env, name = "Env. variable", guide = "none")+
  theme_void()+
  theme(strip.text =  element_text(size = 14))+
  theme(axis.text.x  = element_text(size = 12), 
        legend.text = element_text(size = 14), 
        legend.position = "bottom", 
        legend.direction  = "vertical", 
        legend.title = element_text(size = 14, face="bold"));gano

#View(list_anova_ori[list_anova_ori$NS != "Non-significant",])
check<-merge(list_anova_ori, tax.func, by.x =  "asv", by.y = "Row.names")
gano.ori<-ggplot(data=list_anova_ori,aes(y=asv, x = car,  fill = NS ))+
  facet_grid(.~lab)+
  geom_barh(stat = "identity", position = "stack")+
  scale_x_continuous(limits = c(0,1.1), breaks = c(0,0.5,1),labels = c("0", "0.5", "1"))+
  scale_fill_manual(values=pal.env, name = "Env. variable")+
  theme_void()+
  theme(strip.text =  element_text(size = 12))+
  theme(axis.text.x  = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.position = "bottom", 
        legend.direction  = "vertical", 
        legend.title = element_text(size = 12, face="bold"));gano.ori

ggano<-gano+gano.ori+ plot_layout(widths = c(1,1), guides = "collect")
meta<-gd+ggano+plot_layout(widths = c(5,2))
gt+meta+plot_layout(nrow = 1)

gt+gd

list_lag<-bind_rows(list_lag, .id = "asv")
list_lag$asv<-factor(list_anova_ori$asv,levels = rev(tree$tip.label), ordered = TRUE) 
list_lag$var<-factor(list_lag$var, 
                       levels = (c("time.t.x","light.ac","temp", "sal", "tide", "Winddirection", "Windspeed.avgph")),
                       labels = (c("Time","Accumulated Daylight","Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed")), ordered = TRUE)
list_lag$NS.lag<-ifelse(list_lag$lag.pv>0.05, "Non-significant", as.character(list_lag$var))
list_lag$NS.lag<-ifelse(abs(list_lag$lag)>2, "Non-significant", list_lag$NS.lag)
list_lag$NS.lag<-factor(list_lag$NS.lag, levels = c("Non-significant","Time","Accumulated Daylight","Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed"), ordered = TRUE)

pal.env<-c("gray75","black","gold3","#f0624f", "#913e33", "#224282","#71bf9b","#46705d")
names(pal.env)<-c("Non-significant","Time","Accumulated Daylight","Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed")

ggplot(data=list_lag,aes(y=asv, x = lag.r2,  fill = NS.lag ))+
  geom_bar(stat = "identity", position = "stack")+
  #scale_x_continuous(limits = c(0,1.1), breaks = c(0,0.5,1),labels = c("0", "0.5", "1"))+
  scale_fill_manual(values=pal.env, name = "Env. variable")+
  theme(strip.text =  element_text(size = 14))+
  theme(axis.text.y = element_text(size = 16, angle = 0, hjust = 1, vjust = 0.5))+
  theme(axis.text.x  = element_text(size = 12), 
        legend.text = element_text(size = 14), 
        legend.position = "bottom", 
        legend.direction  = "vertical", 
        legend.title = element_text(size = 14, face="bold"))



