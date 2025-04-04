################################################################################################################
# R script for studying cascabels outputs
################################################################################################################

###############################
# Packages and useful objects
###############################

# Load required packages
library(assertthat)
library(biom)
library(Biostrings)
library(caper)
library(castor)
library(data.table)
library(DECIPHER)
library(DescTools)
library(ggplot2)
library(ggpubr)
library(ggnewscale)
library(ggrepel)
library(ggraph)
library(ggtree)
library(ggtreeExtra)
library(graph)
library(iNEXT)
library(igraph)
library(mapdata)
library(marmap)
library(microseq)
library(phyloseq)
library(phangorn)
library(picante)
library(pspearman)
library(patchwork)
library(reshape2)
library(readxl)
library(scales)
library(ShortRead)
library(stringdist)
library(stringi)
library(stringr)
library(tidyr)
library(tidygraph)
library(tidyverse)
library(treeio)
library(usedist)
library(vegan)

# vector with all the colors available in R, except gray
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

# location of the script
rstudioapi::getSourceEditorContext()$path
options(getClass.msg=FALSE)

###############################
# Data import and format
###############################

# import data filtered by Julia
asv<-as.data.frame(readRDS("~/Desktop/WSNS/asv3_minSkin_20221115.rds"))

# make phyloseq object
seq<-readDNAStringSet("~/Desktop/WSNS/representative_seq_set.fasta", format = "fasta")
seq<-seq[names(seq) %in% rownames(asv)]

# import tree build with MAFT and Fasttree
tree<-read.tree(file = "~/Desktop/WSNS/representative_seq_set.tree")
tree<-get_subtree_with_tips(tree, only_tips=rownames(asv))$subtree

# extract species-taxonomy matrix
taxo<-as.data.frame(observation_metadata(read_biom("~/Desktop/WSNS/asvTable.biom")))
colnames(taxo)<-c("Domain", "Phylum","Class", "Order", "Family", "Genus", "Species")
taxo[taxo == "NA"]<-NA 
taxo[is.na(taxo$Domain),"Domain"]<-"Unclassified"
taxo$level<-colnames(taxo)[apply(taxo, 1, function(x) length(which(!is.na(x))))] # get taxonomic level at which each sequence is annotated
taxo[grep("Proteobacteria", taxo[,2]),2]<-taxo[grep("Proteobacteria", taxo[,2]),3]
taxo<-taxo[rownames(taxo) %in% rownames(asv),]
taxo$asv<-rownames(taxo)

# metadata file
info<-read.csv('~/Desktop/WSNS/metadata_ws_ns.csv', sep = ';')
info<-distinct(info)
colnames(asv) %in% info$sampleID
info$label<-ifelse(info$type == "sediment", "NS\nSed",
                   ifelse(info$type == "sediments", "WS\nSed",
                          ifelse(info$type %in% c("reference" ,"watercolumn") & info$location == "north sea", "NS\nWat",
                                 ifelse(info$type %in% c("reference" ,"watercolumn") & info$location == "wadden sea", "WS\nWat","Lab"))))
info$label<-factor(info$label, levels= c("NS\nWat","NS\nSed", "WS\nWat", "WS\nSed","Lab"), ordered = TRUE)
info[info$type == "sediments","type"]<-"sediment"
info$time<-as.POSIXct(paste(info$date, info$utc.time), format="%d/%m/%Y %H:%M:%S", tz = "CET")
info$hour<-factor(format(info$time, "%H:%M") )
info$DEPTH<-ifelse(grepl("1", info$depth), "surface", "depth")
info$DEPTH<-ifelse(info$type %in% c("sediment", "control"),info$type, info$DEPTH)
info<-info[!info$sampleID %in% c("NIOZ153.027.028","NIOZ309.067.060"),] # samples with weird signal
info$season<-factor(info$season, levels = c("winter", "spring", "summer", "autumn", "lab"), ordered = TRUE)
info$location<-factor(info$location, levels = c("north sea", "wadden sea", "lab"), ordered = TRUE)
info$date<-as.Date(info$date, format = "%d/%m/%Y")

## sunlight info
# calculate sunrise and sunset according to date and coordinates
library("suncalc")
cbind(info$time,getSunlightTimes(data = info, tz = "CET", keep = c("sunrise", "sunset")))
dates<-cbind(info$time,getSunlightTimes(data = info, tz = "CET", keep = c("sunrise", "sunset")))
dates<-na.omit(dates)
colnames(dates)[1]<-"time"
dates<-dates[!duplicated(dates),];rownames(dates)<-NULL

# compute light accumulation (day: +Xh since sunrise/-Yh since sunset)
for (i in 1:nrow(dates)){
  if(dates$time[i] > dates$sunrise[i] & dates$time[i] < dates$sunset[i]){dates$light.ac[i]<-difftime(dates$time[i],dates$sunrise[i], units = "hour" ) }
  if(dates$time[i] < dates$sunrise[i]){dates$light.ac[i]<-difftime(as.POSIXct(dates$sunset[i]-as.difftime(1, unit="days")),dates$time[i], unit = "hour") }
  if(dates$time[i] > dates$sunset[i]){dates$light.ac[i]<- -difftime(dates$time[i],dates$sunset[i], unit = "hour") }
}
# format and merge with info
info<-merge(info, dates[, c("light.ac", "time")], by = "time", all.x = TRUE)

# tide information
tide<-read.csv('~/Desktop/WSNS/tides.csv', sep = ';')
colnames(tide)[6]<-"tide"
tide$time<-as.POSIXct(paste(tide$date, tide$UTC_time), format="%d/%m/%Y %H:%M:%S", tz = "CET")
info<-merge(info, tide[, c("tide", "time")], by = "time", all.x = TRUE)
rownames(info)<-info$sampleID

## wind information and merging with data
# import and format
wind<-read.csv('~/Desktop/WSNS/wind_wsns.csv', sep = ';')
wind$time<-as.POSIXct(paste(wind$Date, wind$Time), format="%d/%m/%Y %H:%M:%S",tz = "GMT")
wind$time.t<-as.numeric(wind$time)
# subset watercolumn data and format for merging with wind data
info2<-info[info$type=="watercolumn",]
info2$time.t<-as.numeric(info2$time)
for (i in 1:nrow(info2)){info2$time.wind[i]<-Closest(a = info2$time.t[i], x = wind$time.t)}
info2$time.wind<-as.POSIXct(info2$time.wind, format="%d/%m/%Y %H:%M:%S",tz = "GMT", origin = "01/01/1970 00:00:00")
# new info file with wind data
info2<-merge(info2, wind, by.x = "time.wind", by.y = "time")
info.w<-info2[, c("sampleID", "Winddirection", "Windspeed.avgph")]

## outputs of FAPROTAX
fapro<-read.table("~/Desktop/WSNS/WS_NS_FAPROTAX_func_annot.list", header = TRUE, row.names = 1)
tax.func<-merge(taxo, fapro, by = "row.names")
ab<-data.frame(sed = rowSums(asv[,colnames(asv) %in% info[info$type == "sediment","sampleID"]]),
           sw = rowSums(asv[, colnames(asv) %in% info[info$type == "watercolumn","sampleID"]]))
ab.tax.func<-merge(ab, tax.func, by.x = "row.names", by.y = "Row.names")

# make phyloseq
ps<-phyloseq(otu_table(as.matrix(asv), taxa_are_rows = TRUE),
             sample_data(info),
             tax_table(as.matrix(taxo)),
             refseq(seq),
             phy_tree(tree))
saveRDS(ps, "~/Desktop/WSNS/phyloseq_clean.RData")

###############################
# Env. Analysis
###############################

# Map of locations
geo<-reshape2::melt(table(info[,c("lon", "lat","label")]))
geo<-geo[geo$value > 0,]
geo$Sampling<-gsub("\n", " ", geo$label)
geo[c(2,5),"Sampling"]<-c("NS Ref", "WS Ref")
geo$Sampling<-factor(geo$Sampling, levels= c("NS Wat", "NS Ref", "NS Sed" ,"WS Wat", "WS Ref", "WS Sed"), ordered = TRUE)
nd.map <- readRDS("~/Desktop/WSNS/gadm36_NLD_2_sp.rds")
nd.map.f <- subset(nd.map, !nd.map$NAME_1  %in% c("Zeeuwse meren", "IJsselmeer"))
nd.map.f <- fortify(nd.map.f)
country_shapes <- geom_polygon(data = nd.map.f, aes(x = long, y = lat, group = group),fill = "#CECECE", color = "#CECECE", size = 0.15)

# formating the labels for geographic coordinates
ewbrks <- seq(3,6,1)
nsbrks <- seq(52,56,1)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(abs(x), "°W"), ifelse(x > 0, paste(abs(x), "°E"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(x, "°S"), ifelse(x > 0, paste(x, "°N"),x))))

#plot map
ggplot()+ country_shapes +
  geom_jitter(data = geo, aes(x = lon, y = lat, color = Sampling, size = value), alpha = 0.90)+
  scale_size(name = "# samples",range = c(5,50), breaks = c(5,50,100), limits = c(0,100))+
  scale_color_manual(values= c("#ebe417", "#99952e", "#4d4b23", "#f22929","#ad2f2f","#4d2323" ))+
  coord_fixed(xlim = c(3.25,6.5), ylim = c(52.5,56.5))+
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0))+
  labs(title = "Sampling Strategy (2018/19)", y = '', x = '')+
  annotate("text",x=5, y=55.5, label = 'bolditalic("North \n     Sea")',color = "grey15", parse = TRUE, size = 7)+
  annotate("text",x=5.85, y=53.2, label = 'bolditalic("Wadden \n     Sea")',color = "grey15", parse = TRUE, size = 7)+
  theme(panel.background = element_rect(fill = "#596673"),
        panel.grid = element_blank(),
        legend.key  = element_blank(),
        legend.text = element_text(size = 16, color = "gray5"),
        legend.title = element_text(size = 18, color = "gray5", face = "bold"),
        axis.text = element_text(size = 16, color = "gray5"),
        title = element_text(size = 18, color = "gray5", face = "bold"))
  #guides(size = FALSE)
  
# seasonality summary
year<-melt(table(format(info[info$location != "lab","time"],"%b-%Y") )); colnames(year)<-c("month","value")
year<-as.data.frame(cbind(str_split_fixed(year$month,"-", 2),year$value)); colnames(year)<-c("month", "year", "value")
year<-rbind(year,data.frame(month = substr(month.name,1,3)[!substr(month.name,1,3) %in% year$month], year = 2018, value = 0))
year$month<-factor(year$month, levels = substr(month.name,1,3) )
year$Year<-factor(year$year,levels = c(2019,2018), ordered = TRUE ); year$value<-as.numeric(year$value); 

# plot number of samples per month
aggregate(value~Year,year, sum)
gy<-ggplot(year ,aes(x = month, y = value, fill = Year) )+
  geom_bar(stat = "identity", alpha = 0.85)+
  scale_fill_manual(values = c("#03396c", "seagreen"))+
  labs(title = "Temporal Span", y = '# of samples', x = '')+
  scale_y_continuous(expand = c(0,0), limits = c(0,30))+
  theme(panel.background = element_rect(fill = "grey88"),
        panel.grid.major.y = element_line(color = "grey78"),
        panel.grid.minor.y = element_line(color = "grey78"),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.key  = element_blank(),
        legend.text = element_text(size = 16, color = "gray5"),
        legend.title = element_text(size = 18, color = "gray5", face = "bold"),
        axis.text = element_text(size = 16, color = "gray5"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1) ,
        title = element_text(size = 16, color = "gray5", face = "bold"));gy

## Envroinmental summary
box.env<-info2[,c("location", "season","DEPTH", "temp", "sal","ox", "light.ac", "tide", "Winddirection", "Windspeed.avgph")]
#box.env[, 4:ncol(box.env)]<-scale(box.env[, 4:ncol(box.env)])
box.env<-melt(box.env)
box.env$variable<-factor(box.env$variable, 
       levels = (c("light.ac","temp", "sal","ox", "tide", "Winddirection", "Windspeed.avgph")),
       labels = (c("Accumulated Daylight (h)","Temperature (ºC)", "Salinity (PSU)","Oxygen (µ mol/kg) ", "Tides (m)", "Wind Direction (0-360º)", "Wind Speed (m/s)")), ordered = TRUE)
box.env$location<-factor(box.env$location,levels = c("north sea", "wadden sea"), labels =c("NS", "WS") , ordered = TRUE )
box.env$season<-factor(box.env$season,levels = c("winter", "spring", "summer","autumn"), labels =c("Winter", "Spring", "Summer", "Autumn") , ordered = TRUE )
box.env$Depth<-factor(box.env$DEPTH , levels = c("surface", 'depth', ""), labels = c("Surface", "Bottom", "") , ordered = TRUE)

aggregate(sal~location+season, data = info2, sd)

##
# Figure S2
##
ggplot(box.env, aes(y = value, x =Depth , fill = variable ))+
  facet_grid(variable~location+season, scales = "free")+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c("gold3","#f0624f", "#913e33", "purple","#224282","#71bf9b","#46705d"), guide = "none")+
  theme(axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text =  element_text(size = 12),
        strip.text.y =  element_text(size = 12, angle = 0, hjust = 0),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face="bold"))

##############################################################
# Rarefaction curves and diversity extrapolation
##############################################################

## Diversity estimates with iNEXT
all<-as.data.frame(asv)
all<-all[,colnames(all) %in% info[info$location != "lab","sampleID"]]

lsf<-data.frame(All = rowSums(all), 
                ns_sum=rowSums(all[, colnames(all)%in%info[info$season =="summer" & info$location == "north sea", "sampleID"]]),
                ns_aut=rowSums(all[, colnames(all)%in%info[info$season =="autumn" & info$location == "north sea", "sampleID"]]),
                ws_win=rowSums(all[, colnames(all)%in%info[info$season =="winter" & info$location == "wadden sea", "sampleID"]]),
                ws_spr=rowSums(all[, colnames(all)%in%info[info$season =="spring" & info$location == "wadden sea", "sampleID"]]),
                ws_sum=rowSums(all[, colnames(all)%in%info[info$season =="summer" & info$location == "wadden sea", "sampleID"]]),
                ws_aut=rowSums(all[, colnames(all)%in%info[info$season =="autumn" & info$location == "wadden sea", "sampleID"]]),
                sed_ns=rowSums(all[, colnames(all)%in%info[info$type =="sediment" & info$location == "north sea", "sampleID"]]),
                sed_ws=rowSums(all[,colnames(all)%in%info[info$type =="sediment" & info$location == "wadden sea", "sampleID"]]),
                sea_ns=rowSums(all[,colnames(all)%in%info[info$type =="watercolumn" & info$location == "north sea", "sampleID"]]),
                sea_ws=rowSums(all[,colnames(all)%in%info[info$type =="watercolumn" & info$location == "wadden sea", "sampleID"]]) )

#div_esti<-iNEXT(x = lsf, datatype="abundance", endpoint = sum(lsf[,"All"])*3)
#saveRDS(div_esti,"/export/lv4/projects/NIOZ153_309/CASCABEL/NIOZ153_309_WS_and_NS/Analysis/iNEXT.RData")
div_esti<-readRDS("~/Desktop/WSNS/iNEXT.RData")

# study output
div_esti$AsyEst
div_esti$iNextEst

# data for plt
df <- bind_rows(div_esti$iNextEst, .id = "site")
df$Assemblage<-factor(df$Assemblage)
df$type<-ifelse(grepl("sea|sed",df$Assemblage), "eco",ifelse(grepl("aut|win|spr|sum",df$Assemblage), "seas", as.character(df$Assemblage)) )
curves_labels<-aggregate( cbind(m,qD) ~ Assemblage+type,  data = df, FUN = max)

diff<-dcast( Assemblage~Method ,data = aggregate(qD~Method+Assemblage,df,max))
diff$coverage<-((diff$Extrapolation-(diff$Extrapolation-diff$Rarefaction))/diff$Extrapolation)*100

# plot
head(df)
grac_seas<-ggplot(df[df$type != "eco",], aes(x=m, y=qD, colour=Assemblage)) + 
  geom_line(mapping = aes(linetype = Method),size = 1)+
  geom_point(data = df[df$type != "eco" & df$Method == "Observed",], size = 2.5)+
  geom_label_repel(data = curves_labels[curves_labels$type != "eco",], aes(label = Assemblage, x = Inf, y = qD, color = Assemblage),alpha = 0.60, size = 4, fontface = "bold")+
  scale_linetype_manual(values = c("dashed","blank", "solid"), guide = FALSE)+
  #rev(c("#1d3c7a", "#727a1d", "#1d7a4d", "#7a421d"))
  scale_color_manual(guide ="none" , values = rev(c("#1d3c7a", "#727a1d", "#1d7a4d", "#7a421d", "#727a1d", "#7a421d", "coral3")))+
  scale_y_continuous(limits = c(0,20000), expand = c(0,0) )+
  scale_x_continuous(limits = c(0,2.1e7), expand = c(0,0) )+
  xlab("# of reads")+ylab("# of OTUs")+
  ggtitle("Rarefaction curves and diversity extrapolates", subtitle = "Across seasons")+xlab("# of reads")+ylab("# of OTUs")+
  #theme(aspect.ratio = 1)+
  theme(title = element_text(size = 18),
        subtitle = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, hjust = 0.5),
        legend.position="none",
        legend.text=element_text(size=16), 
        legend.title = element_text(size=18, face = "bold"),
        axis.line = element_line(size = 1),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = "transparent"), # get rid of legend panel bg
        axis.ticks = element_line(size = 1));grac_seas

grac_eco<-ggplot(df[df$type != "seas",], aes(x=m, y=qD, colour=Assemblage)) + 
  geom_line(aes(linetype = Method), size = 1)+
  geom_point(data = df[df$type != "seas" & df$Method == "Observed",], size = 2.5)+
  geom_label(data = curves_labels[curves_labels$type != "seas",], aes(label = Assemblage, x = Inf, y = qD, color = Assemblage),alpha = 0.60, hjust = 1, size = 4, fontface = "bold")+
  scale_linetype_manual(values = c("dashed", "blank", "solid"), guide = FALSE)+
  #c("#ebe417", "#99952e", "#4d4b23", "#f22929","#ad2f2f","#4d2323" ))+
  scale_color_manual(guide =FALSE , values = rev(c("#4d2323" ,"#4d4b23" ,"#f22929", "#ebe417" , "coral3")))+
  scale_y_continuous(limits = c(0,20000), expand = c(0,0) )+
  scale_x_continuous(limits = c(0,2.1e7), expand = c(0,0) )+
  ggtitle("", subtitle = "Across substrates")+xlab("# of reads")+ylab("# of OTUs")+
  #theme(aspect.ratio = 1)+
  theme(title = element_text(size = 18),
        subtitle = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, hjust = 0.5),
        legend.position="none",
        legend.text=element_text(size=16), 
        legend.title = element_text(size=18, face = "bold"),
        axis.line = element_line(size = 1),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = "transparent"), # get rid of legend panel bg
        axis.ticks = element_line(size = 1));grac_eco

##
# Figure S1
##
(grac_seas/grac_eco)

##############################################################
# Taxonomic plots
##############################################################

## level of taxonomic annotation across samples
glevel<-aggregate(asv, by = list(taxo[,"level"]), FUN = sum)
glevel<-reshape2::melt(glevel);colnames(glevel)<-c("Level", "sample","n_reads")
glevel$Level<-factor(glevel$Level, levels = colnames(taxo)[1:7], ordered = TRUE)
glevel<-merge(glevel, info, by.x = 'sample', by.y = "sampleID")

#plot
gl<-ggplot(glevel, aes(x = sample, y= n_reads, fill = Level))+
  facet_grid(.~label, scale = "free_x", space ="free_x")+
  geom_bar(stat = "identity", position = "fill")+
  scale_y_continuous(expand =c(0,0))+
  xlab("Samples (179)")+ylab("Proportions")+
  scale_fill_viridis_d(direction = -1, drop = FALSE)+
  ggtitle(paste("16S Taxonomic Annotation Level", sep = ""))+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 24,hjust = 0.5))+
  theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5))+
  #theme(axis.text.x = element_text(size = 16, angle = 30, hjust = 1, vjust = 1))+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size = 16))+
  theme(strip.text =  element_text(size = 16, vjust = 0))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.background = element_rect(fill = "transparent",colour = NA))+
  #theme(aspect.ratio = 1)+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gl

## taxonomic annotation across samples at Phylum level
gtaxa<-aggregate(asv, by = list(taxo[,"Phylum"]), FUN = sum)
gtaxa<-reshape2::melt(gtaxa);colnames(gtaxa)<-c("Phylum", "sample","n_reads")
gtaxa$per<-round((gtaxa$n_reads/sum(gtaxa$n_reads))*100,2)
gts<-aggregate(per~ Phylum, data =gtaxa, FUN = sum); gts[order(gts$per),]
gtaxa$taxa<- ifelse(gtaxa$Phylum %in% gts[gts$per >0.005,"Phylum"] ,as.character(gtaxa$Phylum), "Other" )
gtaxa$taxa<-factor(gtaxa$taxa, levels = c(sort(unique(gtaxa$taxa))[-match(x = c("Other", "Unclassified"),sort(unique(gtaxa$taxa)), nomatch = 1000)], "Unclassified", "Other"),ordered = TRUE )
for (i in 1:nrow(gtaxa)){ gtaxa$perc[i]<-(gtaxa$n_reads[i]/sum(gtaxa[gtaxa$sample %in% gtaxa$sample[i], 'n_reads' ]))*100   }
gtaxa<-merge(gtaxa, info, by.x = 'sample', by.y = "sampleID")

# Organize taxonomy by supergroup 
unique(gtaxa$taxa)
labtax<-merge(data.frame(taxa=unique(gtaxa$taxa)), unique(taxo[,1:2]), by.x = "taxa", by.y = "Phylum", all.x = TRUE);rownames(labtax)<-NULL
labtax[labtax$taxa=="Other", "Domain"]<-"Other"
labtax<-labtax[order( factor(labtax$Domain, levels = c("Bacteria", "Archaea", "Other"), ordered = TRUE) ,labtax$taxa),];rownames(labtax)<-NULL
labtax<-unique(labtax);rownames(labtax)=NULL
labtax<-rbind(labtax[-c(match(c("Alphaproteobacteria", "Gammaproteobacteria","Spirochaetota", "Sva0485","Verrucomicrobiota"), labtax$taxa ),  24:28), ],
              labtax[c(match(c("Alphaproteobacteria", "Gammaproteobacteria","Spirochaetota",  "Sva0485","Verrucomicrobiota"), labtax$taxa ),  24:28), ] )

# attribute a color palette (gradient) for each Phylum
labtax$colors<-c(colorRampPalette(c("darkorange4", "darkorange4"))(nrow(labtax[labtax$taxa == "Acidobacteriota",])),
                 colorRampPalette(c("orange", "orange"))(nrow(labtax[labtax$taxa == "Actinobacteriota",])),
                 colorRampPalette(c("pink", "pink"))(nrow(labtax[labtax$taxa == "Bacteroidota",])),
                 colorRampPalette(c("purple", "purple"))(nrow(labtax[labtax$taxa == "Bdellovibrionota",])),
                 colorRampPalette(c("palegreen", "palegreen1"))(nrow(labtax[labtax$taxa == "Calditrichota",])),
                 colorRampPalette(c("palegreen", "palegreen1"))(nrow(labtax[labtax$taxa == "Campylobacterota",])),
                 colorRampPalette(c("#c93e71", "#c93e71"))(nrow(labtax[labtax$taxa == "Chlamydiota",])),
                 colorRampPalette(c("palegreen3", "palegreen3"))(nrow(labtax[labtax$taxa == "Chloroflexi",])),
                 colorRampPalette(c("seagreen", "seagreen"))(nrow(labtax[labtax$taxa == "Cyanobacteria",])),
                 colorRampPalette(c("#940606", "#940606"))(nrow(labtax[labtax$taxa == "Desulfobacterota",])),
                 colorRampPalette(c("slateblue4", "slateblue3"))(nrow(labtax[labtax$taxa == "Firmicutes",])),
                 colorRampPalette(c("pink2", "pink2"))(nrow(labtax[labtax$taxa == "Fusobacteriota",])),
                 colorRampPalette(c("#00B09B", "#96C93D"))(nrow(labtax[labtax$taxa == "Gemmatimonadota",])),
                 colorRampPalette(c("#A1FFCE", "#FAFFD1"))(nrow(labtax[labtax$taxa == "Latescibacterota",])),
                 colorRampPalette(c("royalblue", "royalblue"))(nrow(labtax[labtax$taxa == "Marinimicrobia (SAR406 clade)",])),
                 colorRampPalette(c("yellow", "yellow"))(nrow(labtax[labtax$taxa == "Myxococcota",])),
                 colorRampPalette(c("seagreen2", "seagreen2"))(nrow(labtax[labtax$taxa == "NB1-j",])),
                 colorRampPalette(c("#391a61" , "#391a61"))(nrow(labtax[labtax$taxa == "Nitrospirota",])),
                 colorRampPalette(c("gray75" , "gray75"))(nrow(labtax[labtax$taxa == "Patescibacteria",])),
                 colorRampPalette(c("#61591a", "#61591a"))(nrow(labtax[labtax$taxa == "Planctomycetota",])),
                 colorRampPalette(c("#ffd89b", "#19547b"))(nrow(labtax[labtax$taxa == "Poribacteria",])),
                 colorRampPalette(c("#000046","#000046"))(nrow(labtax[labtax$taxa == "Alphaproteobacteria",])),
                 colorRampPalette(c("#28446b","#28446b"))(nrow(labtax[labtax$taxa == "Gammaproteobacteria",])),
                 colorRampPalette(c("coral","coral"))(nrow(labtax[labtax$taxa == "Spirochaetota",])),
                 colorRampPalette(c("#ff512f","#dd2476"))(nrow(labtax[labtax$taxa == "Sva0485",])),
                 colorRampPalette(c("#06beb6","#48b1bf"))(nrow(labtax[labtax$taxa == "Verrucomicrobiota",])),
                 colorRampPalette(c("gray15", "gray75"))(nrow(labtax[labtax$Domain == "Archaea",])),
                 colorRampPalette(c("gray1", "gray1"))(nrow(labtax[labtax$Domain == "Other",])))

# format color table and create the named palette vector
pal.t<-labtax$colors;names(pal.t)<-labtax$taxa
gtaxa$taxa<-factor(gtaxa$taxa, levels = labtax$taxa, ordered = TRUE)

#plot filters and format
gtaxa<-gtaxa[gtaxa$type != "control",]
gtaxa$type2<-factor(ifelse(gtaxa$type %in% c("reference", "watercolumn"), "WC", "SED"), levels = c("WC", "SED"), ordered = TRUE)
gtaxa$location2<-ifelse(gtaxa$location == "north sea", "NS", "WS") 
gtaxa$season2<-factor(str_to_title(gtaxa$season ), levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)
gtaxa$DEPTH2<-factor(gtaxa$DEPTH, levels = c("surface", "depth", "sediment"), labels = c("Surface", "Bottom", "Sediments"), ordered = TRUE)

# explore %
gtax<-gtaxa
gtax<-gtaxa[gtaxa$type != "reference",]
gtax<-gtax[gtax$type == "sediment",]

gt<-ggplot(gtaxa[gtaxa$type2 != "SED" & gtaxa$type != "reference" ,], aes(x = time, y= n_reads, fill = taxa))+
  facet_grid(DEPTH2~location2+season2, scale = "free_x")+
  geom_area(stat = "identity", position = "fill")+
  scale_y_continuous(expand =c(0,0), breaks = c(0,0.5,1), labels = c("0", "0.5", "1"))+
  scale_x_datetime(expand =c(0,0), date_labels = "%k:%M", date_breaks = "6 hour")+
  xlab("Water-column temporal series\n(Consecutive days, hour)")+ylab("Proportions")+
  scale_fill_manual(name = "Phylum",values = c(pal.t))+
  ggtitle(paste("Microbial community composition", sep = ""))+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0.5))+
  #theme(axis.text.x = element_text(size = 16, angle = 30, hjust = 1, vjust = 1))+
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(strip.text =  element_text(size = 16, vjust = 0))+
  theme(strip.text.y =  element_text(size = 16, angle = 0, vjust = 0.5, hjust = 0))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.background = element_rect(fill = "transparent",colour = NA))+
  #theme(aspect.ratio = 1)+
  theme(panel.spacing = unit(1.25, "lines"))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gt

gs<-ggplot(gtaxa[gtaxa$type2 == "SED" ,], aes(x = depth, y= n_reads, fill = taxa))+
  facet_grid(.~location2+season2, scale = "free_x")+
  geom_bar(stat = "identity", position = "fill")+
  scale_y_continuous(expand =c(0,0), breaks = c(0,0.5,1), labels = c("0", "0.5", "1"))+
  scale_x_continuous(expand =c(0,0))+
  ylab("Proportions")+xlab("Sediment\nDepth (cm)")+
  scale_fill_manual(name = "Phylum",values = c(pal.t))+
  ggtitle(paste("", sep = ""))+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0.5))+
  #theme(axis.text.x = element_text(size = 16, angle = 30, hjust = 1, vjust = 1))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(strip.text =  element_text(size = 16, vjust = 0))+
  theme(strip.text.y =  element_text(size = 16, angle = 0, vjust = 0.5, hjust = 0))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.background = element_rect(fill = "transparent",colour = NA))+
  #theme(aspect.ratio = 1)+
  coord_flip()+
  theme(panel.spacing = unit(1.25, "lines"))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gs


##
# Figure 1
##

# 15 x 10
combined <- gt / gs +theme(strip.text.x =element_blank()) + plot_annotation(tag_levels = 'A')
combined+ plot_layout(guides = "collect", heights = c(1, 0.25)) & theme(legend.position = "right", legend.justification = "top" ) & guides(fill = guide_legend(ncol = 1))

###############################
# Bray-Curtis distance analysis
###############################

# bray curtis distance across ecosystems
# filter out lab controls and sediments for now
info3<-info[ !info$type %in%  c("reference", "control"),]
info3$Depth<-ifelse(info3$depth %in% c(-3, -35), "Deep", "Surface")
info3$Depth<-ifelse(info3$type %in% c("sediment") , "Sediment", info3$Depth)
#info3<-info[info$label %in% c("WS\nSed", "NS\nSed") & info$type != "reference",]
#info3$Depth<-info3$depth
asv3<-asv[,colnames(asv) %in% info3$sampleID]

# Bray curtis distance
tr.asv<-decostand(asv3, MARGIN = 2, method = 'total') # transform as proportions instead of abundances
bray.mat<-vegdist(t(tr.asv), method = "bray") # computes Bray-Curtis distance
mbc<-melt(as.matrix(bray.mat))
mbc$eco_Var1<-NA
mbc$eco_Var2<-NA

# loop to characterize and format the distance matrix for plotting
for (i in 1:nrow(mbc)){ # takes 15-20s
  mbc$eco_Var1[i]<-paste(info3[info3$sampleID == mbc$Var1[i], c("type")],info3[info3$sampleID == mbc$Var1[i], c("location")], info3[info3$sampleID == mbc$Var1[i], c("season")],info3[info3$sampleID == mbc$Var1[i], c("Depth")] , sep = "/")
  mbc$eco_Var2[i]<-paste(info3[info3$sampleID == mbc$Var2[i], c("type")], info3[info3$sampleID == mbc$Var2[i], c("location")], info3[info3$sampleID == mbc$Var2[i], c("season")],info3[info3$sampleID == mbc$Var2[i], c("Depth")] , sep = "/")
}
mbc<-mbc[mbc$Var1 != mbc$Var2,]
mbc$substrate<-ifelse(str_split_fixed(mbc$eco_Var1,  '/', n = 4)[,1] ==str_split_fixed(mbc$eco_Var2,  '/', n = 4)[,1], str_split_fixed(mbc$eco_Var1,  '/', n = 4)[,1], "Inter-Substrate")
mbc$sea<-ifelse(str_split_fixed(mbc$eco_Var1,  '/', n = 4)[,2] ==str_split_fixed(mbc$eco_Var2,  '/', n = 4)[,2], str_split_fixed(mbc$eco_Var1,  '/', n = 4)[,2], "Inter-Sea")
mbc$season<-ifelse(str_split_fixed(mbc$eco_Var1,  '/', n = 4)[,3] ==str_split_fixed(mbc$eco_Var2,  '/', n = 4)[,3], str_split_fixed(mbc$eco_Var1,  '/', n = 4)[,3], "Inter-Season")
mbc$depth<-ifelse(str_split_fixed(mbc$eco_Var1,  '/', n = 4)[,4] ==str_split_fixed(mbc$eco_Var2,  '/', n = 4)[,4], str_split_fixed(mbc$eco_Var1,  '/', n = 4)[,4], "Inter-Depth")
head(mbc)

# format
mbc$Sea<-factor(mbc$sea, levels = c("Inter-Sea","north sea", "wadden sea"), labels = c("Inter-Sea","North Sea", "Wadden Sea"), ordered = TRUE)
mbc$Season<-factor(mbc$season, levels = c("Inter-Season","winter", "spring","summer", "autumn"), labels = c("Inter-Season","Winter", "Spring","Summer", "Autumn"), ordered = TRUE)
mbc$Substrate<-factor(mbc$substrate, levels = c("Inter-Substrate","watercolumn", "sediment"), labels = c("Inter-Substrate","Watercolumn", "Sediment")  , ordered = TRUE)
mbc$Depth<-factor(mbc$depth, levels = c("Surface", "Deep", "Inter-Depth", "Sediment"),labels = c("Surface", "Bottom", "Inter-Depth", "Sediment"), ordered = TRUE)
mbc$Depth2<-NA
mbc$Depth2<-ifelse(mbc$Substrate == "Inter-Substrate", "Inter-Substrate", as.character(mbc$Depth))
mbc$Depth2<-factor(mbc$Depth2, levels = c("Surface", "Bottom", "Inter-Depth", "Sediment", "Inter-Substrate"), ordered = TRUE)

# plot distance
ggturn<-ggplot(data = mbc, aes(x = Season, y = value))+
  facet_grid(Substrate~Sea, space = "free",scales = "free" ,  drop=TRUE)+
  geom_hline(data =na.omit(mbc), aes(yintercept = 1), linetype = "longdash", size = 0.5, color = "gray90")+
  geom_hline(data =na.omit(mbc), aes(yintercept = 0), linetype = "longdash", size = 0.5, color = "gray90")+
  geom_hline(data =na.omit(mbc), aes(yintercept = 0.25), linetype = "longdash", size = 0.5, color = "gray90")+
  geom_hline(data =na.omit(mbc), aes(yintercept = 0.5), linetype = "longdash", size = 0.5, color = "gray90")+
  geom_hline(data =na.omit(mbc), aes(yintercept = 0.75), linetype = "longdash", size = 0.5, color = "gray90")+
  geom_jitter(aes(color = Depth2), size = 2, alpha = 0.5)+
  scale_color_manual(values = c("#7aadff","#063175","gray25", "orange4", "gray75"), name = "Sample Comparison")+
  scale_fill_manual(values = c("#7aadff","#063175","gray25", "orange4", "gray75"), name = "Sample Comparison")+
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1))+
  labs(x = " ", y = "Microbial Community Turnover\n(Bray-Curtis Dissimilarity)", title = "Community Comparisons\n(in Water Column)")+
  ggtitle("")+
  theme(axis.line.x = element_line(colour = 'gray15', size=1, lineend = "square"),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5))+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0.5))+
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1))+
  theme(axis.text.y = element_text(size = 16))+
  theme(strip.text.x =  element_text(size = 16, face = "bold", angle =0, hjust = 0.5))+
  theme(strip.text.y =  element_text(size = 16, face = "bold", hjust = 0))+
  theme(strip.background = element_blank())+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"),
        legend.key = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));ggturn

# Functional distance
fapro<-fapro[rownames(fapro) %in% rownames(asv3),]
cwm<-FD::functcomp(fapro, t(decostand(asv3, MARGIN = 2, method = 'total')), bin.num = TRUE, CWM.type = "all")
bray.cwm<-vegdist(cwm, method = "bray") # computes Bray-Curtis distance
mbc<-merge(mbc, melt(as.matrix(bray.cwm)), by = c("Var1", "Var2"))
mbc$Season<-ifelse(mbc$Season == "Inter-Season", "Inter-Season", "Same-Season")

##
# Figure 4
##

gcomp<-ggplot(data = mbc, aes(x = value.y, y= value.x, color = Depth2, shape = Season))+
  #geom_hex(data = bc.func, aes(x = func, y= comm),size = 4)+
  #scale_fill_viridis(name = "Count")+
  facet_grid(Substrate~Sea,  drop=TRUE)+
  geom_point(size = 2, alpha = 0.25)+
  scale_color_manual(values = c("#7aadff","#063175","gray25", "orange4", "gray75"), name = "Sample Comparison")+
  scale_y_continuous(limits = c(0,1), labels = c("0", "0.25", "0.5", "0.75","1"))+
  scale_x_continuous(limits = c(0,1), labels = c("0", "0.25", "0.5", "0.75","1"))+
  #theme_void()+
  labs(title = "",
       x = "Functional Turnover\n(Bray-Curtis Dissimilarity & CWM of FAPROTAX)",
       y = "Community Composition Turnover\n(Bray-Curtis Dissimilarity based on ASV Distribution)")+
  theme(strip.text.x =  element_text(size = 16, face = "bold", angle =0, hjust = 0.5))+
  theme(strip.text.y =  element_text(size = 16, face = "bold", angle =0, hjust = 0))+
  theme(strip.background = element_blank())+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0))+
  theme(axis.text = element_text(size = 16))+
  theme(axis.title = element_text(size = 16))+
  theme(panel.grid.major= element_line(colour = "white", size = 1))+
  theme(panel.grid.minor= element_line(colour = "white", size = 0.5))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  guides(color = guide_legend(override.aes = list(alpha = 1) ))+
  theme(panel.background = element_rect(fill = "#ebebfa",colour = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(aspect.ratio = 1);gcomp

####################################################
# Bray-Curtis and unifrac distance analysis throughout timeseries
####################################################
library(ecodist)

## Sediment:
# make subsets and compute distance or depth between points within samples of subsets
list_depth<-list()
list_mat_sed<-list()
for (i in c("wadden sea", "north sea")){
  for (j in c("winter", "spring", "summer", "autumn")){
      
      #subset data and info
      sed<-as.data.frame(asv[,colnames(asv) %in% info[info$location == i & info$season == j & info$type == "sediment" , "sampleID"]]) ; sed<-sed[rowSums(sed)>0, ]
      sed<-decostand(sed,MARGIN =2 ,method = "total")

      if(nrow(sed)==0) next
      info.sed<-info[info$sampleID %in% colnames(sed),];rownames(info.sed)<-info.sed$sampleID
      info.sed<-info.sed[info.sed$sampleID %in% colnames(sed),]
      nps<-prune_samples(colnames(otu_table(ps)) %in% rownames(info.sed), ps)
      nps<-prune_taxa(rowSums(nps@otu_table) > 0, nps)
      
      # Info on ASVs abundance
      rsed<-data.frame(nreads = rowSums(sed), cid = names(rowSums(sed)) )
      rsed<-rsed[order(rsed$nreads, decreasing = TRUE),];rsed$rank<-c(1:nrow(rsed))
      rsed$re_rank<-(rsed$nreads/sum(rsed$nreads))*100 
      rsed$rwat<-NA;rsed$cumul[1]<-rsed$re_rank[1]
      for (u in 2:nrow(rsed)) { rsed[u,"cumul"]<-rsed[u,"re_rank"] + rsed[u-1,"cumul"] }
      
      # distances
      bc.sed<-bcdist(t(sed))
      bc.sed.ab<-bcdist(t(sed[rownames(sed) %in% rsed[rsed$re_rank >0.1 ,"cid" ] ,]))
      uni.sed<-UniFrac(nps, weighted = TRUE, normalized = TRUE, FALSE)
      uni.sed.ab<-UniFrac(prune_taxa(rownames(nps@otu_table) %in% rsed[rsed$re_rank >0.1 ,"cid" ], nps), weighted = TRUE, normalized = TRUE, FALSE)
      func<-as.matrix(bray.cwm)[rownames(info.sed),rownames(info.sed)]
      euc.depth<-dist(info.sed[,c("depth")],method = "euclidean");euc.depth<-dist_setNames(euc.depth, rownames(info.sed))
      
      # save bray curtis and other distance
      mat.sed<-Reduce(function(x, y) merge(x, y,by = c("Var1", "Var2"), all=TRUE),list(melt(as.matrix(bc.sed)),melt(as.matrix(bc.sed.ab)),melt(as.matrix(uni.sed)),melt(as.matrix(uni.sed.ab)), melt(as.matrix(euc.depth)), melt(as.matrix(func)) ) )
      colnames(mat.sed)<-c("sp1", "sp2", "bc", "bc.ab","uni", "uni.ab", "depth.df", "func")
      
      # correlograms
      list_depth[[paste(i,j, sep = "/")]]<-as.data.frame(mgram(bc.sed, euc.depth,breaks = seq(1,7,2), nperm = 999  )$mgram )
      list_mat_sed[[paste(i,j, sep = "/")]]<-as.data.frame(mat.sed)
  }
}

# format list
df_mat_sed<-bind_rows(list_mat_sed, .id = "eco")
df_mat_sed$location<-sub("\\/.*","",df_mat_sed$eco)
df_mat_sed$location<-factor(df_mat_sed$location, levels = c("north sea", "wadden sea"), labels = c("NS", "WS"), ordered = TRUE)
df_mat_sed$depth<-sub(".*\\/", "", df_mat_sed$eco)
df_mat_sed$season<-sub("\\/.*","",sub(".*sea/", "", df_mat_sed$eco))
df_mat_sed$season<-factor(df_mat_sed$season, levels = c("winter", "spring", "summer", "autumn"), labels = str_to_title(c("winter", "spring", "summer", "autumn")) , ordered = TRUE)
df_mat_sed<-df_mat_sed[df_mat_sed$sp1 != df_mat_sed$sp2,]
df_mat_sed$depth

mat.summ<-df_mat_sed %>%
  group_by(depth.df, location, season) %>% 
  dplyr::summarize(mean = mean(bc),sd =sd(bc) )

mat.pv<-df_mat_sed %>%
  group_by(location, season) %>% 
  dplyr::summarize(depth.pv = cor.test(bc,depth.df)$p.value)
mat.pv<-mat.pv[mat.pv$depth.pv < 0.05,]

gsed<-ggplot(df_mat_sed, aes(x = depth.df, y= bc, color = location))+
  facet_grid( season~ unname(location) )+
  geom_hline(data =na.omit(df_mat_sed), aes(yintercept = 0.2), linetype = "longdash", size = 0.5, color = "gray90")+
  geom_hline(data =na.omit(df_mat_sed), aes(yintercept = 0.4), linetype = "longdash", size = 0.5, color = "gray90")+
  geom_hline(data =na.omit(df_mat_sed), aes(yintercept = 0.6), linetype = "longdash", size = 0.5, color = "gray90")+
  geom_hline(data =na.omit(df_mat_sed), aes(yintercept = 0.8), linetype = "longdash", size = 0.5, color = "gray90")+
  #geom_hline(data =na.omit(df_mat_sed), aes(yintercept = 0.8), linetype = "longdash", size = 0.5, color = "gray90")+
  geom_point(alpha = 0.1, shape = 19, size = 2)+
  #geom_pointrange(data = mat.summ[paste(mat.summ$location,mat.summ$season) %in% paste(mat.pv$location,mat.pv$season) ,] ,aes(ymin=mean-sd, ymax=mean+sd, y = mean, x = depth.df, color=location))+
  geom_line(data = mat.summ[paste(mat.summ$location,mat.summ$season) %in% paste(mat.pv$location,mat.pv$season) ,],aes( y = mean, x = depth.df, color=location), size= 2)+
  scale_color_manual(values = c("#4d2323", "#4d4b23"),guide = "none")+
  xlab("Sediment Depth Difference (cm)")+
  #ylab("Bray-Curtis Dissimilarity")+
  ylab("Microbial community turnover\n(Bray-Curtis dissimilarity)")+
  scale_x_continuous(breaks = seq(0,8, 2))+
  theme(legend.position = "right")+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 18,hjust = 0.5))+
  theme(axis.text = element_text(size = 16))+
  theme(axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 1))+
  theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5))+
  theme(panel.spacing = unit(0.5, "lines"))+
  theme(strip.background = element_blank())+
  theme(strip.text =  element_text(size = 16, face = "bold"))+
  #theme(strip.text.y  =  element_text(angle = 0))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(legend.key =  element_rect(fill = "#00000000"))+
  theme(panel.grid.minor = element_line(color = "#ffffff", size = 0.5))+
  theme(panel.grid.major = element_line(color = "#ffffff", size = 1))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gsed

## Water Column:
# make subsets and compute distance or time between points within samples of subsets
list_time<-list()
list_tide<-list()
list_mat<-list()
list_perm<-list()
for (i in c("wadden sea", "north sea")){
  for (j in c("winter", "spring", "summer", "autumn")){
    for (z in  c("surface", "depth")){
      
      #subset data and info
      wat<-as.data.frame(asv[,colnames(asv) %in% info[info$location == i & info$season == j & info$DEPTH == z & info$type == "watercolumn" , "sampleID"]]) ; wat<-wat[rowSums(wat)>0, ]
      wat<-decostand(wat, MARGIN = 2, "total")
      #wat<-as.data.frame(asv[,colnames(asv) %in% info[info$location == "wadden sea" & info$season == "winter" & info$DEPTH == "surface" & info$type == "watercolumn" , "sampleID"]]) ; wat<-wat[rowSums(wat)>0, ]
      
      # Info on time and tidal difference of sample selection
      if(nrow(wat)==0) next
      info.wat<-info2[info2$sampleID %in% colnames(wat),];rownames(info.wat)<-info.wat$sampleID
      info.wat<-info.wat[info.wat$sampleID %in% colnames(wat),]
      nps<-prune_samples(colnames(otu_table(ps)) %in% rownames(info.wat), ps)
      nps<-prune_taxa(rowSums(nps@otu_table) > 0, nps)
      
      # Info on ASVs abundance
      rwat<-data.frame(nreads = rowSums(wat), cid = names(rowSums(wat)) )
      rwat<-rwat[order(rwat$nreads, decreasing = TRUE),];rwat$rank<-c(1:nrow(rwat))
      rwat$re_rank<-(rwat$nreads/sum(rwat$nreads))*100 
      rwat$rwat<-NA;rwat$cumul[1]<-rwat$re_rank[1]
      for (u in 2:nrow(rwat)) { rwat[u,"cumul"]<-rwat[u,"re_rank"] + rwat[u-1,"cumul"] }
      
      # distance computation
      bc.wat<-ecodist::bcdist(t(wat))
      bc.wat.ab<-ecodist::bcdist(t(wat[rownames(wat) %in% rwat[rwat$re_rank >1 ,"cid" ] ,]))
      uni.wat<-UniFrac(nps, weighted = TRUE, normalized = TRUE,FALSE)
      uni.wat.ab<-UniFrac(prune_taxa(rownames(nps@otu_table) %in% rwat[rwat$re_rank >0.1 ,"cid" ], nps), weighted = TRUE, normalized = TRUE, FALSE)
      euc.time<-dist(info.wat[,c("time")],method = "euclidean");euc.time<-dist_setNames(euc.time, rownames(info.wat))
      euc.tide<-dist(info.wat[,c("tide")],method = "euclidean");euc.tide<-dist_setNames(euc.tide, rownames(info.wat))
      euc.temp<-dist(info.wat[,c("temp")],method = "euclidean");euc.temp<-dist_setNames(euc.temp, rownames(info.wat))
      euc.sal<-dist(info.wat[,c("sal")],method = "euclidean");euc.sal<-dist_setNames(euc.sal, rownames(info.wat))
      euc.light<-dist(info.wat[,c("light.ac")],method = "euclidean");euc.light<-dist_setNames(euc.light, rownames(info.wat))
      euc.wind.dir<-sapply(info.wat$Winddirection, function(x) sapply(info.wat$Winddirection, function(y) 180 - abs(abs(x - y) - 180)));euc.wind.dir<-dist_setNames(euc.wind.dir, rownames(info.wat))
      euc.wind.sp<-dist(info.wat[,c("Windspeed.avgph")],method = "euclidean");euc.wind.sp<-dist_setNames(euc.wind.sp, rownames(info.wat))
      func<-as.matrix(bray.cwm)[rownames(info.wat),rownames(info.wat)]
      
      # save bray curtis and other distance
      euc.depth<-dist(info.wat[,c("depth")],method = "euclidean");euc.depth<-dist_setNames(euc.depth, rownames(info.wat))
      
      # save bray curtis and other distance
      mat.wat<-Reduce(function(x, y)S4Vectors::merge(x, y,by = c("Var1", "Var2"), all=TRUE),list(melt(as.matrix(bc.wat)),melt(as.matrix(bc.wat.ab)),melt(as.matrix(uni.wat)),melt(as.matrix(uni.wat.ab)),  melt(as.matrix(euc.time)), melt(as.matrix(euc.tide)), melt(as.matrix(euc.light)),melt(as.matrix(euc.wind.dir)),melt(as.matrix(euc.wind.sp)),melt(as.matrix(func)), melt(as.matrix(euc.temp)),melt(as.matrix(euc.sal))) )
      colnames(mat.wat)<-c("sp1", "sp2", "bc", "bc.ab","uni", "uni.ab",  "time.df", "tide.df", "light.df", "wind.dir.df", "wind.sp.df", "func", "temp", "sal")
      mat.wat$hour.df<-mat.wat$time.df/60/60
      
      # Compute PERMANOVA
      perm<-adonis2(bc.wat ~ time + light.ac + temp + sal + tide + Winddirection + Windspeed.avgph,data = info.wat,by = "margin")
      list_perm[[paste(i,j,z, sep = "/")]]<-data.frame(env.var = rownames(perm), R2 = perm$R2, p.value = perm$`Pr(>F)`)
      
      # save correlograms
      list_time[[paste(i,j,z, sep = "/")]]<-as.data.frame(mgram(bc.wat, (euc.time/60)/60,breaks = seq(-2,48,4), nperm = 999  )$mgram )
      list_tide[[paste(i,j,z, sep = "/")]]<-as.data.frame(mgram(bc.wat, euc.tide,breaks = seq(-0.1,2.5,0.2), nperm = 999  )$mgram )
      list_mat[[paste(i,j,z, sep = "/")]]<-as.data.frame(mat.wat)
    }
  }
}

# format list
df_mat<-bind_rows(list_mat, .id = "eco")
df_mat$location<-sub("\\/.*","",df_mat$eco)
df_mat$location<-factor(df_mat$location, levels = c("north sea", "wadden sea"), labels  = c("NS", "WS"), ordered = TRUE)
df_mat$depth<-sub(".*\\/", "", df_mat$eco)
df_mat$depth<-factor(df_mat$depth, levels = c("surface", "depth"),labels = c("Surface", "Bottom"), ordered = TRUE)
df_mat$season<-sub("\\/.*","",sub(".*sea/", "", df_mat$eco))
df_mat$season<-factor(df_mat$season, levels = c("winter", "spring", "summer", "autumn"), labels = str_to_title(c("winter", "spring", "summer", "autumn")), ordered = TRUE)
df_mat<-df_mat[df_mat$sp1 != df_mat$sp2,]
df_mat$hour.df<-plyr::round_any(df_mat$hour.df, 2)

mat.summ<-df_mat %>%
  group_by(hour.df, location, depth, season) %>% 
  dplyr::summarize(mean = mean(bc),sd =sd(bc) )

mat.pv<-df_mat %>%
  group_by(location, season, depth) %>% 
  dplyr::summarize(time.pv = cor.test(bc,time.df)$p.value, tide.pv = cor.test(bc,tide.df)$p.value)
mat.pv<-mat.pv[mat.pv$time.pv < 0.01,]

gth<-ggplot(df_mat, aes(x = hour.df, y= bc, color = depth))+
  facet_grid( season~ location)+
  geom_hline(data =na.omit(df_mat), aes(yintercept = 0.2), linetype = "longdash", size = 0.5, color = "gray90")+
  geom_hline(data =na.omit(df_mat), aes(yintercept = 0.4), linetype = "longdash", size = 0.5, color = "gray90")+
  geom_hline(data =na.omit(df_mat), aes(yintercept = 0.6), linetype = "longdash", size = 0.5, color = "gray90")+
  geom_point(alpha = 0.25, shape = 19, size = 2)+
  #geom_pointrange(data = mat.summ[paste(mat.summ$location,mat.summ$season) %in% paste(mat.pv$location,mat.pv$season) ,] ,aes(ymin=mean-sd, ymax=mean+sd, y = mean, x = hour.df))+
  geom_line(data = mat.summ[paste(mat.summ$location,mat.summ$season) %in% paste(mat.pv$location,mat.pv$season) ,],aes( y = mean, x = hour.df), size= 2)+
  xlab("Hour Difference (h)")+
  ylab("Microbial community turnover\n(Bray-Curtis dissimilarity)")+
  #ylab("Unifrac Distance")+
  scale_color_manual(values = c("#7aadff","#063175"), name = "Depth")+
  scale_x_continuous(breaks = seq(0,48, 12), limits = c(0, 50))+
  theme(legend.position = "bottom", legend.direction = "vertical")+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 16,hjust = 0.5))+
  theme(axis.text = element_text(size = 16))+
  theme(axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 1))+
  theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5))+
  theme(panel.spacing = unit(0.5, "lines"))+
  theme(strip.background = element_blank())+
  theme(strip.text =  element_text(size = 16, face = "bold"))+
  theme(strip.text.y  =  element_text(angle = 0))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(legend.key =  element_rect(fill = "#00000000"))+
  theme(panel.grid.minor = element_line(color = "#ffffff", size = 0.5))+
  theme(panel.grid.major = element_line(color = "#ffffff", size = 1))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gth

#View(df_mat[df_mat$location %in% c("WS") & df_mat$season %in% c("Spring") ,])

# general influence
list_rel<-bind_rows(list_perm, .id = "system")
list_rel$system<-gsub("depth","bottom", list_rel$system)
list_rel$system<-factor(str_to_title(list_rel$system), levels = rev(c("North Sea/Summer/Surface","North Sea/Summer/Bottom", "North Sea/Autumn/Surface","North Sea/Autumn/Bottom","Wadden Sea/Winter/Surface","Wadden Sea/Winter/Bottom", "Wadden Sea/Spring/Surface","Wadden Sea/Spring/Bottom","Wadden Sea/Summer/Surface", "Wadden Sea/Autumn/Surface","Wadden Sea/Autumn/Bottom")), ordered = TRUE)
list_rel$R_name<-ifelse(list_rel$p.value < 0.05, list_rel$env.var, "N.S.")
list_rel<-list_rel[!list_rel$env.var %in% c("Residual", "Total"),]
list_rel$R_name<-factor(list_rel$R_name, 
                              levels = (c("N.S.","time","light.ac","temp", "sal", "tide", "Winddirection", "Windspeed.avgph")),
                              labels = (c( "N.S.","Time","Accumulated Daylight","Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed")), ordered = TRUE)

pal_cond=c("#bd5757", "#6e3232","#9863c9","#37234a","#3b79f5","#1d3c7a", "#e5f53b","#727a1d", "#3bed97", "#ed8139","#7a421d")
names(pal_cond)<-rev(levels(list_rel$system))

pal.env<-c("gray75","black","gold3","#f0624f", "#913e33", "#224282","#71bf9b","#46705d")
names(pal.env)<-c("N.S.","Time","Accumulated Daylight","Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed")

gcar<-ggplot(data = list_rel,aes(x= R2, y=system , fill = R_name))+
  geom_bar(position = "stack",stat ="identity")+
  scale_fill_manual(name = "Env. Variable",values= pal.env)+
  scale_x_continuous(expand = c(0,0), limits = c(0,2))+
  ylab("Water column turnover")+xlab(expression ("R"^2))+ 
  theme(axis.line.x = element_line(colour = 'gray15', linewidth=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', linewidth=1, lineend = "square"))+
  theme(axis.title = element_text(size = 18,hjust = 0.5))+
  theme(axis.text = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16, angle = 0, hjust = 1, vjust = 0.5, color = rev(pal_cond)))+
  theme(panel.spacing = unit(0.5, "lines"))+
  theme(strip.background = element_blank())+
  theme(strip.text =  element_text(size = 16, face = "bold"))+
  theme(strip.text.y  =  element_text(angle = 0, hjust = 0))+
  theme(legend.key =  element_rect(fill = "#00000000"))+
  #theme(legend.position = "bottom", legend.direction = "vertical")+
  theme(panel.grid.major= element_line(colour = "white", linewidth = 1))+
  theme(panel.grid.minor= element_line(colour = "white", linewidth = 0.5))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.grid.minor.y = element_blank())+
  theme(panel.grid.major.y = element_blank())+
  theme(panel.grid.major.x = element_line(color = "gray90", linewidth = 1, linetype ="longdash" ))+
  theme(panel.grid.minor.x  = element_line(color = "gray90", linewidth = 0.5, linetype ="longdash"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0),
        plot.subtitle = element_text(size = 16, face = "plain",hjust = 0) );gcar

gcorel<-gth+theme(strip.text.y = element_blank())+ylab("Microbial Community Turnover\n(Bray-Curtis Dissimilarity)")+   #+ylab("Functional Turnover\n(Bray-Curtis Dissimilarity &\nCWM of FAPROTAX)")+
  gsed+theme(strip.text.y = element_text(hjust =0))  +ylab("")+
  gcar+
  plot_layout(nrow = 1)

##
# Figure 3
##

# 16 x 18
(ggturn/gcorel)+plot_annotation(tag_levels = 'A') & theme(legend.justification = 'left')

# explore specific conditions
exp<-melt(info2[info2$location %in% c("north sea") & info2$season %in% c("summer") & info2$DEPTH %in% c("surface"),c(2, 20,21, 26,27,33,34)], id.vars ="time" )
ggplot(data = exp, aes(x = time, y = value, color = variable))+
  facet_grid(variable~., scales = "free")+
  geom_line()

##################################################
# Connectivity: Networks of average Shared OTUs on a map
##################################################

## Data formating
# order/set dataset
info3<-info[info$type == "watercolumn",]
#info3<-info[info$location != "lab" & info$type != "control",]
asv3<-asv[,colnames(asv) %in% info3$sampleID]
tasv<-t(asv3)
tasv<-tasv[order(rownames(tasv)),] # OG dataset
info4<-info[info$sampleID %in% rownames(tasv),]
info4<-info4[order(info4$sampleID),]
info4$depth<-factor(ifelse(info4$depth %in% c(-1, -10), "Surface", "Bottom" ), levels = c("Surface", "Bottom"), ordered = TRUE)
info4$id<-ifelse(info4$location %in% c("north sea") & info4$season %in% c("summer") ,paste(str_to_title(info4[,4]),str_to_title(info4[,5]),str_to_title(info4[,8]), sep = '/'),paste(str_to_title(info4[,4]),str_to_title(info4[,5]), sep = '/'))

## create data for networking
# Number of shared OTUs across ecosystems / links/edges
beta<-as.matrix(betadiver(tasv)$a);beta[upper.tri(beta)] <- NA
beta<-reshape2::melt(as.matrix(beta))
beta<-na.omit(beta);beta<-beta[beta$Var1 != beta$Var2 ,]
link<-as_tibble(beta);colnames(link)<-c("from", "to","weight")
link<-merge(link, info4[,c(2,ncol(info4))],by.x = 'from', by.y = "sampleID" ); colnames(link)[ncol(link)]<-"from_id"
link<-merge(link, info4[,c(2,ncol(info4))],by.x = 'to', by.y = "sampleID" ); colnames(link)[ncol(link)]<-"to_id"
link<-aggregate(weight ~ ., data = link[,3:5], FUN = mean)
lev<-c("north sea/summer/Bottom", "north sea/summer/Surface","north sea/autumn","wadden sea/winter", "wadden sea/spring", "wadden sea/summer","wadden sea/autumn")
link$from_id<-factor(link$from_id, levels = str_to_title(lev), ordered = TRUE)
link$to_id<-factor(link$to_id, levels = str_to_title(lev), ordered = TRUE)
link<-link[link$from_id != link$to_id,]
link<-link[-4,]

betaa<-merge(beta,rowSums(tasv != 0), by.x = "Var1", by.y = "row.names")
betaa<-merge(betaa,rowSums(tasv != 0), by.x = "Var2", by.y = "row.names")
colnames(betaa)<-c("Var2","Var1", "Connectivity", "Richness.1",  "Richness.2")
ggplot(betaa, aes(x =Richness.1,y=Richness.2,color = Connectivity))+
  geom_point()

link2<-link[link$from_id != "North Sea/Summer/Surface" & link$to_id != "North Sea/Summer/Surface",]
link2$sea_from<-str_split_fixed(link2$from_id, "/", 3)[,1]
link2$sea_to<-str_split_fixed(link2$to_id, "/", 3)[,1]
link2[link2$sea_to == "North Sea" & link2$sea_from == "North Sea", ]

aggregate(weight~sea_from+sea_to,link2, mean)

# format for networks, number of OTUs per ecosystem / Nodes
#notus<-as.data.frame(rowSums(cum_asv != 0));colnames(notus)<-"weight" # total number of OTUs per ecosystem
notus<-aggregate(as.data.frame(rowSums(tasv != 0)), by = list(info4$id), FUN = mean) ;colnames(notus)<-c('id',"weight"); rownames(notus)<-notus[,1] # average number of OTUs per ecosystem
notus$id<-factor(rownames(notus), levels = str_to_title(lev), ordered = TRUE);rownames(notus)<-NULL
notus<-notus[order(notus$id),]
notus$Long<-as.double(c(4.089-0.5, 4.089-0.1,4.089+0.5,5.220-1,5.220-0.5,5.220,5.220))
notus$Lat<-as.double(c(55.306-0.75, 55.306-0.25,55.306,53.180-0.5,53.180-0.5,53.180,53.180+0.5))
notus$Season<-toupper(substr( gsub("/.*","" ,gsub("^[^/]*/","" ,notus$id)), start = 0,stop = 3))
notus$Season<-c("SUM B", "SUM S", "AUT", "WIN", "SPR", "SUM", "AUT")
node<-notus
node[grep("wadden", node$id),"Lat"]<-node[grep("wadden", node$id),"Lat"]+0.4
node[grep("north", node$id),"Lat"]<-node[grep("north", node$id),"Lat"]-0.25
node<-merge(node, aggregate(data.frame(S = vegan::diversity(tasv, index = "simpson"), H = vegan::diversity(tasv, index = "shannon")), by = list(info4$id), FUN = mean), by.y = "Group.1", by.x = "id")
plot(node$weight, node$H)

# Objects for networking plots
net <- graph_from_data_frame(link, directed = F, vertices = node)
net_for_plot <- link %>%
  inner_join(node %>% dplyr::select(id, Long, Lat), by = c('from_id' = 'id')) %>%
  inner_join(node %>% dplyr::select(id, Long, Lat), by = c('to_id' = 'id'))
colnames(net_for_plot)[4:7]<-c("x","y","xend","yend")
node_geo <- node ; colnames(node_geo)[c(3:4)]<-c("x","y")
layout <- create_layout(graph = net, layout = node_geo[,c(1,3,4)])
attributes(layout)

## Graphs preparation
# theme of the plot
maptheme <- 
  theme_bw()+
  theme(panel.border = element_rect(size = 1.5, color = "gray5", fill = FALSE)) +
  theme(panel.grid = element_blank()) +
  theme(title = element_text(size = 16, color = "gray5", face = "bold"))+
  theme(axis.text = element_text(size = 16, color = "gray5"))+
  theme(axis.title = element_blank()) +
  theme(legend.text = element_text(size = 16, color = "gray5"))+
  theme(legend.title = element_text(size = 16, color = "gray5", face = "bold"))+
  theme(legend.position = "right") +
  theme(panel.grid = element_blank()) +
  theme(legend.key = element_rect(colour = NA, fill = NA))+
  theme(panel.background = element_rect(fill = "#596673"))

# Geographic objects
nd.map <- readRDS("~/Desktop/WSNS/gadm36_NLD_2_sp.rds")
nd.map.f <- subset(nd.map, !nd.map$NAME_1  %in% c("Zeeuwse meren", "IJsselmeer"))
nd.map.f <- fortify(nd.map.f)
country_shapes <- geom_polygon(data = nd.map.f, aes(x = long, y = lat, group = group),fill = "#CECECE", color = "#CECECE", size = 0.15)
mapcoords <- coord_fixed(xlim = c(3.25,5.75), ylim = c(52.5,55.5))
sp_coord<-unique(info4[info4$type != "reference" ,c(4,16,17)]); sp_coord$Name<-c("NS", "WS")

# formating the labels for geographic coordinates
ewbrks <- seq(3,6,1)
nsbrks <- seq(52,56,1)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(abs(x), "°W"), ifelse(x > 0, paste(abs(x), "°E"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(x, "°S"), ifelse(x > 0, paste(x, "°N"),x))))

link[order(link$weight),]
layout$face<-"ASVs Connectivity\nin the NS & WS"
# plot
ggnet<-ggraph(graph = layout) + country_shapes +
  facet_grid(.~face)+
  ## edges/links
  geom_edge_arc(aes(edge_width = weight, circular = FALSE), strength = 0, alpha = 0.2) +
  scale_edge_width_continuous(name = "Average # of Shared ASVs",range = c(1, 6), limits = c(95, 801), breaks = seq(100, 800, 200)) +
  #scale_edge_width_continuous(name = "Average #ASVs Shared",range = c(1, 6), limits = c(100, 850), breaks = seq(100, 850, 200)) +
  ## nodes/ecosystems
  geom_node_point(aes(size = weight, fill = id), shape = 21,stroke = 0.15,) +
  scale_size_continuous( name = "Average # of ASVs",range = c(8, 20), limits = c(300, 1700), breaks = seq(300, 1700, 400)) +
  #scale_size_continuous( name = "Average #ASVs",range = c(6, 18), limits = c(600, 1700), breaks = seq(600, 1700, 400)) +
  scale_fill_manual(values = c("#6e3232", "#bd5757","#9863c9","#3b79f5", "#e5f53b", "#3bed97", "#ed8139"), guide = "none")+
  geom_node_text(aes(label = Season), repel = FALSE, size = 5,color = "grey5", fontface = "bold" )+
  ## Place Stations
  #geom_point(data = sp_coord, aes(y= lat,x = lon), shape = 15 , size = 3, color = "gray5")+
  #geom_text(data = sp_coord, aes(y= lat,x = lon, label = Name), hjust=-0.5, vjust=0.5, fontface = "bold")+
  ## plot details
  mapcoords + maptheme+
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0))+
  #ggtitle("Water-column ASVs Connectivity\nIn The Wadden and North Seas")+
  annotate("text",x=5.25, y=55, label = 'italic("North \n     Sea")',color = "grey80", parse = TRUE, size = 7)+
  theme(strip.background = element_blank())+
  theme(strip.text =  element_text(size = 16, face = "bold", hjust = 0))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));ggnet

# add value as a weight attribute
connect.rich<-merge(melt(as_adjacency_matrix(graph.data.frame(rbind(link, link[,c(2,1,3)]), directed=FALSE),type="both",names=TRUE,sparse=FALSE,attr="weight")), layout, by.x = "Var1", by.y = "id")
connect.rich<-connect.rich[connect.rich$Var1 != connect.rich$Var2,]
connect.rich$sea_from<-str_split_fixed(connect.rich$Var1, "/", 3)[,1]
connect.rich$sea_to<-str_split_fixed(connect.rich$Var2, "/", 3)[,1]
connect.rich$cross<-ifelse(connect.rich$sea_from == connect.rich$sea_to, "Same-Sea", "Inter-Sea")

ggplot(connect.rich,aes(x=weight, y=value, color = name, group=cross, shape = cross))+
  geom_point(size=3)+
  geom_smooth(method="lm", color = "gray25")+
  scale_color_manual(values = c("#6e3232", "#bd5757","#9863c9","#3b79f5", "#e5f53b", "#3bed97", "#ed8139"))+
  labs(x = "Average # of ASVs per sample", y= "Average # of ASVs shared with\nsamples from other conditions")+
  theme(plot.title = element_text(size = 24, face = "bold",hjust = 0))+
  scale_y_continuous(limits = c(0,850), expand = c(0,0))+
  #scale_x_continuous(limits = c(200,1700), expand = c(0,0))+
  theme(axis.text = element_text(size = 16))+
  theme(axis.title =  element_text(size = 16, face = "bold"))+
  theme(panel.grid.major= element_line(colour = "white", size = 1))+
  theme(panel.grid.minor= element_line(colour = "white", size = 0.5))+
  theme(panel.background = element_rect(fill = "#ebebfa",colour = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(aspect.ratio = 1)

wc<-link;colnames(wc)[3]<-"wc"
sediment<-link;colnames(sediment)[3]<-"sed"

comp<-merge(wc, sediment, by = c("from_id", "to_id"))
hull_data <- comp %>% slice(chull(wc, sed))

ggplot()+
  geom_point(data = comp, aes(x = wc, y= sed), size = 6, alpha = 0.75)+
  geom_polygon(data = hull_data,aes(x = wc, y= sed),alpha = 0.3,show.legend = FALSE)+
  labs(x = "#ASVs shared in Water Column", y= "#ASVs shared in Sediments")+
  theme(plot.title = element_text(size = 24, face = "bold",hjust = 0))+
  theme(axis.text = element_text(size = 16))+
  theme(axis.title =  element_text(size = 16, face = "bold"))+
  theme(panel.grid.major= element_line(colour = "white", size = 1))+
  theme(panel.grid.minor= element_line(colour = "white", size = 0.5))+
  theme(panel.background = element_rect(fill = "#ebebfa",colour = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(aspect.ratio = 1)

## Testing the connectivity
# generate data to perform testing
betab<-merge(beta, info[, c("sampleID", "DEPTH", "season", "location")], by.x ="Var1", by.y = "sampleID"  )
betab<-merge(betab, info[, c("sampleID", "DEPTH", "season", "location")], by.x ="Var2", by.y = "sampleID"  )
betab$sea<-ifelse(betab$location.x != betab$location.y , "Cross-Sea", str_to_title(betab$location.x))
betab$season<-ifelse(betab$season.x != betab$season.y , "Cross-Season", str_to_title(betab$season.y))
betab$season<-factor(betab$season, levels = c("Winter", "Spring", "Summer", "Autumn", "Cross-Season"), ordered = TRUE )

## test 
o.betab<-betab[betab$sea == "Cross-Sea",]
o.betab[order(o.betab$value, decreasing = TRUE),]

kruskal.test(betab$value,factor(betab$sea))

aggregate(value~sea, betab, mean)
kruskalmc(betab$value,betab$sea, probs = 0.001)
ggplot(data = betab, aes(x = sea ,y = value))+
  theme_bw()+
  #scale_y_continuous(expand = c(0,0), limits = c(0, 3000), breaks = seq(0,3000, 1000))+
  xlab("")+ylab("#ASV shared by\npairs of samples")+ 
  geom_boxplot(fill = "gray5", color = "gray25")+
  theme(legend.position = "none")+
  theme(axis.title = element_text(size = 18),plot.title = element_text(size = 18),legend.text=element_text(size=18))+
  theme(axis.text.y = element_text(size = 16))+theme(axis.text.x = element_text(size = 16, face = "bold", colour = "Black"))+
  theme(strip.text = element_text(size = 18))+
  theme(strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"))+
  theme(panel.grid.major.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(aspect.ratio = 1)

## test 
kruskal.test(betab$value,factor(betab$sea))
kruskalmc(betab$value,betab$sea, probs = 0.001)
ggplot(data = betab, aes(x = season ,y = value))+
  facet_grid(sea~.)+
  theme_bw()+
  #scale_y_continuous(expand = c(0,0), limits = c(0, 3000), breaks = seq(0,3000, 1000))+
  xlab("")+ylab("#ASV shared by\npairs of samples")+ 
  geom_boxplot(fill = "gray5", color = "gray25")+
  theme(legend.position = "none")+
  theme(axis.title = element_text(size = 18),plot.title = element_text(size = 18),legend.text=element_text(size=18))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.text.x = element_text(size = 16, face = "bold", colour = "Black", angle = 45, hjust = 1, vjust=1))+
  theme(strip.text = element_text(size = 18))+
  theme(strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid"))+
  theme(panel.grid.major.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(aspect.ratio = 1)


########################@########################@########################
# Heatmap functions
########################@########################@########################

## Create sub-datasets
# sediment North Sea
sed_ns_aut<-asv[,colnames(asv) %in% info[info$type == "sediment" & info$location == "north sea" & info$season == "autumn","sampleID"] ]; sed_ns_aut<-sed_ns_aut[rowSums(sed_ns_aut)>0,]
sed_ns_sum<-asv[,colnames(asv) %in% info[info$type == "sediment" & info$location == "north sea" & info$season == "summer","sampleID"] ]; sed_ns_sum<-sed_ns_sum[rowSums(sed_ns_sum)>0,]

# sediment Wadden Sea
sed_ws_win<-asv[,colnames(asv) %in% info[info$type == "sediment" & info$location == "wadden sea" & info$season == "winter","sampleID"]]; sed_ws_win<-sed_ws_win[rowSums(sed_ws_win)>0,]
sed_ws_spr<-asv[,colnames(asv) %in% info[info$type == "sediment" & info$location == "wadden sea" & info$season == "spring","sampleID"]]; sed_ws_spr<-sed_ws_spr[rowSums(sed_ws_spr)>0,]
sed_ws_sum<-asv[,colnames(asv) %in% info[info$type == "sediment" & info$location == "wadden sea" & info$season == "summer","sampleID"]]; sed_ws_sum<-sed_ws_sum[rowSums(sed_ws_sum)>0,]
sed_ws_aut<-asv[,colnames(asv) %in% info[info$type == "sediment" & info$location == "wadden sea" & info$season == "autumn","sampleID"]]; sed_ws_aut<-sed_ws_aut[rowSums(sed_ws_aut)>0,]

# watercolumn North Sea
wc_ns_aut<-asv[,colnames(asv) %in% info[info$type == "watercolumn" & info$location == "north sea" & info$season == "autumn","sampleID"] ]; wc_ns_aut<-wc_ns_aut[rowSums(wc_ns_aut)>0,]
wc_ns_sum_surf<-asv[,colnames(asv) %in% info[info$type == "watercolumn" & info$location == "north sea" & info$season == "summer"& info$DEPTH == "surface","sampleID"] ]; wc_ns_sum_surf<-wc_ns_sum_surf[rowSums(wc_ns_sum_surf)>0,]
wc_ns_sum_bott<-asv[,colnames(asv) %in% info[info$type == "watercolumn" & info$location == "north sea" & info$season == "summer"& info$DEPTH == "depth","sampleID"] ]; wc_ns_sum_bott<-wc_ns_sum_bott[rowSums(wc_ns_sum_bott)>0,]

# watercolumn Wadden Sea
wc_ws_win<-asv[,colnames(asv) %in% info[info$type == "watercolumn" & info$location == "wadden sea" & info$season == "winter","sampleID"]]; wc_ws_win<-wc_ws_win[rowSums(wc_ws_win)>0,]
wc_ws_spr<-asv[,colnames(asv) %in% info[info$type == "watercolumn" & info$location == "wadden sea" & info$season == "spring","sampleID"]]; wc_ws_spr<-wc_ws_spr[rowSums(wc_ws_spr)>0,]
wc_ws_sum<-asv[,colnames(asv) %in% info[info$type == "watercolumn" & info$location == "wadden sea" & info$season == "summer","sampleID"]]; wc_ws_sum<-wc_ws_sum[rowSums(wc_ws_sum)>0,]
wc_ws_aut<-asv[,colnames(asv) %in% info[info$type == "watercolumn" & info$location == "wadden sea" & info$season == "autumn","sampleID"]]; wc_ws_aut<-wc_ws_aut[rowSums(wc_ws_aut)>0,]

## create list of datasets
dfList <- list(sed_ns_aut,sed_ns_sum,sed_ws_win,sed_ws_spr,sed_ws_sum, sed_ws_aut,
               wc_ns_aut,wc_ns_sum_surf, wc_ns_sum_bott,wc_ws_win,wc_ws_spr,wc_ws_sum, wc_ws_aut)
names(dfList)<-c("sed_ns_aut","sed_ns_sum","sed_ws_win","sed_ws_spr","sed_ws_sum", "sed_ws_aut",
                    "wc_ns_aut","wc_ns_sum_surf","wc_ns_sum_bott","wc_ws_win","wc_ws_spr","wc_ws_sum", "wc_ws_aut")

## compute ASVs total number of reads per dataset and merge with functional information
# loop across list of datasets
list_funct<-list()
nreads<-list()
for (i in 1:length(dfList)){
  mat<-dfList[[i]]
  df<-data.frame(nreads = rowSums(mat), perc = rowSums(mat)/sum(rowSums(mat)),  cid = rownames(mat), DS = names(dfList)[i] ) 
  dff<-merge(df, fapro, by.x = "cid", by.y = "row.names" )
  list_funct[[names(dfList)[i] ]]<-dff
  nreads[[names(dfList)[i]]]<-data.frame(nreads = sum(rowSums(mat)))
}

# Study amount of annotated reads
list_funct<-bind_rows(list_funct, .id = "DS")
list_funct2<-list_funct
list_funct2$ntraits<-ifelse(rowSums(list_funct2[,5:ncol(list_funct)]) == 0 , "No Functional Annotation","Functional Annotation" )
a.funct2<-merge(aggregate(nreads~DS+ntraits, list_funct2, FUN = sum), bind_rows(nreads, .id = "DS"), by = "DS")
a.funct2$perc<-(a.funct2$nreads.x/a.funct2$nreads.y)*100

# study functions distribution
m.funct<-melt(list_funct, id.vars = c("cid", "nreads","perc", "DS"))
m.funct<-m.funct[m.funct$value >0,]
m.funct$season<-str_split_fixed(m.funct$DS, "_", 4)[,3]
m.funct$eco<-ifelse(str_split_fixed(m.funct$DS, "_", 4)[,1] =="sed", "Sediment", "Watercolumn")
m.funct$eco<-factor(m.funct$eco, levels = c("Watercolumn", "Sediment"), ordered = TRUE)
m.funct$sea<-ifelse(str_split_fixed(m.funct$DS, "_", 4)[,2]=='ns', "North Sea", "Wadden Sea")
m.funct$depth<-ifelse(str_split_fixed(m.funct$DS, "_", 4)[,4]=='bott', "Bottom", 
                      ifelse(str_split_fixed(m.funct$DS, "_", 4)[,4]=='surf', "Surface", ""))

# compute functions total abundance and percentage across all datasets
a.reads<-aggregate(nreads ~ variable, data = m.funct, FUN = sum)
a.reads<-a.reads[order(a.reads$nreads),]
a.reads$perca<-(a.reads$nreads/sum(as.data.frame(lapply(dfList, sum))))*100
a.reads$percca<-factor(cut(a.reads$perca, breaks = c(-Inf,0.1,1, 10, 20), labels = c("< 0.1 %", "0.1-1 %","1-10 %", "10-20 %") ), ordered = TRUE)
a.reads$percb<-(a.reads$nreads/sum(a.funct2[a.funct2$ntraits == "Functional Annotation", "nreads.x"]))*100
a.reads$perccb<-factor(cut(a.reads$percb, breaks = c(-Inf,0.1,1,5, 30, 60, 65), labels = c("< 0.1 %", "0.1-1 %","1-5 %", "5-30 %",  "30-60 %",  "60-65 %" ) ), ordered = TRUE)

# compute functions total abundance and percentage in each dataset
a.funct<-aggregate( cbind(nreads,perc) ~ DS+variable+season+eco+sea+depth, data = m.funct, FUN = sum)
a.funct$percc<-factor(cut(a.funct$perc*100, breaks = c(-Inf,1,2, 5, 20, 36), labels = c("< 1 %", "1-2 %","2-5 %", "5-20 %", "20-35 %") ), ordered = TRUE)
a.funct$variable<-factor(a.funct$variable, a.reads[order(a.reads$nreads),"variable"], ordered =TRUE)
a.funct$season<-factor(a.funct$season, levels = c("win", "spr", "sum", "aut"), labels = c("Winter", "Spring", "Summer", "Autumn") ,ordered =TRUE)
a.funct$conditions<-paste(a.funct$season, a.funct$depth, sep = " ")
a.funct$conditions<-factor(a.funct$conditions, levels = c("Winter ", "Spring ", "Summer ", "Summer Surface",  "Summer Bottom", "Autumn "),ordered =TRUE)
a.funct<-merge(a.funct, a.reads[,c(1,4,6)],by = "variable" )
a.funct$percca<-factor(a.funct$percca, levels = rev(levels(a.funct$percca)),ordered = TRUE)
a.funct$perccb<-factor(a.funct$perccb, levels = rev(levels(a.funct$perccb)),ordered = TRUE)
a.funct<-merge(a.funct, a.funct2[a.funct2$ntraits == "Functional Annotation",], by = "DS")
a.funct$perccannot<-(a.funct$nreads/a.funct$nreads.x)*100
a.funct$perccannotb<-factor(cut(a.funct$perccannot, breaks = c(-Inf,1,5, 20, 50, 100), labels = c("< 1 %", "1-5 %", "5-20 %", "20-50 %", "50-100 %") ), ordered = TRUE)

# plot heatmap of functions across datasets (15 x12)
ggh<-ggplot(data = a.funct, aes(y = variable, x = conditions, fill = perccannotb))+
  facet_grid(perccb~eco+sea, scales = "free", space = "free")+
  geom_tile()+
  scale_fill_brewer(palette="YlOrRd", name = '% per dataset', )+
  labs(x = "Conditions", y = "Functions inferred with FAPROTAX\n(sorted by total read representativity)")+
  theme_bw()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(strip.text =  element_text(size = 16, face = "bold", angle =0, hjust = 0.5),
        strip.text.y =  element_text(size = 16, face = "bold", angle =0, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1))+
  theme(strip.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(panel.background = element_rect(color = NA, fill = "gray90"),
        panel.border = element_rect(colour = "transparent", fill = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA));ggh

a.funct[
        a.funct$variable == "cellulolysis" 
        & a.funct$eco == "Sediment"
        #& a.funct$season == "Summer" 
        #& a.funct$sea == "North Sea"
        #& a.funct$depth == "Surface"
        #& a.funct$perccannotb == "1-5 %"
        ,]

## Compute donuts of number of annotated reads per conditions 
a.funct2$eco<-ifelse(str_split_fixed(a.funct2$DS, "_", 4)[,1] =="sed", "Sediment", "Watercolumn")
a.funct2$eco<-factor(a.funct2$eco, levels = c("Watercolumn", "Sediment"), ordered = TRUE)
a.funct2$sea<-ifelse(str_split_fixed(a.funct2$DS, "_", 4)[,2]=='ns', "North Sea", "Wadden Sea")

don<-aggregate(cbind(nreads.x,nreads.y)~ eco+sea+ntraits, data = a.funct2, FUN = sum )
colnames(don)[4:5]<-c("nreads", "tot")
don$perc<-(don$nreads/don$tot)*100

don<-don[order(don$eco, don$sea),];rownames(don)<-NULL
don$ymax.p<-don$perc;for (i in seq(2,8, 2)){don[i,"ymax.p" ]<-100}
don$ymin.p<-0;for (i in seq(2,8, 2)){don[i,"ymin.p" ]<-don[i-1,"perc" ]}
don$label.p <- (don$ymax.p + don$ymin.p) / 2
don$ntraits<-factor(don$ntraits, levels = c("Functional Annotation", "No Functional Annotation"))

# donuts reads
gdr<-ggplot(don, aes(ymax=ymax.p, ymin=ymin.p, xmax=4, xmin=3, fill=ntraits)) +
  facet_grid(.~eco+sea)+
  geom_rect() +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  scale_fill_manual(name = "",values = c("#bd0026","gray90") )+
  #ggtitle(paste("",lev_tax, " Level", sep = ""))+
  geom_text( x=3.5, aes(y=label.p, label=paste(round(perc,1),"%")), size=4, fontface = "bold", ) +
  #geom_text( x=2, y=2, label= "% reads",size=6, fontface = "bold" ) +
  theme(axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank())+
  theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5))+
  theme(axis.title= element_blank())+
  theme(axis.text = element_blank())+
  theme(strip.text =  element_text(size = 16, face = "bold"))+
  theme(strip.text.y =  element_text(angle = 0))+
  theme(strip.background = element_blank())+
  theme(legend.justification=c(0, 1))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(aspect.ratio = 1)+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0))+
  guides(fill = guide_legend(ncol = 1));gdr

layout <- "
A
B
B
B
"

##
# Figure 5
##

gdr/ggh+theme(strip.text.x = element_blank()) + plot_layout(design = layout) & theme(plot.margin = unit(c(0,0,0,0), "cm"))

#######################################################################
# Characterization of the Origin of ASVs in the WC of the WS in Autumn
#######################################################################

## total abundancance and diversity of origins across samples
# get the total abundance of OTUs across seasons, system and type
infoO<-info
infoO$depth<-factor(ifelse(infoO$depth %in% c(-1, -10), "Surface", "Bottom" ), levels = c("Surface", "Bottom"), ordered = TRUE)
list_mab_ab<-list()
list_mab_per<-list()
for (i in c("wadden sea", "north sea")){
  for (j in c("winter", "spring", "summer","autumn")){
    for (z in  c("sediment", "watercolumn")){
      if(nrow(info[infoO$location == i & infoO$season == j & infoO$type == z ,])==0) next
      if(i == "north sea" & j == "summer" & z == "watercolumn"){
        list_mab_ab[[paste(i,j,z,"Surface", sep = "/")]]<-rowSums(asv[,colnames(asv) %in% infoO[infoO$location == i & infoO$season == j & infoO$type == z & infoO$depth == "Surface" , "sampleID"] ])
        list_mab_per[[paste(i,j,z,"Surface", sep = "/")]]<-rowSums(asv[,colnames(asv) %in% infoO[infoO$location == i & infoO$season == j & infoO$type == z & infoO$depth == "Surface" , "sampleID"] ])/sum(rowSums(asv[,colnames(asv) %in% infoO[infoO$location == i & infoO$season == j & infoO$type == z & infoO$depth == "Surface" , "sampleID"] ]))
        list_mab_ab[[paste(i,j,z,"Bottom", sep = "/")]]<-rowSums(asv[,colnames(asv) %in% infoO[infoO$location == i & infoO$season == j & infoO$type == z & infoO$depth == "Bottom" , "sampleID"] ])
        list_mab_per[[paste(i,j,z,"Bottom", sep = "/")]]<-rowSums(asv[,colnames(asv) %in% infoO[infoO$location == i & infoO$season == j & infoO$type == z & infoO$depth == "Bottom" , "sampleID"] ])/sum(rowSums(asv[,colnames(asv) %in% infoO[infoO$location == i & infoO$season == j & infoO$type == z & infoO$depth == "Bottom" , "sampleID"] ]))
      }else{
      list_mab_ab[[paste(i,j,z, sep = "/")]]<-rowSums(asv[,colnames(asv) %in% info[info$location == i & info$season == j & info$type == z , "sampleID"] ])
      list_mab_per[[paste(i,j,z, sep = "/")]]<-rowSums(asv[,colnames(asv) %in% info[info$location == i & info$season == j & info$type == z , "sampleID"] ])/sum(rowSums(asv[,colnames(asv) %in% info[info$location == i & info$season == j & info$type == z , "sampleID"] ]))
    }}}} 
list_mab<-as.data.frame(bind_cols(list_mab_per));rownames(list_mab)<-rownames(asv)
list_mab$origin<-colnames(list_mab)[apply(list_mab, 1,  which.max)]
list_mab$asv<-rownames(list_mab)
#saveRDS(list_mab, "/export/lv4/projects/NIOZ153_309/Data/Operational/origin_asv.RData")

# ASV origin dataset
asv.ori<-melt(as.matrix(asv[,colnames(asv) %in% info[info$type %in% c("watercolumn", "sediment"), "sampleID"] ]))
asv.ori<-asv.ori[asv.ori$value>0,];colnames(asv.ori)<-c("asv", "sample", "reads")
asv.ori<-merge(asv.ori, list_mab[, c("asv", "origin")], by ="asv")
asv.ori$div<-1
  
# merge with metadata
ori.data<-aggregate(cbind(reads,div) ~ origin+sample ,data= asv.ori,FUN = sum )
ori.data<-merge(ori.data, infoO, by.x = "sample", by.y = "sampleID")
for (i in 1:nrow(ori.data)){ori.data$perc[i]<-100*(ori.data$reads[i]/sum(ori.data[ori.data$sample %in% ori.data$sample[i],"reads"])) }
#ori.data<-melt(origin.ws, id.vars = c("tide","time", "origin"))
#ori.data$origin<-factor(ori.data$origin, levels = c("Water Column", "Sediment","North Sea", "Ubiquitous"), ordered = TRUE)
ori.data$ori.loc<-str_to_title(sub("\\/.*","",ori.data$origin))
ori.data$ori.type<-str_to_title(str_split_fixed(ori.data$origin, "/", 4)[,3] )
ori.data$ori.seas<-str_to_title(str_split_fixed(ori.data$origin, "/", 4)[,2])
ori.data$ori.loc.type<-paste(ori.data$ori.loc, ori.data$ori.type, sep = " / ")

# mise en forme
unique(ori.data$origin)
ori.data$origin<-factor(str_to_title(ori.data$origin),levels =  str_to_title(c("north sea/summer/watercolumn/Bottom","north sea/summer/watercolumn/Surface", "north sea/summer/sediment",
                                           "north sea/autumn/watercolumn", "north sea/autumn/sediment",
                                           "wadden sea/winter/watercolumn","wadden sea/winter/sediment",
                                           "wadden sea/spring/watercolumn","wadden sea/spring/sediment",
                                           "wadden sea/summer/watercolumn","wadden sea/summer/sediment",
                                           "wadden sea/autumn/watercolumn","wadden sea/autumn/sediment")), ordered = TRUE)
pal=c("#6e3232", "#bd5757", "#4d2323","#9863c9","#37234a","#3b79f5","#1d3c7a", "#e5f53b","#727a1d", "#3bed97","#1d7a4d", "#ed8139","#7a421d")
names(pal)<-levels(ori.data$origin)

# aggregation
ori.data.sp<-aggregate(cbind(div, reads, perc)~time+type+depth+location+season+ori.loc.type,ori.data, FUN = sum)
ori.data$location<-factor(ori.data$location,levels = c("north sea", "wadden sea"), labels =c("NS", "WS") , ordered = TRUE )
ori.data$season<-factor(ori.data$season,levels = c("winter", "spring", "summer","autumn"), labels =c("Winter", "Spring", "Summer", "Autumn") , ordered = TRUE )

for (i in 1:nrow(ori.data)){ori.data$per[i]<-(ori.data$reads[i]/sum(ori.data[ori.data$sample %in% ori.data$sample[i],"reads"  ]))*100  }

surf<-ori.data[ori.data$DEPTH %in% c("surface"),]
ag<-aggregate( per~origin, surf[surf$location==  "WS" & surf$season==  "Winter",], mean)
sum(ag[grep("Wadden",ag$origin),"per"])
sum(ag[grep("North",ag$origin),"per"])
sum(ag[grep("Sediment",ag$origin),"per"])
sum(ag[grep("Watercolumn",ag$origin),"per"])

ori<-ggplot(data = ori.data[!ori.data$DEPTH %in% c("sediment"),], 
       aes(x = time, y = per/100))+
  facet_grid(depth~location+season, scales = "free", space = "free")+
  #facet_grid(location ~season, scales = "free")+
  geom_area(aes(fill = origin),stat = "identity", position = "stack")+
  scale_fill_manual(name = "Origin",values= pal)+
  #geom_bar(aes(fill =ori.loc.type ),stat = "identity", position = "stack")+
  #scale_fill_manual(name = "Origin",values= c("#4d4b23", "#ebe417", "#4d2323","#f22929"))+
  #geom_line(aes(color =origin), size = 1)+
  #scale_color_manual(name = "Origin",values= c("#4d4b23", "#ebe417", "#4d2323","#f22929"))+
  ylab("Proportions")+ 
  scale_x_datetime(date_breaks = "6 hour",labels = scales::date_format("%H:%M"), expand = c(0,0) )+
  #coord_flip()+
  scale_y_continuous(expand =c(0,0), breaks = c(0,0.5,1), labels = c("0", "0.5", "1"))+
  theme(panel.spacing.y = unit(1.25, "lines"))+
  theme(axis.line.x = element_line(colour = 'gray15', linewidth=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', linewidth=1, lineend = "square"))+
  theme(axis.title = element_text(size = 18,hjust = 0.5))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 10, angle = 60, hjust = 1, vjust = 1))+
  theme(axis.text.y = element_text(size = 14))+
  theme(strip.text =  element_text(size = 16, vjust = 0))+
  theme(strip.text.y =  element_text(size = 16, angle = 0, vjust = 0.5, hjust = 0))+
  theme(panel.spacing = unit(0.5, "lines"))+
  theme(strip.background = element_blank())+
  theme(legend.text = element_text(size = 14))+
  theme(legend.title = element_text(size = 14, face = "bold"))+
  theme(legend.key =  element_rect(fill = "#00000000"))+
  theme(panel.grid = element_blank())+
  theme(panel.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0),
        plot.subtitle = element_text(size = 16, face = "plain",hjust = 0) ); ori

##
# Figure 6
##

# 20 x 6
(ggnet+ori)+plot_annotation(tag_levels="A") & theme(plot.tag = element_text(size = 16, face = "plain"))

# Further Exploration
#######################################################################

## total abundance and diversity of origins across samples
# Donuts total reads
list_ab<-data.frame(asv = rownames(asv), tot= rowSums(bind_cols(list_mab_ab)), origin = list_mab$origin )
list_ab$div<-1
ori.tot<-aggregate(cbind(tot, div)~ origin,data = list_ab, FUN = sum)
ori.tot$perc<-(ori.tot$tot/sum(ori.tot$tot))*100
ori.tot$origin<-factor(str_to_title(ori.tot$origin), str_to_title(c("north sea/summer/watercolumn/Bottom","north sea/summer/watercolumn/Surface", "north sea/summer/sediment",
                                         "north sea/autumn/watercolumn", "north sea/autumn/sediment",
                                         "wadden sea/winter/watercolumn","wadden sea/winter/sediment",
                                         "wadden sea/spring/watercolumn","wadden sea/spring/sediment",
                                         "wadden sea/summer/watercolumn","wadden sea/summer/sediment",
                                         "wadden sea/autumn/watercolumn","wadden sea/autumn/sediment")), ordered = TRUE)
pal2=c("#6e3232", "#bd5757", "#4d2323","#9863c9","#37234a","#3b79f5","#1d3c7a", "#e5f53b","#727a1d", "#3bed97","#1d7a4d", "#ed8139","#7a421d")
names(pal2)<-levels(ori.tot$origin)

# format for donuts aggregate per season and depth
ori.tot<-ori.tot[order(ori.tot$origin),]
ori.tot$ymax.p<-cumsum(ori.tot$perc)
ori.tot$ymin.p<-c(0, head( ori.tot$ymax.p,n=-1))
ori.tot$ymax.d<-cumsum(ori.tot$div)
ori.tot$ymin.d<-c(0, head( ori.tot$ymax.d,n=-1))
ori.tot$label.p <- (ori.tot$ymax.p + ori.tot$ymin.p) / 2
ori.tot$label.d <- (ori.tot$ymax.d + ori.tot$ymin.d) / 2

# donuts reads
gdr<-ggplot(ori.tot, aes(ymax=ymax.p, ymin=ymin.p, xmax=4, xmin=3, fill=origin)) +
  geom_rect() +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  scale_fill_manual(name = "Origin",values = pal2)+
  #ggtitle(paste("",lev_tax, " Level", sep = ""))+
  geom_text( x=3.5, aes(y=label.p, label=paste(round(perc,1),"%")), size=4, fontface = "bold", ) +
  geom_text( x=2, y=2, label= "% reads",size=6, fontface = "bold" ) +
  theme(axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank())+
  theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5))+
  theme(axis.title= element_blank())+
  theme(axis.text = element_blank())+
  theme(strip.text =  element_text(size = 16, face = "bold"))+
  theme(strip.text.y =  element_text(angle = 0))+
  theme(strip.background = element_blank())+
  theme(legend.justification=c(0, 1))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(aspect.ratio = 1)+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0))+
  guides(fill = guide_legend(ncol = 1));gdr

# donuts div
gdd<-ggplot(ori.tot, aes(ymax=ymax.d, ymin=ymin.d, xmax=4, xmin=3, fill=origin)) +
  geom_rect() +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  scale_fill_manual(name = "Origin",values = pal2)+
  #ggtitle(paste("",lev_tax, " Level", sep = ""))+
  geom_text( x=3.5, aes(y=label.d, label=round(div,0)), size=4, fontface = "bold" ) +
  geom_text( x=2, y=2, label= "# asv",size=6, fontface = "bold" ) +
  theme(axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank())+
  theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5))+
  theme(axis.title= element_blank())+
  theme(axis.text = element_blank())+
  theme(strip.text =  element_text(size = 16, face = "bold"))+
  theme(strip.text.y =  element_text(angle = 0))+
  theme(strip.background = element_blank())+
  theme(legend.justification=c(0, 1))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(aspect.ratio = 1)+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0))+
  guides(fill = guide_legend(ncol = 1));gdd

donuts<-(gdd/gdr) +  plot_layout(guides = 'collect')
donuts+ori+plot_annotation(tag_levels="A") & theme(plot.tag = element_text(size = 16, face = "plain"))

## test %variance of %Origins explained by env. variables in each ecosystem (surface watercolumn) 
# subset of data
surf.ori<-ori.data[ori.data$DEPTH == "surface",]
surf.ori<-merge(surf.ori, info.w, by.x = "sample",  by="sampleID")
# loop anova across ecosystems
library("relaimpo")
list_anova<-list()
for (i in c("WS", "NS")){
  for (j in c("Winter", "Spring", "Summer","Autumn")){
    # make subset of data
    sub<-surf.ori[surf.ori$location %in% c(i) & surf.ori$season %in% c(j),]
    if(nrow(sub)==0) next
    sub.d<-dcast(sample+tide+temp+sal+ox+Winddirection+Windspeed.avgph~ origin, data = sub, value.var = "perc", fill = 0)
    
    # loop all origin perc in ecosystem 
    response_vars<-c(levels(sub$origin)) # 
    res<-rbind()
    for (rvar in response_vars) {
      # linear model
      if(is.null(sub.d[[rvar]])) next
      lm_fit<-lm(sub.d[[rvar]] ~ temp+sal+tide+Winddirection+Windspeed.avgph , data=sub.d)
      # anova of the linea model
      anova_lm_fit <- anova(lm_fit)
      anova_lm_fit_Estimate_ss<- anova_lm_fit$"Sum Sq"
      # save output of anova
      df_anova_test<-cbind(anova_lm_fit,PctExp=(anova_lm_fit_Estimate_ss/sum(anova_lm_fit_Estimate_ss)*100))
      df_anova_test$R_name<-rownames(df_anova_test)
      # relative importance of each variable
      car_relimpo_metrics <- calc.relimp(lm_fit, type = c("car"))
      car_relimpo_metricsDF<-data.frame(cbind(car_relimpo_metrics$car,rownames(car_relimpo_metrics)))
      car_relimpo_metricsDF$R_name<-rownames(car_relimpo_metricsDF)
      colnames(car_relimpo_metricsDF)<-c("car","R_name")
      car_relimpo_metricsDF$ori<-rvar
      res<-rbind(res, car_relimpo_metricsDF)
    }
    list_anova[[paste(i,j, sep = "/")]]<-res
    }
  }
list_anova<-bind_rows(list_anova, .id = "eco")      
list_anova$loc<-factor(sub("\\/.*","",list_anova$eco), levels = c("NS", "WS"),labels = c("North Sea", "Wadden Sea"), ordered = TRUE)
list_anova$season<-factor(str_to_title(str_extract(list_anova$eco, "\\w+$")),  levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)
list_anova$ori<-factor(list_anova$ori, levels = names(pal), ordered = TRUE)
list_anova$R_name<-factor(list_anova$R_name, 
                       levels = (c("temp", "sal", "tide", "Winddirection", "Windspeed.avgph")),
                       labels = (c("Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed")), ordered = TRUE)

ggplot(data = list_anova,aes(x =ori, y=car*100 , fill = R_name))+
  facet_wrap(loc ~season, scales = "free")+
  geom_bar(position = "stack",stat ="identity")+
  scale_fill_manual(name = "Env. Variable",values= c("#f0624f", "#913e33", "#224282","#71bf9b","#46705d"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,100))+
  xlab("Origin Category")+ylab(expression ("R"^2))+ 
  theme(axis.line.x = element_line(colour = 'gray15', linewidth=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', linewidth=1, lineend = "square"))+
  theme(axis.title = element_text(size = 18,hjust = 0.5))+
  theme(axis.text = element_text(size = 16))+
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1, color = pal))+
  theme(panel.spacing = unit(0.5, "lines"))+
  theme(strip.background = element_blank())+
  theme(strip.text =  element_text(size = 16, face = "bold"))+
  theme(strip.text.y  =  element_text(angle = 0, hjust = 0))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(legend.key =  element_rect(fill = "#00000000"))+
  theme(panel.grid.major= element_line(colour = "white", linewidth = 1))+
  theme(panel.grid.minor= element_line(colour = "white", linewidth = 0.5))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.grid.minor = element_line(color = "#ffffff", linewidth = 0.5))+
  theme(panel.grid.major = element_line(color = "#ffffff", linewidth = 1))+
  theme(panel.background = element_rect(fill = "#ebebfa",colour = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0),
        plot.subtitle = element_text(size = 16, face = "plain",hjust = 0) )

## explore taxonomic diversity of the origin categories
ori.tax<-merge(list_ab, taxo, by.x = "asv", by.y = "row.names")
ori.tax$origin<-str_to_title(ori.tax$origin)
tab.ori.tax<-dcast(origin ~ Phylum,data = ori.tax, value.var = "div", fill = 0);rownames(tab.ori.tax)<-tab.ori.tax[,1]; tab.ori.tax<-tab.ori.tax[,-1]

hc=hclust(vegdist(tab.ori.tax, "bray"),method="average")
hc=reorder(hc,wts=-as.matrix(tab.ori.tax)%*%seq(ncol(tab.ori.tax))^2) # vegan::reorder.hclust
tree=ggdendro::dendro_data(as.dendrogram(hc),type="rectangle")

p1=ggplot(ggdendro::segment(tree))+
  geom_segment(aes(x=y,y=x,xend=yend,yend=xend),lineend="round",size=.4)+
  geom_point(data = tree$labels,aes(x=y,y=x, color = label), size = 4 )+
  scale_color_manual(values = pal, guide = "none")+
  #scale_x_continuous(expand=expansion(add=c(0,.01)))+ # don't crop half of line between top-level nodes
  #scale_y_continuous(limits=.5+c(0,nrow(t)),expand=c(0,0))+
  xlab("Bray-Curtis Distance")+
  theme(
    axis.text.y=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_text(size = 11),
    axis.title.x=element_text(size = 11),
    axis.ticks.y=element_blank(),
    #axis.ticks.length=unit(0,"pt"), # remove extra space occupied by ticks
    panel.background=element_rect(fill="white"),
    panel.grid=element_blank(),
    plot.margin=margin(5,5,5,0)
  )

tab.ori.tax=tab.ori.tax[hc$labels[hc$order],]
t2=data.frame(V1=rownames(tab.ori.tax)[row(tab.ori.tax)],V2=colnames(tab.ori.tax)[col(tab.ori.tax)],V3=unname(do.call(c,tab.ori.tax)))
lab=round(100*t2$V3)
lab[lab==0]=""

p2=ggplot(t2,aes(x=factor(V1,level=rownames(tab.ori.tax)),y=V3,fill=V2))+
  #geom_bar(stat="identity",width=1,position=position_fill(reverse=T))+
  geom_bar(stat="identity",width=1,position="stack")+
  #geom_text(aes(label=lab),position=position_stack(vjust=.5,reverse=T),size=3.5)+
  coord_flip()+
  ylab("ASV #")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), limits = c(0,4500 ))+
  scale_fill_manual(name = "Phylum",values= pal.t)+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text=element_text(size=11),
    axis.text.y=element_text(size=11, color = "black"),
    axis.title=element_text(size=11),
    axis.title.y=element_blank(),
    #axis.text.x=element_blank(),
    legend.position = "bottom",
    plot.margin=margin(5,0,5,5)
  )
(gg.taxo=((p2+p1)+  plot_layout(widths = c(3, 1))))

## explore functional diversity of the origin categories
ori.func<-merge(list_ab, fapro, by.x = "asv", by.y = "row.names")
ori.func$origin<-str_to_title(ori.func$origin)
m.ori.func<-melt(ori.func, id.vars = c("origin", "asv", "div", "tot"))
m.ori.func<-m.ori.func[m.ori.func$value>0,]
ori.func.ag<-aggregate(value~origin+variable ,m.ori.func, FUN = sum)
for (i in 1:nrow(ori.func.ag)){ori.func.ag$perc[i]<-(ori.func.ag$value[i]/sum(ori.func.ag[ori.func.ag$origin %in% ori.func.ag$origin[i], "value" ]))*100 }
cont.func<-aggregate(value ~ variable,ori.func.ag, FUN = sum)
cont.func$perc.f<-(cont.func$value/sum(cont.func$value))*100
ori.func.ag<-merge(ori.func.ag,cont.func[, c(1,3)], by = "variable" )

tab.ori<-dcast(origin ~variable,data = ori.func.ag, value.var = "perc",fun.aggregate = sum);rownames(tab.ori)<-tab.ori[,1]; tab.ori<-tab.ori[,-1]

hc1=hclust(vegdist(tab.ori, "bray"),method="average")
hc1=reorder(hc1,wts=-as.matrix(tab.ori)%*%seq(ncol(tab.ori))^2) # vegan::reorder.hclust
tree1=ggdendro::dendro_data(as.dendrogram(hc1),type="rectangle")

p3=ggplot(ggdendro::segment(tree1))+
  geom_segment(aes(x=y,y=x,xend=yend,yend=xend),lineend="round",size=.4)+
  geom_point(data = tree1$labels,aes(x=y,y=x, color = label), size = 4 )+
  scale_color_manual(values = pal, guide = "none")+
  #scale_x_continuous(expand=expansion(add=c(0,.01)))+ # don't crop half of line between top-level nodes
  #scale_y_continuous(limits=.5+c(0,nrow(t)),expand=c(0,0))+
  xlab("Bray-Curtis Distance")+
  theme(
    axis.text.y=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_text(size = 11),
    axis.title.x=element_text(size = 11),
    axis.ticks.y=element_blank(),
    #axis.ticks.length=unit(0,"pt"), # remove extra space occupied by ticks
    panel.background=element_rect(fill="white"),
    panel.grid=element_blank(),
    plot.margin=margin(5,5,5,0)
  )

ori.func.ag$origin=factor(ori.func.ag$origin, levels =  hc1$labels[hc1$order], ordered = TRUE)
ori.func.ag$variable<-factor(ori.func.ag$variable, levels = cont.func[order(cont.func$perc.f, decreasing = TRUE),"variable"],ordered = TRUE)
ori.func.ag$percc<-factor(cut(ori.func.ag$perc, breaks = c(-Inf,1,2, 5, 20, 40), labels = c("< 1 %", "1-2 %","2-5 %", "5-20 %", "20-40 %") ), ordered = TRUE)

str(ori.func.ag)
p4=ggplot(ori.func.ag,aes(x=origin,y=variable,fill=percc))+
  #geom_bar(stat="identity",width=1,position=position_fill(reverse=T))+
  geom_tile(width=1)+
  coord_flip()+
  ylab("Functions (%)")+
  #scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0), position = "right")+
  scale_fill_brewer(palette="YlOrRd", name = '% per dataset')+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(size=11, angle = 90, hjust = 0, vjust = 0),
    axis.text.y=element_text(size=11, color = "black"),
    axis.title=element_text(size=11),
    axis.title.y=element_blank(),
    #axis.text.x=element_blank(),
    legend.position = "bottom",
    plot.margin=margin(5,0,5,5)
  )
(gg.func=((p4+p3)+  plot_layout(widths = c(3, 1))))

(gg.func/gg.taxo)

