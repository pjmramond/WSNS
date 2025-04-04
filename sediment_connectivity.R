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
library(biom)
library(patchwork)


overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# vector with all the colors available in R, except gray
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

# location of the script
rstudioapi::getSourceEditorContext()$path
options(getClass.msg=FALSE)

###############################
# Data import
###############################

# import data filtered by Julia
asv<-as.data.frame(readRDS("Desktop/WSNS/asv3_minSkin_20221115.rds"))

# extract species-taxonomy matrix
taxo<-as.data.frame(observation_metadata(read_biom("Desktop/WSNS/asvTable.biom")))
colnames(taxo)<-c("Domain", "Phylum","Class", "Order", "Family", "Genus", "Species")
taxo[taxo == "NA"]<-NA 
taxo[is.na(taxo$Domain),"Domain"]<-"Unclassified"
taxo$level<-colnames(taxo)[apply(taxo, 1, function(x) length(which(!is.na(x))))] # get taxonomic level at which each sequence is annotated
taxo[grep("Proteobacteria", taxo[,2]),2]<-taxo[grep("Proteobacteria", taxo[,2]),3]
taxo<-taxo[rownames(taxo) %in% rownames(asv),]
taxo$asv<-rownames(taxo)

# metadata file
info<-read.csv('Desktop/WSNS/metadata_ws_ns.csv', sep = ';')
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
tide<-read.csv('Desktop/WSNS/tides.csv', sep = ';')
colnames(tide)[6]<-"tide"
tide$time<-as.POSIXct(paste(tide$date, tide$UTC_time), format="%d/%m/%Y %H:%M:%S", tz = "CET")
info<-merge(info, tide[, c("tide", "time")], by = "time", all.x = TRUE)

## wind information and merging with data
# import and format
wind<-read.csv('Desktop/WSNS/wind_wsns.csv', sep = ';')
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

###############################
# Connectivity
###############################

## Data formating
# order/set dataset
info3<-info[info$type %in% c("sediment","watercolumn"),]
asv3<-asv[,colnames(asv) %in% info3$sampleID]
tasv<-t(asv3)
tasv<-tasv[order(rownames(tasv)),] # OG dataset
info4<-info3[info3$sampleID %in% rownames(tasv),]
info4<-info4[order(info4$sampleID),]
info4$Depth<-factor(ifelse(info4$depth %in% c(-1, -10), "Surface", "Bottom" ), levels = c("Surface", "Bottom"), ordered = TRUE)
info4$id<-paste(str_to_title(info4[,4]),str_to_title(info4[,5]), sep = '/')
info4$depth_sed<-info4$sampleID
info4[info4$type == "sediment", "depth_sed"]<-paste("Sed", info4[info4$type == "sediment", "id"],sep = "/")

# we merge sediments community per season and location (4 depth sampled)
cum_asv<-aggregate(tasv, by = list(info4$depth_sed), FUN = sum) ; rownames(cum_asv)<-cum_asv[,1] ;cum_asv<-cum_asv[,-1] # takes 1m

## create data for networking
# Number of shared OTUs across ecosystems / links/edges
beta<-as.matrix(betadiver(cum_asv)$a);beta[upper.tri(beta)] <- NA
beta<-reshape2::melt(as.matrix(beta))
beta<-na.omit(beta);beta<-beta[beta$value >0,]
link<-as_tibble(beta);colnames(link)<-c("from", "to","weight")
link<-merge(link, info4[,c( 2,4:6,27,1)], by.x = "from", by.y = "sampleID", all.x = TRUE)
link<-merge(link, info4[,c(2,4:6,27,1)], by.x = "to", by.y = "sampleID", all.x = TRUE)
link[is.na(link$location.x), "from"]

# format data for plotting
link[is.na(link$location.x),"location.x"]<-tolower(gsub(".*/","",sub("/[^/]+$", "", link[is.na(link$location.x),"from"])))
link[is.na(link$location.y),"location.y"]<-tolower(gsub(".*/","",sub("/[^/]+$", "", link[is.na(link$location.y),"to"])))
link[is.na(link$season.x),"season.x"]<-tolower(gsub(".*/","",link[is.na(link$season.x),"from"]))
link[is.na(link$season.y),"season.y"]<-tolower(gsub(".*/","",link[is.na(link$season.y),"to"]))
link[is.na(link$type.x),"type.x"]<-"sediment"
link[is.na(link$type.y),"type.y"]<-"sediment"
link<-link[link$season.x == link$season.y & link$location.x == link$location.y & link$type.x != link$type.y,]
link$season.y<-factor(link$season.y, levels = c("winter", 'spring',"summer", "autumn"), labels = c("Winter", 'Spring',"Summer", "Autumn")  , ordered = TRUE)
link<-merge(link, info2, by.x = "to", by.y = "sampleID")

# 12 X 5
aggregate(weight~Depth.y+season.y+location.x , data = link, mean)

gbws<-ggplot(link, aes(x = time.y, y= weight, color = Depth.y))+
  facet_grid(.~str_to_title(location.x)+season.y, scale = "free_x")+
  geom_line(size= 1)+
  ggtitle("")+
  xlab("Time (h)")+ylab("Water-column ASVs shared\nwith the sediments")+
  labs(color = "Depth")+
  scale_color_manual(values = c("#7aadff","#063175"))+
  scale_x_datetime(date_breaks = "6 hour",labels = scales::date_format("%H:%M"))+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 18,hjust = 0.5))+
  theme(axis.text = element_text(size = 16))+
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1))+
  theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5))+
  theme(panel.spacing = unit(0.5, "lines"))+
  theme(strip.background = element_blank())+
  theme(strip.text =  element_blank())+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(legend.key =  element_rect(fill = "#00000000"))+
  theme(panel.grid.major= element_line(colour = "white", size = 1))+
  theme(panel.grid.minor= element_line(colour = "white", size = 0.5))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.grid.minor = element_line(color = "#ffffff", size = 0.5))+
  theme(panel.grid.major = element_line(color = "#ffffff", size = 1))+
  theme(panel.background = element_rect(fill = "#ebebfa",colour = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gbws

gwst<-ggplot(link[link$Depth.y == "Surface",], aes(x = time.y, y= tide))+
  facet_grid(.~str_to_title(location.x)+season.y, scale = "free_x")+
  geom_line(size= 1, color = "#224282")+
  ggtitle("")+
  xlab("Time (h)")+ylab("Water\nHeight (m)")+
  labs(color = "Depth")+
  scale_x_datetime(date_breaks = "6 hour",labels = scales::date_format("%H:%M"))+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 18,hjust = 0.5))+
  theme(axis.text = element_text(size = 16))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5))+
  theme(panel.spacing = unit(0.5, "lines"))+
  theme(strip.background = element_blank())+
  theme(strip.text =  element_text(size = 16, face = "bold"))+
  theme(strip.text.y  =  element_blank())+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(legend.key =  element_rect(fill = "#00000000"))+
  theme(panel.grid.major= element_line(colour = "white", size = 1))+
  theme(panel.grid.minor= element_line(colour = "white", size = 0.5))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.grid.minor = element_line(color = "#ffffff", size = 0.5))+
  theme(panel.grid.major = element_line(color = "#ffffff", size = 1))+
  theme(panel.background = element_rect(fill = "#ebebfa",colour = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gwst

gwlight<-ggplot(link[link$Depth.y == "Surface",], aes(x = time.y, y= light.ac))+
  facet_grid(.~str_to_title(location.x)+season.y, scale = "free_x")+
  geom_line(size= 1, color = "gold3")+
  ggtitle("")+
  xlab("Time (h)")+ylab("Accumulated\nDaylight (h)")+
  labs(color = "Depth")+
  scale_x_datetime(date_breaks = "6 hour",labels = scales::date_format("%H:%M"))+
  theme(axis.line.x = element_line(colour = 'gray15', size=1,lineend = "square"),
        axis.line.y = element_line(colour = 'gray15', size=1, lineend = "square"))+
  theme(axis.title = element_text(size = 18,hjust = 0.5))+
  theme(axis.text = element_text(size = 16))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(plot.title = element_text(size = 24, face = "bold",hjust = 0.5))+
  theme(panel.spacing = unit(0.5, "lines"))+
  theme(strip.background = element_blank())+
  theme(strip.text  =  element_blank())+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(legend.key =  element_rect(fill = "#00000000"))+
  theme(panel.grid.major= element_line(colour = "white", size = 1))+
  theme(panel.grid.minor= element_line(colour = "white", size = 0.5))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.grid.minor = element_line(color = "#ffffff", size = 0.5))+
  theme(panel.grid.major = element_line(color = "#ffffff", size = 1))+
  theme(panel.background = element_rect(fill = "#ebebfa",colour = NA))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));gwlight

##
# Figure S3
##

# 15 x 10
gwst/gwlight /gbws +plot_layout(heights = c(1, 1,4))

###################################
# Connectivity variable explaining
###################################

## data prep
ano<-link
ano$location<-factor(ano$location.x,levels = c("north sea", "wadden sea"), labels =c("North Sea", "Wadden Sea") , ordered = TRUE )
ano$season<-factor(ano$season.x,levels = c("winter", "spring", "summer","autumn"), labels =c("Winter", "Spring", "Summer", "Autumn") , ordered = TRUE )
ano$Depth<-factor(ano$DEPTH , levels = c("surface", 'depth'), labels = c("Surface", "Bottom") , ordered = TRUE)
for (i in grep("time.t.x|light.ac|temp|sal|tide|Winddirection|Windspeed.avgph",colnames(ano)) ){ano[,i]<-as.numeric(ano[,i])}

# run anova for all sub-conditions
list_anova<-list()
for (i in c("North Sea", "Wadden Sea")){
  for (j in c("Winter", "Spring", "Summer", "Autumn")){
    for (z in  c("Surface", "Bottom")){
      # linear model
      anov<-ano[ano$location %in% i & ano$season %in% j  & ano$Depth %in% z,]
      if(nrow(anov)==0) next
      res<-as.data.frame(matrix(ncol= 6, nrow = 7));colnames(res)<-c("var", "r2", "pv", "lag.r2", "lag", "lag.pv")
      res$var<-c("time.t.x","light.ac","temp","sal","tide","Winddirection", "Windspeed.avgph")
      for(y in res$var){
        res[res$var == y, "r2"]<-cor.test(anov$weight, anov[,y])$estimate
        res[res$var == y, "pv"]<-cor.test(anov$weight, anov[,y])$p.value
        corLag<-ccf(anov$weight,anov[,y], plot = FALSE)
        res[res$var == y, "lag.r2"]<-max(corLag$acf)
        res[res$var == y, "lag"]<-corLag$lag[which.max(corLag$acf)]
        res[res$var == y, "lag.pv"]<- 2*(1 - pnorm(max(abs(corLag$acf)), mean = 0, sd = 1/sqrt(corLag$n.used)))
      }
      list_anova[[paste(i,j,z, sep = "/")]]<-res
    }
  }
}

list_anova<-bind_rows(list_anova, .id = "eco")
list_anova$var<-factor(list_anova$var, 
                          levels = (c("time.t.x","light.ac","temp", "sal", "tide", "Winddirection", "Windspeed.avgph")),
                          labels = (c("Time","Accumulated Daylight","Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed")), ordered = TRUE)
list_anova$NS<-ifelse(list_anova$pv>0.1, "Non-significant", as.character(list_anova$var))
list_anova$NS<-factor(list_anova$NS, levels = c("Non-significant","Time","Accumulated Daylight","Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed"), ordered = TRUE)
list_anova$NS.lag<-ifelse(list_anova$lag.pv>0.1, "Non-significant", as.character(list_anova$var))
list_anova$NS.lag<-ifelse(abs(list_anova$lag)>1, "Non-significant", list_anova$NS.lag)
list_anova$NS.lag<-factor(list_anova$NS.lag, levels = c("Non-significant","Time","Accumulated Daylight","Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed"), ordered = TRUE)
list_anova$lag.lab<-ifelse(list_anova$NS.lag %in%  c("Non-significant") | list_anova$lag==0, "", list_anova$lag*4)
list_anova$eco<-factor(list_anova$eco, levels=rev(c("North Sea/Summer/Surface","North Sea/Summer/Bottom", "North Sea/Autumn/Surface", "North Sea/Autumn/Bottom","Wadden Sea/Winter/Surface","Wadden Sea/Winter/Bottom","Wadden Sea/Spring/Surface","Wadden Sea/Spring/Bottom","Wadden Sea/Summer/Surface","Wadden Sea/Summer/Bottom","Wadden Sea/Autumn/Surface","Wadden Sea/Autumn/Bottom")), ordered = TRUE)
#View(list_anova[list_anova$var == "Tides",])
#View(list_anova)

pal2=rev(c("#bd5757", "#4d2323","#9863c9","#37234a","#3b79f5","#1d3c7a", "#e5f53b","#727a1d", "#3bed97", "#ed8139","#7a421d"))
pal.env<-c("gray75","black","gold3","#f0624f", "#913e33", "#224282","#71bf9b","#46705d")
names(pal.env)<-c("Non-significant","Time","Accumulated Daylight","Temperature", "Salinity", "Tides", "Wind Direction", "Wind Speed")


##
# Figure S4
##
gano<-ggplot(data=list_anova,aes(y=eco, x = var,  color = NS.lag , size = lag.r2))+
  geom_point()+
  geom_text(aes(label = lag.lab), color = "gray1", hjust = 2, size = 4)+
  theme_void()+
  labs(size = expression(paste("R"^"2")))+
  scale_color_manual(values=pal.env, name = "Env. variable", guide = "none")+
  theme(strip.text =  element_text(size = 14))+
  theme(axis.text.y = element_text(size = 16, angle = 0, hjust = 1, vjust = 0.5, color =pal2))+
  theme(axis.text.x  = element_text(size = 16, angle = 60, hjust = 1, vjust = 1 ), 
        legend.text = element_text(size = 14), 
        legend.direction  = "vertical", 
        legend.title = element_text(size = 14, face="bold"));gano

