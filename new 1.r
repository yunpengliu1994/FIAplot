###################################################
####### match spname in TRY,FIA and phylogeny #####
###################################################
library(data.table)
Try.splist=as.data.frame(fread("TryAccSpecies.csv"))#TRY + nph18031
corrt.sp=function(dat,spcol){
		require(data.table)
		tpl.ac=unique(as.data.frame(fread("TPL14Ac.csv",header=T)))
		tpl.sy=unique(as.data.frame(fread("TPL14Sy.csv",header=T)))
		corrt.1=cbind(dat,Accepted_ID=tpl.ac[match(dat[,spcol],tpl.ac$Acname),"ID"])
		unmat=subset(corrt.1,is.na(corrt.1$Accepted_ID))
		corrt.2=cbind(unmat[,-dim(unmat)[2]],Accepted_ID=tpl.sy[match(unmat[,spcol],tpl.sy$Syname),"Accepted ID"])
		corrt.3=rbind(subset(corrt.1,!is.na(corrt.1$Accepted_ID)),corrt.2)
		corrt.4=cbind(corrt.3,tpl.ac[match(corrt.3$Accepted_ID,tpl.ac$ID),c("Family","Genus","Acname","Taxonomic status in TPL")])
		return(corrt.4)
	}
sp=read.csv("FIAsp_clean.csv")
sp=	corrt.sp(sp,"Spname")
Try.splist=corrt.sp(Try.splist,"AccSpeciesName")
re=cbind(sp,Try.splist[match(sp$Spname,Try.splist$AccSpeciesName),])
unmat=subset(re[,colnames(sp)],is.na(re$AccSpeciesName))
re2=cbind(unmat,Try.splist[match(unmat$Accepted_ID,Try.splist$Accepted_ID),])
re3=rbind(re[!is.na(re$AccSpeciesID),],re2)
write.csv(re3[,c(1:13)],"TRY.splist.csv") #link spname in TRY and in FIA
write.csv(Try.splist,"TryAccSpecies.csv")
#correct spnames in GRooT data base
Try.splist=read.csv("GRooTAggregateSpeciesVersion.csv")
Try.splist=corrt.sp(sp,"Spname")
Try.splist[is.na(Try.splist$Acname),"Acname"]=Try.splist[is.na(Try.splist$Acname),"Spname"]
write.csv(Try.splist,"GRooTAggregateSpeciesVersion.csv")
#match spnames in FIA and in phylogeny
library(ape)
tre=read.tree("BigTreeSmithBrown2018v0.1/ALLOTB.tre")#356305 tips
tre.list=data.frame(tip1=tre$tip.label,tip2=gsub("_", " ", tre$tip.label))
tre.spname=corrt.sp(tre.list,"tip2")
re=cbind(sp,tre.spname[match(sp$Spname,tre.spname$tip2),])
unmat=subset(re[,colnames(sp)],is.na(re$tip2))
re2=cbind(unmat,tre.spname[match(unmat$Accepted_ID,tre.spname$Accepted_ID),])
re3=rbind(subset(re,!is.na(re$tip2)),re2)
write.csv(re3[,c(1:13)],"Tree.splist.csv")#link spname in Trede and in FIA
write.csv(tre.spname,"TreeAccSpecies.csv")

##########################################################################
####### extract community matrix from FIA and select traits in TRY  ######
##########################################################################
get.fia=function(statename,try.sp){
	require(rFIA)
	require(magrittr)
	require(dplyr)
	state <-  readFIA(paste('FIA/',statename,sep=""),tables=c('PLOT', 'TREE','SURVEY','COND')) # state is the FIA database for that state  #example uses MASS
	##  create pltID which can uniquely identify a plot
	## STATECD: State code NUMBER; UNITCD: Survey unit code; COUNTYCD: County code; PLOT: Plot number
	state$TREE$pltID <-  stringr::str_c(state$TREE$UNITCD, state$TREE$STATECD, state$TREE$COUNTYCD, state$TREE$PLOT,state$TREE$INVYR,sep = '_')
	## #### APPLYING FILTERS ####
	## data frames	
	tree_df <- state$TREE %>%
		# records for alive trees only: STATUSCD == 1
		# A tree that has TPA_UNADJ = NA does not contribute to estimates of tree density, abundance, etc. Therefore, these trees should not be included in any diversity index calculations.
		# some NA TPA values  are due to “legacy trees” – those are trees that were measured as part of the older (pre-2000) inventories, which are still monitored to estimate mortality and growth etc. at the population level. Legacy trees have a Y in P2A_GRM_FLG in the TREE table	
		# remove data missing INVYR
		filter(STATUSCD == 1&INVYR>1995&INVYR!=9999&!is.na(TPA_UNADJ))%>% 
		# UNITCD, STATECD,COUNTYCD,PLOT and INVYR will identify a plot record. Most plot records represent measurements. But a few do not. PLOT_STATUS_CD = 1,  will select plot records that correspond to measurements.	
		# ECOSUBCD: classifying the plot by ecoregions
		left_join(select(state$PLOT,MEASYEAR,MEASMON,LAT,LON,ELEV,CN,PLOT_STATUS_CD,SRV_CN,ECOSUBCD), by=c("PLT_CN"="CN")) %>% filter(PLOT_STATUS_CD == 1) %>%
		# select only the data from annual inventories use the SURVEY table: ANN_INVENTORY == Y	
		# the survey is not for a P3 ozone inventory
		left_join(select(state$SURVEY,CN,ANN_INVENTORY,P3_OZONE_IND), by=c("SRV_CN"="CN")) %>% filter(ANN_INVENTORY == "Y"& P3_OZONE_IND == "N") %>% 
		# remove cases of NA's STDAGE columns;STDAGE: continuous variable derived from estimates of tree age from tree cores of two adult trees from the dominant size class of the forest plot 
		# Filter out plantations (by set STDORGCD==0); 
		# Filter out STDSZCD 5 (non-stocked plots)
		# COND_STATUS_CD == 1 is accessible forest land
		# Filter plots with a single condition using a 95% threshold for CONDPROP_UNADJ
		# PHYSCLCD:  coarsely classify plots as to whether they are on xeric (dry), mesic (moist), and hydric (water-logged) soils. 
		left_join(select(state$COND,PLT_CN,STDAGE,STDORGCD,STDSZCD,COND_STATUS_CD,CONDPROP_UNADJ,PHYSCLCD), by="PLT_CN") %>% 
		filter((!is.na(STDAGE))&(STDAGE!= 9999)&(STDAGE>0)&(STDORGCD == 0)&(STDSZCD!=5)&(COND_STATUS_CD == 1)&CONDPROP_UNADJ >= 0.95) %>%
		# match species in FIA with TRY: SPCD,Spname:FIA;AccSpeciesName:TRY
		left_join(select(try.sp,SPCD,Spname,AccSpeciesName),by="SPCD") %>% 
			select(pltID,SPCD,Spname,AccSpeciesName,HT,DIA,TPA_UNADJ,MEASYEAR,MEASMON,LAT,LON,ECOSUBCD,STDAGE,PHYSCLCD,DRYBIO_AG,DRYBIO_BG)
	tree_df$statename=statename
	return(tree_df)
	rm(state)
}
try.sp=read.csv("TRY.splist.csv") #link spname in TRY and in FIA
state.list=as.character(read.csv("state.list.csv")[,1])
require(parallel)
no_cores <- detectCores() - 1
mycl <- makePSOCKcluster(no_cores); 
fia=do.call(rbind,parLapply(cl=mycl,state.list,get.fia,try.sp))
stopCluster(mycl)
fia$pltID.no.year = substr(fia$pltID,1,nchar(fia$pltID)-5)
write.csv(fia,"fia.csv") 

## read the try data ---
require(data.table)
library(magrittr)
library(dplyr)
TRYdata1 <- fread("TRY/21232.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)
TRYdata2 <- fread("TRY/21233.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)
TRYdata3 <- fread("TRY/21236.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)
TRYdata4 <- fread("TRY/21336.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)
TRYdata5 <- fread("TRY/21337.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)
TRYdata6 <- fread("TRY/21338.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)
TRYdata=rbind(TRYdata1,TRYdata2,TRYdata3,TRYdata4,TRYdata5,TRYdata6)
#extract lon/lat for each record
lat = TRYdata %>% filter(DataID==59) %>% select(AccSpeciesID,AccSpeciesName,ObservationID,OrigValueStr,StdValue)%>%distinct()
lon = TRYdata %>% filter(DataID==60) %>% select(AccSpeciesID,AccSpeciesName,ObservationID,OrigValueStr,StdValue)%>%distinct()
coord = lat %>% left_join(lon, by=c("ObservationID","AccSpeciesID","AccSpeciesName"))
write.csv(coord,"TRYcoord.csv")#adj the orig value mannually, add nph18031 and then change the colname to LON and LAT

height=read.csv("fia.csv")[,-1]
trait=TRYdata %>% filter(TraitID %in% names(stat[stat>0.9])) %>% # 52 traits
	filter(SPCD %in% unique(height$SPCD)) %>% distinct()# 354 sp are found in FIA plots, and then correct the std value of each trait mannually
write.csv(trait,"TraitInTRY.csv")

trait.geo.needed=read.csv("Trait.list.final.csv") %>% select(TraitID,Abbr)
trait=as.data.frame(fread("TraitInTRY.clean.csv"))%>% left_join(trait.geo.needed,by="TraitID") %>% filter(!is.na(Abbr))
write.csv(trait,"TraitInTRY.clean.csv")#removing manually NAs and data with muti-standards， e.g., for seed dry mass, there are mean/max/min, only mean are kept, and also for Species tolerance to planting density, Root rooting depth and Leaf length, Leaf porosity(Summer/winter)
		#Then Updated root dep from Tumber-Dávila, S.J., et. al. (2022). New Phytol. and from GRooT database 
#Update WD using FIA’s REF_SPECIES table (the variable is called “WOOD_SPGR_GREENVOL_DRYWT”)
wd=read.csv("FIADB_REFERENCE/REF_SPECIES.csv")%>% select(SPCD,WOOD_SPGR_GREENVOL_DRYWT)%>%na.omit()
wd.fia=as.data.frame(fread("fia.ENA.csv"))[,-1]%>%filter(!is.na(Spname)&(TPA_UNADJ>0))%>%left_join(wd,by="SPCD")%>%
	select(SPCD,Spname,AccSpeciesName,WOOD_SPGR_GREENVOL_DRYWT)%>%distinct()
write.csv(wd.fia,"wd.fia.csv")#then combine this table to TraitInTRY.clean.csv

##########################################################################################
####### combine growth data with FIA plot data and obtain abundance and trait Metrics  ###
##########################################################################################
## annual growth woody biomass (1000 kg/ha/year) data from Aaron ----
load("Biomass Growth/G_Nov2022/G_Nov2022.Rdata")
head(G)# MEASTIME_avg: middle of t1 and t2; INVYR:t2; 
# B_L_prop: proportion of biomass loss due to tree harvest or mortality, plots with values over 0.5 may should be removed, or make the B_L_prop as random factor
# G_MassBal_MgHaYr: method 1; G_obs_TreeInc_MgHaYr: method 2
#change the INVYR (time 2) to INVYR (time 1)
#Filter out plantations
# Filter plots with a single condition using a 95% threshold for CONDPROP_UNADJ
# Filter out SITECLCD 7 (non productive stands - those with less than 20 ft3/ac/yr)
# Filter out non-stocked plots (with STDSZCD == 5)
# Filter out the the non-accessible plots -- COND_STATUS_CD == 1 is accessible
# Filter by survey table --- annual inventories YES
# Filter by survey table --- B3_OZONE NO
get.biomass=function(statename,FIA_G_calc){
	require(rFIA)
	require(magrittr)
	require(dplyr)
	state <-  readFIA(paste('FIA/',statename,sep=""),tables='PLOT') # state is the FIA database for that state  
	FIA_G_calc2=FIA_G_calc %>% filter(STATE==statename) %>%select(-INVYR)%>% left_join(select(state$PLOT,CN,INVYR),by=c("PREV_PLT_CN"="CN"))
	colnames(FIA_G_calc2)[colnames(FIA_G_calc2)%in%"pltID"]="pltID.no.year"
	FIA_G_calc2$pltID <-  stringr::str_c(FIA_G_calc2$pltID.no.year,FIA_G_calc2$INVYR,sep = '_')
	return(FIA_G_calc2)
}
state.list=unique(G$STATE)
require(parallel)
no_cores <- detectCores() - 1
mycl <- makePSOCKcluster(no_cores); 
FIA_G_calc.new=do.call(rbind,parLapply(cl=mycl,state.list,get.biomass,G))
stopCluster(mycl)	
save(FIA_G_calc.new,file="Biomass Growth/G_Nov2022/FIA_G_calc.new.Rdata")

#obtain FIA data in ENA with both G and phylogeny data---
ENA=c("ME","NH","VT","NY","MA","RI","CT","NJ","PA","DE","MD","MI","OH","IN","IL","WI","WV","VA","NC","TN","KY","SC","GA","AL","MS","FL")
fia0=as.data.frame(fread("fia.csv"))[,-1]
G=get(load("Biomass Growth/G_Nov2022/FIA_G_calc.new.Rdata"))
fia0$ECOPROVCD=substr(fia0$ECOSUBCD,1,nchar(fia0$ECOSUBCD)-2)%>%gsub("^\\s|\\s$", "", .)
tre.spname=na.omit(read.csv("Tree.splist.csv")[,c("Spname","tip1","tip2")])#name table link tips in phylogeny and spname in FIA
fia.ena=fia0 %>% left_join(select(G,pltID,pltID.no.year,nTrees_t1,B_plt_t1_MgHa,B_L_prop,G_MassBal_MgHaYr),by=c("pltID","pltID.no.year"))  %>% 
	filter(G_MassBal_MgHaYr>0&statename%in%ENA&G_MassBal_MgHaYr<15&B_L_prop<0.2&STDAGE<150&(!ECOPROVCD%in%c(234,411))&TPA_UNADJ>0)%>% 
	left_join(tre.spname,by="Spname") %>% filter(!is.na(Spname)&(!is.na(tip1)))
fia.ena$Abund_weight <-  (fia.ena$DIA*0.0254)^2*(fia.ena$TPA_UNADJ * 2.4710538147)   ##TPA base level TPA  in m2/ha	
write.csv(fia.ena,"fia.ENA.csv") #plot data used in further analysis

#obtain the phylogeny of FIA ENA----
library(ape);
tre=read.tree("BigTreeSmithBrown2018v0.1/ALLOTB.tre")#356305 tips
tree_df=as.data.frame(fread("fia.ENA.csv"))[,-1]
tree=keep.tip(tre, tip=unique(tree_df$tip1))#188 tips and 171 internal nodes
tree$tip.label=tree_df[match(tree$tip.label,tree_df$tip1),"Spname"]
save(tree,file="FIA.tree.RData")

##caculate species-trait matrix(trait_mat) and species-plot abundance matrix (abund_mat) ----
require(data.table)
library(dplyr)
tree_df=as.data.frame(fread("fia.ENA.csv"))[,-1]
abund_mat <- tree_df %>%  group_by(pltID,Spname)  %>% summarize(abund_w = sum(Abund_weight,na.rm=TRUE),.groups ="keep")%>% 
		tidyr::pivot_wider(names_from = Spname, values_from = abund_w) %>% tibble::column_to_rownames("pltID")	
abund_mat[is.na(abund_mat)]=0
abund_mat=abund_mat[rowSums(abund_mat,na.rm=T)!=0,] 		
sp.in.plt=unique(tree_df$Spname)
h99 <- tree_df %>% group_by(Spname) %>% summarize(Hmax = quantile(HT,0.99,na.rm=TRUE),.groups ="keep")
try.sp=read.csv("TRY.splist.csv")[,-1] %>% select(AccSpeciesName,Spname) %>% distinct() #link spname in TRY and in FIA
trait.name=c("WD","RootDep.","Le.N","LMA","Hmax")
trait=as.data.frame(fread("TraitInTRY.clean.csv",header=T))[,-1]%>%	full_join(try.sp,by="AccSpeciesName",multiple = "all",relationship ="many-to-many")%>%
	filter(Abbr%in%trait.name&!is.na(Spname))
trait_mat= trait %>% filter(Spname%in%sp.in.plt) %>% group_by(Abbr,Spname) %>% summarize(trait.value = mean(as.numeric(StdValue),na.rm=T),.groups ="keep")%>%
			as.data.frame() %>% tidyr::pivot_wider(names_from = Abbr, values_from = trait.value) %>% full_join(h99, by="Spname") %>%
			.[rowSums(is.na(.[,-1]))!=ncol(.[,-1]),]%>% tibble::column_to_rownames("Spname")
abund_mat=abund_mat[,rownames(trait_mat)]
save(trait_mat,file="trait_mat.Rdata")
save(abund_mat,file="abund_mat.Rdata")

library(dplyr)
mf=function(x)strsplit(x," ")[[1]][1]
trait_mat$genus=do.call(c,lapply(rownames(trait_mat),mf))
ge.mat=trait_mat%>%group_by(genus)%>%summarize(Rmax=mean(RootDep.,na.rm=T),LN=mean(Le.N,na.rm=T),SLA=mean(LMA,na.rm=T))
nrow(trait_mat[trait_mat$genus%in%ge.mat[is.na(ge.mat$Rmax),]$genus,])#Rmax,29 genus is na, containning 41 species, 108 na specie in total
nrow(trait_mat[trait_mat$genus%in%ge.mat[is.na(ge.mat$LN),]$genus,])#LN, 11 genus is na, containning 14 species, 54 na specie in total
nrow(trait_mat[trait_mat$genus%in%ge.mat[is.na(ge.mat$SLA),]$genus,])#SLA, 5 genus is na, containning 8 species, 21 na specie in total

##########################################################
####### evaluate trait missind data and trait space  #####
##########################################################
#Table A3 V1: species with missing trait values
trait_mat=get(load("trait_mat.Rdata"));abund_mat=get(load("abund_mat.Rdata"))
mf=function(x) length(x[x!=0&(!is.na(x))])
mf2=function(x) median(x[x!=0&(!is.na(x))])
sp.stat=data.frame(num.trait.sp=apply(trait_mat,1,mf),
	range.size=apply(abund_mat,2,mf),#num of plts that a sp occupied
	abundance=apply(abund_mat,2,mf2))#total abd of a sp across all plts
range.stat=Rmisc::summarySE(sp.stat,measurevar="range.size", groupvars="num.trait.sp")
abd.stat=Rmisc::summarySE(sp.stat,measurevar="abundance", groupvars="num.trait.sp")

# Table A3 V2: the mean occupied number of plots,mean proportion of BA and SP richness in each plot for each species with missing certain trait values
get.stat=function(spname,abund_mat){
	x=abund_mat[,spname]
	names(x)=rownames(abund_mat)
	range.size=length(x[x>0])
	abd.sp=abund_mat[names(x)[x>0],]
	mf3=function(i,abd.sp,na){
		mat=abd.sp[i,]
		sp.all=mat[mat>0]
		sp.na.t=mat[mat>0&colnames(mat)%in%na]
		ba.na=sum(sp.na.t)/sum(sp.all)*100
		rich.na=length(sp.na.t)/length(sp.all)*100
		return(data.frame(BA=ba.na,Richness=rich.na))
	}
	stat=do.call(rbind,lapply(1:nrow(abd.sp),mf3,abd.sp,spname))
	stat$Spname=spname;stat$range.size=length(x[x>0]);stat$pltID=rownames(abd.sp)
	return(stat)
}

require(parallel)
no_cores <- detectCores() - 1
mycl <- makePSOCKcluster(no_cores);	
trait.name=c("RootDep.","Le.N","LMA")
sp.na.stat.all=c()
for (i in trait.name){
	sp.na=rownames(trait_mat[is.na(trait_mat[,i]),])
	sp.na.stat=do.call(rbind,parLapply(cl=mycl,sp.na,get.stat,abund_mat))
	sp.na.stat$trait=i
	sp.na.stat.all=rbind(sp.na.stat.all,sp.na.stat)
}
stopCluster(mycl)
date()
save(sp.na.stat.all,file="sp.na.stat.all.Rdata")
# Table A3 V3:Within each plot,the proportion of BA and SP richness for species with missing trait values
mf4=function(i,abund_mat,na){
	x=abund_mat[i,]
	sp.all=x[x>0]
	sp.na=x[x>0&colnames(x)%in%na]
	ba.na=sum(sp.na)/sum(sp.all)*100
	rich.na=length(sp.na)/length(sp.all)*100
	pltID.no.year = substr(rownames(abund_mat)[i],1,nchar(rownames(abund_mat)[i])-5)				
	return(data.frame(pltID=rownames(abund_mat)[i],pltID.no.year =pltID.no.year ,BA=ba.na,Richness=rich.na))
}
require(parallel)
no_cores <- detectCores() - 1
mycl <- makePSOCKcluster(no_cores);
trait.name=c("RootDep.","Le.N","LMA")
fin=c()
for (i in trait.name){
	na=rownames(trait_mat[is.na(trait_mat[,i]),])
	re=do.call(rbind,parLapply(cl =mycl,1:nrow(abund_mat),mf4,abund_mat,na))
	re$trait=i
	fin=rbind(fin,re)
}
stopCluster(mycl)
save(fin,file="plt.na.stat.Rdata")
library(dplyr)
stat=fin%>%group_by(trait,pltID.no.year)%>%summarize(Richness=mean(Richness),BA=mean(BA))%>%
	group_by(trait)%>%summarize(SP.mean=mean(Richness),SP.se=sd(Richness)/sqrt(length(Richness)),BA.mean=mean(BA),BA.se=sd(BA)/sqrt(length(BA)))
	
#Fig. S5 plot phylogeny and Trait
#ref:https://yulab-smu.top/treedata-book/chapter7.html
#devtools::install_github("YuLab-SMU/ggtree")
library(ggtree);library(data.table);library(dplyr);library(ggplot2)
trait_mat=get(load("trait_mat.Rdata"))
tree=get(load("FIA.tree.RData"))
ang=tree$tip.label[158:188]
ANG=data.frame(tip1=tree$tip.label,Taxa=ifelse(tree$tip.label%in%ang,"Gymnosperms","Angiosperms"))
trait_mat$tip1=rownames(trait_mat)
trait.sp=ANG %>% left_join(trait_mat,by='tip1')
trait.sp=trait.sp[,c("tip1","Taxa","WD","Le.N","LMA","RootDep.","Hmax")]
colnames(trait.sp)=c("id","Taxa","WD","LN","SLA","Rmax","Hmax")

#trait.sp2[is.na(trait.sp2)]=0
trait.sp2=trait.sp[match(rev(trait.sp$id),trait.sp$id),]
p <- ggtree(tree)  %<+%  trait.sp2 + geom_tippoint(aes(color=Taxa))
p + geom_facet(panel = "WD", data = trait.sp2, geom = geom_col, 
                aes(x = trait.sp2$WD, color = trait.sp2$Taxa, 
                fill = trait.sp2$Taxa), orientation = 'y', width = .6,show.legend=FALSE) +	
    geom_facet(panel = "LN", data = trait.sp2, geom = geom_col, 
                aes(x = trait.sp2$LN, color = trait.sp2$Taxa, 
                fill = trait.sp2$Taxa), orientation = 'y', width = .6,show.legend=FALSE) +	
	geom_facet(panel = "SLA", data = trait.sp2, geom = geom_col, 
                aes(x = trait.sp2$SLA, color = trait.sp2$Taxa, 
                fill = trait.sp2$Taxa), orientation = 'y', width = .6,show.legend=FALSE) +	
	geom_facet(panel = "Rmax", data = trait.sp2, geom = geom_col, 
                aes(x = trait.sp2$Rmax, color = trait.sp2$Taxa, 
                fill = trait.sp2$Taxa), orientation = 'y', width = .6,show.legend=FALSE) +	
	geom_facet(panel = "Hmax", data = trait.sp2, geom = geom_col, 
                aes(x = trait.sp2$Hmax, color = trait.sp2$Taxa, 
                fill = trait.sp2$Taxa), orientation = 'y', width = .6,show.legend=FALSE) +					
	theme_tree2(legend.position=c(.05, .9))+
	guides(color=guide_legend(override.aes = list(size = 4)))+	
	theme(axis.text = element_text(size=15,color='black',angle=15),
		legend.background = element_rect(fill = NA),
		legend.text=element_text(face="bold.italic",size=12),
		legend.title=element_text(face="bold",size=15),
		strip.text = element_text(colour = 'black', face = 'italic', size = rel(1.5)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))
	re=c()
	for (i in c("WD","LN","SLA","Rmax","Hmax")){
		tmp=Rmisc::summarySE(trait.sp,measurevar=i, groupvars="Taxa",na.rm=T)
		tmp2=data.frame(Taxa=tmp$Taxa,Trait=i,Mean=round(tmp[,3],2),SE=round(tmp$se,2))
		re=rbind(re,tmp2)
	}

#Fig. S4 plot trait space
#PCA 
library(ggfortify)
trait.sp2[is.na(trait.sp2)]=0
pca_res <- prcomp(trait.sp2[,-c(1:2)], scale. = TRUE)
autoplot(pca_res, data = trait.sp2, colour = 'Taxa',size=3,alpha=0.8,
         loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE, loadings.label.size = 5, loadings.label.color ="black")+
geom_vline(xintercept=0,col="darkgray",size=1,linetype="longdash",alpha=0.5)+
geom_hline(yintercept=0,col="darkgray",size=1,linetype="longdash",alpha=0.5)+
guides(color=guide_legend(override.aes = list(size = 3)))+	
theme_bw()+theme(axis.text = element_text(size=15,color='black'),
	axis.title = element_text(face="bold",size=18,color='black'),
	legend.position=c(.2, .15),
	legend.background = element_rect(fill = NA),
		legend.text=element_text(face="bold.italic",size=12),
		legend.title=element_text(face="bold",size=15))

#Fig. S2a pairwise trait space
library(GGally)
ggpairs(trait.sp[,-1], upper = list(continuous = wrap("cor", color="black",size=6), combo = "box_no_facet", discrete = "count", na = "na"),
  diag=list(continuous = "barDiag", discrete = "barDiag", na = "naDiag",binwidth=10),
  lower = list(continuous = wrap("points", alpha = 0.8,size=0.2), combo = "box_no_facet", na = "na"))+theme_bw()+
	theme(axis.text = element_text(size=10,color='black'),
		#axis.text.x = element_text(angle=15),
		axis.title = element_text(size=12,color='black'),
		strip.text = element_text(colour = 'black', face = 'italic', size = rel(1.5)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))

#################################################################################################
####### Calculate SP, PD, Hill Evenness, weighted standard deviation of traits within each plot #
#################################################################################################
## caculate SP and PD
library(ape);library(data.table);require(dplyr)
abund_mat=get(load("abund_mat.Rdata"))
comm <- as.matrix(abund_mat)    
comm <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/") #abd(i)/sum(abd)
save(comm,file="comm.Rdata")

tree=get(load("FIA.tree.RData"))
pd.cal=function(k,tree){
	require(ape)
	df2 <- k[k > 0]	
	if (length(df2)==1){
		return(data.frame(SP=1,M_p = 0,M_Ap=0))
	}else{
		tre.plt=keep.tip(tree, tip=names(df2))
		dis2=cophenetic(tre.plt)#pairwise branch length
		#Standardized mean pairwise distance, range (0,1)
		M_p = sum(dis2)/(nrow(dis2)*(nrow(dis2)-1)) 
		#Standardized mean pairwise distance weighted by abundance, Q/S(S-1), Q, Rao’s quadratic entropy; S, Species richness
		M_Ap = sum(df2 %*% dis2 %*% matrix(df2, ncol = 1))/(nrow(dis2)*(nrow(dis2)-1))
		re=data.frame(SP=length(df2),M_p = M_p,M_Ap = M_Ap)
		return(re)
	}
}
re=apply(comm,1,pd.cal,tree)
re2=do.call(rbind,re)
library(hillR)
# pd.shannon=hill_phylo(abund_mat, tree, q = 1)
# pd.simpson=hill_phylo(abund_mat, tree, q = 2)
#mean branch length weighted by abundance
pd=hill_phylo(comm, tree, q = 0)#equal to picante::pd()	
pd.re=data.frame(pltID=names(pd),re2,pd)
save(pd.re,file="pd.re.RData")

## Hill Evenness, as defined by equation 8 in Tuomisto 2012,Oikos 121: 1203-1218 
library(chemodiv);
hilleven=calcDiv(sampleData = abund_mat, type = "HillEven")
hilleven$pltID=rownames(abund_mat)
	
## Calculate weighted standard deviation of traits
library(Hmisc);
tree_df=as.data.frame(fread("fia.ENA.csv"))[,-1]%>%select(pltID,Spname,Abund_weight)
trait_mat=get(load("trait_mat.RData"))
trait_mat$Spname=rownames(trait_mat)
#wtd.var Produced invalid variance estimate if the weights suggest at most one effective observation (sum(Abund_weight) <= 1), 
#so use Abund_weight*10
trait.sd = tree_df %>% left_join(trait_mat,by="Spname")%>% group_by(pltID) %>% 
	summarize(Hmax.sd=sqrt(wtd.var(Hmax, Abund_weight*10,na.rm=TRUE)),WD.sd=sqrt(wtd.var(WD, Abund_weight*10,na.rm=TRUE)),
		RootDep.sd=sqrt(wtd.var(RootDep., Abund_weight*10,na.rm=TRUE)),Le.N.sd=sqrt(wtd.var(Le.N, Abund_weight*10,na.rm=TRUE)),
		LMA.sd=sqrt(wtd.var(LMA, Abund_weight*10,na.rm=TRUE)))
trait.sd = trait.sd%>% left_join(hilleven,by="pltID")
save(trait.sd,file="trait.sd.Rdata")

# Fig.S2b pairwise trait.sd Correlation and hill evenness
library(GGally)
dat=trait.sd[,c("HillEven","WD.sd","Le.N.sd","LMA.sd","RootDep.sd","Hmax.sd")]
colnames(dat)=c("HillEven","WD","LN","SLA","Rmax","Hmax")
ggpairs(dat, upper = list(continuous = wrap("cor", color="black",size=6), combo = "box_no_facet", discrete = "count", na = "na"),
  diag=list(continuous = "barDiag", discrete = "barDiag", na = "naDiag",binwidth=10),
  lower = list(continuous = wrap("points", alpha = 0.01,size=0.2), combo = "box_no_facet"))+theme_bw()+
	theme(axis.text = element_text(size=10,color='black'),
		axis.text.x = element_text(angle=15),
		axis.title = element_text(size=12,color='black'),
		strip.text = element_text(colour = 'black', face = 'italic', size = rel(1.5)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))

#################################################################
####### Calculate FD consider georaphic varaition ###############
#################################################################
### overlay try data with FIA plot using Lon/lat and then caculate FD ###
##set buffer for fia plot and extract trait data within the buffer
library(data.table)
require(magrittr)
require(dplyr)
library(sf)
library(ggpubr)
ENA=c("ME","NH","VT","NY","MA","RI","CT","NJ","PA","DE","MD","MI","OH","IN","IL","WI","WV","VA","NC","TN","KY","SC","GA","AL","MS","FL")
trait=as.data.frame(fread("TraitInTRY.clean.csv"))
coord0=read.csv("TRYcoord.csv")
fia0=as.data.frame(fread("fia.csv"))[,-1]
fia.points=as.data.frame(fread("fia.ENA.csv"))[,-1]%>% select(pltID.no.year,LON,LAT)%>%distinct()%>%
	st_as_sf(coords=c("LON","LAT"))%>% st_set_crs("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=km")

plot.stat=function(traitname,coord3,coord.fia.t,fia.points,limt=NULL,histg=FALSE){
		coord3.t=coord3 %>% select(ObservationID,TRY_ID)
		coord.fia.t2=coord.fia.t %>% full_join(coord3.t,by="TRY_ID")		
		stat=tapply(coord.fia.t2$ObservationID,coord.fia.t2$pltID.no.year,length)
		map.stat=cbind(fia.points,trait.t=stat[match(fia.points$pltID.no.year,names(stat))])	
		if(is.null(limt)) limt=max(map.stat$trait.t,na.rm=T);print(max(map.stat$trait.t,na.rm=T))
		p1=ggplot(map.stat,aes(fill = trait.t,color=trait.t)) +
		geom_sf() +scale_color_gradient(low="palegreen3",high="sandybrown",limits=c(0,limt))+guides(color="none")+
		scale_fill_gradient(low="palegreen3",high="sandybrown",name=traitname,limits=c(0,limt))+		
		theme_light()+ theme(panel.background = element_rect(fill = '#619CFF', colour = 'black'),
			axis.text = element_blank(),legend.background = element_rect(fill = NA),
			legend.text=element_text(face="bold.italic",size=10,angle=30),
			legend.title=element_text(face="bold",size=12),
			legend.position="left")
		if (histg){
			p2=ggplot(map.stat, aes(x=trait.t)) + geom_histogram(colour="black", fill="white",bins=30)+
				geom_vline(aes(xintercept=mean(trait.t)),color="blue", linetype="dashed", size=1)+xlab("Number of trait records")+
				theme_light()+ theme(panel.background = element_rect(fill = '#619CFF', colour = 'black'),
				axis.text.x = element_text(size=10,color='black',angle=0),
				axis.text.y = element_text(size=10,color='black',angle=30),
				axis.title = element_text(size=10,face="bold",color='black',angle=0))
			p=ggarrange(p1,p2,nrow=2,ncol = 1,widths=1,heights=c(4,1))
		} else{
			p=p1#+theme(legend.position="none")
		}			
		return(p)
	}

#buffer for Hmax
coord=fia0 %>% filter(LON>-93.86) %>% group_by(Spname,LAT,LON) %>% 
	summarize(H99 = quantile(HT,0.99,na.rm=TRUE),.groups ="keep") %>% na.omit()
coord$ObservationID=paste("fia",seq(1:nrow(coord)),sep="_")
coord$TraitID=1;coord$TraitType="numeric";coord$Abbr="Hmax"	
coord2=coord[,c("LON","LAT")] %>% distinct()
coord2$TRY_ID=paste("Hmax",seq(1:nrow(coord2)),sep="_")
coord3=coord %>% left_join(coord2,by=c("LON","LAT"))%>%
	as.data.frame()%>%select(ObservationID,LAT,LON,TraitID,TraitType,Abbr,TRY_ID)%>%distinct()
try.points = coord2 %>% st_as_sf(coords=c("LON","LAT")) %>% st_set_crs("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=km") 
coord.fia.t=fia.points %>% st_buffer(1,endCapStyle="ROUND")%>% st_intersection(try.points) %>% 
		as.data.frame() %>% select(c("TRY_ID","pltID.no.year"))%>% distinct()
hmax.p=plot.stat("Hmax",coord3,coord.fia.t,fia.points,,TRUE)

add=data.frame(matrix(ncol = ncol(trait), nrow = nrow(coord)));colnames(add)=colnames(trait)
add$AccSpeciesName=coord$Spname
add$StdValue=coord$H99;add$ObservationID=coord$ObservationID
add$Dataset="FIA";add$TraitID=1;add$TraitName="Max Height"
add$TraitType="numeric";add$DataName="Max Height(99 quantile)"
add$UnitName="m";add$Reference="FIA";add$Abbr="Hmax";add$SpeciesName=coord$Spname
write.csv(add,"TraitInTRY.Hmax.csv")#Amelanchier sanguinea occurs in 2 plots and have NA HT value

#buffer for other traits
get.dat=function(buffer.i,traitname,coord0,fia.points,limt=NULL,get.dat=FALSE){
	coord=coord0 %>% left_join(unique(trait[,c("ObservationID","TraitID","TraitType","Abbr")]),by="ObservationID")%>% 
		na.omit() %>% filter(Abbr%in%traitname&LON<=-60&LON>=-98&LAT>=20&LAT<=55)	
	coord2=unique(coord[,c("LON","LAT")])
	coord2$TRY_ID=paste("try",seq(1:nrow(coord2)),sep="_")
	coord3=coord %>% left_join(coord2,by=c("LON","LAT"))
	coord.fia.t = coord2 %>% select(LON,LAT,TRY_ID)%>% distinct()%>%
		st_as_sf(coords=c("LON","LAT")) %>% st_set_crs("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=km") %>% 
		st_buffer(buffer.i,endCapStyle="ROUND") %>% st_union()%>% st_intersection(fia.points) %>% 
		as.data.frame() %>% select(c("TRY_ID","pltID.no.year"))%>% distinct()
	if(get.dat){
		write.csv(coord3,paste("coord.TRY",traitname,buffer.i,"csv",sep="."))
		write.csv(coord.fia.t,paste("coord.fia",traitname,buffer.i,"csv",sep="."))
	}  
	p=plot.stat(traitname,coord3,coord.fia.t,fia.points,limt)
	return(p)		
}

#Fig. S1 plot the spatial scale of the buffer at which the trait data best cover FIA plot species composition
blank=ggplot()+theme(panel.background = element_rect(fill = 'white', colour = 'white'))
p.wd=get.dat(9,"WD",coord0,fia.points,,TRUE)
pltn=lapply(seq(1,8),get.dat,"WD",coord0,fia.points,230)
pltn[[9]]=p.wd
wd.p=ggarrange(pltn[[1]],pltn[[2]],pltn[[3]],pltn[[4]],pltn[[5]],pltn[[6]],pltn[[7]],pltn[[8]],pltn[[9]],
	nrow=1,ncol = 9,widths=1,heights=c(1,1,1,1,1,1,1,1,1),common.legend=TRUE,legend="left")

p.rd=get.dat(6,"RootDep.",coord0,fia.points,,TRUE)
pltn=lapply(seq(1,5),get.dat,"RootDep.",coord0,fia.points,370)
pltn[[6]]=p.rd
root.p=ggarrange(pltn[[1]],pltn[[2]],pltn[[3]],pltn[[4]],pltn[[5]],pltn[[6]],blank,blank,blank,
	nrow=1,ncol = 9,widths=1,heights=c(1,1,1,1,1,1,1,1,1),common.legend=TRUE,legend="left")

p.len=get.dat(4,"Le.N",coord0,fia.points,,TRUE)	
pltn=lapply(seq(1,3),get.dat,"Le.N",coord0,fia.points,3624)
pltn[[4]]=p.len
leN.p=ggarrange(pltn[[1]],pltn[[2]],pltn[[3]],pltn[[4]],blank,blank,blank,blank,blank,
	nrow=1,ncol = 9,widths=1,heights=c(1,1,1,1,1,1,1,1,1),common.legend=TRUE,legend="left")

p.lma=get.dat(5,"LMA",coord0,fia.points,,TRUE)
pltn=lapply(seq(1,4),get.dat,"LMA",coord0,fia.points,5058)
pltn[[5]]=p.lma
lma.p=ggarrange(pltn[[1]],pltn[[2]],pltn[[3]],pltn[[4]],pltn[[5]],blank,blank,blank,blank,
	nrow=1,ncol = 9,widths=1,heights=c(1,1,1,1,1,1,1,1,1),common.legend=TRUE,legend="left")
	
hmax.p2=ggarrange(hmax.p,blank,blank,blank,blank,blank,blank,blank,blank,
	nrow=1,ncol = 9,widths=1,heights=c(1,1,1,1,1,1,1,1,1),common.legend=TRUE,legend="left")
ggarrange(wd.p,root.p,leN.p,lma.p,hmax.p2,nrow=5,ncol = 1,widths=c(1,1,1,1,1),heights=1)

## caculate the trait diversity(consider lon/lat),FAD and FRic
library(data.table)
require(magrittr)
require(dplyr)
fd.cal.geo=function (k,comm,coord.fia,coord3,trait0,trait_mat0){
	require(magrittr);require(dplyr);require(hillR);require(FD)	
	#caculate FD			
	df2 <- comm[k, ][comm[k, ] > 0] 	
	if (length(df2)<=2){
		return(data.frame(data.frame(pltID=rownames(comm)[k],FRic=NA,FAD=0,MTD=0,MTD_non.wt=0)))
	}else{
		plt.id.k = substr(rownames(comm)[k],1,nchar(rownames(comm)[k])-5)
		sp.in.plt=names(comm[k,][comm[k,]>0])
		tryID.in.plt=coord.fia %>% filter(pltID.no.year==plt.id.k)%>%distinct()
		obsID.in.plt=coord3 %>% filter(TRY_ID%in%tryID.in.plt$TRY_ID)%>%distinct()
		if (nrow(obsID.in.plt)==0) trait_mat=trait_mat0[sp.in.plt,]
		if (nrow(obsID.in.plt)>0){
		trait_mat=trait0 %>% filter(Spname%in%sp.in.plt & ObservationID%in%obsID.in.plt$ObservationID) %>%
				group_by(Abbr,Spname) %>% summarize(trait.value = mean(as.numeric(StdValue),na.rm=T),.groups ="keep")%>%
				as.data.frame() %>% tidyr::pivot_wider(names_from = Abbr, values_from = trait.value)	%>%
				.[rowSums(is.na(.[,-1]))!=ncol(.[,-1]),] %>% tibble::column_to_rownames("Spname")		
		#replace NAs with data that do not consider geo
		mat=t(as.matrix(t(trait_mat)[match(colnames(trait_mat0),colnames(trait_mat)),]))
		colnames(mat)=colnames(trait_mat0);rownames(mat)=rownames(trait_mat)
		trait.add=as.matrix(trait_mat0[rownames(mat),colnames(mat)])
		mat[is.na(mat)]=trait.add[which(is.na(mat),arr.ind = T)]
		trait_mat=mat
		}
		#correct data for gowdis
		#remove duplicate if max=min (which lead to NaN in gowdis)
		gowdis_corect=function(n) {
			n[which.max(n)]=ifelse(max(n,na.rm = T)==min(n,na.rm = T)&length(na.omit(n))>1,
				max(n,na.rm = T)+min(n,na.rm = T)/10000,max(n,na.rm = T));
			return(n)
			}
		trait_mat=apply(trait_mat,2,gowdis_corect)
		G.trait.mat=gowdis(trait_mat,ord = "podani")
		#for V2==TRUE,if sp1 has traits a,b,c, sp2 has traits d,e, then gowdis(sp1,sp2)=NA
		#we assume the gowdis(sp1,sp2) is the mean dis of all other pairwise trait dis within the plot
		G.trait.mat[is.na(G.trait.mat)]=mean(G.trait.mat,na.rm=T)
		if(sum(G.trait.mat)>0){
			df2=df2[names(df2)%in%rownames(trait_mat)]
			dis2 <- as.matrix(G.trait.mat)[names(df2), names(df2)]	
			#Standardized mean pairwise distance, range (0,1)
			M_t = sum(dis2)/(nrow(dis2)*(nrow(dis2)-1)) 
			#Standardized mean pairwise distance weighted by abundance, Q/S(S-1), Q, Rao’s quadratic entropy; S, Species richness
			M_At = sum(df2 %*% dis2 %*% matrix(df2, ncol = 1))/(nrow(dis2)*(nrow(dis2)-1))
			#functional attribute diversity (FAD)
			#ref:Plant Attribute Diversity, Resilience,and Ecosystem Function: The Nature and Significance of Dominant and Minor Species
			FAD=sum(G.trait.mat)
			# hull volume of functional space
			# Some sites had less species than traits so returned FRic is 'NA'
			df3=t(as.matrix(df2)[rownames(trait_mat),])
			rownames(df3)=rownames(comm)[k]
			trait_mat[is.na(trait_mat)]=0
			skip_to_next <- FALSE  
			tryCatch({
				suppressWarnings(re=fundiversity::fd_fric(trait_mat,df3))
			}, error=function(e){skip_to_next <<- TRUE})
			if(skip_to_next) re=data.frame(site=rownames(comm)[k],FRic=NA)
			fin=data.frame(re,FAD=FAD,MTD=M_At,MTD_non.wt=M_t)
			colnames(fin)[1]="pltID"
			return(fin)
		}else{
		return(data.frame(pltID=rownames(comm)[k],FRic=NA,FAD=0,MTD=0,MTD_non.wt=0))
		}					
	}
	}

comm=get(load("comm.Rdata"))
trait_mat0=get(load("trait_mat.Rdata"))[,-3]
#imputation missing values
tree=get(load("FIA.tree.RData"))
library(Rphylopars)
trait_data=trait_mat0
trait_data$species=rownames(trait_data)
trait_data=trait_data[,c("species",colnames(trait_mat0))]
trait_mat1=phylopars(trait_data ,tree,
    pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE)
trait_mat2=	trait_mat1$anc_recon[rownames(trait_mat0),]

trait.name=c("WD","RootDep.","Le.N","LMA","Hmax")
buffer=c(9,6,4,5,1)
coord3=c();coord.fia=c()
for (i in 1:5){	
	coord3=rbind(coord3,as.data.frame(fread(paste("coord.TRY",trait.name[i],buffer[i],"csv",sep="."),header=T))[,-1])
	coord.fia=rbind(coord.fia,as.data.frame(fread(paste("coord.fia",trait.name[i],buffer[i],"csv",sep="."),header=T))[,-1])
	}	
try.sp=read.csv("TRY.splist.csv")[,-1] %>% select(AccSpeciesName,Spname) %>% distinct() %>%na.omit()#link spname in TRY and in FIA
hmax=as.data.frame(fread("TraitInTRY.Hmax.csv",header=T)[,-1])
trait0=as.data.frame(fread("TraitInTRY.clean.csv",header=T))%>% filter(Abbr%in%trait.name)%>%
	full_join(try.sp,by="AccSpeciesName",multiple="all",relationship ="many-to-many")%>% filter((!is.na(Abbr))&(!Abbr%in%"RootDep."))%>%rbind(hmax)
re=c()
for (k in 3581:nrow(comm)){
tem=fd.cal.geo(k,comm,coord.fia,coord3,trait0,trait_mat2)
re=rbind(re,tem)
print(k)
}
require(parallel)
no_cores <- detectCores() - 1
mycl <- makePSOCKcluster(no_cores);	
re=do.call(rbind,parLapply(cl=mycl,1:nrow(comm),fd.cal.geo,comm,coord.fia,coord3,trait0,trait_mat2))
stopCluster(mycl)
date()
save(re,file="FD.nonRmax.imu.considerGEO.Rdata")

##################################################
#### Formal analysis since May, 2023 #############
##################################################

#### link FD,PD with growth data -----
library(data.table)
require(dplyr)
pd.re=get(load("pd.re.RData"))
#fd.geo=get(load("FD.considerGEO.Rdata"))
fd.geo=get(load("FD.nonRmax.considerGEO.Rdata"))
#fd.geo=get(load("FD.nonRmax.imu.considerGEO.Rdata"))
trait.sd=get(load("trait.sd.Rdata"))
tree_df=as.data.frame(fread("fia.ENA.csv"))[,-1]
FRic=get(load("FRic.considerGEO.Rdata"));colnames(FRic)[1]='pltID'
fd.geo$FRic=FRic[match(fd.geo$pltID,FRic$pltID),2]
#remove plots with 0 fd/pd value, or biomass loss>10%,standage>=100 and biogrowth >300. 
fd.bio=list(fd.geo,pd.re,trait.sd,tree_df) %>% purrr::reduce(left_join,by='pltID') %>% 
	filter(STDAGE<100&MTD>0&M_Ap>0&B_L_prop<=0.1&B_plt_t1_MgHa<300&Hmax.sd>0&WD.sd>0&LMA.sd>0&Le.N.sd>0&RootDep.sd>0)%>%
	as.data.frame()
mf=function(x){a=table(x);return(names(a)[which(a==max(a))])}
fd.bio2=fd.bio%>%#na.omit()%>%
	group_by(pltID.no.year)%>%summarize(LON=mean(LON),LAT=mean(LAT),G=mean(G_MassBal_MgHaYr,na.rm=TRUE),SP=mean(SP),MTD=mean(MTD),MBL=mean(M_Ap),FAD=mean(FAD,na.rm=TRUE),PD=mean(pd),
	WD=mean(WD.sd),LN=mean(Le.N.sd),SLA=mean(LMA.sd), Rmax=mean(RootDep.sd), Hmax=mean(Hmax.sd),        
	ECOPROVCD=mf(ECOPROVCD),IniBio=mean(B_plt_t1_MgHa),StdAge=mean(STDAGE),
	FRic=mean(FRic,na.rm=TRUE),HillEven=mean(HillEven,na.rm=TRUE))#Initial.biomass
#fd.bio2$ECOPROVCD=as.factor(fd.bio2$ECOPROVCD)
write.csv(fd.bio2,"fig2.data.all.csv");
measure.time=tree_df%>%select(pltID.no.year,pltID)%>%distinct()%>%group_by(pltID.no.year)%>%tally()%>%group_by(n)%>%tally()
library(ggpubr)
#Fig. S6
fd.p=fd.bio2%>%select(MTD,SP,FRic,ECOPROVCD)%>%na.omit()
fd.p$FRic2=pracma::nthroot(fd.p$FRic,5)
p1=ggplot(fd.p, aes(FRic2,MTD))+geom_point(size=1.2,alpha=0.4)+
labs(x="Functional convex hull volume (5th root)",y="MTD",)+
#geom_smooth(aes(x = FRic2, y = MTD),data=fd.p,method = "loess",color="red",size=1.5,show.legend=FALSE,se =TRUE,linetype=1)+
theme_bw()+
theme(axis.text = element_text(size=20,color='black'),		
	axis.title = element_text(size=25,color='black'))
p2=ggplot(fd.p, aes(SP,FRic2))+geom_point(size=1.2,alpha=0.4)+
labs(x="Species Richness",y="Functional convex hull volume (5th root)")+
#geom_smooth(aes(x = SP, y = FRic2),data=fd.p,method = "loess",color="red",size=1.5,show.legend=FALSE,se =TRUE,linetype=1)+
theme_bw()+
theme(axis.text = element_text(size=20,color='black'),		
	axis.title = element_text(size=25,color='black'))
p=ggarrange(p2,p1,labels="AUTO" ,
	font.label = list(size = 20, color = "black", face = "bold"))
	
p.eco.mtd=ggplot(fd.p, aes(FRic2,MTD))+geom_point(size=1.2,alpha=0.4)+
labs(x="Functional convex hull volume",y="MTD",)+
#geom_smooth(aes(x = FRic2, y = MTD),data=fd.p,method = "loess",color="red",size=1.5,show.legend=FALSE,se =TRUE,linetype=1)+
theme_bw()+facet_wrap(~ECOPROVCD, ncol = 2)+
theme(axis.text = element_text(size=15,color='black'),		
	axis.title = element_text(size=20,color='black'),
	strip.text = element_text(colour = 'black', face = 'italic', size = rel(1.5)), 
	strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))
p.eco.sp=ggplot(fd.p, aes(SP,FRic2))+geom_point(size=1.2,alpha=0.4)+
labs(x="Species Richness",y="Functional convex hull volume")+
#geom_smooth(aes(x = SP, y = FRic2),data=fd.p,method = "loess",color="red",size=1.5,show.legend=FALSE,se =TRUE,linetype=1)+
theme_bw()+facet_wrap(~ECOPROVCD, ncol = 2)+
theme(axis.text = element_text(size=15,color='black'),		
	axis.title = element_text(size=20,color='black'),
	strip.text = element_text(colour = 'black', face = 'italic', size = rel(1.5)), 
	strip.background = element_rect(fill = 'white', colour = 'darkgray', linewidth = rel(2), linetype = 2))
ggarrange(p.eco.sp,p.eco.mtd,labels="AUTO" ,
	font.label = list(size = 20, color = "black", face = "bold"))
	
ggplot(fd.p, aes(SP,HillEven))+geom_point(size=1.2,alpha=0.4)+
labs(x="Species Richness",y="Hill evenness")+
#geom_smooth(aes(x = FRic, y = MTD),data=fd.p,method = "loess",color="red",size=1.5,show.legend=FALSE,se =TRUE,linetype=1)+
theme_bw()+
theme(axis.text = element_text(size=20,color='black'),		
	axis.title = element_text(size=25,color='black'))
#Fig. 1 Fia plot ----
library(sf);library(ggplot2);library(dplyr)
#part1 study area
library(usmap)
ENA=c("ME","NH","VT","NY","MA","RI","CT","NJ","PA","DE","MD","MI","OH","IN","IL","WI","WV","VA","NC","TN","KY","SC","GA","AL","MS","FL")
states <- plot_usmap("states", 
                     color = "black",
                     fill = alpha(0.01)) #this parameter is necessary to get counties to show on top of states
ENA.map <- plot_usmap("states", 
                     color = "red",include =ENA)

NAmap=map_data("world", region = c("USA","Canada","Mexico"))%>%filter(long<=-65&long>=-130&lat<70&lat>22)
ENAmap=map_data("state",region= c("maine","new hampshire","vermont","new york","massachusetts","michigan",
	"ohio","indiana","illinois","wisconsin","west virginia","virginia","north carolina","tennessee","kentucky", 		
	"south carolina","georgia","alabama","mississippi","florida","rhode island","connecticut","new jersey",
	"pennsylvania","delaware","maryland","district of columbia"))
p=ggplot() +  
	geom_polygon(data=NAmap,aes(x=long,y=lat,group=group),color = NA,fill="lightgray",size = 1) + 
	geom_polygon(data=ENAmap,aes(x=long,y=lat,group=group),color = NA,fill="black",size = 1) +
	coord_equal()+theme_bw()+ theme(axis.text = element_blank(),axis.title = element_blank(),
	panel.grid.major =element_blank(),axis.ticks=element_blank(),
	panel.grid.minor =element_blank())
#part 2 main map
leg=data.frame(MAP_UNIT_S=c("211","M211","212","222","221","M221","251","223","231","232"),
	lab=c("Northeastern Mixed Forest","Adirondack-New England Mixed Forest","Laurentian Mixed Forest","Midwest Broadleaf Forest",
	"Eastern Broadleaf Forest",	"Central Appalachian Broadleaf Forest","Temperate Prairie Parkland","Central Interior Broadleaf Forest",
	"Southeastern Mixed Forest","Outer Coastal Plain Mixed Forest"))
fia.points=fd.bio2%>%select(pltID.no.year,LON,LAT,ECOPROVCD)%>%distinct()%>%
		left_join(leg,by=c("ECOPROVCD"="MAP_UNIT_S"))%>%usmap_transform(input_names = c("LON","LAT"))	
the=theme_bw()+theme(panel.background = element_rect(fill = NA, colour = 'black'),
			axis.text = element_blank(),
			axis.title = element_blank(),
			panel.grid =element_blank(),
			#legend.position = c(0.85,0.2),
			legend.background=element_rect(fill='transparent'),
			legend.text=element_text(face="bold",size=10),
			legend.title=element_text(face="bold",size=10),
			legend.key = element_blank())
plot_usmap("states",color = "black",linewidth = 1,include =ENA,fill = alpha(0.6)) + 
	geom_point(data = fia.points, aes(x = x, y = y,color = lab),size = 0.6,alpha = 1) +
	scale_color_manual(values=c("#00BA38","green","brown","darkgreen","#93AA00","#00C19F","#F8766D","#DB72FB","#FF61C3","darkgray"),
	breaks=leg$lab)+
	guides(color=guide_legend(ncol=1,byrow=T,override.aes = list(size = 3)))+	
	labs(color="Ecoprovince")+the+
	patchwork::inset_element(p, left = 0.01, bottom = 0.01, right = 0.35, top = 0.2)
	
#Fig. S7: map diversity and G across Ecoregion ---
library(ggpubr);library(sf)
ecoregion=st_read("S_USA.EcoMapProvinces/S_USA.EcoMapProvinces.shp") 
gwr.map <- fd.bio2 %>% st_as_sf(coords=c("LON","LAT"))%>% st_set_crs("+proj=longlat +datum=WGS84")
theme=theme_light()+ theme(panel.background = element_rect(fill = NA, colour = 'black'),
			panel.grid =element_blank(),
			axis.text = element_text(size=12,color='black',angle=0),
			legend.background = element_rect(fill = NA),
			legend.position = c(0.85,0.28),
			legend.text=element_text(face="bold.italic",size=10),
			legend.title=element_text(face="bold",size=12))
lab=data.frame(vars=c("G","SP","MTD","MBL","FAD","PD"),
	lab=c("Growth","Species \nRichness","MTD","MBL","FAD","Faith PD"))
p=list()
for (i in 1:6){
	mapda=gwr.map
	colnames(mapda)[colnames(mapda)%in%lab$vars[i]]="var"
	if (i%in%c(3,4,5))mapda$var=log(mapda$var)
	p[[i]]=ggplot() + geom_sf(data=mapda,aes(color=var))+	
		geom_sf(data=ecoregion,color="black", alpha = 0) +xlim(-93,-67)+ ylim(25,49)+
		scale_color_gradientn(colors=colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(10),name=lab[i,2])+
		theme
}
ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],nrow=2,ncol =3,
	labels="AUTO" ,
	font.label = list(size = 20, color = "black", face = "bold"))

#cor matrix
# install.packages("psych")
dat=fd.bio2[,c("G","SP","MTD","MBL","FAD","PD")]
#Fig. 1 V1
library(ggcorrplot)	
corr <- round(cor(dat), 2)
p.mat <- cor_pmat(dat)#sig
ggcorrplot(t(corr), hc.order = FALSE, type = "full",tl.srt=0, p.mat = p.mat,
	outline.color = "gray",insig = "pch",lab = TRUE,colors = c("#6D9EC1", "white", "#E46726"))
write.csv(dat,"Fig.2.data.csv")
#Fig. S5 cor among ecogregions
library(ggpubr)
ecoregion=c("211","M211","212","222","221","M221","251","223","231","232")
plt=list()
for (i in ecoregion){
	dat.t=fd.bio2[fd.bio2$ECOPROVCD==i,c("G","SP","MTD","MBL","FAD","PD")]
	corr <- round(cor(dat.t), 2)
	p.mat <- cor_pmat(dat.t)#sig
	plt[[i]]=ggcorrplot(corr, hc.order = FALSE, type = "full",tl.srt=0, p.mat = p.mat,outline.color = "gray",insig = "pch",lab = TRUE,colors = c("#6D9EC1", "white", "#E46726"))
}
ggarrange(plt[[1]],plt[[2]],plt[[3]],plt[[4]],plt[[5]],plt[[6]],plt[[7]],plt[[8]],plt[[9]],plt[[10]],nrow=2,ncol =5,
	labels="AUTO" ,
	font.label = list(size = 20, color = "black", face = "bold"),common.legend=TRUE,legend="right")

#Fig. S9 Scatter plots and histogram of growth and the diversity metrics for each plot.
library(GGally)
p=ggpairs(dat, upper = list(continuous = wrap("cor", color="black",size=6), combo = "box_no_facet", discrete = "count", na = "na"),
  diag=list(continuous = "barDiag", discrete = "barDiag", na = "naDiag",binwidth=10,bins=30),
  lower = list(continuous = wrap("points", alpha = 0.01,size=0.2), combo = "box_no_facet"))+theme_bw()+
	theme(axis.text = element_text(size=10,color='black'),
		axis.text.x = element_text(angle=15),
		axis.title = element_text(size=12,color='black'),
		strip.text = element_text(colour = 'black', face = 'italic', size = rel(1.5)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))
p[2,2] <- p[2,2] + geom_histogram(bins = 10)

#### model regressions -----
#fit LME model
library(lme4);
library(afex)#This will automatically add a p-value column to the output of the lmer for the fixed effects.
fd.bio2[,c("SP","MTD","MBL","FAD","PD","IniBio","StdAge","Hmax","WD","Rmax","LN","SLA")]=
	scale(fd.bio2[,c("SP","MTD","MBL","FAD","PD","IniBio","StdAge","Hmax","WD","Rmax","LN","SLA")])

M=list()
M[[1]] <- lmer(G ~ SP+ FAD + PD  + IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)
M[[2]] <- lmer(G ~ SP+ MTD + MBL  + IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)

p=list()
for (i in 1:2){
dat=as.data.frame(summary(M[[i]])$coef[-c(1,5,6),-3])
colnames(dat)[2]="se"
dat$var=rownames(dat)
dat$type=ifelse(dat$Estimate>0,"positive","negative")
if (i==1) dat$var=factor(dat$var,levels=c("SP","PD","FAD"))
if (i==2) dat$var=factor(dat$var,levels=c("SP","MBL","MTD"))
p[[i]]=ggplot(dat, aes(x=Estimate, y=var,fill=type)) + 	
		#geom_text_repel(aes(x=gradient, y=est,label=est.p,fontface = "italic"),size=2.5)+			
		geom_errorbar(aes(xmin=Estimate-1.96*se, xmax=Estimate+1.96*se,color=type),width=0.15,size=1,show.legend=FALSE,alpha=0.8)+
		# scale_fill_manual(values=mycol)+	
 		#scale_y_continuous(labels = scales::label_number(accuracy = 0.01))+		
		geom_point(size=3,shape=21,color="black",show.legend=FALSE) + 
		labs(x="Estimated slope",y="Diversity metrics") +
		theme_bw()+
		geom_vline(xintercept=0,col="darkgray",size=1,linetype="longdash",alpha=0.5)+
		theme(axis.text = element_text(size=12,color='black'),
			axis.title = element_text(size=15,color='black'),
			plot.margin = margin(t = 0.1, r = 0.6, b = 0.1, l = 0.1, "cm"))
}
ggarrange(p[[1]],p[[2]],nrow=2,ncol=1,labels="AUTO")
	
## Section 2:
#LM models
M=list()
M[[1]] <- lmer(G ~ SP + PD + FAD + IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)
M[[2]] <- lmer(G ~ SP + MBL + MTD + IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)

M[[3]] <- lmer(G ~ SP + MBL+ IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)
M[[4]] <- lmer(G ~ SP + MTD+ IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)
M[[5]] <- lmer(G ~ MBL + MTD+ IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)

M[[6]] <- lmer(G ~ SP + PD+ IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)
M[[7]] <- lmer(G ~ SP + FAD+ IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)
M[[8]] <- lmer(G ~ PD + FAD+ IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)

M[[9]] <- lmer(G ~ SP+ IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)
M[[10]] <- lmer(G ~ PD+ IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)
M[[11]] <- lmer(G ~ FAD+ IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)
M[[12]] <- lmer(G ~ MBL+ IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)
M[[13]] <- lmer(G ~ MTD+ IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)

M[[14]] <- lmer(G ~ SP + MBL + Hmax + IniBio + StdAge+ (1|ECOPROVCD),data=fd.bio2)
M[[15]] <- lmer(G ~ SP + MBL + WD + IniBio + StdAge+ (1|ECOPROVCD),data=fd.bio2)
M[[16]] <- lmer(G ~ SP + MBL + Rmax + IniBio + StdAge+ (1|ECOPROVCD),data=fd.bio2)
M[[17]] <- lmer(G ~ SP + MBL + LN + IniBio + StdAge+ (1|ECOPROVCD),data=fd.bio2)
M[[18]] <- lmer(G ~ SP + MBL + SLA + IniBio + StdAge+ (1|ECOPROVCD),data=fd.bio2)

M.nodiv= lmer(G ~ IniBio + StdAge + (1|ECOPROVCD),data=fd.bio2)

re=c();coef.stat=c()
for (i in 1:18){
tmp=data.frame(model=i,AIC=AIC(M[[i]]),MuMIn::r.squaredGLMM(M[[i]])*100,partR=(MuMIn::r.squaredGLMM(M[[i]])-MuMIn::r.squaredGLMM(M.nodiv))*100)
tmp2=data.frame(model=i,summary(M[[i]])$coef[-1,])
coef.stat=rbind(coef.stat,tmp2)
re=rbind(re,tmp)
}

##conduct model along ecoregion, std and inibio ##
## Using moving window along productivity and stand age gradient, 
## calculate R2 of G-Diversity (incl. PD, FD, SP) within each window slide
fd.bio2=fread("fig2.data.all.csv")
library(dplyr);library(lme4)
GD=function(fd.bio2,var.list,var.t,var.other=NULL,bin.wid=NULL){
	fd.bio=as.data.frame(fd.bio2)
	bin=c(seq(from=0,to=max(fd.bio[,var.t]),by=bin.wid),max(fd.bio[,var.t]))
	re=c();
	for(i in 2:length(bin)){
		if (i<length(bin)){
			fd.bio.t=subset(fd.bio[,c("G",var.other,var.list,"ECOPROVCD")],fd.bio[,var.t]>=bin[i-1]&fd.bio[,var.t]<bin[i])
		} else {
			fd.bio.t=subset(fd.bio[,c("G",var.other,var.list,"ECOPROVCD")],fd.bio[,var.t]>=bin[i-1]&fd.bio[,var.t]<=bin[i])
		}
		j=i-1
		while (nrow(fd.bio.t)<=200&i<=length(bin)){
			i=i+1
			if (i<(length(bin)-1)){
				fd.bio.t.t=subset(fd.bio[,c("G",var.other,var.list,"ECOPROVCD")],fd.bio[,var.t]>=bin[i-1]&fd.bio[,var.t]<bin[i])
			} else {
				fd.bio.t.t=subset(fd.bio[,c("G",var.other,var.list,"ECOPROVCD")],fd.bio[,var.t]>=bin[i-1]&fd.bio[,var.t]<=bin[i])
			}
				fd.bio.t=rbind(fd.bio.t,fd.bio.t.t)		
			}
		if (nrow(fd.bio.t)>10){
			var.valu=apply(fd.bio.t[,c(var.list,var.other)],2,function(x){length(unique(x))})				
			fd.bio.t2=fd.bio.t[,c("G",names(var.valu)[var.valu>=3],"ECOPROVCD")]
			fd.bio.t2=data.frame(G=fd.bio.t2$G,apply(fd.bio.t2[,names(var.valu)[var.valu>=3]],2,scale),ECOPROVCD=fd.bio.t2$ECOPROVCD)
			#fd.bio.t2=fd.bio.t;colnames(fd.bio.t2)[1]="G"
			formu=as.formula(paste("G ~",paste(names(var.valu)[var.valu>=3],collapse="+"),"+(1|ECOPROVCD)",sep=" "))
			M0 <- lmer(formu,data=fd.bio.t2)
			formu.nodiv=as.formula(paste("G ~",paste(var.other,collapse="+"),"+(1|ECOPROVCD)",sep=" "))
			M.nodiv= lmer(formu.nodiv,data=fd.bio.t2)
			
			if (i>length(bin)) bin.t=bin[i-1] else bin.t=bin[i]
			if (bin.t-bin[j]>bin.wid) bin.t=bin[j+1]
			tmp=data.frame(gradient=bin.t,nrow=nrow(fd.bio.t),vars=names(var.valu)[var.valu>=3],
					summary(M0)$coef[-1,],AIC=AIC(M0),MuMIn::r.squaredGLMM(M0)*100,partR=(MuMIn::r.squaredGLMM(M0)[1]-MuMIn::r.squaredGLMM(M.nodiv)[1])*100)			
			re=rbind(re,tmp)
		}else {
			next;}		
		if (i>length(bin)) break;
	}
	return(re)
}

plt.GD=function(re,x.lab,var.list,var.other,reg=TRUE,is.gradient.num=TRUE){
	re$vars=factor(re$vars,levels=c(var.list,var.other))
	p1=ggplot(re, aes(x=gradient, y=Estimate,fill=vars)) 	
	if (is.gradient.num) {
		p1=p1+ geom_errorbar(aes(ymin=Estimate-Std..Error, ymax=Estimate+Std..Error,color=vars),size=0.8,width=max(re$gradient)/50,show.legend=FALSE,alpha=0.6)
	} else {
		p1=p1+ geom_errorbar(aes(ymin=Estimate-Std..Error, ymax=Estimate+Std..Error,color=vars),size=0.8,width=0.5,show.legend=FALSE,alpha=0.6)
	}	
	p1=p1+geom_point(size=3,shape=21,color="black",alpha=0.6) + 
		labs(x=x.lab,y="Standardized \nestimate slope",fill="Parameters") +theme_bw()+
		geom_hline(yintercept=0,col="darkgray",size=1,linetype="longdash",alpha=0.5)+
		guides(fill=guide_legend(override.aes = list(size = 5)))+	
		theme(axis.text = element_text(size=12,color='black'),
			#axis.text.y = element_text(angle=0),		
			axis.title.y = element_text(size=15,color='black'),
			axis.title.x= element_blank(),
			#axis.text.x = element_blank(),
			legend.background=element_rect(fill='transparent'),
			legend.text=element_text(face="bold",size=12),
			legend.title=element_text(face="bold",size=15),
			legend.position = c(0.15,0.75))
	if (reg==TRUE) p1=p1+geom_smooth(aes(x=gradient, y=Estimate,color= vars),data=re,method = "loess",size=1,show.legend=FALSE,se =FALSE,linetype=1)
	
	p3= ggplot(unique(re[,c("gradient","nrow")]), aes(x=gradient, y=nrow))+
			geom_bar(stat="identity")+theme_bw()+
			xlab(x.lab) +ylab("Number \nof plots")+
			theme(axis.text = element_text(size=12,color='black'),
				axis.text.y = element_text(angle=0),	
				axis.title = element_text(size=15,color='black')) 
	p.re=ggarrange(p1,p3,nrow=2,ncol = 1,widths=1,heights=c(3,1),common.legend=TRUE)	
	return(p.re)
}
var.list=c("SP","MBL","MTD");#FDis present similar patterns	
re.std=GD(fd.bio2,var.list,var.other="IniBio",var.t="StdAge",bin.wid=10)#var.other if not null, plot diversity indexes plus var.other
re.inibio=GD(fd.bio2,var.list,var.other="StdAge",var.t="IniBio",bin.wid=50)
GD2=function(fd.bio2,var.list,var.t,var.other=NULL,bin.wid=NULL){
	fd.bio=as.data.frame(fd.bio2)
	fd.bio$bin=round(fd.bio[,var.t]/3)*3
	fd.bio[fd.bio$bin>15,"bin"]=15
	bin=sort(unique(fd.bio$bin))
	re=c();
	for(i in 1:length(bin)){
		fd.bio.t=subset(fd.bio[,c("G",var.other,var.list,"ECOPROVCD")],fd.bio$bin==bin[i])
		var.valu=apply(fd.bio.t[,c(var.list,var.other)],2,function(x){length(unique(x))})				
		fd.bio.t2=fd.bio.t[,c("G",names(var.valu)[var.valu>=3],"ECOPROVCD")]
		fd.bio.t2=data.frame(G=fd.bio.t2$G,apply(fd.bio.t2[,names(var.valu)[var.valu>=3]],2,scale),ECOPROVCD=fd.bio.t2$ECOPROVCD)
		formu=as.formula(paste("G ~",paste(names(var.valu)[var.valu>=3],collapse="+"),"+(1|ECOPROVCD)",sep=" "))
		M0 <- lmer(formu,data=fd.bio.t2)
		formu.nodiv=as.formula(paste("G ~",paste(var.other,collapse="+"),"+(1|ECOPROVCD)",sep=" "))
		M.nodiv= lmer(formu.nodiv,data=fd.bio.t2)
		tmp=data.frame(gradient=bin[i],nrow=nrow(fd.bio.t),vars=names(var.valu)[var.valu>=3],
					summary(M0)$coef[-1,],AIC=AIC(M0),MuMIn::r.squaredGLMM(M0)*100,partR=(MuMIn::r.squaredGLMM(M0)[1]-MuMIn::r.squaredGLMM(M.nodiv)[1])*100)			
		re=rbind(re,tmp)
	}
	return(re)
}
re.sp=GD2(fd.bio2,var.list,var.other=c("StdAge","IniBio"),var.t="SP",bin.wid=3)
re.sp[re.sp$vars=="SP","Estimate"]=abs(re.sp[re.sp$vars=="SP","Estimate"])
leg=data.frame(ECOPROVCD=c("211","M211","212","222","221","M221","251","223","231","232"),
	lab=c("Northeastern Mixed Forest","Adirondack-New England Mixed Forest","Laurentian Mixed Forest","Midwest Broadleaf Forest",
	"Eastern Broadleaf Forest",	"Central Appalachian Broadleaf Forest","Temperate Prairie Parkland","Central Interior Broadleaf Forest",
	"Southeastern Mixed Forest","Outer Coastal Plain Mixed Forest"))
re.eco=c()
for (i in leg$ECOPROVCD){
	fd.bio.t=fd.bio2 %>% filter(ECOPROVCD==i)
	fd.bio.t=data.frame(G=fd.bio.t$G,apply(fd.bio.t[,c("SP","MTD","MBL","FAD","PD","IniBio","StdAge")],2,scale))
	#fd.bio.t[,c("SP","MTD","MBL","FAD","PD","IniBio","StdAge")]=scale(fd.bio.t[,c("SP","MTD","MBL","FAD","PD","IniBio","StdAge")])
	## fit OLS model
	formu=as.formula(paste("G ~",paste(c(var.list,"IniBio","StdAge"),collapse="+"),sep=" "))			
	M0 <- lm(formu ,data=fd.bio.t)
	tmp=data.frame(gradient=i,nrow=nrow(fd.bio.t),vars=c(var.list,"IniBio","StdAge"),
					summary(M0)$coef[-1,],AIC=AIC(M0),Adj.r=summary(M0)$adj.r*100, p=summary(M0)$coefficients[1, 4],ecoregion=leg[match(i,leg$ECOPROVCD),"lab"])			
	re.eco=rbind(re.eco,tmp)
}

library(ggpubr)
p1=plt.GD(re.std,"Stand age (year)",var.list,var.other="IniBio")
p2=plt.GD(re.inibio,"Initial biomass (Mg)",var.list,var.other="StdAge")
p3=plt.GD(re.eco,"Ecoregion",var.list,var.other=c("StdAge","IniBio"),reg=FALSE,is.gradient.num=FALSE)
p4=plt.GD(re.sp,"Species richness",var.list,var.other=c("StdAge","IniBio"),reg=FALSE,is.gradient.num=FALSE)
ggarrange(p1,p2,nrow=1,labels="AUTO",common.legend=TRUE)
ggarrange(p3,p4,nrow=1,labels="AUTO",common.legend=TRUE)
write.csv(re.std,"re.std.csv");write.csv(re.inibio,"re.inibio.csv");write.csv(re.eco,"re.eco.csv");write.csv(re.sp,"re.sp.csv")

##null model
get.sim=function (s, t = 4,tree, r.rep, 
    p = 100, trait_mat=NULL, abun.method = c("lnorm", 
        "norm", "unif"), w.abun = TRUE,para=TRUE) {
	nb.sp <- rep(s, r.rep)
    abun.method <- match.arg(abun.method)    
    abun <- list(rep(0, p))
    abun <- rep(abun, length(nb.sp))
    fill.abun <- function(x, y, z) {
        set <- sample(1:length(x), size = y)
        if (z == "lnorm") 
            x[set] <- rlnorm(length(set))
        if (z == "norm") 
            x[set] <- rnorm(length(set))
        if (z == "unif") 
            x[set] <- runif(length(set))
        return(x)
    }
    abun <- mapply(fill.abun, abun, nb.sp, MoreArgs = list(z = abun.method))
    abun <- t(abun)
    no.tr <- 1:t
    tr <- c("tr")
    names.tr <- paste(tr, no.tr, sep = "")
    no.sp <- 1:p
    sp <- c("sp")
    names.sp <- paste(sp, no.sp, sep = "")
	
	if (is.null(trait_mat)){
		traits <- matrix(NA, p, t)
		if (tr.method == "unif") 
			traits <- apply(traits, 2, function(p) runif(p))
		if (tr.method == "norm") 
			traits <- apply(traits, 2, function(p) rnorm(p))
		if (tr.method == "lnorm") 
			traits <- apply(traits, 2, function(p) rlnorm(p))
	}else{
		#traits=as.matrix(gowdis(trait_mat,ord = "podani"))
		traits=trait_mat
	}  
	dimnames(traits) <- list(names.sp, names.tr)			
    no.com <- 1:(length(nb.sp))
    names.com <- paste("com", no.com, sep = "")
    dimnames(abun) <- list(names.com, names.sp)
	tree$tip.label =sample(colnames(abun)) 
	get.mtd=function(k,abun,traits,tree){
		require(ape);require(FD)    
		df2 <- abun[k, ][abun[k, ] > 0]		
		#MBL
		if (length(df2)==1){
			MBL=0
		}else{
			tre.plt=keep.tip(tree, tip=names(df2))
			dis2=cophenetic(tre.plt)#pairwise branch length
			MBL = sum(df2 %*% dis2 %*% matrix(df2, ncol = 1))/(nrow(dis2)*(nrow(dis2)-1))
		}
		
		#MTD
		if (length(df2)<=2){
			MTD=0
		}else{
			sp.in.plt=names(abun[k,])[abun[k,]>0]	
			trait_mat =	subset(traits,rownames(traits)%in%sp.in.plt)
			gowdis_corect=function(n) {
					n[which.max(n)]=ifelse(max(n,na.rm = T)==min(n,na.rm = T)&length(na.omit(n))>1,
						max(n,na.rm = T)+min(n,na.rm = T)/10000,max(n,na.rm = T));
					return(n)
					}
			trait_mat=apply(trait_mat,2,gowdis_corect)
			G.trait.mat=gowdis(trait_mat,ord = "podani")
			G.trait.mat[is.na(G.trait.mat)]=mean(G.trait.mat,na.rm=T)
			if(sum(G.trait.mat)>0){
				df2=df2[names(df2)%in%rownames(trait_mat)]
				dis2 <- as.matrix(G.trait.mat)[names(df2), names(df2)]	
				MTD = sum(df2 %*% dis2 %*% matrix(df2, ncol = 1))/(nrow(dis2)*(nrow(dis2)-1))				
			}else{
				MTD=0
			}
		}
		return(data.frame(pltID=rownames(abun)[k],SP.cat=length(df2),MTD=MTD,MBL=MBL))
	}
	if (para==TRUE){
		require(parallel)
		no_cores <- detectCores() - 1
		mycl <- makePSOCKcluster(no_cores); 
		results <- do.call(rbind,parLapply(cl=mycl,1:nrow(abun),get.mtd,abun,traits,tree))
		stopCluster(mycl)
	}else{
		results <- do.call(rbind,lapply(1:nrow(abun),get.mtd,abun,traits,tree))
    }
	# values <- FD::dbFD(traits, abun, calc.FDiv = F, w.abun = w.abun, 
        # messages = F, calc.CWM = F)
    # results <- cbind(nb.sp, values$RaoQ)
    # names.var <- c("sp",  "RaoQ")
    # dimnames(results) <- list(names.com, names.var)
	# results=as.data.frame(results)
	# results$MTD=results$RaoQ/(results$sp*(results$sp-1))
    return(results)
}

library(FD);library(data.table);library(dplyr);library(ggcorrplot)	
#tree_df=as.data.frame(fread("fia.ENA.csv"))[,-1]
sp.all=188 #length(unique(tree_df$AccSpeciesName))
fd.bio=fread("fig2.data.all.csv")[,-1]
fd.bio$sp.cat=round(fd.bio$SP/3)*3
r.rep=fd.bio%>%group_by(sp.cat)%>%tally()
trait_mat0=get(load("trait_mat.Rdata"))[,-3]
#imputation missing values
tree=get(load("FIA.tree.RData"))
library(Rphylopars)
trait_data=trait_mat0
trait_data$species=rownames(trait_data)
trait_data=trait_data[,c("species",colnames(trait_mat0))]
trait_mat1=phylopars(trait_data ,tree,
    pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE)
trait_mat2=	trait_mat1$anc_recon[rownames(trait_mat0),]

#sim1=get.sim(s=r.rep$sp.cat,t=4,r.rep=r.rep$n,p=sp.all,tr.method = "norm",abun.method ="lnorm")
sim2=get.sim(s=r.rep$sp.cat,t=4,tree=tree, r.rep=10000,p=sp.all,trait_mat=trait_mat2,abun.method ="lnorm")
#sim3=get.sim(s=r.rep$sp.cat,t=4,r.rep=rep(100,nrow(r.rep)),p=sp.all,tr.method = "norm",abun.method ="lnorm")
load("simulation.rda")
#caculste SES: standardized effect size = (obs value - mean of the expected value)/standard deviation of the expected values. 
sim.stat=sim2%>%group_by(SP.cat)%>%summarize(MTD.m=mean(MTD),MTD.sd=sd(MTD),MBL.m=mean(MBL),MBL.sd=sd(MBL))
fd.bio2=fd.bio%>%left_join(sim.stat,by=c("sp.cat"="SP.cat"))
fd.bio2$MTD.ses =(fd.bio2$MTD-fd.bio2$MTD.m)/fd.bio2$MTD.sd
fd.bio2$MBL.ses =(fd.bio2$MBL-fd.bio2$MBL.m)/fd.bio2$MBL.sd

library(lme4);
library(afex)#This will automatically add a p-value column to the output of the lmer for the fixed effects.
fd.bio3=cbind(fd.bio2[,c("G","ECOPROVCD")],scale(fd.bio2[,c("SP","MTD.ses","MBL.ses","IniBio","StdAge")]))
M <- lmer(G ~ SP+ MTD.ses + MBL.ses  + IniBio + StdAge + (1|ECOPROVCD),data=fd.bio3)
M.nodiv <- lmer(G ~ IniBio + StdAge + (1|ECOPROVCD),data=fd.bio3)
tmp=data.frame(AIC=AIC(M),MuMIn::r.squaredGLMM(M)*100,partR=(MuMIn::r.squaredGLMM(M)-MuMIn::r.squaredGLMM(M.nodiv))*100)
tmp2=data.frame(summary(M)$coef[-1,])

dat=as.data.frame(summary(M)$coef[-c(1,5,6),-3])
colnames(dat)[2]="se"
dat$var=rownames(dat)
dat$type=ifelse(dat$Estimate>0,"positive","negative")
dat$var=factor(dat$var,levels=c("SP","MBL.ses","MTD.ses"))
ggplot(dat, aes(x=Estimate, y=var,fill=type)) + 	
		#geom_text_repel(aes(x=gradient, y=est,label=est.p,fontface = "italic"),size=2.5)+			
		geom_errorbar(aes(xmin=Estimate-se, xmax=Estimate+se,color=type),width=0.15,size=1,show.legend=FALSE,alpha=0.8)+
		# scale_fill_manual(values=mycol)+	
 		#scale_y_continuous(labels = scales::label_number(accuracy = 0.01))+		
		geom_point(size=3,shape=21,color="black",show.legend=FALSE) + 
		labs(x="Estimated slope",y="Diversity metrics") +
		theme_bw()+
		geom_vline(xintercept=0,col="darkgray",size=1,linetype="longdash",alpha=0.5)+
		theme(axis.text = element_text(size=12,color='black'),
			axis.title = element_text(size=15,color='black'),
			plot.margin = margin(t = 0.1, r = 0.6, b = 0.1, l = 0.1, "cm"))


##code graveyard
round(cor(sim[,"MTD"],fd.bio$MTD),2)
corr <- round(cor(sim[,c("sp","MTD")]), 2)#-0.64
p.mat <- cor_pmat(dat)#sig

library(ggpubr)
sim.p1=sim1
sim.p1$sp=factor(sim.p1$sp,levels=r.rep$sp.cat)

sim.p2=sim2
sim.p2$sp=factor(sim.p2$sp,levels=r.rep$sp.cat)
theme=theme_bw()+
	theme(axis.text = element_text(size=20,color='black'),		
	axis.title = element_text(size=25,color='black'))

p1=ggplot(fd.bio, aes(SP,MTD))+geom_point(size=1.2,alpha=0.4)+
labs(x="Species Richness",y="MTD")+ylim(0,2.7)+theme
#geom_smooth(aes(x = SP, y = FRic2),data=fd.p,method = "loess",color="red",size=1.5,show.legend=FALSE,se =TRUE,linetype=1)+

p2=ggplot(sim.p1, aes(sp,MTD))+geom_boxplot(outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE)+
labs(x="Species richness",y="Simulated MTD")+ylim(0,2.7)+theme

p3=ggplot(sim.p2, aes(sp,MTD))+geom_boxplot(outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE)+
labs(x="Species richness",y="Simulated MTD")+ylim(0,2.7)+theme
		
ggarrange(p1,p2,p3,labels="AUTO" ,
	font.label = list(size = 20, color = "black", face = "bold"))
	