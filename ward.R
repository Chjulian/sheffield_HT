#   /.)
#  /)\|   
# // /     
# /'" " 
#################################
#Outbreaker2 for ward studies
#################################

#################################
#Load libraries and outbreaker extra functions
#################################
if(!require(pacman)) install.packages(pacman)
pacman::p_load(ape, linelist, epitrix, fitdistrplus,
               outbreaker2, tidyverse, rio, 
               magrittr, remotes, epicontacts, 
               visNetwork, lubridate, pairsnp,
               ggplot2) #devtools::install_github("gtonkinhill/pairsnp-r")
#https://github.com/reconhub/epicontacts to info on how to install epicontacts
source("functions_maxdist.R")
'%!in%' <- function(x, y)! ('%in%'(x, y))
vector.is.empty <- function(x) return(length(x) ==0 )

#################################
#load linelist data
#################################
# load data, clean dates and filter-out observations
mydata <- read.csv("data/dataset_imputed_2020-12-08.csv")
mydata<- clean_data(mydata, guess_dates = which(names(mydata)=='DateOfOnset_forAnalysis'))
mydata <- mydata[mydata$category%in%c('inpatient', 'outpatient','staff'),]
mydata <- mydata[mydata$loc1%!in%c('com'),] #rm community associated
mydata <- mydata[!is.na(mydata$dateofonset_foranalysis),]

#################################
#Filter but subtype if needed
#################################
#table(mydata$lineage)
# mydata <- mydata[mydata$lineage=='b_2',] 

#################################
#load dna data in DNAbin format, 
#################################
#one sequence per individual
dna <- read.FASTA("data/hoci_phylo-2020-09-03.aln")
mylabels <- clean_data(as.data.frame(names(dna)))
names(dna) <- mylabels$names_dna
dna <- dna[mydata$barcode]; mylabels <- names(dna)
write.FASTA(dna, "data/hoci_phylo-2020-09-03.ord.aln")
dna.SNP <- import_fasta_sparse("data/hoci_phylo-2020-09-03.ord.aln") %>% snp_dist(.)
rownames(dna.SNP) <- colnames(dna.SNP) <- mylabels; rm(mylabels)


#################################
#process onset dates
#################################
date_onset<- mydata$dateofonset_foranalysis
names(date_onset) <- rep("local", length(date_onset)) 

#################################
#Consolidate ward data
#################################
# wards: dataframe with four columns, 
#indicating the case ID, ward ID, start date (i.e. admission to ward) 
#and exit date (i.e. discharge from ward). 
#You can have multiple rows per case for admissions to different wards. 

#create wards object for staff
mydata$ward <- mydata$ward <- paste(mydata$loc1, mydata$loc2,
                     mydata$loc3, mydata$loc4,
                     mydata$loc5, mydata$loc6,
                     mydata$loc7, mydata$loc8)
mydata$ward.count <- as.integer(sapply(mydata$ward, function(rw) length(strsplit(trimws(rw)," ")[[1]])))
mydata$ward <- trimws(mydata$ward)
for(i in 1:nrow(mydata)){
        if(mydata$ward.count[i]==0){
                mydata$ward[i] <- NA
        } else mydata$ward[i] <- ifelse(mydata$ward.count[i]>1 | mydata$loc1[i]=='multiple', 'multiple', mydata$ward[i])
}; rm(i)
wards.staff <- mydata[mydata$category=='staff',c('barcode','ward','dateofcollection', 'dateofcollection')]
wards.staff <- clean_data(wards.staff, guess_dates = c(3,4))
names(wards.staff)<-c('id','ward','adm', 'dis')
wards.staff$adm <- wards.staff$adm-14

#load wards object for patients
wards.patients<- read.csv("data/patient-ward-movements-20201207.csv")
wards.patients<- wards.patients[,c('Barcode','Ward','In','Out')]
names(wards.patients) <- c('id','ward','adm','dis') #rename as examples files for consistency
wards.patients<- clean_data(wards.patients, guess_dates = c(3,4))
wards.patients<- wards.patients[wards.patients$id%in%unique(mydata$barcode),]

#combine wards objects
pre.wards <- rbind(wards.staff,wards.patients)
pre.wards<- pre.wards[order(pre.wards["id"],pre.wards[,"adm"], pre.wards[,"dis"] ),]

#check ward occupation data do not overlap
#subset observations with multiple ward occupation
multiple <- multiple.names <- c()
for (i in unique(pre.wards$id)){
        myind <- pre.wards[pre.wards$id==i,] 
        if(dim(myind)[1]>1) {
                multiple<- c(multiple, dim(myind)[1])
                multiple.names <- c(multiple.names, i)
        }
}
names(multiple) <- multiple.names; rm(i,myind,multiple.names)
wards <- pre.wards[pre.wards$id%!in%names(multiple),]
multiples <- pre.wards[pre.wards$id%in%names(multiple),]
#clean observations with multiple ward occupation
multiples.fixed <- data.frame()
for (i in unique(multiples$id)){
        myind.multiple <- multiples[multiples$id==i,]
        myind.multiple <- myind.multiple %>% distinct() #remove duplicated rows
        while (nrow(myind.multiple)>1) {
                myind <- myind.multiple[1:2,]; myind.multiple <- myind.multiple[-c(1,2),]
                interval1<-interval(myind$adm[1],myind$dis[1])
                interval2<-interval(myind$adm[2],myind$dis[2])
                myintersect <- intersect(interval1, interval2)
                if(!int_overlaps(interval1,interval2)) { #check if times do not overlap
                        multiples.fixed <- rbind(multiples.fixed, myind[1,])
                        myind.multiple <- rbind(myind[2,],myind.multiple)
                } else { #if overlap, then adjust
                        myperiod <- as.period(myintersect)
                        myperiod@day <- myperiod@day+1 #count same day
                        loc <- unique(myind$ward)
                        if(myind$adm[1]!=myind$dis[1] & myind$adm[2]!=myind$dis[2]){ #both intervals are multiple days
                                myind$dis[1] <- myind$dis[1]-myperiod@day
                                myind$adm[2] <- myind$adm[2]+myperiod@day
                                myind[nrow(myind)+1,] <- NA; myind$id[3]<- myind$id[1]
                                myind$adm[3] <- myind$dis[3] <- as.Date(myintersect@start)
                                myind$ward[3] <- ifelse(length(loc)==1, loc, 'multiple')   
                                myind<-myind[c(1,3,2),]
                                multiples.fixed <- rbind(multiples.fixed, myind[1:2,])
                                myind.multiple <- rbind(myind[3,],myind.multiple)
                        } else {#one intervals is a single day
                                if(myind$adm[1]!=myind$dis[1] & myind$adm[2]==myind$dis[2]){ #second interval is one single day
                                        myind$ward[2] <- ifelse(length(loc)==1, loc, 'multiple') #set uncertainty in the location
                                        myind$dis[1] <- myind$dis[1]-myperiod@day
                                        multiples.fixed <- rbind(multiples.fixed, myind[1,])
                                        myind.multiple <- rbind(myind[2,],myind.multiple)
                                }else{
                                        if(myind$adm[1]==myind$dis[1] & myind$adm[2]!=myind$dis[2]) { #first interval is one single day
                                                myind$ward[1] <- ifelse(length(loc)==1, loc, 'multiple') #set uncertainty in the location
                                                myind$adm[2] <- myind$adm[2]+myperiod@day
                                                multiples.fixed <- rbind(multiples.fixed, myind[1,])
                                                myind.multiple <- rbind(myind[2,],myind.multiple)
                                        } else{ #both intervals are one single day
                                                myind$ward <- ifelse(length(loc)==1, loc, 'multiple') #uncertainty in the location
                                                myind.multiple <- rbind(myind[2,],myind.multiple)
                                        }
                                }
                        }
                        rm(myperiod, loc)
                        
                }
                
                rm(myind,interval1,interval2,myintersect)
                
        }
        multiples.fixed <- rbind(multiples.fixed,myind.multiple)
}; rm(i, multiple, multiples, myind.multiple)
#Join ward data
wards <- rbind(wards, multiples.fixed); rm(pre.wards, multiples.fixed, wards.patients, wards.staff)         
wards <- wards.m <- wards[wards$id%in%mydata$barcode,]
#set nncl and multiple as NA
wards$ward[wards$ward%in%c('nncl','multiple')] <- NA 

#Ben, Please double-check this final 'wards' to make sure the results make sense! 

#################################
#set incubation_period: a vector indexed at day = 1 
#################################
n <- rlnorm(10000, meanlog = 1.621, sdlog = 0.418)
fit.gamma <- fitdist(n, distr = "gamma", method = "mle")
# build discretized gamma
incubation_period <- distcrete::distcrete("gamma",
                              shape = fit.gamma$estimate['shape'],
                              scale = fit.gamma$estimate['rate'],
                              w = 0.5, interval = 1)
rm(n,fit.gamma)

#################################
#set generation_time: a vector indexed at day = 1 
#################################
# build discretized gamma
generation_time <- distcrete::distcrete("gamma",
                                          shape = 4.83,
                                          scale = 1.051851,
                                          w = 0.5, interval = 1)
#We can specify multiple incubation periods but we using only one
f <- matrix(
        c(incubation_period$d(1:50), incubation_period$d(1:50)), ncol = 2,
        dimnames = list(NULL, c("local", "import"))
)

#################################
#prepare for run
#################################
## define a random set of transition matrices between wards
unq <- unique(wards$ward)[!is.na(unique(wards$ward))]
n <- length(unq)
transition <- matrix(runif(n*n), n, n, dimnames = list(unq, unq))
diag(transition) <- 0
transition <- t(apply(transition, 1, function(x) x/sum(x)))
rm(unq,n)

#create meta object
meta <- mydata[,c("barcode","dateofonset_foranalysis")]
names(meta) <- c("id", "date_onset")
row.names(meta) <- NULL 

## input data
data <- outbreaker_data(
        dates = date_onset,
        dna = dna,
        ctd_timed = wards,
        ids = meta$id,
        w_dens = generation_time$d(1:50),
        f_dens = f,
        p_trans = list(transition)
)
rm(transition)

## set up outbreaker config file
config <- create_config(
        ## define iterations and sampling (this is a short exploratory run)
        n_iter = 5e3, sample_every = 50,
        ## prior on the reporting probability pi
        move_pi = TRUE, prior_pi = c(5, 5),
        ## eps is the probability of infections occuring between cases registered on the same ward
        move_eps = TRUE, prior_eps = c(1, 1),
        ## tau is the probability of an unobserved cases being moved to a different ward
        move_tau = TRUE, prior_tau = c(1, 1),
        ## leave these options as is for now
        move_joint = TRUE, move_model = TRUE
)


#################################
#find a starting point
#################################
## Identify initial tree (run this for now, it will identify a possible starting point)
#from time to time it fails but usually work after first try
config <- get_initial_tree(data, config, n_iter = 5e3, max_dist = 0)
table(data$D[cbind(seq_along(config$init_alpha), config$init_alpha)],useNA = 'always')

#################################
#run outbreaker
#################################
res <- outbreaker(data, config)
plot(res, burn=0)



#################################
#summaries
#################################

mydf<- summary(res, burnin=1000)
#match contacts
dictionary.tree <- setNames(meta$id, 1:nrow(meta)) #dictionary
mydf$tree$from <- as.vector(dictionary.tree[mydf$tree$from])
mydf$tree$to <- as.vector(dictionary.tree[mydf$tree$to])
tree <- mydf$tree
table(data$D[as.matrix(tree[c("from", "to")])], useNA = 'always')
min.date.i<- which(tree$time==min(tree$time))
if(length(min.date.i)>1) min.date.i<- min.date.i[1]
min.date <- mydata$dateofonset_foranalysis[mydata$barcode== tree$to[min.date.i]]
min.date.i<- abs(tree$time[min.date.i])
tree$time <- tree$time+min.date.i
tree$date <- tree$time+min.date;rm(min.date,min.date.i)

#filtering
tree<-tree[complete.cases(tree), ]
# tree <- tree[tree$support>0.5,]
for(i in 1:nrow(tree)){
        if(any(is.na(c(tree$from[i], tree$to[i])))){
                tree$SNP[i] <- NA
        } else {
                tree$SNP[i] <- dna.SNP[tree$from[i],tree$to[i]]
                }
}
# tree <- tree[tree$SNP<4,]
mydf$tree<- tree
#extract contacts
contacts <- mydf$tree[!is.na(mydf$tree$from),] #c("from", "to")
mydata.contacts <- mydata[mydata$barcode%in%unique(c(mydf$tree$from, mydf$tree$to)),]
#use infection date from tree to sort out ward occupation
contacts$ward.from <- contacts$ward.to <- NA
#add most likely ward to tree based on infection date
for(i in 1:nrow(contacts)){ 
        x<- contacts$from[i]
        if(is.na(mydata.contacts$ward[mydata.contacts$barcode==x])){
                if(x %in% wards$id){ 
                        wards.i <- wards.m[wards.m$id==x,]
                        wards.i$intervals <- interval(wards.i$adm, wards.i$dis)
                        if(nrow(wards.i)!=1){ 
                                j.time <- interval(contacts$date[i],contacts$date[i])
                                j.pos <- which(int_overlaps(wards.i$intervals, j.time))
                                j.ward <- ifelse(!vector.is.empty(j.pos), wards.i$ward[j.pos], NA)
                                contacts$ward.from[i] <- j.ward; rm(j.ward,j.pos,j.time)
                        } else contacts$ward.from[i] <- (wards.i$ward)
                }; rm(wards.i)
        } else {
                if(mydata.contacts$ward[mydata.contacts$barcode==x] != "multiple") {
                        contacts$ward.from[i] <- mydata.contacts$ward[mydata$barcode==x]
                } else { #if so, try to figure out the ward using the occupation ward data and infection date
                        if(x %in% wards$id){ 
                                wards.i <- wards.m[wards.m$id==x,]
                                wards.i$intervals <- interval(wards.i$adm, wards.i$dis)
                                if(nrow(wards.i)!=1){ 
                                        j.time <- interval(contacts$date[i],contacts$date[i])
                                        j.pos <- which(int_overlaps(wards.i$intervals, j.time))
                                        j.ward <- ifelse(!vector.is.empty(j.pos), wards.i$ward[j.pos], NA)
                                        contacts$ward.from[i] <- j.ward; rm(j.ward,j.pos,j.time)
                                } else contacts$ward.from[i] <- (wards.i$ward)
                        }; rm(wards.i)
                }
        }; rm(x)
};rm(i)
for(i in 1:nrow(contacts)){ 
        x<- contacts$to[i]
        if(is.na(mydata.contacts$ward[mydata.contacts$barcode==x])){
                if(x %in% wards$id){ 
                        wards.i <- wards.m[wards.m$id==x,]
                        wards.i$intervals <- interval(wards.i$adm, wards.i$dis)
                        if(nrow(wards.i)!=1){ 
                                j.time <- interval(contacts$date[i],contacts$date[i])
                                j.pos <- which(int_overlaps(wards.i$intervals, j.time))
                                j.ward <- ifelse(!vector.is.empty(j.pos), wards.i$ward[j.pos], NA)
                                contacts$ward.to[i] <- j.ward; rm(j.ward,j.pos,j.time)
                        } else contacts$ward.to[i] <- (wards.i$ward)
                }; rm(wards.i)
        } else {
                if(mydata.contacts$ward[mydata.contacts$barcode==x] != "multiple") {
                        contacts$ward.to[i] <- mydata.contacts$ward[mydata$barcode==x]
                } else { #if so, try to figure out the ward using the occupation ward data and infection date
                        if(x %in% wards$id){ 
                                wards.i <- wards.m[wards.m$id==x,]
                                wards.i$intervals <- interval(wards.i$adm, wards.i$dis)
                                if(nrow(wards.i)!=1){ 
                                        j.time <- interval(contacts$date[i],contacts$date[i])
                                        j.pos <- which(int_overlaps(wards.i$intervals, j.time))
                                        j.ward <- ifelse(!vector.is.empty(j.pos), wards.i$ward[j.pos], NA)
                                        contacts$ward.to[i] <- j.ward; rm(j.ward,j.pos,j.time)
                                } else contacts$ward.to[i] <- (wards.i$ward)
                        }; rm(wards.i)
                }
        }; rm(x)
};rm(i)

#get cluster info
cons_tree <- make_epicontacts(mydata.contacts, contacts, id = "barcode",
                              from = 1, to = 2, directed = TRUE)
# summary(cons_tree)
myclusters<-get_clusters(cons_tree, output = c("data.frame"))
myclusters.counts <- myclusters[,c("cluster_member","cluster_size")] %>% distinct()
#df to save data
dictionary.cat <- setNames(mydata$category, mydata$barcode) #dictionary
tChains.df <- tibble("id"=integer(), "size"=numeric(),"n.wards"=numeric(),
                     "n.nncl"=numeric(),"n.multiple"=numeric(),
                     "NA.wards"=numeric(), 
                     "staff_to_staff"=integer(), 'inpatient_to_inpatient'=integer(),
                     "staff_to_inpatient"=integer(), 'inpatient_to_staff'=integer(),
                     "ini"=character(), "end"=character(), "duration"=integer(),
                     "lineage"=character(), "n.lineages"=integer(),"index.category"=character(),
                     "index.association"=character(), "aveSNPpair"=integer(),"maxSNPpair"=integer(),
                     "aveGenerations"=numeric(), "firstOnset"=numeric(),
)
tChains.id <- which(myclusters.counts$cluster_size!=1)
for(tChain in tChains.id){
        i<- myclusters$id[myclusters$cluster_member==tChain] #get ids in cluster
        j<-tree[tree$from%in%i | tree$to%in%i,]
        contacts.cat  <- contacts[contacts$from%in%i | contacts$to%in%i,]
        contacts.cat$from <- as.vector(dictionary.cat[contacts.cat$from])
        contacts.cat$to <- as.vector(dictionary.cat[contacts.cat$to])
        contacts.cat$from_to <- paste(contacts.cat$from, contacts.cat$to, sep="_to_")
        #try to get wards info
        w<-c()
        for(x in i){ #check if the individual is observed in multiple wards
                if(is.na(mydata$ward[mydata$barcode==x])){
                        if(x %in% wards$id){ #"shef_cba03"
                                wards.i <- wards.m[wards.m$id==x,]
                                wards.i$intervals <- interval(wards.i$adm, wards.i$dis)
                                if(nrow(wards.i)!=1){ #check from(s) and to(s)
                                        for(myrow in 1:nrow(j)){
                                                if(j$from[myrow] == x | j$to[myrow] == x){
                                                        j.time <- interval(j$date[myrow],j$date[myrow])
                                                        j.pos <- which(int_overlaps(wards.i$intervals, j.time))
                                                        j.ward <- ifelse(!vector.is.empty(j.pos), wards.i$ward[j.pos], NA)
                                                        w <- c(w, j.ward); rm(j.ward,j.pos,j.time)
                                                }
                                        }
                                } else w <- c(w, wards.i$ward)
                        } else w <- c(w, NA); rm(wards.i)
                } else {
                        if(mydata$ward[mydata$barcode==x] != "multiple") {
                                w <- c(w, mydata$ward[mydata$barcode==x])
                        } else { #if so, try to figure out the ward using the occupation ward data and infection date
                                if(x %in% wards$id){ #"shef_cba03"
                                        wards.i <- wards.m[wards.m$id==x,]
                                        wards.i$intervals <- interval(wards.i$adm, wards.i$dis)
                                        if(nrow(wards.i)!=1){ #check from(s) and to(s)
                                                for(myrow in 1:nrow(j)){
                                                        if(j$from[myrow] == x | j$to[myrow] == x){
                                                                j.time <- interval(j$date[myrow],j$date[myrow])
                                                                j.pos <- which(int_overlaps(wards.i$intervals, j.time))
                                                                j.ward <- ifelse(!vector.is.empty(j.pos), wards.i$ward[j.pos], NA)
                                                                w <- c(w, j.ward); rm(j.ward,j.pos,j.time)
                                                        }
                                                }
                                        } else w <- c(w, wards.i$ward)
                                } else w <- c(w, NA); rm(wards.i)
                        }
                        
                }
                
                
        };rm(x)
        while(length(which(j$time==min(j$time)))>1)  j<- j[(sample(which(j$time==min(j$time)),1)*-1),]  
        tChains.line <- tibble("id"=tChain, "size"=length(mydata$ward[mydata$barcode%in%i]),
                               "n.wards"=sum(unique(w)%!in%c(NA, "nncl", "multiple")),
                               "n.nncl"=sum(w=='nncl', na.rm = T),
                               "n.multiple"=sum(w=='multiple', na.rm = T),
                               "NA.wards"=sum(is.na(w)), 
                               "staff_to_staff"=sum(contacts.cat$from_to=='staff_to_staff'), 
                               'inpatient_to_inpatient'=sum(contacts.cat$from_to=='inpatient_to_inpatient'),
                               "staff_to_inpatient"=sum(contacts.cat$from_to=='staff_to_inpatient'), 
                               'inpatient_to_staff'=sum(contacts.cat$from_to=='inpatient_to_staff'),
                               'ini'= min(wards$adm[wards$id%in%i]),
                               'end'=max(wards$adm[wards$id%in%i]),
                               'duration'=ifelse(length(i)!=2, max(wards$adm[wards$id%in%i])-min(wards$adm[wards$id%in%i]), 0),
                               'lineage'=paste(unique(mydata$lineage[mydata$barcode%in%i]), collapse = " "),
                               'n.lineages'=length(unique(mydata$lineage[mydata$barcode%in%i])),
                               "index.category"=mydata$category[mydata$barcode==j$from[which(j$time==min(j$time))]],
                               "index.association"=mydata$healthcareassociation[mydata$barcode==j$from[which(j$time==min(j$time))]],
                               "aveSNPpair"=mean(j$SNP), "maxSNPpair"=max(j$SNP), "aveGenerations"=mean(j$generations),
                               "firstOnset"=min(j$date)
                                   ); rm(i,contacts.cat, j)
        tChains.df <- rbind(tChains.df, tChains.line); rm(tChains.line)
}; tChains.df$aveSNPpair <- round (tChains.df$aveSNPpair,2)






