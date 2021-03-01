pacman::p_load(ape, linelist, epitrix, fitdistrplus,
               outbreaker2, tidyverse, rio, 
               magrittr, remotes, epicontacts, 
               visNetwork, lubridate, pairsnp,
               ggplot2, chron, tidyverse, glue) 
source("src/functions_wards.R")
source("src/functions.R")
source("src/onset_distribution.R")

get_tpairs.2 <- function(res, data) {
        lapply(
                seq_len(sum(grepl("alpha", names(res)))),
                function(i) {
                        tibble(
                                from = res[[glue("alpha_{i}")]],
                                to = i,
                                step = res[[glue("step")]],
                                t_inf = res[[glue("t_inf_{i}")]],
                                t_onw = res[[glue("t_onw_{i}")]],
                                kappa = res[[glue("kappa_{i}")]],
                                ward_from = mapply(get_ward, from, ifelse(is.na(kappa) | kappa == 1, t_inf, t_onw), MoreArgs = list(data = data)),
                                ward_to = mapply(get_ward, to, t_inf, MoreArgs = list(data = data)),
                                dist = data$D[cbind(from, to)]
                        )
                }
        ) %>%
                bind_rows() %>%
                drop_na(from)
}




## get data
data <- d$data
config <- d$config
data <- add_convolutions(data, config)
res <- d$res
raw <- d$raw

##create dictionaries
dic.id <- setNames(data$ids, 1:length(data$ids)) 
dic.ward <- setNames(unique(data$ctd_timed$ward), 1:length(unique(data$ctd_timed$ward))) 
dic.cat <- setNames(raw$category, raw$barcode) 



## check proper burning by vis
plot(res, burn=2500)
res <- res[51:nrow(res),]

## get param object
param <- create_param(data, config)$current

## update param object
param <- update_par(param, nrow(res), res, data, config)

## extract tpair data with burnin
tpairs <- get_tpairs.2(res, data)

tpairs$pair <- paste(tpairs$from, tpairs$to, sep="_")

##translate using dictionaries
tpairs$from.id <- as.vector(dic.id[tpairs$from])
tpairs$to.id <- as.vector(dic.id[tpairs$to])
tpairs$ward_from.id <- as.vector(dic.ward[tpairs$ward_from])
tpairs$ward_to.id <- as.vector(dic.ward[tpairs$ward_to])
tpairs$from.cat <- as.vector(dic.cat[tpairs$from.id])
tpairs$to.cat <- as.vector(dic.cat[tpairs$to.id])

tpairs$pair.cat <- paste(tpairs$from.cat, tpairs$to.cat, sep="_")

#summarise data
#short summary by pair
df <- tpairs %>%
        group_by(pair) %>%
        summarise(
                'n' = n(),
                'p' = n()/nrow(res),
                'n.t_inf'=length(unique(t_inf)),
                'n.kappa'=length(unique(kappa)),
                'n.ward_from'=length(unique(ward_from)),
                'n.ward_to'=length(unique(ward_to)),
                'n.dist'=length(unique(dist)),
                'ward_from'= as.numeric(names(sort(table(ward_from), decreasing = TRUE)[1])),
                'ward_to'= as.numeric(names(sort(table(ward_to), decreasing = TRUE)[1]))
                )
df

myseq <- seq(0, 1, 0.01)
sameWard <- numeric()
UniquePairs <- numeric()
for(i in myseq){
        df.t<- df[df$p>=i,]
        sameWard<- c(sameWard, round(as.numeric(table(df.t$ward_from==df.t$ward_to)[2])*100/nrow(df.t),2))
        UniquePairs <- c(UniquePairs, nrow(df.t))
}; rm(i,df.t)

df.t <- data.frame('p'=myseq,'sameWard'=sameWard, 'UniquePairs'=UniquePairs)

require(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, max(df.t$UniquePairs)))

p1<-ggplot(df.t)+geom_point(aes(x=p, y=sameWard, color=UniquePairs)) + theme_minimal() +
        ylab('% linked pairs--same ward') + xlab('pp of linkage in outbreaker') +
        ylim(0,100) + sc

p2<-ggplot(df.t[df.t$p>0.25,])+geom_point(aes(x=p, y=UniquePairs, color=UniquePairs)) + theme_minimal() +
        ylab('# linked pairs') + xlab('pp of linkage in outbreaker') + xlim(0,1) +
        sc
ggpubr::ggarrange(p1,p2, ncol = 1, common.legend = T, legend = c('right'))
        
#summary of category by mcmc step
df.cat <- tpairs %>%
        group_by(step) %>% 
        summarise(
                'inpatient_inpatient' = sum(pair.cat=='inpatient_inpatient'),
                'inpatient_staff' = sum(pair.cat=='inpatient_staff'),
                'staff_inpatient' = sum(pair.cat=='staff_inpatient'),
                'staff_staff' = sum(pair.cat=='staff_staff') 
                ) %>%
        mutate('total' = select(., 2:5) %>% rowSums(na.rm = TRUE)) %>%
        mutate('inpatient_inpatient%'=inpatient_inpatient*100/total,
               'inpatient_staff%'=inpatient_staff*100/total,
               'staff_inpatient%'=staff_inpatient*100/total,
               'staff_staff%'=staff_staff*100/total)

df.cat.long <- df.cat %>%
        select('inpatient_inpatient%', 'inpatient_staff%',
                'staff_inpatient%', 'staff_staff%') %>%
        pivot_longer(cols = c('inpatient_inpatient%', 'inpatient_staff%',
                                             'staff_inpatient%', 'staff_staff%'))

require(ggplot2)
ggplot(df.cat.long) + geom_boxplot(aes(x=name, y=value, color=name)) + 
        theme_minimal()

