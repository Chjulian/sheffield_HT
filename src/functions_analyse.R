add_meta <- function(run, mydata, data){
        run$from <- mydata$barcode[run$from]
        run$to <- mydata$barcode[run$to]
        run$ward_from <- ifelse(run$ward_from==-1, NA, run$ward_from)
        run$ward_to <- ifelse(run$ward_from==-1, NA, run$ward_to)
        run$ward_from <- unique(data$ctd_timed$ward)[run$ward_from]
        run$ward_to <- unique(data$ctd_timed$ward)[run$ward_to]
        run$same_ward <- ifelse(run$ward_from==run$ward_to, TRUE, FALSE)
        run$from.cat <- mydata$category[match(run$from, mydata$barcode)]
        run$to.cat <- mydata$category[match(run$to, mydata$barcode)]
        run$trans.cat <- paste0(run$from.cat, "-", run$to.cat)
        run$from.hosp.assoc <- mydata$healthcareassociation[match(run$from, mydata$barcode)]
        run$to.hosp.assoc <- mydata$healthcareassociation[match(run$to, mydata$barcode)]
        return(run)
}


check_pairs <- function(run){
        fun.df <- run
        fun.df$ward_snp <- paste0(fun.df$same_ward, "-SNP", fun.df$dist)
        sum.table <- prop.table(table(fun.df$ward_snp))
        out.tib <- tibble(
                "SameWard-0SNP" = round(sum.table["TRUE-SNP0"], 2),
                "SameWard-1SNP" = round(sum.table["TRUE-SNP1"], 2),
                "SameWard-2SNP" = round(sum.table["TRUE-SNP2"], 2),
                "SameWard-3SNP" = round(sum.table["TRUE-SNP3"], 2),
                "Other" = round(1 - (sum.table["TRUE-SNP0"] +
                                       sum.table["TRUE-SNP1"]+
                                       sum.table["TRUE-SNP2"] +
                                       sum.table["TRUE-SNP3"]),2)
        ) 
        return(out.tib)
}



analyse_run <- function(run, mydata, baylink.pairs){
        network.all <- run[!is.na(run$from),]
        total.cases <- length(unique(network.all$to))
        network.direct <- network.all[network.all$kappa==1,]
        direct.pairs <- length(unique(network.direct$to))
        pair.cat.prop <- prop.table(table(network.direct$trans.cat))
        pair.cat.n <- table(network.direct$trans.cat)
        contacts <- network.all[,c("from", "to")]
        known.infectors <- unique(contacts$to)
        all_infector_infectees <- unique(c(as.character(contacts$from), as.character(contacts$to)))
        tree <- make_epicontacts(mydata[mydata$barcode %in% all_infector_infectees,], contacts, id = "barcode",
                                 from = 1, to = 2, directed = TRUE)
        tree <- get_clusters(tree, output = "epicontacts")
        network.all$cluster <- tree$linelist$cluster_member[match(network.all$to, tree$linelist$id)]
        cluster.num <- length(unique(tree$linelist$cluster_member))
        cluster.index <- sapply(unique(tree$linelist$cluster_member),
                                function(cluster){
                                        cluster.df <- network.all[network.all$cluster==cluster,]
                                        cluster.df <- cluster.df %>% arrange(t_inf)
                                        cluster.index <- cluster.df$from.cat[1]
                                        return(cluster.index)
                                }
        ) 
        
        
        cluster.index.prop <- prop.table(table(cluster.index))
        cluster.index.n <- table(cluster.index)
        cluster.duration <- sapply(unique(tree$linelist$cluster_member),
                                   function(cluster){
                                           cluster.df <- tree$linelist[tree$linelist$cluster_member==cluster,]
                                           cluster.df <- cluster.df %>% arrange(dateofonset_foranalysis)
                                           cluster.time <- cluster.df$dateofonset_foranalysis[nrow(cluster.df)] -
                                                   cluster.df$dateofonset_foranalysis[1]
                                           return(cluster.time)
                                   }
        )
        cluster.wards.n <- sapply(unique(network.all$cluster),
                                   function(cluster){
                                           cluster.df <- network.all[network.all$cluster==cluster,]
                                           wards <- unique(c(cluster.df$ward_from, cluster.df$ward_to))
                                           ward.n <- length(wards)
                                           return(ward.n)
                                   }
        )
        cluster.size <- sapply(unique(tree$linelist$cluster_member),
                               function(cluster){
                                       size <- length(tree$linelist$cluster_member[
                                               tree$linelist$cluster_member==cluster
                                       ])
                               })
        
        cluster.tib <- tibble(
                cluster = unique(tree$linelist$cluster_member),
                size = cluster.size,
                duration = cluster.duration
        )
        secondary.infections <- lapply(unique(network.all$from),
                                       function(infector){
                                               sec.infections <- length(network.all$from[network.all$from==infector])
                                               infector.tib <- tibble(barcode = infector,
                                                                      secondary.infections = sec.infections)
                                               return(infector.tib)
                                       }
        ) %>% bind_rows()
        no.infectee <- mydata$barcode[mydata$barcode %!in% secondary.infections$barcode]
        secondary.infections <- rbind(
                secondary.infections,
                tibble(
                        barcode = no.infectee,
                        secondary.infections = 0
                )
        )
        secondary.infections <- left_join(secondary.infections,
                                          select(mydata, "barcode", "category",
                                                 "dateofonset_foranalysis"))
        
        network.direct$pair <- paste0(network.direct$from, "_", network.direct$to)
        network.direct$bay.link <- ifelse(network.direct$pair %in% baylink.pairs, TRUE, FALSE)
        network.direct.inpt.links <- network.direct[network.direct$trans.cat=="inpatient-inpatient",]
        baylink.table <- prop.table(table(network.direct.inpt.links$bay.link))

        out.tib <-  tibble(
                total.cases = length(mydata$barcode),
                linked.infections.n = length(known.infectors),
                linked.infections.prop = length(known.infectors)/length(mydata$barcode),
                linked.infections.staff.n = length(mydata$barcode[mydata$barcode %in% known.infectors & mydata$category=="staff"]),
                linked.infections.pt.n = length(mydata$barcode[mydata$barcode %in% known.infectors & mydata$category=="inpatient"]),
                direct.links.n = direct.pairs,
                inpatient.inpatient.prop = pair.cat.prop["inpatient-inpatient"],
                inpatient.staff.prop = pair.cat.prop["inpatient-staff"],
                staff.inpatient.prop = pair.cat.prop["staff-inpatient"],
                staff.staff.prop = pair.cat.prop["staff-staff"],
                inpatient.inpatient.n = pair.cat.n["inpatient-inpatient"],
                inpatient.staff.n = pair.cat.n["inpatient-staff"],
                staff.inpatient.n = pair.cat.n["staff-inpatient"],
                staff.staff.n = pair.cat.n["staff-staff"],
                total.clusters = cluster.num,
                max.cluster.size = max(cluster.size),
                mean.cluster.size = mean(cluster.size),
                median.cluster.size = median(cluster.size), 
                longest.cluster = max(cluster.duration),
                mean.cluster.duration = mean(cluster.duration),
                median.cluster.duration = median(cluster.duration),
                wards.per.cluster.mean = mean(cluster.wards.n),
                wards.per.cluster.median = median(cluster.wards.n),
                wards.per.cluster.max = max(cluster.wards.n),
                index.patient.prop = cluster.index.prop["inpatient"],
                index.staff.prop = cluster.index.prop["staff"],
                index.patient.n = cluster.index.n["inpatient"],
                index.staff.n = cluster.index.n["staff"],
                sec.inf.all.mean = mean(secondary.infections$secondary.infections),
                sec.inf.staff.mean = mean(secondary.infections$secondary.infections[secondary.infections$category=="staff"]),
                sec.inf.inpatients.mean = mean(secondary.infections$secondary.infections[secondary.infections$category=="inpatient"]),
                sec.inf.all.median = median(secondary.infections$secondary.infections),
                sec.inf.staff.median = median(secondary.infections$secondary.infections[secondary.infections$category=="staff"]),
                sec.inf.inpatients.median = median(secondary.infections$secondary.infections[secondary.infections$category=="inpatient"]),
                sec.inf.all.max = max(secondary.infections$secondary.infections),
                sec.inf.staff.max = max(secondary.infections$secondary.infections[secondary.infections$category=="staff"]),
                sec.inf.inpatients.max = max(secondary.infections$secondary.infections[secondary.infections$category=="inpatient"]),
                inpatient.baylink.prop = baylink.table["TRUE"]
                )
        
        return(out.tib)
}
