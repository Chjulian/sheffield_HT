#convert indexed and scale data back to original metadata and unscaled dates
add_meta <-
        function(run,
                 mydata,
                 data,
                 scaled_date_onset,
                 origin,
                 scale) {
                run$from <- mydata$barcode[run$from]
                run$to <- mydata$barcode[run$to]
                run$ward_from <-
                        ifelse(run$ward_from == -1, NA, run$ward_from)
                run$ward_to <- ifelse(run$ward_from == -1, NA, run$ward_to)
                run$ward_from <- unique(data$ctd_timed$ward)[run$ward_from]
                run$ward_to <- unique(data$ctd_timed$ward)[run$ward_to]
                run$same_ward <-
                        ifelse(run$ward_from == run$ward_to, TRUE, FALSE)
                run$from.cat <-
                        mydata$category[match(run$from, mydata$barcode)]
                run$to.cat <- mydata$category[match(run$to, mydata$barcode)]
                run$trans.cat <- paste0(run$from.cat, "-", run$to.cat)
                run$from.hosp.assoc <-
                        mydata$healthcareassociation[match(run$from, mydata$barcode)]
                run$to.hosp.assoc <-
                        mydata$healthcareassociation[match(run$to, mydata$barcode)]
                run$t_inf <-
                        as.Date(run$t_inf, origin = min(scaled_date_onset))
                run$t_inf <-
                        scale_onset_dates(run$t_inf, origin, to = "down", scale = scale)
                run$from.onset <-
                        mydata$dateofonset_foranalysis[match(run$from, mydata$barcode)]
                run$to.onset <-
                        mydata$dateofonset_foranalysis[match(run$to, mydata$barcode)]
                run$incubation <-
                        difftime(run$to.onset, run$t_inf, units = "days")
                run$serial.interval <-
                        difftime(run$to.onset, run$from.onset,  units = "days")
                return(run)
        }

#generate outcomes of interest from outbreaker run output
analyse_run <- function(run, mydata, baylink.pairs, wave) {
        network.all <- run[!is.na(run$from), ]
        network.all$unsampled.case <-
                ifelse(!is.na(network.all$kappa),
                       network.all$kappa - 1,
                       network.all$kappa)
        total.cases <- length(unique(network.all$to))
        network.direct <- network.all[network.all$kappa == 1, ]
        direct.pairs <- length(unique(network.direct$to))
        pair.cat.prop <- prop.table(table(network.direct$trans.cat))
        pair.cat.n <- table(network.direct$trans.cat)
        contacts <- network.all[, c("from", "to")]
        known.infectors <- unique(contacts$to)
        all_infector_infectees <-
                unique(c(
                        as.character(contacts$from),
                        as.character(contacts$to)
                ))
        tree <-
                make_epicontacts(
                        mydata[mydata$barcode %in% all_infector_infectees, ],
                        contacts,
                        id = "barcode",
                        from = 1,
                        to = 2,
                        directed = TRUE
                )
        tree <- get_clusters(tree, output = "epicontacts")
        network.all$cluster <-
                tree$linelist$cluster_member[match(network.all$to, tree$linelist$id)]
        cluster.num <- length(unique(tree$linelist$cluster_member))
        cluster.index <-
                sapply(unique(tree$linelist$cluster_member),
                       function(cluster) {
                               cluster.df <- network.all[network.all$cluster == cluster, ]
                               cluster.df <-
                                       cluster.df %>% arrange(t_inf)
                               cluster.index <-
                                       cluster.df$from.cat[1]
                               return(cluster.index)
                       })
        
        cluster.index.hoci <-
                sapply(unique(tree$linelist$cluster_member),
                       function(cluster) {
                               cluster.df <- network.all[network.all$cluster == cluster, ]
                               cluster.df <-
                                       cluster.df %>% arrange(t_inf)
                               cluster.index <-
                                       cluster.df$from.hosp.assoc[1]
                               return(cluster.index)
                       })
        
        cluster.index.prop <- prop.table(table(cluster.index))
        cluster.index.n <- table(cluster.index)
        cluster.index.hoci.prop <-
                prop.table(table(cluster.index.hoci))
        cluster.index.hoci.n <- table(cluster.index.hoci)
        cluster.duration <-
                sapply(unique(tree$linelist$cluster_member),
                       function(cluster) {
                               cluster.df <- tree$linelist[tree$linelist$cluster_member == cluster, ]
                               cluster.df <-
                                       cluster.df %>% arrange(dateofonset_foranalysis)
                               cluster.time <-
                                       cluster.df$dateofonset_foranalysis[nrow(cluster.df)] -
                                       cluster.df$dateofonset_foranalysis[1]
                               return(cluster.time)
                       })
        cluster.wards.n <- sapply(unique(network.all$cluster),
                                  function(cluster) {
                                          cluster.df <- network.all[network.all$cluster == cluster, ]
                                          wards <-
                                                  unique(c(cluster.df$ward_from, cluster.df$ward_to))
                                          ward.n <- length(wards)
                                          return(ward.n)
                                  })
        cluster.size <- sapply(unique(tree$linelist$cluster_member),
                               function(cluster) {
                                       size <- length(tree$linelist$cluster_member[tree$linelist$cluster_member ==
                                                                                           cluster]) +
                                               sum(network.all$unsampled.case[network.all$cluster ==
                                                                                      cluster])
                                       
                                       
                               })
        cluster.quantile <- quantile(cluster.size)
        cluster.duration.quantile <- quantile(cluster.duration)
        wards.cluster.quantile <- quantile(cluster.wards.n)
        
        secondary.infections <- lapply(unique(network.all$from),
                                       function(infector) {
                                               sec.infections <-
                                                       length(network.all$from[network.all$from == infector])
                                               infector.tib <-
                                                       tibble(barcode = infector,
                                                              secondary.infections = sec.infections)
                                               return(infector.tib)
                                       }) %>% bind_rows()
        no.infectee <-
                mydata$barcode[mydata$barcode %!in% secondary.infections$barcode]
        secondary.infections <- rbind(
                secondary.infections,
                tibble(barcode = no.infectee,
                       secondary.infections = 0)
        )
        secondary.infections <- left_join(
                secondary.infections,
                dplyr:::select(
                        mydata,
                        "barcode",
                        "category",
                        "healthcareassociation",
                        "dateofonset_foranalysis"
                ),
                by = "barcode"
        )
        sec.inf.all.quantile <-
                quantile(secondary.infections$secondary.infections)
        sec.inf.staff.quantile <-
                quantile(secondary.infections$secondary.infections[secondary.infections$category ==
                                                                           "staff"])
        sec.inf.inpatients.quantile <-
                quantile(secondary.infections$secondary.infections[secondary.infections$category ==
                                                                           "inpatient"])
        sec.inf.coca.quantile <-
                quantile(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                           "community_onset_community_associated"])
        sec.inf.cosha.quantile <-
                quantile(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                           "community_onset_suspected_healthcare_associated"])
        sec.inf.hoha.quantile <-
                quantile(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                           "hospital_onset_healthcare_associated"])
        sec.inf.hoiha.quantile <-
                quantile(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                           "hospital_onset_indetermite_healthcare_associated"])
        sec.inf.hosha.quantile <-
                quantile(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                           "hospital_onset_suspected_healthcare_associated"])
        
        network.direct$pair <-
                paste0(network.direct$from, "_", network.direct$to)
        network.direct$bay.link <-
                ifelse(network.direct$pair %in% baylink.pairs, TRUE, FALSE)
        network.direct.inpt.links <-
                network.direct[network.direct$trans.cat == "inpatient-inpatient", ]
        baylink.table <-
                prop.table(table(network.direct.inpt.links$bay.link))
        cosha.prop <-
                length(network.all$to[network.all$to.hosp.assoc == "community_onset_suspected_healthcare_associated"]) /
                length(mydata$barcode[mydata$healthcareassociation == "community_onset_suspected_healthcare_associated"])
        hoha.prop <-
                length(network.all$to[network.all$to.hosp.assoc == "hospital_onset_healthcare_associated"]) /
                length(mydata$barcode[mydata$healthcareassociation == "hospital_onset_healthcare_associated"])
        hosha.prop <-
                length(network.all$to[network.all$to.hosp.assoc == "hospital_onset_suspected_healthcare_associated"]) /
                length(mydata$barcode[mydata$healthcareassociation == "hospital_onset_suspected_healthcare_associated"])
        hoiha.prop <-
                length(network.all$to[network.all$to.hosp.assoc == "hospital_onset_indetermite_healthcare_associated"]) /
                length(mydata$barcode[mydata$healthcareassociation == "hospital_onset_indetermite_healthcare_associated"])
        
        top10.wards <-
                c(
                        "loc0124",
                        "loc0121",
                        "loc0111",
                        "loc0147",
                        "loc0113",
                        "loc0085",
                        "loc0157",
                        "loc0182",
                        "loc0112",
                        "loc0082"
                )
        
        if (wave == "wave2") {
                top10.wards <-
                        c(
                                "loc0080",
                                "loc0123",
                                "loc0200",
                                "loc0116",
                                "loc0088",
                                "loc0113",
                                "loc0032",
                                "loc0124",
                                "loc0112",
                                "loc0033"
                        )
        }
        
        ward.outcomes <- lapply(top10.wards, function(i) {
                loop.df <- network.all[network.all$ward_from == i, ]
                n <- nrow(loop.df)
                n.clusters <- length(unique(loop.df$cluster))
                n.staff <-
                        length(loop.df$to[loop.df$to.cat == "staff"])
                n.pt <-
                        length(loop.df$to[loop.df$to.cat == "inpatient"])
                ward.tib <- tibble(
                        ward =  i,
                        ward.transmission.n = n,
                        ward.clusters.n = n.clusters,
                        staff.n = n.staff,
                        patient.n = n.pt,
                        staff.prop = n.staff / n,
                        patient.prop = n.pt / n
                )
                
                return(ward.tib)
        }) %>% bind_rows()
        
        wards.transmissions <-
                spread(
                        select(ward.outcomes, "ward", "ward.transmission.n"),
                        ward,
                        ward.transmission.n
                )
        names(wards.transmissions) <-
                paste0(names(wards.transmissions), ".transmissions.n")
        wards.clusters <-
                spread(select(ward.outcomes, "ward", "ward.clusters.n"),
                       ward,
                       ward.clusters.n)
        names(wards.clusters) <-
                paste0(names(wards.clusters), ".clusters.n")
        wards.staff.n <-
                spread(select(ward.outcomes, "ward", "staff.n"),
                       ward,
                       staff.n)
        names(wards.staff.n) <-
                paste0(names(wards.staff.n), ".staff.n")
        wards.pt.n <-
                spread(select(ward.outcomes, "ward", "patient.n"),
                       ward,
                       patient.n)
        names(wards.pt.n) <- paste0(names(wards.pt.n), ".patient.n")
        wards.staff.prop <-
                spread(select(ward.outcomes, "ward", "staff.prop"),
                       ward,
                       staff.prop)
        names(wards.staff.prop) <-
                paste0(names(wards.staff.prop), ".staff.prop")
        wards.pt.prop <-
                spread(select(ward.outcomes, "ward", "patient.prop"),
                       ward,
                       patient.prop)
        names(wards.pt.prop) <-
                paste0(names(wards.pt.prop), ".patient.prop")
        
        
        out.tib <-  tibble(
                total.cases = length(mydata$barcode),
                linked.infections.n = length(known.infectors),
                unsampled.linked.cases = sum(network.all$unsampled.case),
                linked.infections.prop = length(known.infectors) / length(mydata$barcode),
                linked.infections.staff.n = length(mydata$barcode[mydata$barcode %in% known.infectors &
                                                                          mydata$category == "staff"]),
                linked.infections.pt.n = length(mydata$barcode[mydata$barcode %in% known.infectors &
                                                                       mydata$category == "inpatient"]),
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
                min.cluster.size = min(cluster.size),
                mean.cluster.size = mean(cluster.size),
                sd.cluster.size = sd(cluster.size),
                median.cluster.size = median(cluster.size),
                cluster.quantile.25 = cluster.quantile["25%"],
                cluster.quantile.75 = cluster.quantile["75%"],
                longest.cluster = max(cluster.duration),
                shortest.cluster = min(cluster.duration),
                mean.cluster.duration = mean(cluster.duration),
                sd.cluster.duration = sd(cluster.duration),
                median.cluster.duration = median(cluster.duration),
                cluster.duration.quantile.25 = cluster.duration.quantile["25%"],
                cluster.duration.quantile.75 = cluster.duration.quantile["75%"],
                wards.per.cluster.max = max(cluster.wards.n),
                wards.per.cluster.min = min(cluster.wards.n),
                wards.per.cluster.mean = mean(cluster.wards.n),
                wards.per.cluster.sd = sd(cluster.wards.n),
                wards.per.cluster.median = median(cluster.wards.n),
                wards.cluster.quantile.25 = wards.cluster.quantile["25%"],
                wards.cluster.quantile.75 = wards.cluster.quantile["75%"],
                index.patient.prop = cluster.index.prop["inpatient"],
                index.staff.prop = cluster.index.prop["staff"],
                index.patient.n = cluster.index.n["inpatient"],
                index.staff.n = cluster.index.n["staff"],
                index.coca.prop = cluster.index.hoci.prop["community_onset_community_associated"],
                index.coca.n = cluster.index.hoci.n["community_onset_community_associated"],
                index.cosha.prop = cluster.index.hoci.prop["community_onset_suspected_healthcare_associated"],
                index.cosha.n = cluster.index.hoci.n["community_onset_suspected_healthcare_associated"],
                index.hoha.prop = cluster.index.hoci.prop["hospital_onset_healthcare_associated"],
                index.hoha.n = cluster.index.hoci.n["hospital_onset_healthcare_associated"],
                index.hoiha.prop = cluster.index.hoci.prop["hospital_onset_indetermite_healthcare_associated"],
                index.hoiha.n = cluster.index.hoci.n["hospital_onset_indetermite_healthcare_associated"],
                index.hosha.prop = cluster.index.hoci.prop["hospital_onset_suspected_healthcare_associated"],
                index.hosha.n = cluster.index.hoci.n["hospital_onset_suspected_healthcare_associated"],
                sec.inf.all.mean = mean(secondary.infections$secondary.infections),
                sec.inf.staff.mean = mean(secondary.infections$secondary.infections[secondary.infections$category ==
                                                                                            "staff"]),
                sec.inf.inpatients.mean = mean(secondary.infections$secondary.infections[secondary.infections$category ==
                                                                                                 "inpatient"]),
                sec.inf.coca.mean = mean(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                           "community_onset_community_associated"]),
                sec.inf.cosha.mean = mean(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                            "community_onset_suspected_healthcare_associated"]),
                sec.inf.hoha.mean = mean(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                           "hospital_onset_healthcare_associated"]),
                sec.inf.hoiha.mean = mean(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                            "hospital_onset_indetermite_healthcare_associated"]),
                sec.inf.hosha.mean = mean(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                            "hospital_onset_suspected_healthcare_associated"]),
                sec.inf.all.sd = sd(secondary.infections$secondary.infections),
                sec.inf.staff.sd = sd(secondary.infections$secondary.infections[secondary.infections$category ==
                                                                                        "staff"]),
                sec.inf.inpatients.sd = sd(secondary.infections$secondary.infections[secondary.infections$category ==
                                                                                             "inpatient"]),
                sec.inf.coca.sd = sd(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                       "community_onset_community_associated"]),
                sec.inf.cosha.sd = sd(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                        "community_onset_suspected_healthcare_associated"]),
                sec.inf.hoha.sd = sd(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                       "hospital_onset_healthcare_associated"]),
                sec.inf.hoiha.sd = sd(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                        "hospital_onset_indetermite_healthcare_associated"]),
                sec.inf.hosha.sd = sd(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                        "hospital_onset_suspected_healthcare_associated"]),
                sec.inf.all.median = median(secondary.infections$secondary.infections),
                sec.inf.staff.median = median(secondary.infections$secondary.infections[secondary.infections$category ==
                                                                                                "staff"]),
                sec.inf.inpatients.median = median(secondary.infections$secondary.infections[secondary.infections$category ==
                                                                                                     "inpatient"]),
                sec.inf.coca.median = median(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                               "community_onset_community_associated"]),
                sec.inf.cosha.median = median(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                                "community_onset_suspected_healthcare_associated"]),
                sec.inf.hoha.median = median(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                               "hospital_onset_healthcare_associated"]),
                sec.inf.hoiha.median = median(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                                "hospital_onset_indetermite_healthcare_associated"]),
                sec.inf.hosha.median = median(secondary.infections$secondary.infections[secondary.infections$healthcareassociation ==
                                                                                                "hospital_onset_suspected_healthcare_associated"]),
                sec.inf.all.quantile.25 = sec.inf.all.quantile["25%"],
                sec.inf.all.quantile.75 = sec.inf.all.quantile["75%"],
                sec.inf.staff.quantile.25 = sec.inf.staff.quantile["25%"],
                sec.inf.staff.quantile.75 = sec.inf.staff.quantile["75%"],
                sec.inf.inpatients.quantile.25 = sec.inf.inpatients.quantile["25%"],
                sec.inf.inpatients.quantile.75 = sec.inf.inpatients.quantile["75%"],
                sec.inf.coca.quantile.25 = sec.inf.coca.quantile["25%"],
                sec.inf.coca.quantile.75 = sec.inf.coca.quantile["75%"],
                sec.inf.cosha.quantile.25 = sec.inf.cosha.quantile["25%"],
                sec.inf.cosha.quantile.75 = sec.inf.cosha.quantile["75%"],
                sec.inf.hoha.quantile.25 = sec.inf.hoha.quantile["25%"],
                sec.inf.hoha.quantile.75 = sec.inf.hoha.quantile["75%"],
                sec.inf.hoiha.quantile.25 = sec.inf.hoiha.quantile["25%"],
                sec.inf.hoiha.quantile.75 = sec.inf.hoiha.quantile["75%"],
                sec.inf.hosha.quantile.25 = sec.inf.hosha.quantile["25%"],
                sec.inf.hosha.quantile.75 = sec.inf.hosha.quantile["75%"],
                sec.inf.all.max = max(secondary.infections$secondary.infections),
                sec.inf.staff.max = max(secondary.infections$secondary.infections[secondary.infections$category ==
                                                                                          "staff"]),
                sec.inf.inpatients.max = max(secondary.infections$secondary.infections[secondary.infections$category ==
                                                                                               "inpatient"]),
                #sec.inf.coca.max = max(secondary.infections$secondary.infections[secondary.infections$healthcareassociation=="community_onset_community_associated"]),
                #sec.inf.cosha.max = max(secondary.infections$secondary.infections[secondary.infections$healthcareassociation=="community_onset_suspected_healthcare_associated"]),
                #sec.inf.hoha.max = max(secondary.infections$secondary.infections[secondary.infections$healthcareassociation=="hospital_onset_healthcare_associated"]),
                #sec.inf.hoiha.max = max(secondary.infections$secondary.infections[secondary.infections$healthcareassociation=="hospital_onset_indetermite_healthcare_associated"]),
                #sec.inf.hosha.max = max(secondary.infections$secondary.infections[secondary.infections$healthcareassociation=="hospital_onset_suspected_healthcare_associated"]),
                sec.inf.all.min = min(secondary.infections$secondary.infections),
                sec.inf.staff.min = min(secondary.infections$secondary.infections[secondary.infections$category ==
                                                                                          "staff"]),
                sec.inf.inpatients.min = min(secondary.infections$secondary.infections[secondary.infections$category ==
                                                                                               "inpatient"]),
                inpatient.baylink.prop = baylink.table["TRUE"],
                cosha.prop = cosha.prop,
                hoha.prop = hoha.prop,
                hosha.prop = hosha.prop,
                hoiha.prop = hoiha.prop,
                sec.inf.0.n = length(secondary.infections$barcode[secondary.infections$secondary.infections ==
                                                                          0]),
                sec.inf.1.n = length(secondary.infections$barcode[secondary.infections$secondary.infections ==
                                                                          1]),
                sec.inf.2.n = length(secondary.infections$barcode[secondary.infections$secondary.infections ==
                                                                          2]),
                sec.inf.3.n = length(secondary.infections$barcode[secondary.infections$secondary.infections ==
                                                                          3]),
                sec.inf.4.n = length(secondary.infections$barcode[secondary.infections$secondary.infections ==
                                                                          4]),
                sec.inf.5.n = length(secondary.infections$barcode[secondary.infections$secondary.infections ==
                                                                          5]),
                sec.inf.6.n = length(secondary.infections$barcode[secondary.infections$secondary.infections ==
                                                                          6]),
                sec.inf.7.n = length(secondary.infections$barcode[secondary.infections$secondary.infections ==
                                                                          7]),
                sec.inf.8.n = length(secondary.infections$barcode[secondary.infections$secondary.infections ==
                                                                          8]),
                sec.inf.9.n = length(secondary.infections$barcode[secondary.infections$secondary.infections ==
                                                                          9]),
                sec.inf.10.n = length(secondary.infections$barcode[secondary.infections$secondary.infections ==
                                                                           10])
        )
        
        out.tib <- cbind(
                out.tib,
                wards.transmissions,
                wards.clusters,
                wards.staff.n,
                wards.staff.prop,
                wards.pt.n,
                wards.pt.prop
        )
        return(out.tib)
}
