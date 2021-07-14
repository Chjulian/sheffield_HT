'%!in%' <- function(x, y)
        ! ('%in%'(x, y))

vector.is.empty <- function(x)
        return(length(x) == 0)

##convert locations variables in the input metadata file into the correct format for outbreaker
ward_occupation <- function(wards.patients = wards.patients) {
        #create wards object for staff
        mydata$wards <-
                mydata$ward <- paste(
                        mydata$loc1,
                        mydata$loc2,
                        mydata$loc3,
                        mydata$loc4,
                        mydata$loc5,
                        mydata$loc6,
                        mydata$loc7,
                        mydata$loc8
                )
        mydata$ward.count <-
                as.integer(sapply(mydata$ward, function(rw)
                        length(strsplit(
                                trimws(rw), " "
                        )[[1]])))
        mydata$ward <- trimws(mydata$ward)
        staff.i <- which(mydata$category == 'staff')
        for (i in 1:nrow(mydata)) {
                if (mydata$ward.count[i] == 0) {
                        mydata$ward[i] <- NA
                } else {
                        if (i %in% staff.i) {
                                #if staff
                                if (mydata$ward.count[i] == 1) {
                                        #if one entry in ward occupation
                                        mydata$ward[i] <-
                                                ifelse(
                                                        mydata$loc1[i] == 'multiple',
                                                        'multiple',
                                                        mydata$ward[i]
                                                )
                                } else{
                                        #multiples entries in ward occupation
                                        occupation <-
                                                unlist(strsplit(
                                                        str_replace_all(
                                                                mydata$ward[i],
                                                                "multiple",
                                                                ""
                                                        ),
                                                        " "
                                                ))
                                        cap <-
                                                ifelse(
                                                        length(occupation) < 3,
                                                        length(occupation),
                                                        3
                                                )
                                        occupation <-
                                                occupation[1:cap]
                                        mydata$ward[i] <-
                                                sample(occupation, 1)
                                }
                                
                        } else
                                mydata$ward[i] <-
                                        ifelse(
                                                mydata$ward.count[i] > 1 |
                                                        mydata$loc1[i] == 'multiple',
                                                'multiple',
                                                mydata$ward[i]
                                        )
                        
                        
                }
        }
        rm(i)
        wards.staff <-
                mydata[mydata$category == 'staff', c('barcode',
                                                     'ward',
                                                     'dateofcollection',
                                                     'dateofonset_foranalysis')]
        
        wards.staff <-
                clean_data(wards.staff, guess_dates = c(3, 4))
        wards.staff$adm <- wards.staff$dateofonset_foranalysis - 14
        wards.staff$dis <- wards.staff$dateofcollection - 1
        wards.staff$dateofcollection <- NULL
        wards.staff$dateofonset_foranalysis <- NULL
        names(wards.staff) <- c('id', 'ward', 'adm', 'dis')
        wards.staff$dis <- paste(wards.staff$dis, "12:00:00")
        wards.staff$adm <- paste(wards.staff$adm, "12:00:00")
        
        wards.patients.adm <- wards.patients$adm
        wards.patients.dis <- wards.patients$dis
        
        wards.patients <-
                clean_data(wards.patients, guess_dates = c(3, 4))
        wards.patients$adm <- wards.patients.adm
        wards.patients$dis <- wards.patients.dis
        wards.patients <-
                wards.patients[wards.patients$id %in% unique(mydata$barcode),]
        
        wards <- rbind(wards.staff, wards.patients)
        wards$adm <- as.POSIXct(wards$adm)
        wards$dis <- as.POSIXct(wards$dis)
        
        
        mydata <<- mydata
        return(as.data.frame(wards))
}

log_sum <- function(u, v) {
        return(max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v))))
}

log_sum_vec <- function(w) {
        total <- w[1]
        if (length(w) < 2)
                return(total)
        
        for (i in 2:length(w)) {
                total <- log_sum(total, w[i])
                
        }
        return(total)
}

convolve_log <- function(x, y) {
        n <- length(x)
        m <- length(y)
        
        r <- lapply(1:(n + m - 1), function(k) {
                i <- 1:max(m, n)
                i <- i[((i <= m) & ((k - m + i) <= n)) &
                               ((k - m + i) > 0)]
                log_sum_vec(x[k - m + i] + y[i])
        })
        return(unlist(r))
}
                                  
## add convolutions to data$log_w_dens and data$log_a_dens
## rows = kappa value
## columns = time interval
add_convolutions <- function(data, config) {
        ## COMPUTE CONVOLUTIONS IF NEEDED ##
        if (config$max_kappa > 1) {
                ## first compute convolutions on natural scale
                for (i in 2:config$max_kappa) {
                        if (!is.null(data$log_a_dens)) {
                                data$log_a_dens[[i]] <- log(exp(data$log_a_dens[[i - 1]]) %*%
                                                                    exp(data$log_a_dens[[1]]))
                        }
                        if (!is.null(data$log_w_dens)) {
                                data$log_w_dens <- rbind(
                                        data$log_w_dens,
                                        convolve_log(
                                                data$log_w_dens[i - 1,],
                                                rev(data$log_w_dens)
                                        )[seq_len(ncol(data$log_w_dens))]
                                )
                        }
                }
        }
        if (any(is.infinite(data$log_w_dens)))
                data$log_w_dens[is.na(data$log_w_dens)] <-
                        min(data$log_w_dens[is.finite(data$log_w_dens)])
        
        ## name rows/columns (useful if internal debugging needed)
        if (!is.null(data$log_w_dens)) {
                rownames(data$log_w_dens) <-
                        paste("kappa", seq_len(nrow(data$log_w_dens)),
                              sep = "=")
                colnames(data$log_w_dens) <-
                        seq_len(ncol(data$log_w_dens))
        }
        if (!is.null(data$log_a_dens)) {
                names(data$log_a_dens) <- 1:config$max_kappa
        }
        return(data)
}


## function which returns individual summaries of each MCMC
summary.all.steps <- function(res, burnin) {
        burnt_res <- res[res$step > burnin,]
        lapply(seq_len(length(burnt_res$step)),
               function(i) {
                       out_df <- summary.res(res = burnt_res[i,],
                                             burnin = 0,
                                             support = 0)
                       out_df$step <- burnt_res[i, "step"]
                       return(out_df)
               }) %>% bind_rows()
}

## convert outbreaker res object to wide dataframe which each cases ancestor and associated variables
get_steps <- function(res, data) {
        out_df <- lapply(seq_len(sum(grepl(
                "alpha", names(res)
        ))),
        function(i) {
                out.tib <- tibble(
                        from = res[[glue("alpha_{i}")]],
                        to = i,
                        t_inf = res[[glue("t_inf_{i}")]],
                        t_onw = res[[glue("t_onw_{i}")]],
                        kappa = res[[glue("kappa_{i}")]],
                        ward_from = mapply(
                                get_ward,
                                from,
                                ifelse(is.na(kappa) |
                                               kappa == 1, t_inf, t_onw),
                                MoreArgs = list(data = data)
                        ),
                        ward_to = mapply(get_ward, to, t_inf, MoreArgs = list(data = data)),
                        dist = data$D[cbind(from, to)]
                )
                if (is.na(out.tib$from[1])) {
                        if (typeof(out.tib$ward_from) == "list") {
                                out.tib$ward_from <- NA
                        }
                        if (typeof(out.tib$ward_to) == "list") {
                                out.tib$ward_to <- NA
                        }
                }
                return(out.tib)
        }) %>% bind_rows()
        out_df$step <-
                rep(res$step, sum(grepl("alpha", names(res))))
        return(out_df)
}
