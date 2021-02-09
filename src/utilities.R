d<- d_fix0.999
raw <- d$raw
contacts <- d$cons_tree$contacts

#check if two linkeds individual in outbreaker2 share at any time same ward
vector.is.empty <- function(x) return(length(x) ==0 )
raw$wards <- raw$ward <- paste(raw$loc1, raw$loc2,
                               raw$loc3, raw$loc4,
                               raw$loc5, raw$loc6,
                               raw$loc7, raw$loc8)
raw$wards <- trimws(raw$wards)
contacts$overlap <- NA
for(i in 1:nrow(contacts)){
        x<- contacts[i,'from']
        y<- contacts[i,'to']
        x.w <- raw$wards[raw$barcode==x]
        y.w <- raw$wards[raw$barcode==y]
        i.w <- intersect(x.w,y.w)
        if(!vector.is.empty(i.w)) contacts$overlap[i]<-i.w
        rm(x,y,x.w,y.w,i.w)
}
contacts$overlap <- stringr::str_replace_na(contacts$overlap)
contacts$overlap<-stringr::str_replace_all(contacts$overlap, "multiple", "NA")
contacts$overlap<-stringr::str_replace_all(contacts$overlap, "nncl", "NA")
sum(contacts$overlap!="NA")
