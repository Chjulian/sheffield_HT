#visualization of consensus tree
vis_epicontacts(cons_tree,
                node_shape="category",
                node_color ="ward", 
                shapes=c(inpatient="user", outpatient="user-plus", staff="user-md"),
                edge_label="SNP",
                edge_color = "ward.from")

#visualize trace
plot(res, burn=1000)
# 
#ancestries
plot(res, type='alpha', burnin=1000)
# 
# infection dates
plot(res,  type = "t_inf", burnin = 1000) 

##generation between cases
plot(res, type = "kappa", burnin = 1000)
# 
#mutation rate
plot(res, "mu", burn = 1000, type = "density") 

# transmission trees
mynetwork<-plot(res,  type = "network", burnin = 1000, 
                min_support = .05, label=meta$id)
mynetwork
# fancier transmission trees
dictionary.icons <-setNames(c('f007','f234','f0f0'),c('inpatient','outpatient','staff')) 
mynetwork$x$nodes$shape <-  rep("icon",length(mynetwork$x$nodes$shape))
mynetwork$x$nodes$icon.code <-  as.vector(dictionary.icons[mydata$category])
mynetwork$x$nodes$icon.color <-  mynetwork$x$nodes$color
mynetwork$x$nodes$icon.size <-  (mynetwork$x$nodes$value+1)*25
mynetwork%>%
        addFontAwesome() #https://fontawesome.com/cheatsheet    


ggplot(tChains.df) + geom_histogram(aes(x=size), fill="#3182bd") + theme_minimal()
ggplot(tChains.df) + geom_bar(aes(x=n.wards), fill="#3182bd") + theme_minimal()
ggplot(tChains.df) + geom_bar(aes(x=duration), fill="#3182bd") + theme_minimal()
ggplot(tChains.df) + geom_bar(aes(x=lineage), fill="#3182bd") + theme_minimal() +  theme(axis.text.x = element_text(angle = 90))
ggplot(tChains.df) + geom_bar(aes(x=n.lineages), fill="#3182bd") + theme_minimal() 
ggplot(tChains.df) + geom_bar(aes(x=index.category), fill="#3182bd") + theme_minimal()
ggplot(tChains.df) + geom_bar(aes(x=index.association), fill="#3182bd") + theme_minimal() +  theme(axis.text.x = element_text(angle = 90))
ggplot(tChains.df) + geom_histogram(aes(x=aveSNPpair), fill="#3182bd") + theme_minimal() 
ggplot(tChains.df) + geom_bar(aes(x=firstOnset), fill="#3182bd") + theme_minimal() 
ggplot(x) + geom_bar(aes(x=from_to, y=Freq), stat='identity', fill="#3182bd") + theme_minimal() +  theme(axis.text.x = element_text(angle = 90)) 

# #one Run (quick) summary
# contacts.cat <- mydf$tree[!is.na(mydf$tree$from),c("from", "to")]
# dictionary.cat <- setNames(mydata$category, mydata$barcode) #dictionary
# contacts.cat$from <- as.vector(dictionary.cat[contacts.cat$from])
# contacts.cat$to <- as.vector(dictionary.cat[contacts.cat$to])
# contacts.cat$from_to <- paste(contacts.cat$from, contacts.cat$to, sep="_to_")
# x<-table(contacts.cat$from_to); x <- as.data.frame(x); names(x) <- c("from_to","Freq")