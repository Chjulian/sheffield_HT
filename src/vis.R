d<- d0

#visualization of consensus tree
cons_tree <- vis_epicontacts(d$cons_tree,
                             node_shape="category",
                             node_color ="ward", 
                             shapes=c(inpatient="user", outpatient="user-plus", staff="user-md"),
                             edge_label="SNP",
                             edge_color = "ward.from")
htmlwidgets::saveWidget(cons_tree, "~/Desktop/network_d0.html",   selfcontained = TRUE)

res<- d$res
#visualize trace
plot(res, burn=0)
# 
#ancestries
plot(res, type='alpha', burnin=1000)
# 
# infection dates
plot(res,  type = "t_inf", burnin = 1000) 

##generation between cases
plot(res, type = "kappa", burnin = 1000)
# 
#mutation rate and other params
plot(res, "mu", burn = 1000, type = "density") 
hist(res$eps)
hist(res$tau)
hist(res$pi)



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


# #one Run (quick) summary
contacts.cat <- d$res.sum$tree[!is.na(d$res.sum$tree$from),c("from", "to")]
dictionary.cat <- setNames(mydata$category, mydata$barcode) #dictionary
contacts.cat$from <- as.vector(dictionary.cat[contacts.cat$from])
contacts.cat$to <- as.vector(dictionary.cat[contacts.cat$to])
contacts.cat$from_to <- paste(contacts.cat$from, contacts.cat$to, sep="_to_")
x<-table(contacts.cat$from_to); x <- as.data.frame(x); names(x) <- c("from_to","Freq")


tChains.df <- d$df
mylabel<-"1snps"

p1<-ggplot(tChains.df) + geom_bar(aes(x=size), fill="#D1362F") + theme_minimal() + ggtitle(paste("size_",mylabel, sep=))
p2<-ggplot(tChains.df) + geom_bar(aes(x=n.wards), fill="#D1362F") + theme_minimal() + ggtitle(paste("n.wards_",mylabel, sep=))
p2b<- ggplot(x) + geom_bar(aes(x=from_to, y=Freq), stat='identity', fill="#3182bd") + theme_minimal() +  theme(axis.text.x = element_text(angle = 90)) + ggtitle(paste("size_",mylabel, sep=))
p3<-ggplot(tChains.df) + geom_bar(aes(x=duration), fill="#D1362F") + theme_minimal() + ggtitle(paste("duration_",mylabel, sep=))
p4<-ggplot(tChains.df) + geom_bar(aes(x=lineage), fill="#D1362F") + theme_minimal() +  theme(axis.text.x = element_text(angle = 90)) + ggtitle(paste("lineage_",mylabel, sep=))
p5<-ggplot(tChains.df) + geom_bar(aes(x=n.lineages), fill="#D1362F") + theme_minimal() + ggtitle(paste("n.lineages_",mylabel, sep=))
p6<-ggplot(tChains.df) + geom_bar(aes(x=index.category), fill="#D1362F") + theme_minimal() + ggtitle(paste("index.category_",mylabel, sep=))
p7<-ggplot(tChains.df) + geom_bar(aes(x=index.association), fill="#D1362F") + theme_minimal() +  theme(axis.text.x = element_text(angle = 90)) + ggtitle(paste("index.association_",mylabel, sep=))
p8<-ggplot(tChains.df) + geom_histogram(aes(x=aveSNPpair), fill="#D1362F") + theme_minimal()  + ggtitle(paste("aveSNPpair_",mylabel, sep=))
p9<-ggplot(tChains.df) + geom_bar(aes(x=firstOnset), fill="#D1362F") + theme_minimal() + ggtitle(paste("firstOnset_",mylabel, sep=))

pdf('~/Desktop/0snps.pdf',w=5,h=5,bg='white')
print(p1)
print(p2)
print(p2b)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)
dev.off()


#links per staff
staff <-mydata$barcode[mydata$category=="staff"]
staff.tree <- d$res.sum$tree[d$res.sum$tree$from%in%staff,]
length(unique(staff.tree$from))
sort(table(staff.tree$from), decreasing = T)
summary(as.vector(table(staff.tree$from)))
#wards
sort(table(c(d$cons_tree$contacts$ward.to, d$cons_tree$contacts$ward.from)), decreasing = T)

