if(!require("limma")){install.packages("limma")}
if(!require("Biobase")){install.packages("Biobase")}
if(!require("ggplot2")){install.packages("ggplot2")}
if(!require("pamr")){install.packages("pamr")}
if(!require("gplots")){install.packages("gplots")}
if(!require("fpc")){install.packages("fpc")}
if(!require("cluster")){install.packages("cluster")}

#nascondo i warning, presenti a causa degli NA nel dataset di partenza
options(warn = -1)

dir.create("./immagini", showWarnings = FALSE)
dir.create("./immagini/plot base", showWarnings = FALSE)
dir.create("./immagini/biomarcatori/", showWarnings = FALSE)
dir.create("./immagini/clustering/", showWarnings = FALSE)


dati = read.table("./dati.csv", header = T, sep =",")
concentrazioni = read.table("./concentrazioni.csv", header = T, sep =",", check.names = F, dec=",")
colnames(concentrazioni)[1] = "sample.id"

#escludo dal dati quelli col campo dati.non.pronti =1
dati=dati[dati$dati.non.pronti==0,]
dati$dati.non.pronti = NULL

#Manipolazione dataset
dati$sample.id = levels(dati$sample.id)[dati$sample.id] #from factor to string
concentrazioni$sample.id = levels(concentrazioni$sample.id)[concentrazioni$sample.id] #from factor to string
dati$dob = as.Date(dati$dob, format="%d/%m/%Y")
dati$css = as.factor(dati$css)
dati$qit = as.numeric(as.character(dati$qit))
rownames(dati) = as.vector(dati[,1])
rownames(concentrazioni) = as.vector(concentrazioni[,1])


cat("\n Summary di dati: \n")
print(summary(dati[1:6]))
print(summary(dati[7:12]))
print(summary(dati[13:15]))

cat("\n Summary di concentrazioni: \n")
print(summary(concentrazioni[1:6]))
print(summary(concentrazioni[7:12]))
print(summary(concentrazioni[13:18]))
print(summary(concentrazioni[19:24]))
print(summary(concentrazioni[25:30]))
print(summary(concentrazioni[31:36]))
print(summary(concentrazioni[37:42]))
print(summary(concentrazioni[42:45]))


#creazione dataset "completo" (per analisi e clustering)
db = merge(dati[,c(1,3,4)], concentrazioni)


#----------------------------ANALISI DEL DATASET
dati$eta2 = cut(dati$eta, breaks = c(0,5,10,17), labels = c("1-5 anni","6-10 anni", "11-17 anni"))
dati_patologico = dati[dati$tipo=="patologico",]


ggplot(data = dati, aes(tipo))+geom_bar() + theme(axis.title.y = element_text(size = rel(1.8), angle = 90),axis.text=element_text(size=20)) + labs(x="")
ggsave(filename="./immagini/plot base/barplot tipo.jpg", width=8.5, dpi=300)

ggplot(data = dati, aes(dob, fill = tipo, group = tipo))+geom_histogram(position = "stack")+scale_x_date(date_breaks = "1 year", date_labels = "%Y")+theme(axis.text=element_text(size=8),axis.title.y = element_text(size = rel(1.8), angle = 90),legend.text=element_text(size=15), legend.position="bottom")+labs(x="",fill="")
ggsave(filename="./immagini/plot base/Istogramma date tipo.jpg", width=8.5, dpi=300)

ggplot(data = dati, aes(dob, fill = factor(eta2), group = eta2))+geom_histogram(position = "stack")+scale_x_date(date_breaks = "1 year", date_labels = "%Y")+theme(axis.text=element_text(size=8),axis.title.y = element_text(size = rel(1.8), angle = 90),legend.text=element_text(size=15), legend.position="bottom")+labs(x="",fill="")
ggsave(filename="./immagini/plot base/Istogramma date eta.jpg", width=10,height = 6, dpi=300)

ggplot(data = dati, aes(eta2, fill = tipo, group = tipo))+geom_bar(position = "dodge")+ theme(axis.title.y = element_text(size = rel(1.8), angle = 90),axis.text=element_text(size=20),legend.text=element_text(size=15), legend.position="bottom")+labs(x="", fill="")
ggsave(filename="./immagini/plot base/Barplot eta per tipo.jpg", width=8.5, dpi=300)
  
ggplot(data = dati_patologico, aes(eta2 ,fill = factor(css), group = css))+geom_bar(position = "stack")+ theme(axis.title.y = element_text(size = rel(1.8), angle = 90),axis.text=element_text(size=20),legend.text=element_text(size=15), legend.position="bottom")+labs(x="", fill="CSS")
ggsave(filename="./immagini/plot base/Barplot patologici eta per css.jpg", width=8.5, dpi=300)

dati$eta2 = NULL

#----------------------------RICERCA BIOMARCATORI CON LIMMA

#trasposta di concentrazioni (necessaria per il match della matrice di design)
concentrazioni = t(concentrazioni[,2:ncol(concentrazioni)])

#seleziono sample comuni in concentrazioni con dati
new_concentrazioni= concentrazioni[,which(colnames(concentrazioni) %in% rownames(dati))]

#Matrice di design
design = model.matrix(~0+tipo,data=dati)
colnames(design)= c("sano", "patologico")

#I dati non provengono da NGS quindi non bisogna normalizzarli
#normalized.concentrazioni = voom(new_concentrazioni, design)

initial.model = lmFit(new_concentrazioni, design)
limma.model = eBayes(initial.model)

# campiono i dati per le heatmap
sani = dati$sample.id[dati$tipo=="controllo sano"]
patologici = dati$sample.id[dati$tipo=="patologico"]
campioni = as.vector(na.omit(c(sani,patologici)))
colori = character(length (campioni))
names (colori) = campioni
colori[patologici] = "red"
colori[sani] = "blue"


#ANALISI BASE

result1 = toptable( limma.model, number = nrow(new_concentrazioni), adjust.method = "BH",coef ="patologico", p.value=0.01, lfc=1)

#heatmap completa
geni_diff=rownames(result1)
espressioni_diff= new_concentrazioni[geni_diff,]
jpeg("./immagini/biomarcatori/heatmapLimma_base.jpg", width = 800,  quality=100)
heatmap.2(espressioni_diff, col=redgreen(100), scale="row",key= T, symkey=F, density.info="none", trace="none", ColSideColors = colori, cexRow = 0.8)
dev.off()

#heatmap sani patologici
espressioni_diff=apply(espressioni_diff[,2:3],1:2, as.numeric)
colnames(espressioni_diff) = c("sano", "patologico")
png("./immagini/biomarcatori/heatmapLimma_base.png")
heatmap.2(espressioni_diff, col=redgreen(100), scale="col",key= T, symkey=F, density.info="none", trace="none", srtCol=45, cexCol = 1.2,cexRow = 0.8)
dev.off()


#ANALISI AVANZATA
contrasts = makeContrasts(patologico-sano, levels = design)
contrasts.model =  eBayes(contrasts.fit(limma.model,contrasts))
result = toptable( contrasts.model, number = nrow(new_concentrazioni), adjust.method = "BH", p.value=0.1, lfc=1.5)

#heatmap completa
geni_diff=rownames(result)
espressioni_diff= new_concentrazioni[geni_diff,]
jpeg("./immagini/biomarcatori/heatmapLimmacomplete_adv.png", width = 800,  quality=100)
heatmap.2(espressioni_diff, col=redgreen(100), scale="row",key= T, symkey=F, density.info="none", trace="none", ColSideColors = colori)
dev.off()


#heatmap sani patologici
espressioni_diff=apply(espressioni_diff[,2:3],1:2, as.numeric)
colnames(espressioni_diff) = c("sano", "patologico")
png("./immagini/biomarcatori/heatmapLimma.png")
heatmap.2(espressioni_diff, col=redgreen(100), scale="col",key= T, symkey=F, density.info="none", trace="none", srtCol=45, cexCol = 1.2,cexRow = 2)
dev.off()



#----------------------------RICERCA BIOMARCATORI CON PAMR
all_geni=rownames(new_concentrazioni)
names(all_geni) = rownames(new_concentrazioni)

exp.data= list(
  "x"= new_concentrazioni,
  "y" = dati$tipo,
  "genenames"= all_geni,
  "geneid" = all_geni
)


pamr.trained = pamr.train(exp.data, geni_diff)
pamr.cv= pamr.cv(pamr.trained,exp.data)
png("./immagini/biomarcatori/pamrCrossValidationTest.png")
pamr.plotcv(pamr.cv)
dev.off()

pamr.threshold= pamr.cv$threshold[which.min(pamr.cv$error)]
cat("\n Pamr soglia scelta:",pamr.threshold, "\n")
pamr.biomarkers = pamr.listgenes(pamr.trained, exp.data, pamr.threshold)

png("./immagini/biomarcatori/pamrPlotCen.png")
pamr.plotcen(pamr.trained, exp.data, pamr.threshold)
dev.off()

png("./immagini/biomarcatori/pamrGenePlot.png", height = 600)
pamr.geneplot(pamr.trained, exp.data, pamr.threshold)
dev.off()

#estrarre nome dei geni differenzialmente espressi
geni_diff_pamr= as.vector(pamr.biomarkers[,"id"])

#heatmap sani patologici
rownames(pamr.biomarkers)=as.vector(pamr.biomarkers[,1])
pamr.biomarkers=apply(pamr.biomarkers[,2:3],1:2, as.numeric)
colnames(pamr.biomarkers) = c("sano", "patologico")
png("./immagini/biomarcatori/heatmapPamr.png")
heatmap.2(pamr.biomarkers, col=redgreen(100), scale="col",key= T, symkey=F, density.info="none", trace="none", srtCol=45, cexCol = 1.2,cexRow = 2)
dev.off()


#----------------------------BOXPLOT BIOMARCATORI

ggplot(data = db, aes(tipo,C2))+geom_boxplot(aes(fill = tipo))+ theme(legend.position="bottom",axis.text=element_text(size=20))+ labs(x="",fill="")
ggsave(filename="./immagini/biomarcatori/C2.jpg", width=8.5, dpi=300)

ggplot(data = db, aes(tipo,ARG))+geom_boxplot(aes(fill = tipo))+ theme(axis.text=element_text(size=20),legend.position="bottom")+ labs(x="",fill="")
ggsave(filename="./immagini/biomarcatori/ARG.jpg", width=8.5, dpi=300)

ggplot(data = db, aes(tipo,CIT))+geom_boxplot(aes(fill = tipo))+ theme(axis.text=element_text(size=20),legend.position="bottom")+ labs(x="",fill="")
ggsave(filename="./immagini/biomarcatori/CIT.jpg", width=8.5, dpi=300)



#----------------------------CLUSTERING
dbclear = db
dbclear$SA = NULL #anche se setto l'unico valore NA a 0 non appera tra i più correlati (anzi tra i meno)
dbclear$`C5:1`= NULL
dbclear$`C10:2`= NULL
dbclear$`C18-OH`= NULL
dbclear =subset(dbclear, is.na(dbclear$eta) != T)
db$tipo = as.numeric(db$tipo) #!!!!!!!!!!!!!!!  
dbclear$tipo = as.numeric(dbclear$tipo) #!!!!!!!!!!!!!!!
rownames(db) = db$sample.id

kmeans = kmeans(dbclear[,c(-1)], 2, iter.max=500, trace=2)
cat ("\nClustering con kmeans:\n")
print(table(dbclear$tipo, kmeans$cluster))
png("./immagini/clustering/cluster con kmeans.png")
plotcluster(dbclear[,c(-1,-3)], kmeans$cluster, main="Clustering con kmeans",pch = 19)
dev.off()


cat ("\n\nClustering gerarchico distanza di manhattan metodo single\n")
d <- dist(db[, c(-1,-2)], method = "manhattan") # distance matrix
fit <- hclust(d, method="single") 
groups <- cutree(fit, k=2) 
png("./immagini/clustering/Distanza manhattan metodo single.png")
plot(as.dendrogram(fit), main="Distanza manhattan metodo single", col=groups, leaflab = "none")
rect.hclust(fit, k=2, border="red")
dev.off()
print(table(db$tipo,groups))
cat("Media diametri cluster: ", mean(cluster.stats(d, groups)$diameter))
cat("\nSeparazione: ", cluster.stats(d, groups)$separation[1])

cat ("\n\nClustering gerarchico distanza minkowski metodo complete:\n")
d <- dist(db[, c(-1,-2)], method = "minkowski")
fit <- hclust(d, method="complete") 
groups <- cutree(fit, k=2)
png("./immagini/clustering/Distanza minkowski metodo complete.png")
plot(as.dendrogram(fit), main="Distanza euclidea metodo complete", col=groups, leaflab = "none")
rect.hclust(fit, k=2, border="red")
dev.off()
print(table(db$tipo,groups))
cat("Media diametri cluster: ", mean(cluster.stats(d, groups)$diameter))
cat("\nSeparazione: ", cluster.stats(d, groups)$separation[1])


cat ("\n\nClustering gerarchico distanza euclidea metodo ward.D:\n")
d <- dist(db[, c(-1)], method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
groups <- cutree(fit, k=2)
png("./immagini/clustering/Distanza euclidea metodo wardD.png")
plot(as.dendrogram(fit), main="Distanza euclidea metodo ward.D",  col=groups, leaflab = "none")
rect.hclust(fit, k=2, border="red")
dev.off()
print(table(db$tipo,groups))
cat("Media diametri cluster: ", mean(cluster.stats(d, groups)$diameter))
cat("\nSeparazione: ", cluster.stats(d, groups)$separation[1])

db$tipo = factor(db$tipo, labels = c("controllo sano", "patologico"))
