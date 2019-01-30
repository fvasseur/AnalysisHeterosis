#===========================================
#===========================================
#
# 1. PREPARE DATAFILEs
#
#===========================================
#===========================================



#===========================================
# Measure pairwise PHENOTYPIC distances
#===========================================

# Euclidian phenotypic distance
# --------------------------------------------------
tabtot <- read.csv("S1_Data.csv", header=T, sep=",")
for(i in (c(8:11)))
{
  tabtot[,i] <- as.numeric(as.character(tabtot[,i])) 
}
for(i in (c(1:7)))
{
  tabtot[,i] <- as.factor(as.character(tabtot[,i])) 
}

tab <- aggregate(tabtot[,c(8:11)], by=list(idExp=tabtot$idExp,
                                           idGenotype=tabtot$idGenotype,
                                           idMother=tabtot$idMother,
                                           idFather=tabtot$idFather,
                                           Type=tabtot$Type), FUN=mean, na.rm=T)


acc <- tab[tab$idExp=="GT01" & tab$Type=="Inbred",]
acc <- droplevels(acc)
geno <- read.csv("1001genomes-accessions.csv", header=T, sep=",")
names(geno)[1] <- "idMother" 
geno <- geno[,c(1,6,7,8,9)]
acc <- merge(geno, acc, all.y=T)
acc <- acc[,c(7,6,1,8,9,4,3,5,2,10:13)]
acc$idMother <- as.factor(as.character(acc$idMother))
acc <- droplevels(acc)

acctemp <- acc
acctemp[,c(10:13)] <- log10(acctemp[,c(10:13)])
matrix_EucDistPhenoTot <- dist(acctemp[,c(10:13)])
matrix_EucDistPhenoTot <- as.matrix(matrix_EucDistPhenoTot)
rownames(matrix_EucDistPhenoTot) <- acc$idMother
colnames(matrix_EucDistPhenoTot) <- acc$idMother





#===========================================
# Measure pairwise GEOGRAPHIC distances
#===========================================

# Import GEOGRAPHIC data
#----------------------------------------
GeoAcc <- acc[,c(1,6,7)]
matgeo <- as.matrix(GeoAcc[,c(2,3)])
matrix_GeoDist <- distm(matgeo, matgeo, fun=distGeo)
rownames(matrix_GeoDist) <- acc$idMother
colnames(matrix_GeoDist) <- acc$idMother





#===========================================
# Measure pairwise GENETIC distances
#===========================================
# use vcftools from terminal to filter SNPs of the whole 1001 genomes data
# vcftools --gzvcf PATH/1001genomes_snp-short-indel_only_ACGTN.vcf.gz --minDP 10 --keep PATH/list_acc.txt --max-missing 0.85 --recode
#
# Then PLINK to calculate ibs genetic distances
# ./plink --vcf PATH/407acc --distance ibs allele-ct square
# (unit = allele count)

# Import Genetic data
#----------------------------------------
matrix_GeneticDist <- read.table("plink.dist", header=F, sep="\t")
names_GeneticDist <- read.table("plink.dist.id", header=F, sep="\t")
names(matrix_GeneticDist) <- names_GeneticDist$V1
rownames(matrix_GeneticDist) <- names_GeneticDist$V1 
matrix_GeneticDist <- as.matrix(matrix_GeneticDist)
str(matrix_GeneticDist)





# Combine distances
#------------------------------------------------
fc_transform_DistMatrix <- function(matrix, name, disp.all=T){
  
  matemp <- as.matrix(matrix)
  tabtemp <- melt(matemp)
  names(tabtemp)[1] <- "Geno1"
  names(tabtemp)[2] <- "Geno2"
  names(tabtemp)[3] <- "Dist"
  tabtemp$Geno1 <- as.factor(tabtemp$Geno1)
  tabtemp$Geno2 <- as.factor(tabtemp$Geno2)
  tabtemp$comb <- tabtemp$Geno1:tabtemp$Geno2
  tabtemp <- tabtemp[,c(4,1:3)]
  names(tabtemp)[4] <- name
  
  if(disp.all==F){
    fac <- expand.grid(rownames(matemp),colnames(matemp))
    fac$key <- apply(fac, 1, function(x)paste(sort(x), collapse=':'))
    str(fac) 
    fac <- subset(fac, !duplicated(fac$key)) 
    names(fac) <- c("Geno1", "Geno2", "comb")
    fac <- fac[,c(3,1,2)]
    tabtemp <- merge(fac, tabtemp[,c(1,4)], by="comb", all.x=T)
  }
  
  tabDist1 <- tabtemp
  return(tabDist1)
}


table_GeneticDist <- fc_transform_DistMatrix(matrix=matrix_GeneticDist, name="GeneticDist", disp.all=F)
table_DistPhenoTot <- fc_transform_DistMatrix(matrix=matrix_EucDistPhenoTot, name="DistPhenoTot", disp.all=F)
table_GeoDist <- fc_transform_DistMatrix(matrix=matrix_GeoDist, name="GeoDist", disp.all=F)

table_GeoDist <- table_GeoDist[,c(1,4)]
tabDist <- merge(table_GeneticDist, table_GeoDist, by="comb", all.y=T, all.x=T)
table_DistPhenoTot <- table_DistPhenoTot[,c(1,4)]
tabDist <- merge(tabDist, table_DistPhenoTot, by="comb", all.x=T, all.y=T)

for(i in c(4:6))
{
  tabDist[,i] <- as.numeric(as.character(tabDist[,i]))
}
tabDist$comb <- as.factor(as.character(tabDist$comb))
tabDist$Geno1 <- as.factor(as.character(tabDist$Geno1))
tabDist$Geno2 <- as.factor(as.character(tabDist$Geno2))
tabDist <- droplevels(tabDist)



# Merge with accessions info about genetic AVegetativeDryMassIXTURE group 
#-------------------------------------------------------------
acc1 <- acc[,c(3,8,7)]
names(acc1)[1] <- "Geno1"
names(acc1)[2] <- "GroupGeno1"
names(acc1)[3] <- "RelictGeno1"
acc1 <- merge(acc1, tabDist, by="Geno1", all.y=T)

acc2 <- acc[,c(3,8,7)]
names(acc2)[1] <- "Geno2"
names(acc2)[2] <- "GroupGeno2"
names(acc2)[3] <- "RelictGeno2"
acc2 <- merge(acc1, acc2, by="Geno2", all.y=T, all.x=T)

acc2 <- acc2[,c(5,2,1,3,4,9,10,6:8)]
acc2$GroupGeno1 <- as.factor(as.character(acc2$GroupGeno1))
acc2$GroupGeno2 <- as.factor(as.character(acc2$GroupGeno2))

acc2$Group <- NA
for(i in 1:dim(acc2)[1])
{
  if(is.na(as.character(acc2[i, "GroupGeno1"]))){acc2[i, "Group"] <- "unknown"} else{
    if(is.na(as.character(acc2[i, "GroupGeno2"]))){acc2[i, "Group"] <- "unknown"} else{
      if(as.character(acc2[i, "GroupGeno1"])!="relict") {if(as.character(acc2[i, "GroupGeno2"])!="relict") {acc2[i, "Group"] <- "nonrelict_vs_nonrelict"}}
      if(as.character(acc2[i, "GroupGeno1"])=="relict") {if(as.character(acc2[i, "GroupGeno2"])=="relict") {acc2[i, "Group"] <- "relict_vs_relict"}}
      if(as.character(acc2[i, "GroupGeno1"])=="relict") {if(as.character(acc2[i, "GroupGeno2"])!="relict") {acc2[i, "Group"] <- "relict_vs_nonrelict"}}
      if(as.character(acc2[i, "GroupGeno1"])!="relict") {if(as.character(acc2[i, "GroupGeno2"])=="relict") {acc2[i, "Group"] <- "relict_vs_nonrelict"}}
    }
  }
  print(i)
} 

acc2$Group <- as.factor(acc2$Group)
tabDist <- acc2

#write.table(tabDist, "Distance_accessions_UniquePairs.csv", sep=",", dec=".", col.names = T, row.names = F)









#===========================================
# Categorization of heterosis effect
#===========================================

tabDist2 <- tabDist[,c(1,4:11)]
names(tabDist2)[1:5] <- c("idGenotype","GroupMother", "RelictMother","GroupFather","RelictFather")

tabacc <- tab[tab$Type=="Inbred",]
tabacc <- aggregate(tabacc[,c(6:10,12)], by=list(idGenotype=tabacc$idGenotype,
                                                 idMother=tabacc$idMother,
                                                 idFather=tabacc$idFather), FUN=mean, na.rm=T)
tabacc <- droplevels(tabacc)

hyb <- tabtot[tabtot$Type=="Hybrid",]
hyb <- droplevels(hyb)

tabhyb <- tab[tab$Type=="Hybrid",]
tabhyb <- droplevels(tabhyb)
tabhyb <- merge(tabhyb, tabDist2, by="idGenotype", all.x=T)
tabhyb <- tabhyb[,c(1,3,4,6,7,10,9,8,12:19)]
names(tabhyb)[4:9] <- c("VegetativeDryMass","Lifespan","VegetativeDryMass","GR","Fruit","Age50")



# Parental trait values
#----------------------------------------------
for(i in levels(tabhyb$idGenotype))
{
  moth <- as.character(tabhyb[tabhyb$idGenotype==i, "idMother"])
  genomoth <- paste(moth, ":", moth, sep="")
  fath <- as.character(tabhyb[tabhyb$idGenotype==i, "idFather"])
  genofath <- paste(fath, ":", fath, sep="")
  
  tabhyb[tabhyb$idGenotype==i, "MotherVegetativeDryMass"] <- tabacc[tabacc$idGenotype==genomoth,"VegetativeDryMass"]
  tabhyb[tabhyb$idGenotype==i, "FatherVegetativeDryMass"] <- tabacc[tabacc$idGenotype==genofath,"VegetativeDryMass"]
  tabhyb[tabhyb$idGenotype==i, "MotherGR"] <- tabacc[tabacc$idGenotype==genomoth,"GrowthRate"]
  tabhyb[tabhyb$idGenotype==i, "FatherGR"] <- tabacc[tabacc$idGenotype==genofath,"GrowthRate"]
  tabhyb[tabhyb$idGenotype==i, "MotherLifespan"] <- tabacc[tabacc$idGenotype==genomoth,"Lifespan"]
  tabhyb[tabhyb$idGenotype==i, "FatherLifespan"] <- tabacc[tabacc$idGenotype==genofath,"Lifespan"]
  tabhyb[tabhyb$idGenotype==i, "MotherFruit"] <- tabacc[tabacc$idGenotype==genomoth,"FruitNumber"]
  tabhyb[tabhyb$idGenotype==i, "FatherFruit"] <- tabacc[tabacc$idGenotype==genofath,"FruitNumber"]
}

tabhyb$distLifespan <- abs(tabhyb$MotherLifespan - tabhyb$FatherLifespan)
tabhyb$distVegetativeDryMass <- abs(tabhyb$MotherVegetativeDryMass - tabhyb$FatherVegetativeDryMass)
tabhyb$distGR <- abs(tabhyb$MotherGR - tabhyb$FatherGR)
tabhyb$distFruit <- abs(tabhyb$MotherFruit - tabhyb$FatherFruit)




# Categorization of hybrids  
#-----------------------------------------------------
tabhyb$heterosis_type_Lifespan <- "null"
tabhyb$heterosis_type_VegetativeDryMass <- "null"
tabhyb$heterosis_type_GR <- "null"
tabhyb$heterosis_type_Fruit <- "null"
threshold <- 0.05/(450*3*4)

for(i in levels(tabhyb$idGenotype))
{
  
  # Age at reproduction (= lifespan)
  maxpar <- max(tabhyb[tabhyb$idGenotype==i, c("MotherLifespan","FatherLifespan")], na.rm=T)
  minpar <- min(tabhyb[tabhyb$idGenotype==i, c("MotherLifespan","FatherLifespan")], na.rm=T)
  meanpar <- mean(as.numeric(tabhyb[tabhyb$idGenotype==i, c("MotherLifespan","FatherLifespan")]), na.rm=T)
  
  if(length(na.omit(hyb[hyb$idGenotype==i,"Lifespan"]))>1){
    boot_hybLifespan <- one.boot(na.omit(hyb[hyb$idGenotype==i,"Lifespan"]), mean, R=1000)
    mod1 <- t.test(boot_hybLifespan$t, mu=minpar, alternative = "less")
    mod2 <- t.test(boot_hybLifespan$t, mu=meanpar, alternative = "two.sided")
    mod3 <- t.test(boot_hybLifespan$t, mu=maxpar, alternative = "greater")
    est <- as.numeric(mod2$estimate)
    tabhyb[tabhyb$idGenotype==i, "CIhybLifespaninf"] <- mod1$conf.int[1]
    tabhyb[tabhyb$idGenotype==i, "CIhybLifespansup"] <- mod1$conf.int[2] 
    
    if(est>meanpar & mod2$p.value < threshold) {if(est>maxpar & mod3$p.value < threshold) {
      tabhyb[tabhyb$idGenotype==i, "heterosis_type_Lifespan"] <- "AboveBestPar"} else{
        tabhyb[tabhyb$idGenotype==i, "heterosis_type_Lifespan"] <- "AboveMeanPar"}} else{
          if(est<minpar & mod1$p.value < threshold) {
            tabhyb[tabhyb$idGenotype==i, "heterosis_type_Lifespan"] <- "BelowWorstPar"} else{ 
              if(est<meanpar & mod2$p.value < threshold) {tabhyb[tabhyb$idGenotype==i, "heterosis_type_Lifespan"] <- "BelowMeanPar"} else{}}}
  } else{
    CIhybLifespaninf <- NA 
    CIhybLifespansup <- NA
    tabhyb[tabhyb$idGenotype==i, "CIhybLifespaninf"] <- NA
    tabhyb[tabhyb$idGenotype==i, "CIhybLifespansup"] <- NA
    tabhyb[tabhyb$idGenotype==i, "heterosis_type_Lifespan"] <-NA
  }
  
  
  # VegetativeDryMass
  maxpar <- max(tabhyb[tabhyb$idGenotype==i, c("MotherVegetativeDryMass","FatherVegetativeDryMass")], na.rm=T)
  minpar <- min(tabhyb[tabhyb$idGenotype==i, c("MotherVegetativeDryMass","FatherVegetativeDryMass")], na.rm=T)
  meanpar <- mean(as.numeric(tabhyb[tabhyb$idGenotype==i, c("MotherVegetativeDryMass","FatherVegetativeDryMass")]), na.rm=T)
  
  if(length(na.omit(hyb[hyb$idGenotype==i,"VegetativeDryMass"]))>1){
    boot_hybLifespan <- one.boot(na.omit(hyb[hyb$idGenotype==i,"VegetativeDryMass"]), mean, R=1000)
    mod1 <- t.test(boot_hybLifespan$t, mu=minpar, alternative = "less")
    mod2 <- t.test(boot_hybLifespan$t, mu=meanpar, alternative = "two.sided")
    mod3 <- t.test(boot_hybLifespan$t, mu=maxpar, alternative = "greater")
    est <- as.numeric(mod2$estimate)
    tabhyb[tabhyb$idGenotype==i, "CIhybVegetativeDryMassinf"] <- mod1$conf.int[1]
    tabhyb[tabhyb$idGenotype==i, "CIhybVegetativeDryMasssup"] <- mod1$conf.int[2] 
    
    if(est>meanpar & mod2$p.value < threshold) {if(est>maxpar & mod3$p.value < threshold) {
      tabhyb[tabhyb$idGenotype==i, "heterosis_type_VegetativeDryMass"] <- "AboveBestPar"} else{
        tabhyb[tabhyb$idGenotype==i, "heterosis_type_VegetativeDryMass"] <- "AboveMeanPar"}} else{
          if(est<minpar & mod1$p.value < threshold) {
            tabhyb[tabhyb$idGenotype==i, "heterosis_type_VegetativeDryMass"] <- "BelowWorstPar"} else{ 
              if(est<meanpar & mod2$p.value < threshold) {tabhyb[tabhyb$idGenotype==i, "heterosis_type_VegetativeDryMass"] <- "BelowMeanPar"} else{}}}
  } else{
    CIhybVegetativeDryMassinf <- NA 
    CIhybVegetativeDryMasssup <- NA
    tabhyb[tabhyb$idGenotype==i, "CIhybVegetativeDryMassinf"] <- NA
    tabhyb[tabhyb$idGenotype==i, "CIhybVegetativeDryMasssup"] <- NA
    tabhyb[tabhyb$idGenotype==i, "heterosis_type_VegetativeDryMass"] <-NA
  }
  
  # GrowthRate
  maxpar <- max(tabhyb[tabhyb$idGenotype==i, c("MotherGR","FatherGR")], na.rm=T)
  minpar <- min(tabhyb[tabhyb$idGenotype==i, c("MotherGR","FatherGR")], na.rm=T)
  meanpar <- mean(as.numeric(tabhyb[tabhyb$idGenotype==i, c("MotherGR","FatherGR")]), na.rm=T)
  
  if(length(na.omit(hyb[hyb$idGenotype==i,"GrowthRate"]))>1){
    boot_hybLifespan <- one.boot(na.omit(hyb[hyb$idGenotype==i,"GrowthRate"]), mean, R=1000)
    mod1 <- t.test(boot_hybLifespan$t, mu=minpar, alternative = "less")
    mod2 <- t.test(boot_hybLifespan$t, mu=meanpar, alternative = "two.sided")
    mod3 <- t.test(boot_hybLifespan$t, mu=maxpar, alternative = "greater")
    est <- as.numeric(mod2$estimate)
    tabhyb[tabhyb$idGenotype==i, "CIhybGRinf"] <- mod1$conf.int[1]
    tabhyb[tabhyb$idGenotype==i, "CIhybGRsup"] <- mod1$conf.int[2] 
    
    if(est>meanpar & mod2$p.value < threshold) {if(est>maxpar & mod3$p.value < threshold) {
      tabhyb[tabhyb$idGenotype==i, "heterosis_type_GR"] <- "AboveBestPar"} else{
        tabhyb[tabhyb$idGenotype==i, "heterosis_type_GR"] <- "AboveMeanPar"}} else{
          if(est<minpar & mod1$p.value < threshold) {
            tabhyb[tabhyb$idGenotype==i, "heterosis_type_GR"] <- "BelowWorstPar"} else{ 
              if(est<meanpar & mod2$p.value < threshold) {tabhyb[tabhyb$idGenotype==i, "heterosis_type_GR"] <- "BelowMeanPar"} else{}}}
  } else{
    CIhybGRinf <- NA 
    CIhybGRsup <- NA
    tabhyb[tabhyb$idGenotype==i, "CIhybGRinf"] <- NA
    tabhyb[tabhyb$idGenotype==i, "CIhybGRsup"] <- NA
    tabhyb[tabhyb$idGenotype==i, "heterosis_type_GR"] <-NA
  }
  
  
  # FruitNumber
  maxpar <- max(tabhyb[tabhyb$idGenotype==i, c("MotherFruit","FatherFruit")], na.rm=T)
  minpar <- min(tabhyb[tabhyb$idGenotype==i, c("MotherFruit","FatherFruit")], na.rm=T)
  meanpar <- mean(as.numeric(tabhyb[tabhyb$idGenotype==i, c("MotherFruit","FatherFruit")]), na.rm=T)
  
  if(length(na.omit(hyb[hyb$idGenotype==i,"FruitNumber"]))>1){
    boot_hybLifespan <- one.boot(na.omit(hyb[hyb$idGenotype==i,"FruitNumber"]), mean, R=1000)
    mod1 <- t.test(boot_hybLifespan$t, mu=minpar, alternative = "less")
    mod2 <- t.test(boot_hybLifespan$t, mu=meanpar, alternative = "two.sided")
    mod3 <- t.test(boot_hybLifespan$t, mu=maxpar, alternative = "greater")
    est <- as.numeric(mod2$estimate)
    tabhyb[tabhyb$idGenotype==i, "CIhybFruitinf"] <- mod1$conf.int[1]
    tabhyb[tabhyb$idGenotype==i, "CIhybFruitsup"] <- mod1$conf.int[2] 
    
    if(est>meanpar & mod2$p.value < threshold) {if(est>maxpar & mod3$p.value < threshold) {
      tabhyb[tabhyb$idGenotype==i, "heterosis_type_Fruit"] <- "AboveBestPar"} else{
        tabhyb[tabhyb$idGenotype==i, "heterosis_type_Fruit"] <- "AboveMeanPar"}} else{
          if(est<minpar & mod1$p.value < threshold) {
            tabhyb[tabhyb$idGenotype==i, "heterosis_type_Fruit"] <- "BelowWorstPar"} else{ 
              if(est<meanpar & mod2$p.value < threshold) {tabhyb[tabhyb$idGenotype==i, "heterosis_type_Fruit"] <- "BelowMeanPar"} else{}}}
  } else{
    CIhybFruitinf <- NA 
    CIhybFruitsup <- NA
    tabhyb[tabhyb$idGenotype==i, "CIhybFruitinf"] <- NA
    tabhyb[tabhyb$idGenotype==i, "CIhybFruitsup"] <- NA
    tabhyb[tabhyb$idGenotype==i, "heterosis_type_Fruit"] <-NA
  }
  
}






#===========================================
# Quantification of MPH and BPH
#===========================================

for(i in levels(tabhyb$idGenotype))
{
  # Lifespan
  maxpar <- max(tabhyb[tabhyb$idGenotype==i, c("MotherLifespan","FatherLifespan")], na.rm=T)
  minpar <- min(tabhyb[tabhyb$idGenotype==i, c("MotherLifespan","FatherLifespan")], na.rm=T)
  meanpar <- mean(as.numeric(tabhyb[tabhyb$idGenotype==i, c("MotherLifespan","FatherLifespan")]), na.rm=T)
  tabhyb[tabhyb$idGenotype==i, "MPH_Lifespan"] <-  (tabhyb[tabhyb$idGenotype==i, "Lifespan"] - meanpar) / meanpar
  tabhyb[tabhyb$idGenotype==i, "BPH_Lifespan"] <-  (tabhyb[tabhyb$idGenotype==i, "Lifespan"] - maxpar) / maxpar
  maxpar <- NA
  minpar <- NA
  meanpar <- NA
  # VegetativeDryMass
  maxpar <- max(tabhyb[tabhyb$idGenotype==i, c("MotherVegetativeDryMass","FatherVegetativeDryMass")], na.rm=T)
  minpar <- min(tabhyb[tabhyb$idGenotype==i, c("MotherVegetativeDryMass","FatherVegetativeDryMass")], na.rm=T)
  meanpar <- mean(as.numeric(tabhyb[tabhyb$idGenotype==i, c("MotherVegetativeDryMass","FatherVegetativeDryMass")]), na.rm=T)
  tabhyb[tabhyb$idGenotype==i, "MPH_VegetativeDryMass"] <-  (tabhyb[tabhyb$idGenotype==i, "VegetativeDryMass"] - meanpar) / meanpar
  tabhyb[tabhyb$idGenotype==i, "BPH_VegetativeDryMass"] <-  (tabhyb[tabhyb$idGenotype==i, "VegetativeDryMass"] - maxpar) / maxpar
  maxpar <- NA
  minpar <- NA
  meanpar <- NA
  # GrowthRate
  maxpar <- max(tabhyb[tabhyb$idGenotype==i, c("MotherGR","FatherGR")], na.rm=T)
  minpar <- min(tabhyb[tabhyb$idGenotype==i, c("MotherGR","FatherGR")], na.rm=T)
  meanpar <- mean(as.numeric(tabhyb[tabhyb$idGenotype==i, c("MotherGR","FatherGR")]), na.rm=T)
  tabhyb[tabhyb$idGenotype==i, "MPH_GR"] <-  (tabhyb[tabhyb$idGenotype==i, "GR"] - meanpar) / meanpar
  tabhyb[tabhyb$idGenotype==i, "BPH_GR"] <-  (tabhyb[tabhyb$idGenotype==i, "GR"] - maxpar) / maxpar
  maxpar <- NA
  minpar <- NA
  meanpar <- NA
  # FruitNumber
  maxpar <- max(tabhyb[tabhyb$idGenotype==i, c("MotherFruit","FatherFruit")], na.rm=T)
  minpar <- min(tabhyb[tabhyb$idGenotype==i, c("MotherFruit","FatherFruit")], na.rm=T)
  meanpar <- mean(as.numeric(tabhyb[tabhyb$idGenotype==i, c("MotherFruit","FatherFruit")]), na.rm=T)
  tabhyb[tabhyb$idGenotype==i, "MPH_Fruit"] <-  (tabhyb[tabhyb$idGenotype==i, "Fruit"] - meanpar) / meanpar
  tabhyb[tabhyb$idGenotype==i, "BPH_Fruit"] <-  (tabhyb[tabhyb$idGenotype==i, "Fruit"] - maxpar) / maxpar
  
}






#===========================================
# Model trait-trait allometric relationships
#===========================================

# Vegetative dry mass vs GrowthRate
#------------------------------------------
tabtot2 <- na.omit(tabtot[,c("idGenotype","Type","VegetativeDryMass","GrowthRate")])
tabtot2 <- aggregate(tabtot2[3:4], by=list(idGenotype=tabtot2$idGenotype, Type=tabtot2$Type), FUN=mean, na.rm=T)
tabtot2 <- droplevels(tabtot2)

# Inbreds
tabt1 <- tabtot2[tabtot2$Type=="Inbred",]
tabt1 <- droplevels(tabt1)
nls1 <- nls(GrowthRate ~ A2*(VegetativeDryMass^(B2+C2*log10(VegetativeDryMass))), data=tabt1, start=list(A2=0.005,B2=1.5,C2=-0.07))
summary(nls1)
nls2 <- nls(GrowthRate ~ A2*(VegetativeDryMass^(B2)), data=tabt1, start=list(A2=0.005,B2=0.5))
summary(nls2)
AIC(nls2) - AIC(nls1)

A2nlsI <- coef(nls1)[1]
B2nlsI <- coef(nls1)[2]
C2nlsI <- coef(nls1)[3]

# Hybrids
tabt1 <- tabtot2[tabtot2$Type=="Hybrid",]
tabt1 <- droplevels(tabt1)
nls1 <- nls(GrowthRate ~ A2*(VegetativeDryMass^(B2+C2*log10(VegetativeDryMass))), data=tabt1, start=list(A2=0.005,B2=1.5,C2=-0.07))
summary(nls1)
nls2 <- nls(GrowthRate ~ A2*(VegetativeDryMass^(B2)), data=tabt1, start=list(A2=0.005,B2=0.5))
summary(nls2)
AIC(nls2) - AIC(nls1)

A2nlsH <- coef(nls1)[1]
B2nlsH <- coef(nls1)[2]
C2nlsH <- coef(nls1)[3]


# VegetativeDryMass vs Fruit number
#------------------------------------------
tabtot2 <- na.omit(tabtot[,c("idGenotype","Type","VegetativeDryMass","FruitNumber")])
tabtot2 <- aggregate(tabtot2[3:4], by=list(idGenotype=tabtot2$idGenotype, Type=tabtot2$Type), FUN=mean, na.rm=T)
tabtot2 <- droplevels(tabtot2)

# Inbreds
tabt1 <- tabtot2[tabtot2$Type=="Inbred",]
tabt1 <- droplevels(tabt1)
nls1 <- nls(FruitNumber ~ VegetativeDryMass/(A3+ B3*VegetativeDryMass+ C3*(VegetativeDryMass^2)), data=tabt1, 
            start=list(A3=2,B3=1.3,C3=0.01),
            nls.control(maxiter=1000, minFactor = 1/100000000000))
summary(nls1)
nls2 <- nls(GrowthRate ~ A2*(VegetativeDryMass^(B2)), data=tabt1, start=list(A2=0.005,B2=0.5))
summary(nls2)
AIC(nls2) - AIC(nls1)

A3nlsI <- coef(nls1)[1]
B3nlsI <- coef(nls1)[2]
C3nlsI <- coef(nls1)[3]

# Hybrids
tabt1 <- tabtot2[tabtot2$Type=="Hybrid",]
tabt1 <- droplevels(tabt1)
nls1 <- nls(FruitNumber ~ VegetativeDryMass/(A3+ B3*VegetativeDryMass+ C3*(VegetativeDryMass^2)), data=tabt1, 
            start=list(A3=2,B3=1.3,C3=0.01),
            nls.control(maxiter=1000, minFactor = 1/100000000000))
summary(nls1)
nls2 <- nls(GrowthRate ~ A2*(VegetativeDryMass^(B2)), data=tabt1, start=list(A2=0.005,B2=0.5))
summary(nls2)
AIC(nls2) - AIC(nls1)

A3nlsH <- coef(nls1)[1]
B3nlsH <- coef(nls1)[2]
C3nlsH <- coef(nls1)[3]









# NLDev from nonlinear equations
#----------------------------------------------

for(i in levels(tabhyb$idGenotype))
{
  meandm <- mean(as.numeric(tabhyb[tabhyb$idGenotype==i, c("MotherVegetativeDryMass","FatherVegetativeDryMass")]), na.rm=T)
  maxdm <- max(as.numeric(tabhyb[tabhyb$idGenotype==i, c("MotherVegetativeDryMass","FatherVegetativeDryMass")]), na.rm=T)
  mindm <- min(as.numeric(tabhyb[tabhyb$idGenotype==i, c("MotherVegetativeDryMass","FatherVegetativeDryMass")]), na.rm=T)
  hybdm <- as.numeric(tabhyb[tabhyb$idGenotype==i, "VegetativeDryMass"])
  tabhyb[tabhyb$idGenotype==i, "VegetativeDryMass_meanPar"] <- meandm
  
  
  # GrowthRate 
  maxpar <- max(tabhyb[tabhyb$idGenotype==i, c("MotherGR","FatherGR")], na.rm=T)
  meanpar <- mean(as.numeric(tabhyb[tabhyb$idGenotype==i, c("MotherGR","FatherGR")]), na.rm=T)
  tabhyb[tabhyb$idGenotype==i, "GR_meanPar"] <- meanpar
  predmeanGR <- mean(c(A2nlsI*(mindm^(B2nlsI+C2nlsI*log10(mindm))),A2nlsI*(maxdm^(B2nlsI+C2nlsI*log10(maxdm)))))
  predmaxGR <- max(c(A2nlsI*(mindm^(B2nlsI+C2nlsI*log10(mindm))),A2nlsI*(maxdm^(B2nlsI+C2nlsI*log10(maxdm)))))
  tabhyb[tabhyb$idGenotype==i, "predHyb_GR_fromMeanParVegetativeDryMass"] <- A2nlsI*(meandm^(B2nlsI+C2nlsI*log10(meandm)))
  tabhyb[tabhyb$idGenotype==i, "predHyb_GR_fromHybVegetativeDryMass"] <- A2nlsI*(hybdm^(B2nlsI+C2nlsI*log10(hybdm)))
  tabhyb[tabhyb$idGenotype==i, "predMPH_GR_fromMeanParVegetativeDryMass"] <- ((A2nlsI*(meandm^(B2nlsI+C2nlsI*log10(meandm)))) - predmeanGR) / predmeanGR
  tabhyb[tabhyb$idGenotype==i, "predBPH_GR_fromMeanParVegetativeDryMass"] <- ((A2nlsI*(meandm^(B2nlsI+C2nlsI*log10(meandm)))) - predmaxGR) / predmaxGR
  tabhyb[tabhyb$idGenotype==i, "predMPH_GR_fromHybVegetativeDryMass"] <- ((A2nlsI*(hybdm^(B2nlsI+C2nlsI*log10(hybdm)))) - predmeanGR) / predmeanGR
  tabhyb[tabhyb$idGenotype==i, "predBPH_GR_fromHybVegetativeDryMass"] <- ((A2nlsI*(hybdm^(B2nlsI+C2nlsI*log10(hybdm)))) - predmaxGR) / predmaxGR
  maxpar <- NA
  meanpar <-NA
  predmaxGR <- NA
  predmeanGR <- NA
  
  # FruitNumber 
  maxpar <- max(tabhyb[tabhyb$idGenotype==i, c("MotherFruit","FatherFruit")], na.rm=T)
  meanpar <- mean(as.numeric(tabhyb[tabhyb$idGenotype==i, c("MotherFruit","FatherFruit")]), na.rm=T)
  tabhyb[tabhyb$idGenotype==i, "Fruit_meanPar"] <- meanpar
  predmeanFruit <- mean(c(mindm/(A3nlsI+ B3nlsI*mindm+ C3nlsI*(mindm^2)),maxdm/(A3nlsI+ B3nlsI*maxdm+ C3nlsI*(maxdm^2))))
  predmaxFruit <- max(c(mindm/(A3nlsI+ B3nlsI*mindm+ C3nlsI*(mindm^2)),maxdm/(A3nlsI+ B3nlsI*maxdm+ C3nlsI*(maxdm^2))))
  tabhyb[tabhyb$idGenotype==i, "predHyb_Fruit_fromMeanParVegetativeDryMass"] <- meandm/(A3nlsI+ B3nlsI*meandm+ C3nlsI*(meandm^2))
  tabhyb[tabhyb$idGenotype==i, "predHyb_Fruit_fromHybVegetativeDryMass"] <- hybdm/(A3nlsI+ B3nlsI*hybdm+ C3nlsI*(hybdm^2))
  tabhyb[tabhyb$idGenotype==i, "predMPH_Fruit_fromMeanParVegetativeDryMass"] <- ((meandm/(A3nlsI+ B3nlsI*meandm+ C3nlsI*(meandm^2))) - predmeanFruit) / predmeanFruit
  tabhyb[tabhyb$idGenotype==i, "predBPH_Fruit_fromMeanParVegetativeDryMass"] <- ((meandm/(A3nlsI+ B3nlsI*meandm+ C3nlsI*(meandm^2))) - predmaxFruit) / predmaxFruit
  tabhyb[tabhyb$idGenotype==i, "predMPH_Fruit_fromHybVegetativeDryMass"] <- ((hybdm/(A3nlsI+ B3nlsI*hybdm+ C3nlsI*(hybdm^2))) - predmeanFruit) / predmeanFruit
  tabhyb[tabhyb$idGenotype==i, "predBPH_Fruit_fromHybVegetativeDryMass"] <- ((hybdm/(A3nlsI+ B3nlsI*hybdm+ C3nlsI*(hybdm^2))) - predmaxFruit) / predmaxFruit
  
  maxpar <- NA
  meanpar <-NA
  meandm <- NA
  maxdm <- NA
  mindm <- NA
  hybdm <- NA
  meanAge50 <- NA
  maxAge50 <- NA
  minAge50 <- NA
  hybAge50 <- NA
  predmaxFruit <- NA
  predmeanFruit <- NA
}



#write.table(tabhyb, "Table_Hybrid_Analyzed.csv", sep=",", dec=".", col.names = T, row.names = F)







#===========================================
#===========================================
#
# 2. Summary of number of hybrids, parents, etc.
#
#===========================================
#===========================================

# hyb with all phenotypes
count1 <- na.omit(tabhyb[,c(1,5:10)])
dim(count1)[1]

# hyb with parental genotypic information
count2 <- na.omit(tabhyb[,c(1,13)])
dim(count2)[1]

# total number of parents
count3 <- levels(as.factor(as.character(c(levels(tabhyb[,2]),levels(tabhyb[,3])))))
length(count3)


# number of genotypes used as both male and female parent
tabtemp <- as.data.frame(cbind(levels(as.factor(as.character(c(levels(tabhyb[,2]),levels(tabhyb[,3]))))), rep(NA,415)))
names(tabtemp) <- c("idGenotype","used_n")
tabtemp$used_n <- as.numeric(tabtemp$used_n)
str(tabtemp)

for(i in levels(tabtemp$idGenotype))
{
  tabtemp[tabtemp$idGenotype==i, "used_Nmother"] <- length(na.omit(tabhyb[tabhyb$idMother==i,"idGenotype"]))
  tabtemp[tabtemp$idGenotype==i, "used_Nfather"] <- length(na.omit(tabhyb[tabhyb$idFather==i,"idGenotype"]))
  
}
tabtemp$used_n <- tabtemp$used_Nmother + tabtemp$used_Nfather

write.table(tabtemp, "Table_Hybrid_Counting.csv", sep=",", dec=".", col.names = T, row.names = F)













#===========================================
#===========================================
#
# 3. MAKE FIGURES
#
#===========================================
#===========================================


colhyb <- rgb(col2rgb("dodgerblue3")[1], col2rgb("dodgerblue3")[2],col2rgb("dodgerblue3")[3], 
              maxColorValue = 255, alpha=70)
colinb <- rgb(col2rgb("firebrick3")[1], col2rgb("firebrick3")[2],col2rgb("firebrick3")[3], 
              maxColorValue = 255, alpha=70)

colMPH <- rgb(col2rgb("lightseagreen")[1], col2rgb("lightseagreen")[2],col2rgb("lightseagreen")[3], 
              maxColorValue = 255, alpha=70)
colBPH <- rgb(col2rgb("navy")[1], col2rgb("navy")[2],col2rgb("navy")[3], 
              maxColorValue = 255, alpha=70)

colmeanpar <- rgb(col2rgb("gold2")[1], col2rgb("gold2")[2],col2rgb("gold2")[3], 
                  maxColorValue = 255, alpha=70)
colhybval <- rgb(col2rgb("chocolate4")[1], col2rgb("chocolate4")[2],col2rgb("chocolate4")[3], 
                 maxColorValue = 255, alpha=70)







#=============================================
#  Box 1
#=============================================

# Import Arabidopsis Data
#----------------------------------------

acc <- read.csv("MeanPheno_accessions_FloSub_GT01.csv", header=T, sep=",") 
for(i in c(1)) {acc[,i] <- as.factor(as.character(acc[,i]))}
for(i in c(2:15)) {acc[,i] <- as.numeric(as.character(acc[,i]))}
acc <- droplevels(acc)
acc$M <- acc$DMmax
acc$G <- acc$GRtot
str(acc)

ril <- read.csv("CVL_AllData_av.csv", header=T, sep=",")   # data accessions collected in MPI
ril <- ril[ril$idCondition=="CTxWW",]
ril$M <- ril$shoot_DM
ril$G <- ril$Growth_rate
ril <- droplevels(ril)
str(ril)


intersp <- read.csv("S1_Data_interSP.csv", header=T, sep=",") 
intersp$M <- intersp$M * 1000000 
intersp$G <- intersp$G * 1000000
str(intersp)



# Figure interspecific
#----------------------------------------

pdf("Box_Fig1_log.pdf", width = 6, height = 6)
plot(log10(intersp$M), log10(intersp$G), ylim=c(-0.6,9),xlim=c(-1.3,11.3),
     pch=16, col="grey60", cex=1.1, axes=F,
     xlab="log[Plant dry mass M (g)]", 
     ylab="log[Growth rate (mg d-1)]")
box()
axis(1, at=seq(-10,20,1), labels=F, tcl=0.5, mgp=c(3,0.5,0))
axis(1, at=seq(0,12,3), labels=10^(seq(0,12,3)), tcl=0.5, mgp=c(3,0.5,0))
axis(1, at=log10(c(seq(0,1,0.1),seq(1,10,1),seq(10,100,10),seq(100,1000,100),
                   seq(1000,10000,1000), seq(10000,100000,10000), seq(100000,1000000,100000),
                   seq(1000000,10000000,1000000),seq(10000000,100000000,10000000),seq(100000000,1000000000,100000000),
                   seq(1000000000,10000000000,1000000000),seq(10000000000,100000000000,10000000000),
                   seq(100000000000,1000000000000,100000000000),seq(1000000000000,10000000000000,1000000000000))),
     labels=F, tcl=0.3, mgp=c(3,0.5,0))
axis(2, at=seq(-10,20,1), labels=F, tcl=0.5, mgp=c(3,0.5,0))
axis(2, at=seq(-2,12,2), labels=10^(seq(-2,12,2)), tcl=0.5, mgp=c(3,0.5,0), las=2)
axis(2, at=log10(c(seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,1,0.1),seq(1,10,1),seq(10,100,10),seq(100,1000,100),
                   seq(1000,10000,1000), seq(10000,100000,10000), seq(100000,1000000,100000),
                   seq(1000000,10000000,1000000),seq(10000000,100000000,10000000),seq(100000000,1000000000,100000000),
                   seq(1000000000,10000000000,1000000000),seq(10000000000,100000000000,10000000000),
                   seq(100000000000,1000000000000,100000000000),seq(1000000000000,10000000000000,1000000000000))),
     labels=F, tcl=0.3, mgp=c(3,0.5,0))

mod1 <- sma(log10(intersp$G) ~ log10(intersp$M))
abline(mod1, col="red", lwd=3)
legend("topleft", legend="log(y) = 0.83 + 0.75log(x)", cex=1, bty="n")
dev.off()




# Figures accessions A. thaliana
#----------------------------------------
acc$M <- log10(acc$M)
acc$G <- log10(acc$G)

pdf("Box_Fig2_log.pdf", width = 6, height = 6)
plot(acc$M, acc$G, xlim=c(1,3.8),ylim=c(-0.5,1.7),
     pch=16, col="grey50", cex=1.1, axes=F,
     xlab="log[Plant dry mass M (g)]", 
     ylab="log[Growth rate (mg d-1)]")
box()
axis(1, at=seq(-10,20,1), labels=F, tcl=0.5, mgp=c(3,0.5,0))
axis(1, at=seq(0,12,1), labels=10^(seq(0,12,1)), tcl=0.5, mgp=c(3,0.5,0))
axis(1, at=log10(c(seq(0,1,0.1),seq(1,10,1),seq(10,100,10),seq(100,1000,100),
                   seq(1000,10000,1000), seq(10000,100000,10000), seq(100000,1000000,100000),
                   seq(1000000,10000000,1000000),seq(10000000,100000000,10000000),seq(100000000,1000000000,100000000),
                   seq(1000000000,10000000000,1000000000),seq(10000000000,100000000000,10000000000),
                   seq(100000000000,1000000000000,100000000000),seq(1000000000000,10000000000000,1000000000000))),
     labels=F, tcl=0.3, mgp=c(3,0.5,0))
axis(2, at=seq(-10,20,1), labels=F, tcl=0.5, mgp=c(3,0.5,0))
axis(2, at=seq(-2,12,1), labels=10^(seq(-2,12,1)), tcl=0.5, mgp=c(3,0.5,0), las=2)
axis(2, at=log10(c(seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,1,0.1),seq(1,10,1),seq(10,100,10),seq(100,1000,100),
                   seq(1000,10000,1000), seq(10000,100000,10000), seq(100000,1000000,100000),
                   seq(1000000,10000000,1000000),seq(10000000,100000000,10000000),seq(100000000,1000000000,100000000),
                   seq(1000000000,10000000000,1000000000),seq(10000000000,100000000000,10000000000),
                   seq(100000000000,1000000000000,100000000000),seq(1000000000000,10000000000000,1000000000000))),
     labels=F, tcl=0.3, mgp=c(3,0.5,0))

mod1 <- sma(acc$G ~ acc$M)
abline(mod1, col="black", lwd=3)

#mod1 <- glm(G ~ poly(M, degree=2, raw=T), family = gaussian, data= acc)
#summary(mod1)
#coef(mod1) 
#curve(coef(mod1)[1] + coef(mod1)[2]*x + coef(mod1)[3]*(x^2), from=1.2, to=3.53, lwd=3, col="red", add=T, lty=1)

dev.off()



# Figures RILs in A. thaliana
#----------------------------------------
ril$M <- log10(ril$M)
ril$G <- log10(ril$G)

pdf("Box_Fig3_log.pdf", width = 6, height = 6)
plot(ril$M, ril$G, xlim=c(0.33,2.8),ylim=c(-0.7,1.55),
     pch=16, col="grey50", cex=1.1, axes=F,
     xlab="log[Plant dry mass M (g)]", 
     ylab="log[Growth rate (mg d-1)]")
box()
axis(1, at=seq(-10,20,1), labels=F, tcl=0.5, mgp=c(3,0.5,0))
axis(1, at=seq(0,12,1), labels=10^(seq(0,12,1)), tcl=0.5, mgp=c(3,0.5,0))
axis(1, at=log10(c(seq(0,1,0.1),seq(1,10,1),seq(10,100,10),seq(100,1000,100),
                   seq(1000,10000,1000), seq(10000,100000,10000), seq(100000,1000000,100000),
                   seq(1000000,10000000,1000000),seq(10000000,100000000,10000000),seq(100000000,1000000000,100000000),
                   seq(1000000000,10000000000,1000000000),seq(10000000000,100000000000,10000000000),
                   seq(100000000000,1000000000000,100000000000),seq(1000000000000,10000000000000,1000000000000))),
     labels=F, tcl=0.3, mgp=c(3,0.5,0))
axis(2, at=seq(-10,20,1), labels=F, tcl=0.5, mgp=c(3,0.5,0))
axis(2, at=seq(-2,12,1), labels=10^(seq(-2,12,1)), tcl=0.5, mgp=c(3,0.5,0), las=2)
axis(2, at=log10(c(seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,1,0.1),seq(1,10,1),seq(10,100,10),seq(100,1000,100),
                   seq(1000,10000,1000), seq(10000,100000,10000), seq(100000,1000000,100000),
                   seq(1000000,10000000,1000000),seq(10000000,100000000,10000000),seq(100000000,1000000000,100000000),
                   seq(1000000000,10000000000,1000000000),seq(10000000000,100000000000,10000000000),
                   seq(100000000000,1000000000000,100000000000),seq(1000000000000,10000000000000,1000000000000))),
     labels=F, tcl=0.3, mgp=c(3,0.5,0))

mod1 <- sma(ril$G ~ ril$M)
abline(mod1, col="black", lwd=3)
#legend("topleft", legend="log(y) = -0.80 + 0.75log(x)", cex=1, bty="n")

#mod1 <- glm(G ~ poly(M, degree=2, raw=T), family = gaussian, data= ril)
#summary(mod1)
#coef(mod1) 
#curve(coef(mod1)[1] + coef(mod1)[2]*x + coef(mod1)[3]*(x^2), from=0.34, to=2.75, lwd=3, col="red", add=T, lty=1)

dev.off()






#=============================================
#  FIG. 1
#=============================================

# Michaelis-Menten curve
#-------------------------------
S <- c(0,1,2,5,8,12,30,50,60,70,80)
v <- c(0,11.1,25.4,44.8,54.5,58.2,72.0,60.1,73, 75, 72)
kinData <- data.frame(S,v)
m1 <- drm(v ~ S, data = kinData, fct = MM.2())

pdf("Fig_MMcurve.pdf", width = 7, height = 5)
par(mar=c(4.5,6,.5,.5), tcl=0)
plot(m1, log = '', type = "n", main = "", xlab="", ylab="", 
     lwd=5, col="darkred")
curve(18+0.77*x, from=1.5, to=70, lwd=1, col="grey30", add=T, lty=1)
axis(2, at=c(-10,69), pos=35.75, col="black", lwd.ticks=0, labels = NA, lty=2)
axis(2, at=c(-10,19.26), pos=1.5, col="black", lwd.ticks=0, labels = NA, lty=2)
axis(2, at=c(-10,72), pos=70, col="black", lwd.ticks=0, labels = NA, lty=2)
axis(1, at=c(-10,1.5), pos=19.3, col="black", lwd.ticks=0, labels = NA, lty=2)
axis(1, at=c(-10,35.75), pos=45.55, col="black", lwd.ticks=0, labels = NA, lty=2)
axis(1, at=c(-10,35.75), pos=68.6, col="black", lwd.ticks=0, labels = NA, lty=2)
axis(1, at=c(-10,70), pos=72.1, col="black", lwd.ticks=0, labels = NA, lty=2)
dev.off()



pdf("Fig_Gaussiancurve.pdf", width = 7, height = 5)
par(mar=c(4.5,6,.5,.5), tcl=0)
curve(7+1.2*x*(exp(x*-0.003)), from=0, to=2000, 
      lty=1, xlab="", ylab="", 
      lwd=5, col="darkred")
curve(55+0.05*x, from=50, to=800, lwd=1, col="grey30", add=T, lty=1)

axis(2, at=c(-10,58), pos=50, col="black", lwd.ticks=0, labels = NA, lty=2)
axis(2, at=c(-10,94), pos=800, col="black", lwd.ticks=0, labels = NA, lty=2)
axis(2, at=c(-10,150), pos=425, col="black", lwd.ticks=0, labels = NA, lty=2)
axis(1, at=c(-100,50), pos=58, col="black", lwd.ticks=0, labels = NA, lty=2)
axis(1, at=c(-100,425), pos=76.2, col="black", lwd.ticks=0, labels = NA, lty=2)
axis(1, at=c(-100,425), pos=150, col="black", lwd.ticks=0, labels = NA, lty=2)
axis(1, at=c(-100,800), pos=95.2, col="black", lwd.ticks=0, labels = NA, lty=2)
dev.off()











#=============================================
#  FIG. 2
#=============================================

# draw crosses on a map
#-------------------------------
hyb <- tab[tab$Type=="Hybrid",]
hyb <- droplevels(hyb)

for(i in 1:dim(hyb)[1])
{
  mother <- as.character(hyb[i,"idMother"])
  father <- as.character(hyb[i,"idFather"])
  hyb[i,"lat_mother"] <- acc[acc$idMother==mother,"latitude"]
  hyb[i,"lon_mother"] <- acc[acc$idMother==mother,"longitude"]
  hyb[i,"lat_father"] <- acc[acc$idFather==father,"latitude"]
  hyb[i,"lon_father"] <- acc[acc$idFather==father,"longitude"]
}
hybcoord <- na.omit(hyb[,c(13:16)])


# function to call
plot_my_connection=function( dep_lon, dep_lat, arr_lon, arr_lat, ...){
  inter <- gcIntermediate(c(dep_lon, dep_lat), c(arr_lon, arr_lat), n=50, addStartEnd=TRUE, breakAtDateLine=F)             
  inter=data.frame(inter)
  diff_of_lon=abs(dep_lon) + abs(arr_lon)
  if(diff_of_lon > 180){
    lines(subset(inter, lon>=0), ...)
    lines(subset(inter, lon<0), ...)
  }else{
    lines(inter, ...)
  }
}



#pdf("map_hybrids.pdf", width = 8, height = 5)
jpeg("map_hybrids.jpg", width = 4000, height = 2400, quality = 100, res=800)
par(mar=c(0,0,0,0))
map('world',col="grey80", fill=TRUE, bg="white", lwd=0.05,mar=rep(0,4),
    border=0, ylim=c(13,67), xlim=c(-30,140))
points(acc$longitude, acc$latitude, pch=16, col="firebrick3", cex=0.5)

for(i in 1:dim(hyb)[1])
{
  plot_my_connection(hybcoord[i,2], hybcoord[i,1], 
                     hybcoord[i,4], hybcoord[i,3], 
                     col="dodgerblue3", lwd=0.35, lty=1)
}
dev.off()













#-----------------------------------------------------------------
#-----------------------------------------------------------------
# FIG. 3
#-----------------------------------------------------------------
#-----------------------------------------------------------------

# VegetativeDryMass
p2_1 <- ggplot(tab, aes(VegetativeDryMass))+ theme_classic() +
  geom_density(aes(fill=Type), alpha=0.55, bw=150) +
  scale_fill_manual(values=c("dodgerblue3", "firebrick3")) +
  theme(axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), plot.title=element_text(size=20)) +
  scale_x_continuous(expression(paste("Vegetative dry mass (mg)")),
                     breaks=seq(0,2000,1000), labels=round(seq(0,2000,1000), digits=2), limits=c(0,2000)) +
  scale_y_continuous(expression(paste("Density")))
pdf("distrib_VegetativeDryMass.pdf", width = 6, height = 4)
p2_1
dev.off()


# model genetic effect
mean(tab[tab$Type=="Inbred", "VegetativeDryMass"], na.rm=T)
sd(tab[tab$Type=="Inbred", "VegetativeDryMass"], na.rm=T)
mean(tab[tab$Type=="Hybrid", "VegetativeDryMass"], na.rm=T)
sd(tab[tab$Type=="Hybrid", "VegetativeDryMass"], na.rm=T)

tab2 <- na.omit(tab[,c("idGenotype","Type","VegetativeDryMass")])
tab2 <- droplevels(tab2)
mod <- lmer(VegetativeDryMass ~ 1 + (1|idGenotype), data = tab2)
summary(mod)

mod <- lm(aov(VegetativeDryMass ~ Type, data=tab))
summary(mod)



# Age at reproduction (= Lifespan)
p2_1 <- ggplot(tab, aes(Lifespan))+ theme_classic() +
  geom_density(aes(fill=Type), alpha=0.55, bw=9) +
  scale_fill_manual(values=c("dodgerblue3", "firebrick3")) +
  theme(axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), plot.title=element_text(size=20)) +
  scale_x_continuous(expression(paste("Plant lifespan (days)")),
                     breaks=seq(0,150,50), labels=round(seq(0,150,50), digits=2), limits=c(0,170)) +
  scale_y_continuous(expression(paste("Density")))
pdf("distrib_Lifespan.pdf", width = 6, height = 4)
p2_1
dev.off()

# model genetic effect
mean(tab[tab$Type=="Inbred", "Lifespan"], na.rm=T)
sd(tab[tab$Type=="Inbred", "Lifespan"], na.rm=T)
mean(tab[tab$Type=="Hybrid", "Lifespan"], na.rm=T)
sd(tab[tab$Type=="Hybrid", "Lifespan"], na.rm=T)

tab2 <- na.omit(tab[,c("idGenotype","Type","Lifespan")])
tab2 <- droplevels(tab2)
mod <- lmer(Lifespan ~ 1 + (1|idGenotype), data = tab2)
summary(mod)

mod <- lm(aov(Lifespan ~ Type, data=tab))
summary(mod)




# GrowthRate
p2_1 <- ggplot(tab, aes(GrowthRate))+ theme_classic() +
  geom_density(aes(fill=Type), alpha=0.55, bw=3) +
  scale_fill_manual(values=c("dodgerblue3", "firebrick3")) +
  theme(axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), plot.title=element_text(size=20)) +
  scale_x_continuous(expression(paste("Growth rate (mg d"^-1,")")),
                     breaks=seq(0,40,10), labels=round(seq(0,40,10), digits=2), limits=c(0,40)) +
  scale_y_continuous(expression(paste("Density")))
pdf("distrib_GrowthRate.pdf", width = 6, height = 4)
p2_1
dev.off()


# model genetic effect
mean(tab[tab$Type=="Inbred", "GrowthRate"], na.rm=T)
sd(tab[tab$Type=="Inbred", "GrowthRate"], na.rm=T)
mean(tab[tab$Type=="Hybrid", "GrowthRate"], na.rm=T)
sd(tab[tab$Type=="Hybrid", "GrowthRate"], na.rm=T)

tab2 <- na.omit(tab[,c("idGenotype","Type","GrowthRate")])
tab2 <- droplevels(tab2)
mod <- lmer(GrowthRate ~ 1 + (1|idGenotype), data = tab2)
summary(mod)

mod <- lm(aov(GrowthRate ~ Type, data=tab))
summary(mod)




# Fruit
p2_1 <- ggplot(tab, aes(FruitNumber))+ theme_classic() +
  geom_density(aes(fill=Type), alpha=0.55, bw=25) +
  scale_fill_manual(values=c("dodgerblue3", "firebrick3")) +
  theme(axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), plot.title=element_text(size=20)) +
  scale_x_continuous(expression(paste("Fruit number")),
                     breaks=seq(0,350,50), labels=round(seq(0,350,50), digits=2), limits=c(0,370)) +
  scale_y_continuous(expression(paste("Density")))
pdf("distrib_Fruit.pdf", width = 6, height = 4)
p2_1
dev.off()

# model genetic effect
mean(tab[tab$Type=="Inbred", "FruitNumber"], na.rm=T)
sd(tab[tab$Type=="Inbred", "FruitNumber"], na.rm=T)
mean(tab[tab$Type=="Hybrid", "FruitNumber"], na.rm=T)
sd(tab[tab$Type=="Hybrid", "FruitNumber"], na.rm=T)

tab2 <- na.omit(tab[,c("idGenotype","Type","FruitNumber")])
tab2 <- droplevels(tab2)
mod <- lmer(FruitNumber ~ 1 + (1|idGenotype), data = tab2)
summary(mod)

mod <- lm(aov(FruitNumber ~ Type, data=tab))
summary(mod)




# Pie chart
tabhyb <- read.csv("Table_Hybrid_Analyzed.csv", header=T, sep=",")
for(i in (c(4:7,12:14,16:27, 32:62)))
{
  tabhyb[,i] <- as.numeric(as.character(tabhyb[,i])) 
}
for(i in (c(1:3,8:11,15, 28:31)))
{
  tabhyb[,i] <- as.factor(as.character(tabhyb[,i])) 
}
str(tabhyb)

tabhyb <- na.omit(tabhyb[,c(1,28:31)])
tabhyb <- droplevels(tabhyb)
str(tabhyb)


# Age at reproduction
percent_negBPH <- 100*(length((tabhyb[tabhyb$heterosis_type_Lifespan=="BelowWorstPar", "heterosis_type_Lifespan"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_Lifespan"])))
percent_negMPH <- 100*(length((tabhyb[tabhyb$heterosis_type_Lifespan=="BelowMeanPar", "heterosis_type_Lifespan"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_Lifespan"])))
percent_posMPH <- 100*(length((tabhyb[tabhyb$heterosis_type_Lifespan=="AboveMeanPar", "heterosis_type_Lifespan"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_Lifespan"])))
percent_posBPH <- 100*(length((tabhyb[tabhyb$heterosis_type_Lifespan=="AboveBestPar", "heterosis_type_Lifespan"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_Lifespan"])))
percent_noHet <- 100*(length((tabhyb[tabhyb$heterosis_type_Lifespan=="null", "heterosis_type_Lifespan"])) / 
                        length(na.omit(tabhyb[, "heterosis_type_Lifespan"])))

pdf("PieChart_FT.pdf", width = 4, height = 4)
slices <- c(percent_noHet, percent_negBPH, percent_negMPH, percent_posMPH, percent_posBPH) 
lbls <- c(paste(as.character(round(percent_noHet, digits = 0)), "%", sep=""), 
          paste(as.character(round(percent_negBPH, digits = 0)), "%", sep=""),
          paste(as.character(round(percent_negMPH, digits = 0)), "%", sep=""),
          paste(as.character(round(percent_posMPH, digits = 0)), "%", sep=""),
          paste(as.character(round(percent_posBPH, digits = 0)), "%", sep=""))
pie3D(slices,labels=lbls,explode=0.1, theta = 1.15, start=0,
      main="", col=c("azure3","coral4","coral","darkolivegreen3","darkolivegreen"))
dev.off()


# DM
percent_negBPH <- 100*(length((tabhyb[tabhyb$heterosis_type_VegetativeDryMass=="BelowWorstPar", "heterosis_type_VegetativeDryMass"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_VegetativeDryMass"])))
percent_negMPH <- 100*(length((tabhyb[tabhyb$heterosis_type_VegetativeDryMass=="BelowMeanPar", "heterosis_type_VegetativeDryMass"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_VegetativeDryMass"])))
percent_posMPH <- 100*(length((tabhyb[tabhyb$heterosis_type_VegetativeDryMass=="AboveMeanPar", "heterosis_type_VegetativeDryMass"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_VegetativeDryMass"])))
percent_posBPH <- 100*(length((tabhyb[tabhyb$heterosis_type_VegetativeDryMass=="AboveBestPar", "heterosis_type_VegetativeDryMass"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_VegetativeDryMass"])))
percent_noHet <- 100*(length((tabhyb[tabhyb$heterosis_type_VegetativeDryMass=="null", "heterosis_type_VegetativeDryMass"])) / 
                        length(na.omit(tabhyb[, "heterosis_type_VegetativeDryMass"])))

pdf("PieChart_DM.pdf", width = 4, height = 4)
slices <- c(percent_noHet, percent_negBPH, percent_negMPH, percent_posMPH, percent_posBPH) 
lbls <- c(paste(as.character(round(percent_noHet, digits = 0)), "%", sep=""), 
          paste(as.character(round(percent_negBPH, digits = 0)), "%", sep=""),
          paste(as.character(round(percent_negMPH, digits = 0)), "%", sep=""),
          paste(as.character(round(percent_posMPH, digits = 0)), "%", sep=""),
          paste(as.character(round(percent_posBPH, digits = 0)), "%", sep=""))
pie3D(slices,labels=lbls,explode=0.1, theta = 1.15, start=0,
      main="", col=c("azure3","coral4","coral","darkolivegreen3","darkolivegreen"))
dev.off()


# GR
percent_negBPH <- 100*(length((tabhyb[tabhyb$heterosis_type_GR=="BelowWorstPar", "heterosis_type_GR"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_GR"])))
percent_negMPH <- 100*(length((tabhyb[tabhyb$heterosis_type_GR=="BelowMeanPar", "heterosis_type_GR"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_GR"])))
percent_posMPH <- 100*(length((tabhyb[tabhyb$heterosis_type_GR=="AboveMeanPar", "heterosis_type_GR"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_GR"])))
percent_posBPH <- 100*(length((tabhyb[tabhyb$heterosis_type_GR=="AboveBestPar", "heterosis_type_GR"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_GR"])))
percent_noHet <- 100*(length((tabhyb[tabhyb$heterosis_type_GR=="null", "heterosis_type_GR"])) / 
                        length(na.omit(tabhyb[, "heterosis_type_GR"])))

pdf("PieChart_GR.pdf", width = 4, height = 4)
slices <- c(percent_noHet, percent_negBPH, percent_negMPH, percent_posMPH, percent_posBPH) 
lbls <- c(paste(as.character(round(percent_noHet, digits = 0)), "%", sep=""), 
          paste(as.character(round(percent_negBPH, digits = 0)), "%", sep=""),
          paste(as.character(round(percent_negMPH, digits = 0)), "%", sep=""),
          paste(as.character(round(percent_posMPH, digits = 0)), "%", sep=""),
          paste(as.character(round(percent_posBPH, digits = 0)), "%", sep=""))
pie3D(slices,labels=lbls,explode=0.1, theta = 1.15, start=0,
      main="", col=c("azure3","coral4","coral","darkolivegreen3","darkolivegreen"))
dev.off()


# FN
percent_negBPH <- 100*(length((tabhyb[tabhyb$heterosis_type_Fruit=="BelowWorstPar", "heterosis_type_Fruit"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_Fruit"])))
percent_negMPH <- 100*(length((tabhyb[tabhyb$heterosis_type_Fruit=="BelowMeanPar", "heterosis_type_Fruit"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_Fruit"])))
percent_posMPH <- 100*(length((tabhyb[tabhyb$heterosis_type_Fruit=="AboveMeanPar", "heterosis_type_Fruit"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_Fruit"])))
percent_posBPH <- 100*(length((tabhyb[tabhyb$heterosis_type_Fruit=="AboveBestPar", "heterosis_type_Fruit"])) / 
                         length(na.omit(tabhyb[, "heterosis_type_Fruit"])))
percent_noHet <- 100*(length((tabhyb[tabhyb$heterosis_type_Fruit=="null", "heterosis_type_Fruit"])) / 
                        length(na.omit(tabhyb[, "heterosis_type_Fruit"])))

pdf("PieChart_FN.pdf", width = 4, height = 4)
slices <- c(percent_noHet, percent_negBPH, percent_negMPH, percent_posMPH, percent_posBPH) 
lbls <- c(paste(as.character(round(percent_noHet, digits = 0)), "%", sep=""), 
          paste(as.character(round(percent_negBPH, digits = 0)), "%", sep=""),
          paste(as.character(round(percent_negMPH, digits = 0)), "%", sep=""),
          paste(as.character(round(percent_posMPH, digits = 0)), "%", sep=""),
          paste(as.character(round(percent_posBPH, digits = 0)), "%", sep=""))
pie3D(slices,labels=lbls,explode=0.1, theta = 1.15, start=0,
      main="", col=c("azure3","coral4","coral","darkolivegreen3","darkolivegreen"))
dev.off()












#-----------------------------------------------------------------
#-----------------------------------------------------------------
# FIG. 4
#-----------------------------------------------------------------
#----------------------------------------------------------------

pdf("VegetativeDryMass_vs_GR.pdf", width = 6, height = 5.5)
#jpeg("GenDist_vs_GeoDist.jpg", width = 4000, height = 3400, quality = 100, res=700)
par(mar=c(4,4,.5,.5), tcl=0.3)
plot(tabtot2$VegetativeDryMass, tabtot2$GrowthRate, 
     xlab="Plant dry mass", ylab="Growth rate", type="n")
points(tabtot2[tabtot2$Type=="Inbred", "VegetativeDryMass"], 
       tabtot2[tabtot2$Type=="Inbred", "GrowthRate"], 
       pch=16, col=colinb, cex=1.8)
points(tabtot2[tabtot2$Type=="Hybrid", "VegetativeDryMass"], 
       tabtot2[tabtot2$Type=="Hybrid", "GrowthRate"], 
       pch=16, col=colhyb, cex=1.8)

# Inbreds
tabt1 <- tabtot2[tabtot2$Type=="Inbred",]
tabt1 <- droplevels(tabt1)
nls1 <- nls(GrowthRate ~ A2*(VegetativeDryMass^(B2+C2*log10(VegetativeDryMass))), data=tabt1, start=list(A2=0.005,B2=1.5,C2=-0.07))
summary(nls1)
A2nlsI <- coef(nls1)[1]
B2nlsI <- coef(nls1)[2]
C2nlsI <- coef(nls1)[3]
curve.S <- function(dm, A2, B2, C2) { A2*(dm^(B2+C2*log10(dm)))}
curve(curve.S(dm = x, A2 = A2nlsI, B2 = B2nlsI, C2=C2nlsI), 0, 
      max(tabt1[,"VegetativeDryMass"]), lwd = 2, add=T, col="black", lty=1)

# Hybrids
tabt1 <- tabtot2[tabtot2$Type=="Hybrid",]
tabt1 <- droplevels(tabt1)
nls1 <- nls(GrowthRate ~ A2*(VegetativeDryMass^(B2+C2*log10(VegetativeDryMass))), data=tabt1, start=list(A2=0.005,B2=1.5,C2=-0.07))
summary(nls1)
A2nlsH <- coef(nls1)[1]
B2nlsH <- coef(nls1)[2]
C2nlsH <- coef(nls1)[3]
curve.S <- function(dm, A2, B2, C2) { A2*(dm^(B2+C2*log10(dm)))}
curve(curve.S(dm = x, A2 = A2nlsH, B2 = B2nlsH, C2=C2nlsH), 0, 
      max(tabt1[,"VegetativeDryMass"]), lwd = 3, add=T, col="black", lty=2)

legend("topleft", legend=c("Inbred","Hybrid"), pch=16,
       col =c(colinb, colhyb), 
       lty=c(1,2), lwd=c(3,2), bty="n", cex=1.4, inset=c(0.05,0))
dev.off()




# VegetativeDryMass vs Fruit number
#------------------------------------------
tabtot2 <- na.omit(tabtot[,c("idGenotype","Type","VegetativeDryMass","FruitNumber")])
tabtot2 <- aggregate(tabtot2[3:4], by=list(idGenotype=tabtot2$idGenotype, Type=tabtot2$Type), FUN=mean, na.rm=T)
tabtot2 <- droplevels(tabtot2)

pdf("VegetativeDryMass_vs_Fruit.pdf", width = 6, height = 5.5)
#jpeg("GenDist_vs_GeoDist.jpg", width = 4000, height = 3400, quality = 100, res=700)
par(mar=c(4,4,.5,.5), tcl=0.3)
plot(tabtot2$VegetativeDryMass, tabtot2$FruitNumber, 
     xlab="Vegetative dry mass", ylab="Fruit number", type="n")
points(tabtot2[tabtot2$Type=="Inbred", "VegetativeDryMass"], 
       tabtot2[tabtot2$Type=="Inbred", "FruitNumber"], 
       pch=16, col=colinb, cex=1.8)
points(tabtot2[tabtot2$Type=="Hybrid", "VegetativeDryMass"], 
       tabtot2[tabtot2$Type=="Hybrid", "FruitNumber"], 
       pch=16, col=colhyb, cex=1.8)

# Inbreds
tabt1 <- tabtot2[tabtot2$Type=="Inbred",]
tabt1 <- droplevels(tabt1)
nls1 <- nls(FruitNumber ~ VegetativeDryMass/(A3+ B3*VegetativeDryMass+ C3*(VegetativeDryMass^2)), data=tabt1, 
            start=list(A3=2,B3=1.3,C3=0.01),
            nls.control(maxiter=1000, minFactor = 1/100000000000))
summary(nls1)
A3nlsI <- coef(nls1)[1]
B3nlsI <- coef(nls1)[2]
C3nlsI <- coef(nls1)[3]
curve.S <- function(dm, A3, B3, C3) {dm/(A3+ B3*dm+ C3*(dm^2))}
curve(curve.S(dm = x, A3 = A3nlsI, B3 = B3nlsI, C3=C3nlsI), 20, 
      max(tabt1[,"VegetativeDryMass"]), lwd = 2, add=T, col="black", lty=1)

# Hybrids
tabt1 <- tabtot2[tabtot2$Type=="Hybrid",]
tabt1 <- droplevels(tabt1)
nls1 <- nls(FruitNumber ~ VegetativeDryMass/(A3+ B3*VegetativeDryMass+ C3*(VegetativeDryMass^2)), data=tabt1, 
            start=list(A3=2,B3=1.3,C3=0.01),
            nls.control(maxiter=1000, minFactor = 1/100000000000))
summary(nls1)
A3nlsH <- coef(nls1)[1]
B3nlsH <- coef(nls1)[2]
C3nlsH <- coef(nls1)[3]
curve.S <- function(dm, A3, B3, C3) {dm/(A3+ B3*dm+ C3*(dm^2))}
curve(curve.S(dm = x, A3 = A3nlsH, B3 = B3nlsH, C3=C3nlsH), 20, 
      max(tabt1[,"VegetativeDryMass"]), lwd = 3, add=T, col="black", lty=2)

legend("topright", legend=c("Inbred","Hybrid"), pch=16,
       col =c(colinb, colhyb), 
       lty=c(1,2), lwd=c(3,2), bty="n", cex=1.4, inset=c(0.05,0))
dev.off()







#-----------------------------------------------------------------
#-----------------------------------------------------------------
# FIG. 5
#-----------------------------------------------------------------
#----------------------------------------------------------------

pdf("Fig5_corrDist.pdf", width=6, height = 5)
par(mfrow=c(2,2), oma=c(2,2,0,0), mar=c(.5,.5,.5,.5))

# A
plot(log10(tabhyb$GeneticDist), tabhyb$MPH_GR, 
     type="n", las=2,
     xlab="Genetic distance", ylab="Heterosis on Growth rate",
     xlim=c(3.85,4.65), ylim=c(-1.05,3))
axis(1, at=c(-100,100), labels=NULL, lwd.ticks = 0, pos=0, lwd=1, lty=2, col="firebrick3")
points(log10(tabhyb$GeneticDist), tabhyb$BPH_GR, 
       pch=16, col=colBPH, cex=1.3)
points(log10(tabhyb$GeneticDist), tabhyb$MPH_GR, 
       pch=16, col=colMPH, cex=1.3)
mod1 <- sma(tabhyb$MPH_GR ~ log10(tabhyb$GeneticDist))
#mod2 <- lm(tabhyb$MPH_GR ~ log10(tabhyb$GeneticDist))
mod2 <- lm(tabhyb$MPH_GR ~ poly(log10(tabhyb$GeneticDist), degree=2, raw=T))
mod3 <- sma(tabhyb$BPH_GR ~ log10(tabhyb$GeneticDist))
#mod4 <- lm(tabhyb$BPH_GR ~ log10(tabhyb$GeneticDist))
mod4 <- lm(tabhyb$BPH_GR ~ poly(log10(tabhyb$GeneticDist), degree=2, raw=T))

mod5 <- slope.test(y=tabhyb$BPH_GR, x=log10(tabhyb$GeneticDist), test.value=4.039655)
mod5

#curve(coef(mod1)[1] + coef(mod1)[2]*x, from=4.17, 
#      to=max(log10(tabhyb$GeneticDist), na.rm=T), lwd=2, add=T, col="lightseagreen")
curve(coef(mod2)[1] + coef(mod2)[2]*x + coef(mod2)[3]*(x^2), from=3.92, 
      to=max(log10(tabhyb$GeneticDist), na.rm=T), lwd=2, add=T, col="lightseagreen")
#curve(coef(mod1)[1] + coef(mod1)[2]*x, from=4.17, 
#      to=max(log10(tabhyb$GeneticDist), na.rm=T), lwd=2, add=T, col="lightseagreen")
curve(coef(mod4)[1] + coef(mod4)[2]*x + coef(mod4)[3]*(x^2), from=3.92, 
      to=max(log10(tabhyb$GeneticDist), na.rm=T), lwd=2, add=T, col="navy")
summary(mod2)
summary(mod4)

legend("topleft",  bty="n", cex=1.5,
       legend=c(expression(paste("MPH: ", italic(r), ""^2, " = 0.07***")), 
                expression(paste("BPH: ", italic(r), ""^2, " = 0.06***"))),
       text.col = c("lightseagreen","navy"), inset=c(-0.09,-0.03))


# B
plot(log10(tabhyb$distVegetativeDryMass), tabhyb$MPH_GR, 
     type="n", las=2,
     xlab="Phenotypic distance", ylab="Heterosis on Growth rate",
     xlim=c(-1,3.2), ylim=c(-1.05,3))
axis(1, at=c(-100,100), labels=NULL, lwd.ticks = 0, pos=0, lwd=1, lty=2, col="firebrick3")
points(log10(tabhyb$distVegetativeDryMass), tabhyb$BPH_GR, 
       pch=16, col=colBPH, cex=1.3)
points(log10(tabhyb$distVegetativeDryMass), tabhyb$MPH_GR, 
       pch=16, col=colMPH, cex=1.3)

mod1 <- sma(tabhyb$MPH_GR ~ log10(tabhyb$distVegetativeDryMass))
#mod2 <- lm(tabhyb$MPH_GR ~ log10(tabhyb$))
mod2 <- lm(tabhyb$MPH_GR ~ poly(log10(tabhyb$distVegetativeDryMass), degree=2, raw=T))
confint(mod2)
mod3 <- sma(tabhyb$BPH_GR ~ log10(tabhyb$distVegetativeDryMass))
#mod4 <- lm(tabhyb$BPH_GR ~ log10(tabhyb$distVegetativeDryMass))
mod4 <- lm(tabhyb$BPH_GR ~ poly(log10(tabhyb$distVegetativeDryMass), degree=2, raw=T))
confint(mod4)
summary(mod2)
summary(mod4)
mod5 <- slope.test(y=tabhyb$BPH_GR, x=log10(tabhyb$distVegetativeDryMass), test.value=0.5619933)
mod5

#curve(coef(mod1)[1] + coef(mod1)[2]*x, from=0.7, 
#      to=max(log10(tabhyb$distVegetativeDryMass), na.rm=T), lwd=2, add=T, col="lightseagreen")
#curve(coef(mod3)[1] + coef(mod3)[2]*x, from=4.19, 
#      to=max(log10(tabhyb$distVegetativeDryMass), na.rm=T), lwd=2, add=T, col="navy")
curve(coef(mod2)[1] + coef(mod2)[2]*x + coef(mod2)[3]*(x^2), from=-0.52, 
      to=max(log10(tabhyb$distVegetativeDryMass), na.rm=T), lwd=2, add=T, col="lightseagreen")
curve(coef(mod4)[1] + coef(mod4)[2]*x + coef(mod4)[3]*(x^2), from=-0.52, 
      to=max(log10(tabhyb$distVegetativeDryMass), na.rm=T), lwd=2, lty=2, add=T, col="navy")

legend("topleft",  bty="n", cex=1.5,
       legend=c(expression(paste("MPH: ", italic(r), ""^2, " = 0.05***")), 
                expression(paste("BPH: ", italic(r), ""^2, " = 0.009"^NS))),
       text.col = c("lightseagreen","navy"), inset=c(-0.09,-0.03))



# C
plot(log10(tabhyb$GeneticDist), tabhyb$MPH_Fruit, 
     type="n", las=2,
     xlab="Genetic distance", ylab="Heterosis on Fruit number",
     xlim=c(3.85,4.7), ylim=c(-1,2.3))
axis(1, at=c(-100,100), labels=NULL, lwd.ticks = 0, pos=0, lwd=1, lty=2, col="firebrick3")
points(log10(tabhyb$GeneticDist), tabhyb$MPH_Fruit, 
       pch=16, col=colMPH, cex=1.3)
points(log10(tabhyb$GeneticDist), tabhyb$BPH_Fruit, 
       pch=16, col=colBPH, cex=1.3)
mod1 <- sma(tabhyb$MPH_Fruit ~ log10(tabhyb$GeneticDist))
#mod2 <- lm(tabhyb$MPH_Fruit ~ log10(tabhyb$GeneticDist))
mod2 <- lm(tabhyb$MPH_Fruit ~ poly(log10(tabhyb$GeneticDist), degree=2, raw=T))
mod3 <- sma(tabhyb$BPH_Fruit ~ log10(tabhyb$GeneticDist))
#mod4 <- lm(tabhyb$BPH_Fruit ~ log10(tabhyb$GeneticDist))
mod4 <- lm(tabhyb$BPH_Fruit ~ poly(log10(tabhyb$GeneticDist), degree=2, raw=T))
mod5 <- slope.test(y=tabhyb$BPH_Fruit, x=log10(tabhyb$GeneticDist), test.value=-3.759696)
mod5


curve(coef(mod2)[1] + coef(mod2)[2]*x + coef(mod2)[3]*(x^2),lty=2, from=3.91, 
      to=4.63, lwd=2, add=T, col="lightseagreen")
curve(coef(mod4)[1] + coef(mod4)[2]*x + coef(mod4)[3]*(x^2), lty=2, from=3.91, 
      to=4.62, lwd=2, add=T, col="navy")
summary(mod2)
summary(mod4)

legend("topleft",  bty="n", cex=1.5,
       legend=c(expression(paste("MPH: ", italic(r), ""^2, " = 0.007"^NS)), 
                expression(paste("BPH: ", italic(r), ""^2, " = 0.01*"))),
       text.col = c("lightseagreen","navy"), inset=c(-0.09,-0.03))


# D
plot(log10(tabhyb$distVegetativeDryMass), tabhyb$MPH_Fruit, 
     type="n", las=2,
     xlab="Phenotypic distance", ylab="Heterosis on Growth rate",
     xlim=c(-1,3.2), ylim=c(-1,2.3))
axis(1, at=c(-100,100), labels=NULL, lwd.ticks = 0, pos=0, lwd=1, lty=2, col="firebrick3")
points(log10(tabhyb$distVegetativeDryMass), tabhyb$MPH_Fruit, 
       pch=16, col=colMPH, cex=1.3)
points(log10(tabhyb$distVegetativeDryMass), tabhyb$BPH_Fruit, 
       pch=16, col=colBPH, cex=1.3)

mod1 <- sma(tabhyb$MPH_Fruit ~ log10(tabhyb$distVegetativeDryMass))
#mod2 <- lm(tabhyb$MPH_Fruit ~ log10(tabhyb$distVegetativeDryMass))
mod2 <- lm(tabhyb$MPH_Fruit ~ poly(log10(tabhyb$distVegetativeDryMass), degree=2, raw=T))
confint(mod2)
mod3 <- sma(tabhyb$BPH_Fruit ~ log10(tabhyb$distVegetativeDryMass))
#mod4 <- lm(tabhyb$BPH_Fruit ~ log10(tabhyb$distVegetativeDryMass))
mod4 <- lm(tabhyb$BPH_Fruit ~ poly(log10(tabhyb$distVegetativeDryMass), degree=2, raw=T))
confint(mod4)
summary(mod2)
summary(mod4)

mod5 <- slope.test(y=tabhyb$BPH_Fruit, x=log10(tabhyb$distVegetativeDryMass), test.value=-0.5022512)
mod5

curve(coef(mod2)[1] + coef(mod2)[2]*x + coef(mod2)[3]*(x^2), lty=2, from=-0.44, 
      to=3.1, lwd=2, add=T, col="lightseagreen")
curve(coef(mod4)[1] + coef(mod4)[2]*x + coef(mod4)[3]*(x^2), from=-0.44, 
      to=3.1, lwd=2, add=T, col="navy")

legend("topleft",  bty="n", cex=1.5,
       legend=c(expression(paste("MPH: ", italic(r), ""^2, " = 0.007"^NS)), 
                expression(paste("BPH: ", italic(r), ""^2, " = 0.09***"))),
       text.col = c("lightseagreen","navy"), inset=c(-0.09,-0.03))

dev.off()











#-----------------------------------------------------------------
#-----------------------------------------------------------------
# FIG. 6
#-----------------------------------------------------------------
#-----------------------------------------------------------------

pdf("Fig6.pdf", width=7, height = 3)
par(mfrow=c(2,2), oma=c(2,0,0,0), mar=c(.5,4.2,.5,.5))

# A
tabtot2 <- na.omit(tabtot[,c("idGenotype","Type","VegetativeDryMass","GrowthRate")])
tabtot2 <- droplevels(tabtot2)
str(tabtot2)

#pdf("Fig6A.pdf", width = 6, height = 5.5)
#jpeg("GenDist_vs_GeoDist.jpg", width = 4000, height = 3400, quality = 100, res=700)
#par(mar=c(4,4,.5,.5), tcl=0.3)
plot(tabtot2$VegetativeDryMass, tabtot2$GrowthRate, 
     xlab="Plant dry mass", ylab="Growth rate", type="n", ylim=c(0,45))
points(tabtot2[tabtot2$Type=="Hybrid", "VegetativeDryMass"], 
       tabtot2[tabtot2$Type=="Hybrid", "GrowthRate"], 
       pch=16, col="grey80", cex=1)
tabt1 <- tabtot2[tabtot2$Type=="Hybrid",]
tabt1 <- droplevels(tabt1)
nls1 <- nls(GrowthRate ~ A2*(VegetativeDryMass^(B2+C2*log10(VegetativeDryMass))), data=tabt1, start=list(A2=0.005,B2=1.5,C2=-0.07))
summary(nls1)
A2nlsH <- coef(nls1)[1]
B2nlsH <- coef(nls1)[2]
C2nlsH <- coef(nls1)[3]
curve.S <- function(dm, A2, B2, C2) { A2*(dm^(B2+C2*log10(dm)))}
curve(curve.S(dm = x, A2 = A2nlsH, B2 = B2nlsH, C2=C2nlsH), 0, 
      max(tabt1[,"VegetativeDryMass"]), lwd = 1, add=T, col="grey50", lty=1)

points(tabtot[tabtot$idExp=="GT01" & tabtot$idPot=="58", "VegetativeDryMass"], #611
       tabtot[tabtot$idExp=="GT01" & tabtot$idPot=="58", "GrowthRate"], 
       pch=16, col="firebrick3", cex=1.5)
points(tabtot[tabtot$idExp=="GT05" & tabtot$idPot=="299", "VegetativeDryMass"], #406
       tabtot[tabtot$idExp=="GT05" & tabtot$idPot=="299", "GrowthRate"], 
       pch=16, col="firebrick3", cex=1.5)
points(tabtot[tabtot$idExp=="GT05" & tabtot$idPot=="476", "VegetativeDryMass"], #1167
       tabtot[tabtot$idExp=="GT05" & tabtot$idPot=="476", "GrowthRate"], 
       pch=16, col="dodgerblue3", cex=1.5)


legend("topleft", legend=c("Phenotypic non-linearity = 0.22", "MPH = 0.18"), 
       bty="n", cex=1.3, inset=c(-0.05,0))
dev.off()





# B
plot(tabhyb$predMPH_GR_fromHybVegetativeDryMass, tabhyb$MPH_GR, 
     type="n", las=2,
     xlab="Non-linear deviation", ylab="Heterosis on Growth rate",
     xlim=c(-1.05,1.6), ylim=c(-1.08,2.4))
axis(1, at=c(-100,100), labels=NULL, lwd.ticks = 0, pos=0, lwd=1, lty=2, col="firebrick3")
abline(0,1, lty=2)
points(tabhyb$predMPH_GR_fromHybVegetativeDryMass, tabhyb$MPH_GR, 
       pch=16, col=colMPH, cex=1.3)
points(tabhyb$predBPH_GR_fromHybVegetativeDryMass, tabhyb$BPH_GR, 
       pch=16, col=colBPH, cex=1.3)
mod2 <- lm(tabhyb$MPH_GR ~ tabhyb$predMPH_GR_fromHybVegetativeDryMass)
mod4 <- lm(tabhyb$BPH_GR ~ tabhyb$predBPH_GR_fromHybVegetativeDryMass)
summary(mod2)
summary(mod4)

legend("topleft",  bty="n", cex=1.5,
       legend=c(expression(paste("MPH: ", italic(r), ""^2, " = 0.75***")), 
                expression(paste("BPH: ", italic(r), ""^2, " = 0.66***"))),
       text.col = c("lightseagreen","navy"), inset=c(-0.12,-0.05))




# C
tabtot2 <- na.omit(tabtot[,c("idGenotype","Type","VegetativeDryMass","FruitNumber")])
tabtot2 <- droplevels(tabtot2)
str(tabtot2)

#pdf("Fig6C.pdf", width = 6, height = 5.5)
#par(mar=c(4,4,.5,.5), tcl=0.3)
plot(tabtot2$VegetativeDryMass, tabtot2$FruitNumber, 
     xlab="Vegetative dry mass", ylab="Fruit number", type="n", ylim=c(-10,360))
points(tabtot2[tabtot2$Type=="Hybrid", "VegetativeDryMass"], 
       tabtot2[tabtot2$Type=="Hybrid", "FruitNumber"], 
       pch=16, col="grey80", cex=1)
tabt1 <- tabtot2[tabtot2$Type=="Hybrid",]
tabt1 <- droplevels(tabt1)
nls1 <- nls(FruitNumber ~ VegetativeDryMass/(A3+ B3*VegetativeDryMass+ C3*(VegetativeDryMass^2)), data=tabt1, 
            start=list(A3=2,B3=1.3,C3=0.01),
            nls.control(maxiter=1000, minFactor = 1/100000000000))
summary(nls1)
A3nlsH <- coef(nls1)[1]
B3nlsH <- coef(nls1)[2]
C3nlsH <- coef(nls1)[3]
curve.S <- function(dm, A3, B3, C3) {dm/(A3+ B3*dm+ C3*(dm^2))}
curve(curve.S(dm = x, A3 = A3nlsH, B3 = B3nlsH, C3=C3nlsH), 20, 
      max(tabt1[,"VegetativeDryMass"]), lwd = 1, add=T, col="grey50", lty=1)

points(tabtot[tabtot$idExp=="GT01" & tabtot$idPot=="77", "VegetativeDryMass"], #611
       tabtot[tabtot$idExp=="GT01" & tabtot$idPot=="77", "FruitNumber"], 
       pch=16, col="firebrick3", cex=1.5)
points(tabtot[tabtot$idExp=="GT01" & tabtot$idPot=="704", "VegetativeDryMass"], #406
       tabtot[tabtot$idExp=="GT01" & tabtot$idPot=="704", "FruitNumber"], 
       pch=16, col="firebrick3", cex=1.5)
points(tabtot[tabtot$idExp=="GT05" & tabtot$idPot=="852", "VegetativeDryMass"], #1167
       tabtot[tabtot$idExp=="GT05" & tabtot$idPot=="852", "FruitNumber"], 
       pch=16, col="dodgerblue3", cex=1.5)


legend("topright", legend=c("Phenotypic non-linearity = 0.20", "MPH = 0.35"), 
       bty="n", cex=1.3, inset=c(0.0,0))
dev.off()




# D
plot(tabhyb$predMPH_Fruit_fromHybVegetativeDryMass, tabhyb$MPH_Fruit, 
     type="n", las=2,
     xlab="Non-linear deviation", ylab="Heterosis on Fruit number",
     xlim=c(-0.76,0.82), ylim=c(-0.94,2.32))
axis(1, at=c(-100,100), labels=NULL, lwd.ticks = 0, pos=0, lwd=1, lty=2, col="firebrick3")
abline(0,1, lty=2)
points(tabhyb$predMPH_Fruit_fromHybVegetativeDryMass, tabhyb$MPH_Fruit, 
       pch=16, col=colMPH, cex=1.3)
points(tabhyb$predMPH_Fruit_fromHybVegetativeDryMass, tabhyb$BPH_Fruit, 
       pch=16, col=colBPH, cex=1.3)
mod2 <- lm(tabhyb$MPH_Fruit ~ tabhyb$predMPH_Fruit_fromHybVegetativeDryMass)
mod4 <- lm(tabhyb$BPH_Fruit ~ tabhyb$predMPH_Fruit_fromHybVegetativeDryMass)
summary(mod2)
summary(mod4)

legend("topleft",  bty="n", cex=1.5,
       legend=c(expression(paste("MPH: ", italic(r), ""^2, " = 0.14***")), 
                expression(paste("BPH: ", italic(r), ""^2, " = 0.10***"))),
       text.col = c("lightseagreen","navy"), inset=c(-0.12,-0.05))

dev.off()

















#-----------------------------------------------------------------
#-----------------------------------------------------------------
# FIG. S1
#-----------------------------------------------------------------
#-----------------------------------------------------------------

dat <- tab
dat$M <- as.numeric(as.character(dat$VegetativeDryMass))
dat$GrowthRate <- as.numeric(as.character(dat$GrowthRate))
dat$GrowthRate <- log10(dat$GrowthRate)
dat$M <- log10(dat$M)

acc <- dat[dat$idExp=="Exp1" & dat$Type=="Accession",] 
acc <- aggregate(acc[,c(8:11)], by=list(idGenotype=acc$idGenotype), FUN=mean, na.rm=T)
acc <- droplevels(acc)

hyb <- dat[dat$idExp=="Exp2" & dat$Type=="Hybrid",]
hyb <- aggregate(hyb[,c(8:11)], by=list(idGenotype=hyb$idGenotype), FUN=mean, na.rm=T)
hyb <- droplevels(hyb)



mod0 <- slope.com(y=GrowthRate, x=M, groups = Type, data=dat)
mod0 <- sma(GrowthRate ~ M*Type, data=dat)

mod1 <- sma(GrowthRate ~ M, data=acc)
mod11 <- slope.test(y=acc$GrowthRate, x=acc$M, test.value = 0.75)
mod1

mod2 <- sma(GrowthRate ~ M, data=hyb)
mod22 <- slope.test(y=hyb$GrowthRate, x=hyb$M, test.value = 0.75)
mod2


pdf("Plot_AllomRelatiomships_log10.pdf", width = 6, height = 6)
#jpeg("GenDist_vs_GeoDist.jpg", width = 4000, height = 3400, quality = 100, res=700)
par(mar=c(4,4,.5,.5), tcl=0.3)
plot(dat[dat$Type=="Accession","M"], dat[dat$Type=="Accession","GrowthRate"], type="n", xlim=c(0.76,3.52),ylim=c(-0.5,1.7), axes=F,
     xlab="log[Vegetative dry mass M (g)]", 
     ylab="log[Growth rate (mg d-1)]")
box()
axis(1, at=seq(-10,20,1), labels=F, tcl=0.5, mgp=c(3,0.5,0))
axis(1, at=seq(0,12,1), labels=10^(seq(0,12,1)), tcl=0.5, mgp=c(3,0.5,0))
axis(1, at=log10(c(seq(0,1,0.1),seq(1,10,1),seq(10,100,10),seq(100,1000,100),
                   seq(1000,10000,1000), seq(10000,100000,10000), seq(100000,1000000,100000),
                   seq(1000000,10000000,1000000),seq(10000000,100000000,10000000),seq(100000000,1000000000,100000000),
                   seq(1000000000,10000000000,1000000000),seq(10000000000,100000000000,10000000000),
                   seq(100000000000,1000000000000,100000000000),seq(1000000000000,10000000000000,1000000000000))),
     labels=F, tcl=0.3, mgp=c(3,0.5,0))
axis(2, at=seq(-10,20,1), labels=F, tcl=0.5, mgp=c(3,0.5,0))
axis(2, at=seq(-2,12,1), labels=10^(seq(-2,12,1)), tcl=0.5, mgp=c(3,0.5,0), las=2)
axis(2, at=log10(c(seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,1,0.1),seq(1,10,1),seq(10,100,10),seq(100,1000,100),
                   seq(1000,10000,1000), seq(10000,100000,10000), seq(100000,1000000,100000),
                   seq(1000000,10000000,1000000),seq(10000000,100000000,10000000),seq(100000000,1000000000,100000000),
                   seq(1000000000,10000000000,1000000000),seq(10000000000,100000000000,10000000000),
                   seq(100000000000,1000000000000,100000000000),seq(1000000000000,10000000000000,1000000000000))),
     labels=F, tcl=0.3, mgp=c(3,0.5,0))


points(dat[dat$Type=="Accession","M"], dat[dat$Type=="Accession","GrowthRate"], 
       pch=16, col=colinb, cex=1.8)
points(dat[dat$Type=="Hybrid","M"], dat[dat$Type=="Hybrid","GrowthRate"], 
       pch=16, col=colhyb, cex=1.8)

abline(mod1, col="firebrick3", lwd=3)
abline(mod2, col="dodgerblue4", lwd=3)

legend("topleft", legend=c("Slope in accessions = 0.74","Slope in hybrids = 0.78"), pch=16,
       col =c(colinb, colhyb), pt.cex=1.8,
       lty=c(1,1), lwd=c(3,3), bty="n", cex=1.1, inset=c(0.05,0))

dev.off()







#-----------------------------------------------------------------
#-----------------------------------------------------------------
# FIG. S2
#-----------------------------------------------------------------
#-----------------------------------------------------------------

pdf("FigS2.pdf", width=7, height = 3)
par(mfrow=c(1,2), oma=c(2,0,0,0), mar=c(.5,4.2,.5,.5))

# A
plot(tabhyb$predHyb_GR_fromHybVegetativeDryMass, tabhyb$GR, 
     type="n", las=2,
     xlab="Predicted hybrid growth rate", ylab="Observed hybrid growth rate",
     xlim=c(-1,27), ylim=c(-1,40))
abline(0,1, lty=2)
points(tabhyb$GR_meanPar, tabhyb$GR, 
       pch=16, col=colmeanpar, cex=1.3)
points(tabhyb$predHyb_GR_fromHybVegetativeDryMass, tabhyb$GR, 
       pch=16, col=colhybval, cex=1.3)
mod2 <- lm(tabhyb$GR ~ tabhyb$GR_meanPar)
mod4 <- lm(tabhyb$GR ~ tabhyb$predHyb_GR_fromHybVegetativeDryMass)
summary(mod2)
summary(mod4)

legend("topleft",  bty="n", cex=0.85,
       legend=c(expression(paste("Genetic additivity: ", italic(r), ""^2, " = 0.78***")), 
                expression(paste("Phenotypic non-linearity: ", italic(r), ""^2, " = 0.95***"))),
       text.col = c("gold2","chocolate4"), inset=c(-0.06,-0.0))


# B
plot(tabhyb$predHyb_Fruit_fromHybVegetativeDryMass, tabhyb$Fruit, 
     type="n", las=2,
     xlab="Predicted hybrid fruit number", ylab="Observed hybrid number",
     xlim=c(20,280), ylim=c(17,380))
abline(0,1, lty=2)
points(tabhyb$Fruit_meanPar, tabhyb$Fruit, 
       pch=16, col=colmeanpar, cex=1.3)
points(tabhyb$predHyb_Fruit_fromHybVegetativeDryMass, tabhyb$Fruit, 
       pch=16, col=colhybval, cex=1.3)
mod2 <- lm(tabhyb$Fruit ~ tabhyb$Fruit_meanPar)
mod4 <- lm(tabhyb$Fruit ~ tabhyb$predHyb_Fruit_fromHybVegetativeDryMass)
summary(mod2)
summary(mod4)

legend("topleft",  bty="n", cex=0.85,
       legend=c(expression(paste("Genetic additivity: ", italic(r), ""^2, " = 0.45***")), 
                expression(paste("Phenotypic non-linearity: ", italic(r), ""^2, " = 0.08***"))),
       text.col = c("gold2","chocolate4"), inset=c(-0.06,-0.0))

dev.off()








#-----------------------------------------------------------------
#-----------------------------------------------------------------
# FIG. S3
#-----------------------------------------------------------------
#-----------------------------------------------------------------

pdf("FigS3.pdf", width=6, height = 5)
par(mfrow=c(2,2), oma=c(0.5,0.5,0,0), mar=c(3,3,.5,.5))


plot(tabhyb$GR, tabhyb$MPH_GR, pch=16, cex=0.8, col="grey20", xlab="Growth rate (mg d-1) in hybrids",
     ylab="MPH of fruit number")
(cor.test(tabhyb$GR, tabhyb$MPH_GR)$estimate)^2

plot(tabhyb$GR, tabhyb$BPH_GR, pch=16, cex=0.8, col="grey20", xlab="Growth rate (mg d-1) in hybrids",
     ylab="BPH of fruit number")
(cor.test(tabhyb$GR, tabhyb$BPH_GR)$estimate)^2

plot(tabhyb$Fruit, tabhyb$MPH_Fruit, pch=16, cex=0.8, col="grey20", xlab="Fruit number in hybrids",
     ylab="MPH of fruit number")
(cor.test(tabhyb$Fruit, tabhyb$MPH_Fruit)$estimate)^2

plot(tabhyb$Fruit, tabhyb$BPH_Fruit, pch=16, cex=0.8, col="grey20", xlab="Fruit number in hybrids",
     ylab="BPH of fruit number")
(cor.test(tabhyb$Fruit, tabhyb$BPH_Fruit)$estimate)^2

dev.off()


