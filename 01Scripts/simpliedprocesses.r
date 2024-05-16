#######################################
## DOWNLOAD AND CLEAN DATA FROM GBIF ##
#######################################
library(rgbif)
library(maps)
library(ggplot2)
library(terra)
library(geodata)

setwd(file.choose())
getwd()

invspp <- data.table::fread('00RawData/taxon.csv')
invspp <- invspp[invspp$taxonomicStatus %in% 'ACCEPTED',]
spp <- invspp$scientificName


results <- read.csv('00RawData/RawData.csv', header = T)
(table(results$taxonRank))
results2 <- results[results$taxonRank %in% 'SPECIES',]
length(table(results$species))
(table(results$species))

# Precompute the count of each species
species_counts <- as.data.frame(table(results$species))
mapp_qc_t <- read.csv('02ProcessedData/qb_invasive_spp.csv', header = T)
mapp_qc_t_70 <- mapp_qc_t[mapp_qc_t$year > 1970,]

length(unique(mapp_qc_t_70$species))
length(unique(mapp_qc_t$species))

(table(mapp_qc_t_70$institutionCode))

# initial plots

g <- ggplot(mapp_qc_t, aes(species))
g + geom_bar()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#names
names(mapp_qc_t)
#### area per spp per year
yr <- unique(mapp_qc_t_70$year)
yr<- yr[-46]

yr<- yr[order(yr)]
sppq <- unique(mapp_qc_t_70$species)

## new dataframe
spp_area <- as.data.frame(matrix(nrow = length(sppq), ncol =  1+length(yr)))
names(spp_area) <- c( yr, 'species')
spp_area$species <- sppq

spp_records <- spp_area

## count number of records per year, per species
for(i in 1:length(sppq)) {
  i.t <- mapp_qc_t[mapp_qc_t$species %in% sppq[i],]
  for(j in 1:length(yr)) {
    i.t.y <- i.t[i.t$year %in% yr[j],]
    if (nrow(i.t.y) > 2) {
      spp_records[i,j] <- nrow(i.t.y)
    } else if (nrow(i.t.y) == 2 | nrow(i.t.y) == 1) {
      spp_records[i,j] <- nrow(i.t.y)
    } else if (nrow(i.t.y) == 0) {
      spp_records[i,j] <- 0
    }
  }
}

m2_to_ha <- function(area_m2) {
  area_ha <- area_m2 / 10000
  return(area_ha)
}
## calculate area per specie per record

for(i in 1:length(sppq)) {
  pb <- txtProgressBar(min = 0, max = 100, style = 3)
  i.t <- mapp_qc_t[mapp_qc_t$species %in% sppq[i],]
  for(j in 1:length(yr)) {
    setTxtProgressBar(pb, i)
    i.t.y <- i.t[i.t$year %in% yr[j],]
    # species with more of 3 records
    if (nrow(i.t.y) > 2) {
      i.t.y.p <- vect(i.t.y, geom =c('x', 'y'))
      i.t.y.ch<- convHull(i.t.y.p)
      spp_area[i,j] <-m2_to_ha(expanse(i.t.y.ch))
      # species with 1 or 2 records
    } else if (nrow(i.t.y) == 2 | nrow(i.t.y) == 1) {
      spp_area[i,j] <- nrow(i.t.y)
    # species without records
    } else if (nrow(i.t.y) == 0) {
      spp_area[i,j] <- NA
    }
    setTxtProgressBar(pb, i)
  }
}

#write.csv(spp_area, '02ProcessedData/_invasive_spp.csv', row.names = F)
spp_area <- read.csv('/Users/ccruzr/Library/Mobile Documents/com~apple~CloudDocs/Cristian/Documents/Estudios/Postgrado/PhD/Courses/Modeling and Indicators Bios2/Alien_Indicator/02ProcessedData/area_invasive_spp.csv', header = T)
spp_area2 <- (log1p(spp_area[,c(1:54)]))
spp_area2$species <- spp_area$species
names(spp_area)


