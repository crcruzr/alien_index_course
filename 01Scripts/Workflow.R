#######################################
## DOWNLOAD AND CLEAN DATA FROM GBIF ##
#######################################
library(rgbif)
library(maps)
library(ggplot2)
library(terra)
library(geodata)
library(tidyverse)

setwd('/Users/ccruzr/Library/Mobile Documents/com~apple~CloudDocs/Cristian/Documents/Estudios/Postgrado/PhD/Courses/Modeling and Indicators Bios2/alien_index_course')
getwd()


qc <- gadm('CAN', level=1, path=tempdir())
qc<- qc[qc$NAME_1 %in% 'Québec',]
plot(qc)


invspp <- data.table::fread('00RawData/taxon.csv')
invspp <- invspp[invspp$taxonomicStatus %in% 'ACCEPTED',]
spp <- invspp$scientificName

#Contarinia baeri.  "Paraphytomyza populicola" don't have records


# download GBIF occurrence data for this species; this takes time if there are many data points!
gbif_data <- occ_data(scientificName = spp, country = 'CA', stateProvince= "Québec", year = '1970,2024', hasCoordinate = TRUE)

# create and fill a list with only the 'data' section for each species:
myspecies_coords_list <- vector("list", length(spp))
names(gbif_data) <- spp
nm<- c( "key",  "scientificName", "kingdom", "phylum",  "family", "genus", "specificEpithet", "taxonRank",  "year", "continent",  "countryCode", "decimalLongitude", "decimalLatitude",  "occurrenceStatus", "institutionCode")

## Centralize data
results<-as.data.frame(matrix(ncol = length(nm)))
names(results) <- nm
  for (s in 1:length(spp)) {
    coords <- gbif_data[[s]]$data
    if (!is.null(coords)) {
    if (!'year' %in% colnames(coords)) {
      # Subset the data frame while maintaining the 'year' column
      coords$year <- NA
      coords <- coords[ ,nm]
    } else {
      # If 'year' column doesn't exist, subset without it
      coords <- gbif_data[[s]]$data[ , nm]
    }
    results <- rbind(results, coords)
  }
}

results <- results[-1,]

results$species <- paste(results$genus, results$specificEpithet)
(table(results$taxonRank))

results2 <- results[results$taxonRank %in% 'SPECIES',]
length(table(results$species))
(table(results$species))

(table(results$year))

write.csv(results, '00RawData/RawData.csv')
results <- read.csv('00RawData/RawData.csv', header = T)


# Precompute the count of each species
species_counts <- as.data.frame(table(results$species))

mapp <- vect(results, geom = c('decimalLongitude', 'decimalLatitude'), crs="+proj=longlat +datum=WGS84")
mapp_qc <- intersect(qc, mapp)

plot(qc)
plot(mapp_qc, add = T)

coor_geo <- as.data.frame(mapp_qc,  geom="XY" )
#project to planar coordiantes
plot(qc, add = T)
mapp2 <- project(mapp_qc, 'EPSG:32198')
plot(mapp2)

qc_pla <- project(qc, 'EPSG:32198' )
plot(qc_pla)
plot(mapp2, add = T)

mapp_qc_t <- as.data.frame(mapp2 , geom="XY")
mapp_qc_t$decimalLongitude <- coor_geo$x
mapp_qc_t$decimalLatitude <- coor_geo$y

rm(results, coor_geo, gbif_data, myspecies_coords_list)

mapp_qc_t$year <- as.numeric(mapp_qc_t$year)
mapp_qc_t<-mapp_qc_t[!is.na(mapp_qc_t$year),]

write.csv(mapp_qc_t, '02ProcessedData/qb_invasive_spp.csv', row.names = F)
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

ggplot(data = mapp_qc_t, aes(x = decimalLongitude,
                         y = decimalLatitude,
                         colour = species)) +
  geom_polygon(data = map_data("world"),
        aes(x = long, y = lat, group = group),
                    fill = "grey95",
                   color = "gray40",
                   size = 0.1) +
  geom_point(size = 0.7, alpha = 0.5) +
  coord_fixed(xlim = range(mapp_qc_t$decimalLongitude, na.rm = TRUE),
              ylim = range(mapp_qc_t$decimalLatitude, na.rm = TRUE)) +
  theme_bw()

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
spp_area <- data.table::fread('02ProcessedData/area_invasive_spp.csv', header = T)
spp_area2 <- (log1p(spp_area[,c(1:54)]))
spp_area2$species <- spp_area$species
names(spp_area)


##Rasterize 
#1000 = 1 Km@
r <- rast(qc_pla, crs= 'EPSG:32198')
res(r) <- 10000  # in meters
qc_r <- rasterize(qc_pla, r)
qc_r <- crop(qc_r, qc_pla)
plot(qc_r)
r
table(is.na(mapp_qc_t_70$year))
mapp_qc_t_70<- mapp_qc_t_70[!is.na(mapp_qc_t_70$year),]
table(is.na(mapp_qc_t_70$year))

##

spp_f <- unique(mapp_qc_t_70$species)

mapp_qc_t_70$id <- 1:nrow(mapp_qc_t_70)

mapp_qc_t_70<- mapp_qc_t_70[-51462,] ## record outside of the area


win_area_ras <- as.data.frame(matrix(nrow = length(spp_f), ncol = 51))
names(win_area_ras) <- c(paste0('win_',1970:2020))
rownames(win_area_ras) <- spp_f

#w<- win_area_ras
win_area_ras <- w

## to calculate the number of pixels per window

for(i in 1:ncol(win_area_ras)){
    pb <- txtProgressBar(min = 0, max = 100, style = 3)
    setTxtProgressBar(pb, (i*100)/ncol(win_area_ras))

window.i <- mapp_qc_t_70[mapp_qc_t_70$year >= (1969+i) & mapp_qc_t_70$year <= (1973+i),] 
  for (j in 1:length(spp_f)) {
    window.i.j <- window.i[window.i$species %in% spp_f[j],]
      if (!nrow(window.i.j) == 0) {
      sp <- vect(window.i.j, geom = c('x', 'y'), crs= 'EPSG:32198')
      spp.j <- rasterize(sp, qc_r, "species", function(x, ...) length(unique(na.omit(x))))
      num_valid_pixels <- sum(freq(spp.j)$count)
      win_area_ras[spp_f[j],i] <- num_valid_pixels
      } else {
       win_area_ras[spp_f[j],i] <- 0 
      }
    }
}

i
j
copy <- win_area_ras
win_area_ras[win_area_ras == 0] <- NA


### Graphics
win_area_ras$species <- rownames(win_area_ras)
# create data with less columns 
win_area_ras_long <- pivot_longer(win_area_ras, cols = -species, names_to = "Window", values_to = "Value")
row.names(win_area_ras)

win_area_ras_long$Valuelog <- (log1p(win_area_ras_long$Value))

## without log
win_area_ras_long %>%
  group_by(Window) %>%
  summarize(median = quantile(Value, p= 0.5, na.rm = TRUE) , 
            q1 = quantile(Value, p = 0.25, na.rm = TRUE),
            q3 = quantile(Value, p = 0.75, na.rm = TRUE)) %>%
ggplot(aes(x = Window, y = median)) +
  geom_point() +
  geom_errorbar(aes(ymin = q1, ymax = q3)) +
    geom_line(aes(x = Window, y = q1))+
     theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

win_area_ras_long$Window <- gsub('win_', '', win_area_ras_long$Window)

## with log
graph1 <-win_area_ras_long %>%
  group_by(Window) %>%
  summarize(mean = mean(Valuelog,na.rm = TRUE) , 
            q1 = quantile(Valuelog, p = 0.25, na.rm = TRUE),
            q3 = quantile(Valuelog, p = 0.75, na.rm = TRUE)) %>%
ggplot( aes(x = Window, y = mean)) +
  geom_errorbar(aes(ymin = q1, ymax = q3), col = '#000099') +
  geom_point(size = 10, col = '#FF8000') +
    #geom_line(aes(x = Window, y = q1)) +
     theme_bw()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 24),
    axis.text.y = element_text(size = 24),
     axis.title.x = element_text(size = 26),  # Change x-axis title size
  axis.title.y = element_text(size = 26),  # Change y-axis title size
    )
    

png("04Graphics/rec_Xwind_5Y_10km2_log.png", width = 2080, height = 1080)
print(graph1)
dev.off()


## to calculate the area per poligon

win_area_vec <- as.data.frame(matrix(nrow = length(spp_f), ncol = 51))
names(win_area_vec) <- c(paste0('win_',1970:2020))
rownames(win_area_vec) <- spp_f

for(i in 1:ncol(win_area_vec)){
    pb <- txtProgressBar(min = 0, max = 100, style = 3)
    setTxtProgressBar(pb, (i*100)/ncol(win_area_vec))

window.i <- mapp_qc_t_70[mapp_qc_t_70$year >= (1969+i) & mapp_qc_t_70$year <= (1973+i),] 
  for (j in 1:length(spp_f)) {
    window.i.j <- window.i[window.i$species %in% spp_f[j],]
      if (nrow(window.i.j) > 2) {
      i.t.y.p <- vect(window.i.j, geom =c('x', 'y'), crs= 'EPSG:32198' )
      i.t.y.ch<- convHull(i.t.y.p)
      num_valid_poligons <- m2_to_ha(expanse(i.t.y.ch))
      win_area_vec[spp_f[j],i] <- num_valid_poligons
      } else if (nrow(window.i.j) == 2 | nrow(window.i.j) == 1) {
      win_area_vec[spp_f[j],i] <- nrow(window.i.j)
    # species without records
    } else if (nrow(window.i.j) == 0) {
      win_area_vec[spp_f[j],i] <- NA
    }
  }
}

### Graphics
win_area_vec$species <- rownames(win_area_vec)
# create data with less columns 
win_area_vec_long <- pivot_longer(win_area_vec, cols = -species, names_to = "Window", values_to = "Value")

win_area_vec_long$Valuelog <- (log1p(win_area_vec_long$Value))

## without log
win_area_vec_long %>%
  group_by(Window) %>%
  summarize(median = quantile(Value, p= 0.5, na.rm = TRUE) , 
            q1 = quantile(Value, p = 0.25, na.rm = TRUE),
            q3 = quantile(Value, p = 0.75, na.rm = TRUE)) %>%
ggplot(aes(x = Window, y = median)) +
  geom_point() +
  geom_errorbar(aes(ymin = q1, ymax = q3)) +
    geom_line(aes(x = Window, y = q1)) 

## with log
win_area_vec_long %>%
  group_by(Window) %>%
  summarize(median = quantile(Valuelog, p= 0.5, na.rm = TRUE) , 
            q1 = quantile(Valuelog, p = 0.25, na.rm = TRUE),
            q3 = quantile(Valuelog, p = 0.75, na.rm = TRUE)) %>%
ggplot(aes(x = Window, y = median)) +
  geom_point() +
  geom_errorbar(aes(ymin = q1, ymax = q3)) +
    geom_line(aes(x = Window, y = q1)) 


write.csv(win_area_vec_long, '02ProcessedData/Windo_5Yr_vector.csv')
write.csv(win_area_ras_long, '02ProcessedData/Windo_5Yr_raster.csv')
write.csv(win_area_ras, '02ProcessedData/Windo_5Yr_raster_ext.csv')

writeRaster(qc_r, '02ProcessedData/Quebec_planar_raster.tif')
qc_r
plot(qc_r)


### plots 

