# Cecilie Kristensen 201506140 
# Thesis

# NNMDS for presence/absence data - both the data with and without the seedmixtures

#-------------------- loading the packages and the data (renaming and structuring) ------------------------
# packages ####
#install.packages("spdep", dependencies=TRUE)
#install.packages("ncf", dependencies=TRUE)

# loading packages ####
library(dplyr) #structur data
library(vegan) # NMDS
library(gtools) # converting NA to 0
library(ggplot2) # plotting


# loading data for each plot seperately (not needed) ####

# before combining data if nessecary 
#presence_33 <- read.csv('presence_33.csv',header = TRUE)
#presence_34 <- read.csv('presence_34.csv',header = TRUE)
#presence_35 <- read.csv('presence_35.csv',header = TRUE)
#presence_40 <- read.csv('presence_40.csv',header = TRUE)
#presence_41 <- read.csv('presence_41.csv',header = TRUE)
#presence_42 <- read.csv('presence_42.csv',header = TRUE)
#presence_43 <- read.csv('presence_43.csv',header = TRUE)
#presence_all <- cbind.data.frame(presence_33, presence_34, presence_35, presence_40,presence_41, presence_42, presence_43)
#write.csv(presence_all,"presence_all.csv", row.names = FALSE)

# loading the combinded data ####

setwd('/Users/ceciliekristensen/Documents/UNI/Kandidat/Speciale_hjem /data/speciale_2020/presence_csv')
presence_all <- read.csv('presence_all.csv',header = TRUE)

setwd('/Users/ceciliekristensen/Documents/UNI/Kandidat/Speciale_hjem /data/data_2020_12_04/')
tot_artslist<- read.csv('total_artsliste.csv',header = TRUE, sep = ';')


name <- c("33_1", "33_2", "33_3", "33_4", "38_1", "34_01", "34_11", "38_2", "38_3", "38_4", "38_5",
          "34_12", "34_02", "35_1", "35_2", "35_3", "40_1", "40_2", "40_3", "41_1", "41_2", "41_3", 
          "41_4", "42_1", "42_2", "42_3", "43_1", "43_2", "43_3", "43_4", "43_5", "43_6", "43_7")

names(presence_all) <- name

rowname <- as.vector(tot_artslist$arter)
rownames(presence_all) <- rowname

# old fiels division ####

pre_markareal <- c("33", "33", "33", "33","38", "34_1", "34_0", "38", "38", "38", "38", 
                   "34_1", "34_0", "35", "35", "35", "40", "40", "40", "41", "41", "41",
                   "41", "42", "42", "42", "43", "43", "43", "43", "43", "43", "43")

Jordpre <- c("ager", "ager", "ager", "ager","ager", "ager", "ager", "ager", "ager", "ager", "ager","ager","ager","ager","ager","ager","MVJ",
          "MVJ", "MVJ","ager","ager","ager","ager","ager","ager","ager", "MVJ", "MVJ","MVJ","MVJ","MVJ","MVJ",
          "MVJ")

# loading crop data ####

setwd('/Users/ceciliekristensen/Documents/UNI/Kandidat/Speciale_hjem /data/miljo_data')
crops <- read.csv('NMDS_afgrode.csv', header = TRUE, sep = ';')

afgrode <- crops %>%
  dplyr::select(afgrode)

names(crops) [1] <- 'grass_score'

# loading environmental data ####

setwd('/Users/ceciliekristensen/Documents/UNI/Kandidat/Speciale_hjem /data/')
miljo_data <- read.csv('scale_miljo.csv')

miljo_data <- cbind(miljo_data, afgrode, pre_markareal)

miljo <- miljo_data %>% 
  dplyr::select(log_TWI,
                SRI,
                dist_aa, 
                log_dist_kant,
                log_pH,
                C,
                N, 
                P, 
                afgrode)

names(miljo) [1] <- 'TWI' 
names(miljo) [2] <- 'SRI'
names(miljo) [3] <- 'Afstand_til_aa'
names(miljo) [4] <- 'Afstand_til_markkant'
names(miljo) [5] <- 'pH'
names(miljo) [6] <- 'C'
names(miljo) [7] <- 'N'
names(miljo) [8] <- 'P' 
names(miljo) [9] <- 'Afgrode'

# creating the data without the herb-mixture ####

presence_all_nonsown <- (presence_all[-c(6,12,23,31,40,43,56,64,69,85), ])

## removale of the herbmixtures and other sown plants 
# 6 alm byg
# 12 alm kællingetand
# 23 bibernelle
# 31 cikorie 
# 40 eng-svingel 
# 43 esparsette 
# 56 hvidkløver
# 64 kommen 
# 69 lancet-vejbred
# 85 omvendt kløver

# der kunne godt overvejes om der var nogle som skulle fjernes eller tilføjes
# der er her valgt at bibeholdt rødkløver, alm rajgræs, 
# da de blev fundet på tværs af plottene og ikke pga en blanding


# This code is remowing the species which have not been observed in any of the plots 
presence_allplot <- presence_all[apply(presence_all, 1, function(x) !all(x==0)),]

#write.csv(presence_allplot,"presence_allplot.csv")
as.matrix(presence_allplot)

#---------------------- NMDS with herb mixture ------------------------

# dissimilarity matrix and plot of ordination ####

#turning the matrix - ready for the analysis
tpresence_allplot <- t(presence_allplot)

# measuring the distance "jaccard" between the plots                           
jaccard <- vegdist(tpresence_allplot, method = "jaccard", binary = TRUE, wascores=TRUE)

set.seed(1009)
NMDSSpecies <- vegan::metaMDS(tpresence_allplot, distance = "jaccard", # Jaccard dissimilarity
                              k = 3, # Dimensions
                              trymax = 200) # Number of tries
                

# plotting the NMDS with the name for the plots
plot(NMDSSpecies)

vegan::ordiplot(NMDSSpecies, type = "n") # ordination space
orditorp(NMDSSpecies, display = "species", col = "blue", air= 0.01) # species inside ordination space
orditorp(NMDSSpecies, display = "sites", col = "red", air = 0.01) # sites 

# stress level and number of dimensions ####

# - see meaning in zuur p. 264 

stressplot(NMDSSpecies)
stress.species <- vegan::goodness(NMDSSpecies) # 
mean(stress.species)
goodness(NMDSSpecies)

NMDSSpecies$stress

# environmental data #

set.seed(1009)
miljo_NMDS <- envfit(NMDSSpecies, miljo, permu = 999,strata = NULL, choices = c(1,2,3), na.rm = TRUE)
plot(miljo_NMDS)
miljo_NMDS

# choices is the axes to be plotted - here 3

# testing for the number of dimensions ###
# https://stackoverflow.com/questions/49223740/cumulative-variance-explained-for-nmds-in-r/49224573

n = 10
stress <- vector(length = n)
for (i in 1:n) {
  stress[i] <- metaMDS(tpresence_allplot, distance = "jaccard", k = i)$stress
}

names(stress) <- paste0(1:n, "Dim")
# x11(width = 10/2.54, height = 7/2.54)
par(mar = c(3.5,3.5,1,1), mgp = c(2, 0.6, 0), cex = 0.8, las = 2)
barplot(stress, ylab = "stress")


#---------------------------- NMDS without herb mixture ------------------------

# This code is remowing the species which have not been observed in any of the plots
presence_allplot_nonsown <- presence_all_nonsown[apply(presence_all_nonsown, 1, function(x) !all(x==0)),]

write.csv(presence_allplot_nonsown,"presence_allplot_nonsown.csv")

#### dissimilarity matrix and plot of ordination ####

#write.csv(presence_allplot,"presence_allplot.csv")
as.matrix(presence_allplot_nonsown)

# dissimilarity matrix and ordination ####

#turning the matrix - ready for the analysis
tpresence_nonsown <- t(presence_allplot_nonsown)

# measuring the distance "jaccard" between the plots - not needed                           
#jaccard <- vegdist(tpresence_nonsown, method = "jaccard", binary = TRUE, wascores=TRUE)

set.seed(1009)
NMDSSpecies_non <- vegan::metaMDS(tpresence_nonsown, distance = "jaccard", #Jaccard dissimilarity
                              k = 3, # Dimensions
                              trymax = 200) # Number of tries


# plotting the NMDS with the name for the plots
plot(NMDSSpecies_non)
vegan::ordiplot(NMDSSpecies_non, type = "n") # drawing ordination space
orditorp(NMDSSpecies_non, display = "species", col = "blue", air= 0.01) #plotting species inside ordination space
orditorp(NMDSSpecies_non, display = "sites", col = "red", air = 0.01)


# stress level and number of dimensions ####

# - see meaning in zuur p. 264 

stressplot(NMDSSpecies_non)
stress.species_non <- vegan::goodness(NMDSSpecies_non) # 
mean(stress.species_non)
goodness(NMDSSpecies_non)

NMDSSpecies_non$stress

# environmental data ###

set.seed(1009)
miljo_NMDS_non <- envfit(NMDSSpecies_non, miljo, permu = 999,strata = NULL, choices = c(1,2,3), na.rm = TRUE)
plot(miljo_NMDS_non)
miljo_NMDS_non

# testing for the number of dimensions ###
# https://stackoverflow.com/questions/49223740/cumulative-variance-explained-for-nmds-in-r/49224573

n = 10
stress <- vector(length = n)
for (i in 1:n) {
  stress[i] <- metaMDS(tpresence_nonsown, distance = "jaccard", k = i)$stress
}

names(stress) <- paste0(1:n, "Dim")
# x11(width = 10/2.54, height = 7/2.54)
par(mar = c(3.5,3.5,1,1), mgp = c(2, 0.6, 0), cex = 0.8, las = 2)
barplot(stress, ylab = "stress")

# - see meaning in zuur p. 264 

#--------------------------- plotting with herbs (ggplot) ---------------------

# first the NMDS with herbs
# Extracting info from nmds to be used in ggplot
species.scores = as.data.frame(scores(NMDSSpecies)) # Creating data frame containing coordinates from nmds
species.scores$pre_markareal = miljo$pre_markareal # Adding field from data
View(species.scores)

# Extracting data for environment vectors
miljo_coord_cont = as.data.frame(scores(miljo_NMDS, "vectors")) * ordiArrowMul(miljo_NMDS) # extracting continuous variables as vectors
View(miljo_coord_cont)

# Adding environmental variables ####

set.seed(1009)
miljo_pre = envfit(NMDSSpecies, miljo, permutations = 999, choices = c(1,2,3), strata = NULL, na.rm = TRUE)
miljo_pre

plot(NMDSSpecies)
plot(miljo_pre)

# Plotting presence ordination with field shapes using ggplot ##

# NMDS1 + 2 ####

pre_gg1 = ggplot(data = species.scores, # using data from species raunkiaer ordination
                       aes(x = NMDS1, y = NMDS2)) + # defining axis
  geom_point(data = species.scores, # Using scatterplot geom
             aes(colour = pre_markareal, shape = jord), size = 5) +
  scale_colour_manual(values = c("darkgoldenrod","darkorange", "firebrick", 
                                 "red4", "red", "darkolivegreen1",
                                 "goldenrod1", "tomato2",
                                 "darkolivegreen4"))+
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"), # adjusting 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "black"), 
        legend.text = element_text(size = 9, colour = "black")) +
  labs(colour = "Markareal")

pre_gg1 # calling plot

# NMDS1 + 3 ####

pre_gg2 = ggplot(data = species.scores, # using data from species presence ordination
                       aes(x = NMDS1, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores, # Using scatterplot geom
             aes(colour = pre_markareal), size = 5) +
  scale_colour_manual(values = c("darkgoldenrod","darkorange", "firebrick", 
                                 "red4", "red", "darkolivegreen1",
                                 "goldenrod1", "tomato2",
                                 "darkolivegreen4"))+
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"), # adjusting apearence of plot
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "black"), 
        legend.text = element_text(size = 9, colour = "black")) + 
  labs(colour = "Markareal")

pre_gg2

# NMDS2 + 3 ####

pre_gg3 = ggplot(data = species.scores, # using data from species presence ordination
                       aes(x = NMDS2, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores, # Using scatterplot geom
             aes(colour = pre_markareal), size = 5) +
  scale_colour_manual(values = c("darkgoldenrod","darkorange", "firebrick", 
                                 "red4", "red", "darkolivegreen1",
                                 "goldenrod1", "tomato2",
                                 "darkolivegreen4"))+
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"), # adjusting apearence of plot
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "black"), 
        legend.text = element_text(size = 9, colour = "black")) +
  labs(colour = "Markareal")

pre_gg3


#---------------------------- plotting without herbs-------------------

# first the NMDS with herbs
# Extracting info from nmds to be used in ggplot
species.scores_non = as.data.frame(scores(NMDSSpecies_non)) # Creating data frame containing coordinates from nmds
species.scores_non$pre_markareal = miljo$pre_markareal # Adding field from data
View(species.scores_non)

# Extracting data for environment vectors
miljo_coord_cont = as.data.frame(scores(miljo_NMDS_non, "vectors")) * ordiArrowMul(miljo_NMDS_non) # extracting continuous variables as vectors
View(miljo_coord_cont)

# Adding environmental variables  ####

set.seed(1009)
miljo_pre_non = envfit(NMDSSpecies_non, miljo, permutations = 999, choices = c(1,2,3), strata = NULL, na.rm = TRUE)
miljo_pre_non

plot(NMDSSpecies_non)
plot(miljo_pre_non)

# Plotting presence ordination with field shapes using ggplot ##

# NMDS1 + 2 ####

pre_non_gg1 = ggplot(data = species.scores_non, # using data from species presence ordination
                           aes(x = NMDS1, y = NMDS2)) + # defining axis in esthetics
  geom_point(data = species.scores_non, # Using scatterplot geom
             aes(colour = pre_markareal), size = 5) +
  scale_colour_manual(values = c("darkgoldenrod","darkorange", "firebrick", 
                                 "red4", "red", "darkolivegreen1",
                                 "goldenrod1", "tomato2",
                                 "darkolivegreen4"))+
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"), # adjusting apearence of plot
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "black"), 
        legend.text = element_text(size = 9, colour = "black")) +
  labs(colour = "Markareal")

pre_non_gg1 # calling plot

# NMDS1 + 3 ####

pre_non_gg2 = ggplot(data = species.scores_non, # using data from species presence ordination
                           aes(x = NMDS1, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores_non, # Using scatterplot geom
             aes(colour = pre_markareal), size = 5) +
  scale_colour_manual(values = c("darkgoldenrod","darkorange", "firebrick", 
                                 "red4", "red", "darkolivegreen1",
                                 "goldenrod1", "tomato2",
                                 "darkolivegreen4"))+
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"), # adjusting apearence of plot
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "black"), 
        legend.text = element_text(size = 9, colour = "black")) +
  labs(colour = "Markareal")

pre_non_gg2

# NMDS2 + 3 ####

pre_non_gg3 = ggplot(data = species.scores_non, # using data from species presence ordination
                           aes(x = NMDS2, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores_non, # Using scatterplot geom
             aes(shape = Jord, colour = pre_markareal), size = 5) +
  scale_colour_manual(values = c("mediumorchid4", "mediumpurple1", 
                                 "lightskyblue", "cornflowerblue", "gray75",
                                 "royalblue3", "navy",
                                 "gray52"))+
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"), # adjusting apearence of plot
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "black"), 
        legend.text = element_text(size = 9, colour = "black")) +
  labs(colour = "Markareal")

pre_non_gg3 

# -------------------------------- plots with environmental vectors (herbs) -------------------

# NMDS 1 + 2 ####

pre_en_gg1 = ggplot(data = species.scores, # using data from species presence ordination
                        aes(x = NMDS1, y = NMDS2)) + # defining axis in esthetics
  geom_point(data = species.scores, # Using scatterplot geom
             aes(shape = Jordpre, colour = pre_markareal), size = 6) +
  scale_colour_manual(values = c("mediumorchid2","mediumorchid4", "mediumpurple1", 
                                 "lightskyblue", "cornflowerblue", "gray75",
                                 "royalblue3", "navy",
                                 "gray52")) + 
  geom_segment(aes(x = 0, y = 0, 
                   xend = NMDS1, 
                   yend = NMDS2),
               data = miljo_coord_cont, size = 0.5, alpha = 0.5, colour = "grey30") +
  geom_text(data = miljo_coord_cont, 
            aes(x = NMDS1, y = NMDS2), 
            colour = "grey30", 
            fontface = "bold",
            size = 5,
            label = row.names(miljo_coord_cont)) +
  theme(axis.title = element_text(size = 12, face = "bold", colour = "black"), # adjusting apearence of plot
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 12, face = "bold", colour = "black"), 
        legend.text = element_text(size = 13, colour = "black")) +
  labs(colour = "Markareal")

pre_en_gg1


# NMDS 1 + 3 ####

pre_en_gg2 = ggplot(data = species.scores, # using data from species presence ordination
                    aes(x = NMDS1, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores, # Using scatterplot geom
             aes(shape = Jordpre, colour = pre_markareal), size = 6) +
  scale_colour_manual(values = c("mediumorchid2","mediumorchid4", "mediumpurple1", 
                                 "lightskyblue", "cornflowerblue", "gray75",
                                 "royalblue3", "navy",
                                 "gray52")) + 
  geom_segment(aes(x = 0, y = 0, 
                   xend = NMDS1, 
                   yend = NMDS3),
               data = miljo_coord_cont, size = 0.5, alpha = 0.5, colour = "grey30") +
  geom_text(data = miljo_coord_cont, 
            aes(x = NMDS1, y = NMDS3), 
            colour = "grey30", 
            fontface = "bold",
            size = 5,
            label = row.names(miljo_coord_cont)) +
  theme(axis.title = element_text(size = 12, face = "bold", colour = "black"), # adjusting apearence of plot
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 12, face = "bold", colour = "black"), 
        legend.text = element_text(size = 13, colour = "black")) +
  labs(colour = "Markareal")

pre_en_gg2

# NMDS 2 + 3 ####

pre_en_gg3 = ggplot(data = species.scores, # using data from species presence ordination
                    aes(x = NMDS2, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores, # Using scatterplot geom
             aes(shape = Jordpre, colour = pre_markareal), size = 6) +
  scale_colour_manual(values = c("mediumorchid2","mediumorchid4", "mediumpurple1", 
                                 "lightskyblue", "cornflowerblue", "gray75",
                                 "royalblue3", "navy",
                                 "gray52")) + 
  geom_segment(aes(x = 0, y = 0, 
                   xend = NMDS2, 
                   yend = NMDS3),
               data = miljo_coord_cont, size = 0.5, alpha = 0.5, colour = "grey30") +
  geom_text(data = miljo_coord_cont, 
            aes(x = NMDS2, y = NMDS3), 
            colour = "grey30", 
            fontface = "bold",
            size = 5,
            label = row.names(miljo_coord_cont)) +
  theme(axis.title = element_text(size = 12, face = "bold", colour = "black"), # adjusting apearence of plot
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 12, face = "bold", colour = "black"), 
        legend.text = element_text(size = 13, colour = "black")) +
  labs(colour = "Markareal")

pre_en_gg3

# -------------------------------- plots with environmental vectors (without herbs) -------------------

# NMDS 1 + 2 ####

  
pre_en_non_gg1 = ggplot(data = species.scores_non, # using data from species presence ordination
                        aes(x = NMDS1, y = NMDS2)) + # defining axis in esthetics
  geom_point(data = species.scores_non, # Using scatterplot geom
             aes(shape = Jordpre, colour = pre_markareal), size = 6) +
  scale_colour_manual(values = c("mediumorchid2","mediumorchid4", "mediumpurple1", 
                                 "lightskyblue", "cornflowerblue", "gray75",
                                 "royalblue3", "navy",
                                 "gray52")) + 
  geom_segment(aes(x = 0, y = 0, 
                   xend = NMDS1, 
                   yend = NMDS2),
               data = miljo_coord_cont, size = 0.5, alpha = 0.5, colour = "grey30") +
  geom_text(data = miljo_coord_cont, 
            aes(x = NMDS1, y = NMDS2), 
            colour = "grey30", 
            fontface = "bold",
            size = 5,
            label = row.names(miljo_coord_cont)) +
  theme(axis.title = element_text(size = 12, face = "bold", colour = "black"), # adjusting apearence of plot
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 12, face = "bold", colour = "black"), 
        legend.text = element_text(size = 13, colour = "black")) +
  labs(colour = "Markareal")

pre_en_non_gg1

# NMDS 1 + 3 ####
pre_en_non_gg2 = ggplot(data = species.scores_non, # using data from species presence ordination
                     aes(x = NMDS1, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores_non, # Using scatterplot geom
             aes(shape = Jordpre, colour = pre_markareal), size = 6) +
  scale_colour_manual(values = c("mediumorchid2","mediumorchid4", "mediumpurple1", 
                                 "lightskyblue", "cornflowerblue", "gray75",
                                 "royalblue3", "navy",
                                 "gray52")) + 
  geom_segment(aes(x = 0, y = 0, 
                xend = NMDS1, 
                yend = NMDS3),
  data = miljo_coord_cont, size = 0.5, alpha = 0.5, colour = "grey30") +
  geom_text(data = miljo_coord_cont, 
            aes(x = NMDS1, y = NMDS3), 
            colour = "grey30", 
            fontface = "bold",
            size = 5,
            label = row.names(miljo_coord_cont)) +
  theme(axis.title = element_text(size = 12, face = "bold", colour = "black"), # adjusting apearence of plot
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 12, face = "bold", colour = "black"), 
        legend.text = element_text(size = 13, colour = "black")) +
  labs(colour = "Markareal")

pre_en_non_gg2


# NMDS 2 + 3 ####

pre_en_non_gg3 = ggplot(data = species.scores_non, # using data from species presence ordination
                        aes(x = NMDS2, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores_non, # Using scatterplot geom
             aes(shape = Jordpre, colour = pre_markareal), size = 6) +
  scale_colour_manual(values = c("mediumorchid2","mediumorchid4", "mediumpurple1", 
                                 "lightskyblue", "cornflowerblue", "gray75",
                                 "royalblue3", "navy",
                                 "gray52")) + 
  geom_segment(aes(x = 0, y = 0, 
                   xend = NMDS2, 
                   yend = NMDS3),
               data = miljo_coord_cont, size = 0.5, alpha = 0.5, colour = "grey30") +
  geom_text(data = miljo_coord_cont, 
            aes(x = NMDS2, y = NMDS3), 
            colour = "grey30", 
            fontface = "bold",
            size = 5,
            label = row.names(miljo_coord_cont)) +
  theme(axis.title = element_text(size = 12, face = "bold", colour = "black"), # adjusting apearence of plot
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 12, face = "bold", colour = "black"), 
        legend.text = element_text(size = 13, colour = "black")) +
  labs(colour = "Markareal")

pre_en_non_gg3
