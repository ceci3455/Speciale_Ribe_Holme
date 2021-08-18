# Cecilie Kristensen 201506140 
# Thesis

# NNMDS for raunkiaer data - both the data with and without the seedmixtures 

#--------------------------- Loading packages and data (renaming and combinding) ------------

# loading packages ####
library(dplyr) #structur data
library(vegan) # NMDS
library(gtools) # converting NA to 0
library(ggplot2)

# loading raunkiaer data and renaming ####
setwd('/Users/ceciliekristensen/Documents/UNI/Kandidat/Speciale_hjem /data/speciale_2020/raunkiaer')
raunkiaer <- read.csv('p_raunkiaer_all.csv', header = TRUE)

setwd('/Users/ceciliekristensen/Documents/UNI/Kandidat/Speciale_hjem /data/data_2020_12_04/')
tot_artslist<- read.csv('total_artsliste.csv',header = TRUE, sep = ';')

raunkiaer <- raunkiaer[-1]

name <- c("38_1", "34_01", "34_11", "38_2", "38_3", "38_4", "38_5",
          "34_12", "34_02", "35_1", "35_2", "35_3", "40_1", "40_2", "40_3", "41_1", "41_2", "41_3", 
          "41_4", "42_1", "42_2", "42_3", "43_1", "43_2", "43_3", "43_4", "43_5", "43_6", "43_7")

names(raunkiaer) <- name

rowname <- as.vector(tot_artslist$arter)
rownames(raunkiaer) <- rowname

pre_markareal <- c("38", "34_1", "34_0", "38", "38", "38", "38", 
  "34_1", "34_0", "35", "35", "35", "40", "40", "40", "41", "41", "41",
  "41", "42", "42", "42", "43", "43", "43", "43", "43", "43", "43")

Jord <- c("ager", "ager", "ager", "ager", "ager", "ager", "ager","ager","ager","ager","ager","ager","MVJ",
          "MVJ", "MVJ","ager","ager","ager","ager","ager","ager","ager", "MVJ", "MVJ","MVJ","MVJ","MVJ","MVJ",
          "MVJ")

# downloading environmental data ####

setwd('/Users/ceciliekristensen/Documents/UNI/Kandidat/Speciale_hjem /data/')
miljo_data <- read.csv('scale_miljo.csv')
# not including the crops data - which is loaded seperately 

# removing the data with no raunkiaer data
miljo_data <- miljo_data[5:33,]

# loading the crops data ####
setwd('/Users/ceciliekristensen/Documents/UNI/Kandidat/Speciale_hjem /data/miljo_data')
crops <- read.csv('NMDS_afgrode.csv', header = TRUE, sep = ';')

crops <- crops %>%
  dplyr::select(afgrode)

names(crops) [1] <- 'grass_score'

afgrode <- crops[5:33,]

# combinding the data and renaming ####
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

# creating the data - without the herb mixtures ####

raunkiaer_nonsown <- (raunkiaer[-c(6,12,23,31,40,43,56,64,69,85), ])

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

#--------------------------- NMDS with herb-mixtures -----------
# Dissimilarity coeffincient and plot ####
#This code is remowing the species which have not been observed in any of the plots 
#Therefor only the presence (p) of species


set.seed(1009)

p_raunkiaer <- raunkiaer[apply(raunkiaer,1 , function(x) !all(x==0)),]

tp_raunkiaer <- t(p_raunkiaer)

# measuring the distance "bray" between the plots                           
bray <- vegdist(tp_raunkiaer, method = "bray", binary = TRUE, wascores=TRUE)

# this code is not necessary

set.seed(1009)
NMDSraunkiaer <- vegan::metaMDS(tp_raunkiaer, distance = "bray", # bray curtis 
                              k = 3, # Number of dimensions
                              trymax = 200) # Number of tries

# plotting the NMDS with the name for the plots
plot(NMDSraunkiaer)

vegan::ordiplot(NMDSraunkiaer, type = "n") # drawing ordination space
orditorp(NMDSraunkiaer, display = "species", col = "blue", air= 0.01) #plotting species inside ordination space
orditorp(NMDSraunkiaer, display = "sites", col = "red", air = 0.01) #plotting the 5-m circles

# - see meaning in zuur p. 264 

# testing the stress level and dimensions ####
stressplot(NMDSraunkiaer)
stress.species <- vegan::goodness(NMDSraunkiaer) # 
mean(stress.species)
goodness(NMDSraunkiaer)

NMDSraunkiaer$stress

# quite good (interval 0.1 - 0.015 - which is the next best categori)

# environmental data ##

miljo_NMDS <- envfit(NMDSraunkiaer, miljo, permu = 999,strata = NULL, choices = c(1,2,3), na.rm = TRUE)
plot(miljo_NMDS)

# choices is the axes to be plotted - here 3
# testing for the number of dimensions ##

# https://stackoverflow.com/questions/49223740/cumulative-variance-explained-for-nmds-in-r/49224573

n = 10
stress <- vector(length = n)
for (i in 1:n) {
  stress[i] <- metaMDS(tp_raunkiaer, distance = "bray", k = i)$stress
}

names(stress) <- paste0(1:n, "Dim")
# x11(width = 10/2.54, height = 7/2.54)
par(mar = c(3.5,3.5,1,1), mgp = c(2, 0.6, 0), cex = 0.8, las = 2)
barplot(stress, ylab = "stress")

#--------------------------- NMDS without the herb-mixtures -------------

#This code is remowing the species which have not been observed in any of the plots
p_nonsown <- raunkiaer_nonsown[apply(raunkiaer_nonsown, 1, function(x) !all(x==0)),]

# dissimilarity matrix and plot (without herbs) ####
#turning the matrix - ready for the analysis

tp_nonsown <- t(p_nonsown)

# measuring the distance "jaccard" between the plots                           
jaccard <- vegdist(tp_nonsown, method = "bray", binary = TRUE, wascores=TRUE)

set.seed(1009)
NMDSraunkiaer_non <- vegan::metaMDS(tp_nonsown, distance = "bray", # Using bray dissimilarity
                              k = 3, # Number of dimensions
                              trymax = 200) # Number of tries


# plotting the NMDS with the name for the plots (without herbs )
plot(NMDSraunkiaer_non)
vegan::ordiplot(NMDSraunkiaer_non, type = "n") # ordination space
orditorp(NMDSraunkiaer_non, display = "species", col = "blue", air= 0.01) #s pecies inside ordination space
orditorp(NMDSraunkiaer_non, display = "sites", col = "red", air = 0.01) # sites

# - see meaning in zuur p. 264 

stressplot(NMDSraunkiaer_non)
stress.species <- vegan::goodness(NMDSraunkiaer_non) # 
mean(stress.species)
goodness(NMDSraunkiaer_non)

NMDSraunkiaer_non$stress

# environmental data (NDMS - without herbs) ##

miljo_NMDS_non <- envfit(NMDSraunkiaer_non, miljo, permu = 999,choices = c(1,2,3), strata = NULL, na.rm = TRUE)

plot(miljo_NMDS_non)

#testing the stress level and the number of dimensions (NMDS - without herbs) ####

# https://stackoverflow.com/questions/49223740/cumulative-variance-explained-for-nmds-in-r/49224573

n = 10

stress <- vector(length = n)
for (i in 1:n) {
  stress[i] <- metaMDS(tp_nonsown, distance = "bray", k = i)$stress
}

names(stress) <- paste0(1:n, "Dim")
# x11(width = 10/2.54, height = 7/2.54)
par(mar = c(3.5,3.5,1,1), mgp = c(2, 0.6, 0), cex = 0.8, las = 2)
barplot(stress, ylab = "stress")


#--------------------------- plotting with herbs (ggplot) ---------------------

# first the NMDS with herbs
# Extracting info from nmds to be used in ggplot
species.scores = as.data.frame(scores(NMDSraunkiaer)) # Creating data frame containing coordinates from nmds
species.scores$pre_markareal = miljo$pre_markareal # Adding field from data
View(species.scores)

# Extracting data for environment vectors
miljo_coord_cont = as.data.frame(scores(miljo_NMDS, "vectors")) * ordiArrowMul(miljo_NMDS) # extracting continuous variables as vectors
View(miljo_coord_cont)

# Adding environmental variables ####

set.seed(1234)
miljo_raun = envfit(NMDSraunkiaer, miljo, permutations = 999, choices = c(1,2,3), strata = NULL, na.rm = TRUE)
miljo_raun

plot(NMDSraunkiaer)
plot(miljo_raun)

# Plotting raunkiaer ordination with field shapes using ggplot ##

# NMDS1 + 2 ####

raunkiaer_gg1 = ggplot(data = species.scores, # using data from species raunkiaer ordination
                    aes(x = NMDS1, y = NMDS2)) + # defining axis
  geom_point(data = species.scores, # Using scatterplot geom
             aes(shape = jord, colour = pre_markareal,), size = 5) +
  scale_colour_manual(values = c("darkorange", "firebrick", 
                                 "red4", "red", "darkolivegreen1",
                                 "goldenrod1", "tomato2",
                                 "darkolivegreen4")) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "black"), # adjusting 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "black"), 
        legend.text = element_text(size = 9, colour = "black")) +
  labs(colour = "Markareal")

raunkiaer_gg1 # calling plot

# NMDS1 + 3 ####

raunkiaer_gg2 = ggplot(data = species.scores, # using data from species presence ordination
                       aes(x = NMDS1, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores, # Using scatterplot geom
             aes(colour = pre_markareal), size = 5) +
  scale_colour_manual(values = c("darkorange", "firebrick", 
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

raunkiaer_gg2

# NMDS2 + 3 ####

raunkiaer_gg3 = ggplot(data = species.scores, # using data from species presence ordination
                       aes(x = NMDS2, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores, # Using scatterplot geom
             aes(colour = pre_markareal), size = 5) +
  scale_colour_manual(values = c("darkorange", "firebrick", 
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

raunkiaer_gg3

#---------------------------- plotting without herbs-------------------

# first the NMDS with herbs
# Extracting info from nmds to be used in ggplot
species.scores_non = as.data.frame(scores(NMDSraunkiaer_non)) # Creating data frame containing coordinates from nmds
species.scores_non$pre_markareal = miljo$pre_markareal # Adding field from data
View(species.scores_non)

# Extracting data for environment vectors
miljo_coord_cont = as.data.frame(scores(miljo_NMDS_non, "vectors")) * ordiArrowMul(miljo_NMDS_non) # extracting continuous variables as vectors
View(miljo_coord_cont)

# Adding environmental variables  ####

set.seed(1009)
miljo_raun_non = envfit(NMDSraunkiaer_non, miljo, permutations = 999, choices = c(1,2,3), strata = NULL, na.rm = TRUE)
miljo_raun_non

plot(NMDSraunkiaer_non)
plot(miljo_raun_non)

# Plotting raunkiaer ordination with field shapes using ggplot ##

# NMDS1 + 2 ####

raunkiaer_non_gg1 = ggplot(data = species.scores_non, # using data from species presence ordination
                       aes(x = NMDS1, y = NMDS2)) + # defining axis in esthetics
  geom_point(data = species.scores_non, # Using scatterplot geom
             aes(colour = pre_markareal), size = 5) +
  scale_colour_manual(values = c("darkorange", "firebrick", 
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

raunkiaer_non_gg1 # calling plot

# NMDS1 + 3 ####

raunkiaer_non_gg2 = ggplot(data = species.scores_non, # using data from species presence ordination
                       aes(x = NMDS1, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores_non, # Using scatterplot geom
             aes(colour = pre_markareal), size = 5) +
  scale_colour_manual(values = c("darkorange", "firebrick", 
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

raunkiaer_non_gg2

# NMDS2 + 3 ####

raunkiaer_non_gg3 = ggplot(data = species.scores_non, # using data from species presence ordination
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

raunkiaer_non_gg3 


#--------------------------plotting with environmental data (herbs)-----------

# NMDS1 + 2 ####

raun_en_gg1 = ggplot(data = species.scores, # using data from species presence ordination
                    aes(x = NMDS1, y = NMDS2)) + # defining axis in esthetics
  geom_point(data = species.scores, # Using scatterplot geom
             aes(shape = Jord, colour = pre_markareal), size = 6) +
  scale_colour_manual(values = c("mediumorchid4", "mediumpurple1", 
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

raun_en_gg1


# NMDS1 + 3 ####

raun_en_gg2 = ggplot(data = species.scores, # using data from species presence ordination
                    aes(x = NMDS1, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores, # Using scatterplot geom
             aes(shape = Jord, colour = pre_markareal), size = 6) +
  scale_colour_manual(values = c("mediumorchid4", "mediumpurple1", 
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

raun_en_gg2

# NMDS2 + 3 ####

raun_en_gg3 = ggplot(data = species.scores, # using data from species presence ordination
                     aes(x = NMDS2, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores, # Using scatterplot geom
             aes(shape = Jord, colour = pre_markareal), size = 6) +
  scale_colour_manual(values = c("mediumorchid4", "mediumpurple1", 
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

raun_en_gg3

#------------------------ plotting with environmental data (without herbs) ----

# NMDS1 + 2 ####
raun_non_en_gg1 = ggplot(data = species.scores_non, # using data from species presence ordination
                     aes(x = NMDS1, y = NMDS2)) + # defining axis in esthetics
  geom_point(data = species.scores_non, # Using scatterplot geom
             aes(shape = Jord, colour = pre_markareal), size = 6) +
  scale_colour_manual(values = c("mediumorchid4", "mediumpurple1", 
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

raun_non_en_gg1

# NMDS1 + 3 ####

raun_non_en_gg2 = ggplot(data = species.scores_non, # using data from species presence ordination
                         aes(x = NMDS1, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores_non, # Using scatterplot geom
             aes(shape = Jord, colour = pre_markareal), size = 6) +
  scale_colour_manual(values = c("mediumorchid4", "mediumpurple1", 
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

raun_non_en_gg2

# NMDS2 + 3 ####

raun_non_en_gg3 = ggplot(data = species.scores_non, # using data from species presence ordination
                         aes(x = NMDS2, y = NMDS3)) + # defining axis in esthetics
  geom_point(data = species.scores_non, # Using scatterplot geom
             aes(shape = Jord, colour = pre_markareal), size = 6) +
  scale_colour_manual(values = c("mediumorchid4", "mediumpurple1", 
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

raun_non_en_gg3


#----------------------------- end --------------------------


