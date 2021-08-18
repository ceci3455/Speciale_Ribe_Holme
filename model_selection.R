#------------ Exploration - Chapter 4 Zuur et al.--------------
# loading packages and files ####

library(ggplot2) # plotting
library(dplyr) # sorting and selection
library(nlme) # model - for gls
library(lme4) # model - for glmer
library(tidyr) # sorting
library(MASS) # model selection
library(stargazer) # p-values and model values
library(multcompView) # plotting
library(car) # VIF
library(ggpubr)
library(tidyverse)
library(rstatix)
library(ggplot2)


setwd('/Users/ceciliekristensen/Documents/UNI/Kandidat/Speciale_hjem /data/miljo_data')
plotdata <- read.csv('miljo_data_indeks_sum.csv', header = TRUE, sep =';')


setwd('/Users/ceciliekristensen/Documents/UNI/Kandidat/Speciale_hjem /data/miljo_data')
crops <- read.csv('NMDS_afgrode.csv', header = TRUE, sep = ';')

crops <- crops %>%
  dplyr::select(afgrode)

names(crops) [1] <- 'afgrode'

#loooking at data - converting into factors ####

str(plotdata)

# some of the varibles are defined as numeric and therefor the have to be defined as factors

plotdata$Graesning <- factor(plotdata$Graesning, labels = c("0", "1"))
plotdata$plot_id <- factor(plotdata$plot_id)
plotdata$markareal <- factor(plotdata$markareal)
plotdata$plot <- factor(plotdata$plot)

str(plotdata)

plotdata <- cbind(plotdata, crops)

# boxplots and histograms ####

#predictor
boxplot(plotdata$SRI)
boxplot(plotdata$hojde)
boxplot(plotdata$veg_hojde)
boxplot(plotdata$dist_kant)
boxplot(plotdata$dist_aa)
boxplot(plotdata$TWI)
boxplot(plotdata$pH)
boxplot(plotdata$P)
boxplot(plotdata$N)
boxplot(plotdata$C)
boxplot(plotdata$afgrode)


#predictor
hist(plotdata$SRI, breaks = 30)
hist(plotdata$hojde, breaks = 30)
hist(plotdata$veg_hojde, breaks = 30)
hist(plotdata$SRI, breaks = 30)
hist(plotdata$TWI, breaks = 30)
hist(plotdata$dist_kant, breaks = 30)
hist(plotdata$dist_aa, breaks = 30)
hist(plotdata$pH, breaks = 30)
hist(plotdata$P, breaks = 30)
hist(plotdata$N, breaks = 30)
hist(plotdata$C, breaks = 30)
hist(plotdata$afgrode, breaks = 15)


#response
boxplot(plotdata$NOVANA)
boxplot(plotdata$Naturlig)
boxplot(plotdata$Mark)

hist(plotdata$NOVANA, breaks = 30)
hist(plotdata$Naturlig, breaks = 30)
hist(plotdata$Mark, breaks = 30)

# transforming the data which is skewed ####

hist(log(plotdata$hojde), breaks = 30)
hist(log(plotdata$veg_hojde), breaks = 30)
hist(log(plotdata$TWI), breaks = 30)
hist(log(plotdata$dist_kant), breaks = 30)
hist(log(plotdata$pH), breaks = 30)
hist(log(plotdata$afgrode), breaks = 15) 


hist(log(plotdata$Naturlig), breaks = 30)

# if i don't change  the value it is not possible to make the shapiro-test
log_hojde <- log(plotdata$hojde)
log_hojde[log_hojde == '-Inf'] <- 0

# I'm making a new column for each of the new transformed data

plotdata_trans <- mutate(plotdata, log_hojde = log_hojde, log_TWI = log(plotdata$TWI)
                         , log_veghojde = log(plotdata$veg_hojde), 
                         log_dist_kant = log(plotdata$dist_kant), log_pH = log(plotdata$pH), 
                         log_Naturlig = log(plotdata$Naturlig), log_NOVANA = log(plotdata$NOVANA), 
                         log_afgrode = log(plotdata$afgrode))

#testing for normality with the shapiro-test 
# collinearity ####

# I've already remove the colinerarity between the slop and aspect by using SRI instead

trans_data <- data.frame(plotdata_trans$log_NOVANA, 
                         plotdata_trans$log_Naturlig, 
                         plotdata_trans$Mark, 
                         plotdata_trans$SRI, 
                         plotdata_trans$dist_aa, 
                         plotdata_trans$log_hojde, 
                         plotdata_trans$log_TWI, 
                         plotdata_trans$log_veghojde,
                         plotdata_trans$log_dist_kant,
                         plotdata_trans$log_pH, 
                         plotdata_trans$N, 
                         plotdata_trans$P,
                         plotdata_trans$C, 
                         plotdata_trans$log_afgrode)

# renaming ####

names(trans_data) [1] <- "NOVANA"
names(trans_data) [2] <- "Naturlig"
names(trans_data) [3] <- "Mark"
names(trans_data) [4] <- "SRI"
names(trans_data) [5] <- "dist_aa"
names(trans_data) [6] <- "log_hojde"
names(trans_data) [7] <- "log_TWI"
names(trans_data) [8] <- "log_veghojde"
names(trans_data) [9] <- "log_dist_kant"
names(trans_data) [10] <- "log_pH"
names(trans_data) [11] <- "N"
names(trans_data) [12] <- "P"
names(trans_data) [13] <- "C"
names(trans_data) [14] <- "afgrode"

####

matrix_cor <- cor(trans_data[,4:14], method = "spearman")
View(matrix_cor)

#The different concentrations of nutrients are correlated - 
#I will therefore choose one for the model 

#### testing P, N, C and pH for themselves ###

soil <- trans_data %>%
  dplyr::select(log_pH, 
                N, 
                P, 
                C)
# removing the NA value and the other for that point and testing it without

soil_p <- soil[-18,]

soil_matrix <- cor(soil_p, method = "spearman") 
View(soil_matrix)

# by testing it seperate I can see that there is only a correlation between C and N. 
# I therefor have to choose one 

# scaling, renaming and saving (both scaled and non scaled) ####

modeldata <- data.frame(plotdata_trans$NOVANA, plotdata_trans$log_Naturlig, plotdata_trans$Mark, plotdata_trans$SRI, plotdata_trans$log_hojde, 
                        plotdata_trans$log_TWI, plotdata_trans$log_veghojde, plotdata_trans$dist_aa,
                        plotdata_trans$log_dist_kant, plotdata_trans$log_pH, plotdata_trans$C, plotdata_trans$N, plotdata_trans$P, 
                        plotdata_trans$log_afgrode)

names(modeldata) [1] <- "NOVANA"
names(modeldata) [2] <- "Naturlig"
names(modeldata) [3] <- "Mark"
names(modeldata) [4] <- "SRI"
names(modeldata) [5] <- "log_hojde"
names(modeldata) [6] <- "log_TWI"
names(modeldata) [7] <- "log_veghojde"
names(modeldata) [8] <- "dist_aa"
names(modeldata) [9] <- "log_dist_kant"
names(modeldata) [10] <- "log_pH"
names(modeldata) [11] <- "C"
names(modeldata) [12] <- "N"
names(modeldata) [13] <- "P"
names(modeldata) [14] <- "log_afgrode"


modeldata$P[18] <- 11.235

# har testet hvilke modeller der giver den største forklaringsgrad: 
# modellen med den nye gennemsnitsværdi #R2c 0.6032145
# modellen uden P #R2c 0.600748
# modellen med P med uden datapunkt 18 #R2c 0.589497
# alle tests er udført på NOVANA - mixed modellen 


scale_variables <- scale(modeldata)

scale_data <- cbind(plotdata_trans$plot_id, scale_variables, plotdata_trans$jord, plotdata_trans$Graesning,
                    plotdata_trans$markareal)

model_data1 <- data.frame(scale_data)

names(model_data1) [1] <- "id"
names(model_data1) [16] <- "jord"
names(model_data1) [17] <-"graesning"
names(model_data1) [18] <- "markareal"

str(model_data1)

model_data1$graesning <- factor(model_data1$graesning, labels = c("0","1"))
model_data1$id <- factor(model_data1$id)
model_data1$markareal <- factor(model_data1$markareal, labels = c("33", "34", "35", "40", 
                                                                  "41", "42", "43"))
model_data1$jord <- factor(model_data1$jord, labels = c("ager", "mvj"))

str(model_data1)

# without scaling 
modeldata_nonscale <- cbind(plotdata_trans$plot_id, modeldata, plotdata_trans$jord, plotdata_trans$Graesning,
                            plotdata_trans$markareal)

names(modeldata_nonscale) [1] <- "id"
names(modeldata_nonscale) [16] <- "jord"
names(modeldata_nonscale) [17] <-"graesning"
names(modeldata_nonscale) [18] <- "markareal"


setwd('/Users/ceciliekristensen/Documents/UNI/Kandidat/Speciale_hjem /data/')
write.csv(modeldata_nonscale,"no_scale_miljo_afgrode.csv", row.names = FALSE)
write.csv(model_data1,"scale_miljo_afgrode.csv", row.names = FALSE)

# ---------------- Model building ----------------------------------
# Linear model ####

# Can't use these models - because the assumptions is continues variables

lm_NOVANA <-lm(NOVANA ~ SRI 
              + dist_aa 
              + graesning 
              + log_hojde 
              + log_TWI 
              + log_veghojde 
              + log_dist_kant
              + jord 
              +  markareal
              + log_pH
              + N
              + P
              + log_afgrode, 
              data = model_data1)

anova(lm_NOVANA)
# The variable jord is significant and dist aa

plot(lm_NOVANA)
summary(lm_model_NOVANA)

# The points 14, 27, and 29 have a tendensy 
# YOU HAVE TO CHECK THEM AGAIN

lm_Naturlig <-lm(Naturlig ~ SRI 
                     + dist_aa 
                     + graesning 
                     + log_hojde 
                     + log_TWI 
                     + log_veghojde 
                     + log_dist_kant
                     + jord 
                     +  markareal
                     + log_pH
                     + N
                     + P
                 + log_afgrode, 
                     data = model_data1)

anova(lm_Naturlig)

# The variable jord is significant samt dist aa

plot(lm_Naturlig)
summary(lm_Naturlig)

#The points 14, 27, 29 have a tendensy

lm_Mark <-lm(Mark ~ SRI 
                       + dist_aa 
                       + graesning 
                       + log_hojde 
                       + log_TWI 
                       + log_veghojde 
                       + log_dist_kant
                       + jord 
                       +  markareal
                       + log_pH
                       + N
                       + P
             + log_afgrode, 
                       data = model_data1)

anova(lm_Mark)
# The variable jord is significant

plot(lm_Mark)
summary(lm_Mark)

# the points 14, 27, 31 have a tendensy 

hist(residuals(lm_NOVANA),breaks = 50)
hist(residuals(lm_Naturlig),breaks = 50)
hist(residuals(lm_Mark),breaks = 50)

shapiro.test(residuals(lm_NOVANA))
shapiro.test(residuals(lm_Naturlig))
shapiro.test(residuals(lm_Mark))

# GLS ####

#Generalized least squares fitted linear model 
#The assumption of GLS is that the errors are independent and identically distributed. 
#Furthermore, other assumptions include:
  
# The error variances are homoscedastic
# Errors are uncorrelated
# Normally distributed

GLS_NOVANA <- gls(NOVANA ~ SRI 
                  + dist_aa 
                  + graesning 
                  + log_hojde 
                  + log_TWI 
                  + log_veghojde 
                  + log_dist_kant
                  + jord 
                  + log_pH
                  + N
                  + P
                  + log_afgrode, 
                  data = model_data1)

anova(GLS_NOVANA)
summary(GLS_NOVANA)

stargazer(GLS_NOVANA, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

#----------------Naturlig---------------

GLS_Naturlig <- gls(Naturlig ~ SRI 
                    + dist_aa 
                    + graesning 
                    + log_hojde 
                    + log_TWI 
                    + log_veghojde 
                    + log_dist_kant
                    + jord 
                    + log_pH
                    + N
                    + P
                    + log_afgrode, 
                    data = model_data1,
                    method = "ML",
                    na.action = "na.omit") 

anova(GLS_Naturlig)
summary(GLS_Naturlig)

stargazer(GLS_Naturlig, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

#-----------------Mark-----------
GLS_Mark <- gls(Mark ~ SRI 
                + dist_aa 
                + graesning 
                + log_hojde 
                + log_TWI 
                + log_veghojde 
                + log_dist_kant
                + jord 
                + log_pH
                + N
                + P
                + log_afgrode, 
                data = model_data1,
                method = "ML",
                na.action = "na.omit") 

anova(GLS_Mark)
summary(GLS_Mark)

stargazer(GLS_Mark, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")


# GLM ####

# Generalized linear model 
#(Generalized) Linear models make some strong assumptions concerning the data structure:
# Independance of each data points.
# Correct distribution of the residuals.
# Correct specification of the variance structure.
# Linear relationship between the response and the linear predictor.


#----------------NOVANA-----------------

GLM_NOVANA <- glm(NOVANA ~ SRI 
                    + dist_aa 
                    + graesning 
                    + log_hojde 
                    + log_TWI 
                    + log_veghojde 
                    + log_dist_kant
                    + jord 
                    + log_pH
                    + N
                    + P
                  + log_afgrode, 
                    data = model_data1,
                    na.action = "na.omit") 

anova(GLM_NOVANA)
summary(GLM_NOVANA)

MuMIn::r.squaredGLMM(GLM_NOVANA)

stargazer(GLM_NOVANA, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")


plot(GLM_NOVANA)

# punkterne som udgøre et problem er 14 og 27, som er to vådepunkter

#----------------Naturlig---------------

GLM_Naturlig <- glm(Naturlig ~ SRI 
                  + dist_aa 
                  + graesning 
                  + markareal
                  + log_hojde 
                  + log_TWI 
                  + log_veghojde 
                  + log_dist_kant
                  + jord 
                  + log_pH
                  + N
                  + P
                  + log_afgrode, 
                  data = model_data1) 

anova(GLM_Naturlig)
summary(GLM_Naturlig)

MuMIn::r.squaredGLMM(GLM_Naturlig)

stargazer(GLM_Naturlig, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

plot(GLM_Naturlig)

# samme gældende


#-----------------Mark-----------
GLM_Mark <- glm(Mark ~ SRI 
                  + markareal
                  + dist_aa 
                  + graesning 
                  + log_hojde 
                  + log_TWI 
                  + markareal
                  + log_veghojde 
                  + log_dist_kant
                  + jord 
                  + log_pH
                  + N
                  + P
                + log_afgrode, 
                  data = model_data1) 

anova(GLM_Mark)
summary(GLM_Mark)

MuMIn::r.squaredGLMM(GLM_Mark)

stargazer(GLM_Mark, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")


plot(GLM_Mark)

# her er det samme våde punkter, men så også 31 (43-5) - 
# hvilket er svære at forklare end at det må været et artsfattigt spot. 

# Mixed effect model ####
#-------------NOVANA----------------

Mixed_NOVANA <- nlme::lme(NOVANA ~ SRI 
                        + dist_aa 
                        + graesning 
                        + log_hojde 
                        + log_TWI 
                        + log_veghojde 
                        + log_dist_kant
                        + jord 
                        + log_pH
                        + N
                        + P
                        + log_afgrode, 
                        (random = ~ 1|markareal), 
                        data = model_data1, method = "ML")

anova(Mixed_NOVANA)
summary(Mixed_NOVANA)
MuMIn::r.squaredGLMM(Mixed_NOVANA)

plot(Mixed_NOVANA)

library(stargazer)

stargazer(Mixed_NOVANA, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

vif(Mixed_NOVANA)

#-----------------Naturlig-----------------------------

Mixed_Naturlig <- nlme::lme(Naturlig ~ SRI 
                          + dist_aa 
                          + graesning 
                          + log_hojde 
                          + log_TWI 
                          + log_veghojde 
                          + log_dist_kant
                          + jord 
                          + log_pH
                          + N
                          + P
                          + log_afgrode, 
                          (random = ~ 1|markareal), 
                          data = model_data1, method = "ML") 


anova(Mixed_Naturlig)
summary(Mixed_Naturlig)
MuMIn::r.squaredGLMM(Mixed_Naturlig)

plot(Mixed_Naturlig)

stargazer(Mixed_Naturlig, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")


vif(Mixed_Naturlig)

#-----------------------Mark-----------------------

Mixed_Mark <- nlme::lme(Mark ~ SRI 
                          + dist_aa 
                          + graesning 
                          + log_hojde 
                          + log_TWI 
                          + jord
                          + log_veghojde 
                          + log_dist_kant
                          + log_pH
                          + N
                          + P
                        + log_afgrode, 
                          (random = ~ 1|markareal), 
                          data = model_data1, method = "ML")

anova(Mixed_Mark)
summary(Mixed_Mark)
MuMIn::r.squaredGLMM(Mixed_Mark)

plot(Mixed_Mark)

stargazer(Mixed_Mark, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")


vif(Mixed_Mark)

# Model selection ####

stepAIC(Mixed_NOVANA, method=ML, direction = "backward")
stepAIC(Mixed_Naturlig, method=ML, direction = "backward")
stepAIC(Mixed_Mark, method=ML, direction = "backward")

#------------------- reduced models -------------------

reduced_NOVANA <- nlme::lme(NOVANA ~ graesning + jord, (random = ~ 1|markareal), 
                            data = model_data1, method = "ML")


reduced_Naturlig <- nlme::lme(Naturlig ~ graesning + log_dist_kant + jord, 
                              (random =  ~ 1 |markareal), data = model_data1, method = "ML")


reduced_Mark <- nlme::lme(Mark ~ graesning + jord + log_afgrode + N, 
                              (random =  ~ 1 |markareal), data = model_data1, method = "ML")



reduced_Mark

AIC(Mixed_NOVANA)
AIC(reduced_NOVANA)
AIC(Mixed_Naturlig)
AIC(reduced_Naturlig)
AIC(Mixed_Mark)
AIC(reduced_Mark)

MuMIn::r.squaredGLMM(reduced_NOVANA)
MuMIn::r.squaredGLMM(reduced_Naturlig)
MuMIn::r.squaredGLMM(reduced_Mark)

#so the best model is the reduced "mark model"


# testing if there is a multiple collinearity problem in the models - all value is below 5 
# there is therfor not any problems 

vif(reduced_NOVANA)
vif(reduced_Naturlig)
vif(reduced_Mark)

# looking at the models (r2) the mixed model is best 
# (from the GLM-models)

#---------------- comparing all the models --------------

# NOVANA 
MuMIn::r.squaredGLMM(Mixed_NOVANA)
MuMIn::r.squaredGLMM(GLM_NOVANA)
MuMIn::r.squaredGLMM(reduced_NOVANA)

# Naturlig 
MuMIn::r.squaredGLMM(Mixed_Naturlig)
MuMIn::r.squaredGLMM(GLM_Naturlig)
MuMIn::r.squaredGLMM(reduced_Naturlig)

# MARK
MuMIn::r.squaredGLMM(Mixed_Mark)
MuMIn::r.squaredGLMM(GLM_Mark)
MuMIn::r.squaredGLMM(reduced_Mark)

stargazer(reduced_NOVANA, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

stargazer(reduced_Naturlig, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

stargazer(reduced_NOVANA, 
          reduced_Naturlig,
          reduced_Mark, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")


?`stargazer`

summary(reduced_NOVANA)
summary(reduced_Naturlig)
summary(reduced_Mark)

str(reduced_Mark)

#---------------- plotting ------------
boxplot(NOVANA ~ markareal, data = model_data1)
boxplot(Naturlig ~ markareal, data = model_data1)
boxplot(Mark ~ markareal, data = model_data1)

#---------------- testing for normality (field-level) -----------

mark33 <- model_data1 %>% filter(markareal == 33)

shapiro.test(mark33$NOVANA) #normally
shapiro.test(mark33$Naturlig) #normally
shapiro.test(mark33$Mark) #normally

mark34 <- model_data1 %>% filter(markareal == 34)

shapiro.test(mark34$NOVANA) # not normally
shapiro.test(mark34$Naturlig) #normally
shapiro.test(mark34$Mark) #normally

mark35 <- model_data1 %>% filter(markareal == 35)

shapiro.test(mark35$NOVANA) #normally
shapiro.test(mark35$Naturlig) #normally
shapiro.test(mark35$Mark) #normally

mark40 <- model_data1 %>% filter(markareal == 40)

shapiro.test(mark40$NOVANA) #normally
shapiro.test(mark40$Naturlig) #normally
shapiro.test(mark40$Mark) #normally

mark41 <- model_data1 %>% filter(markareal == 41)

shapiro.test(mark41$NOVANA) #not normally
shapiro.test(mark41$Naturlig) #normally
shapiro.test(mark41$Mark) #normally

mark42 <- model_data1 %>% filter(markareal == 42)

shapiro.test(mark42$NOVANA) #normally
shapiro.test(mark42$Naturlig) #normally
shapiro.test(mark42$Mark) #normally

mark43 <- model_data1 %>% filter(markareal == 43)

shapiro.test(mark43$NOVANA) #not normally
shapiro.test(mark43$Naturlig) #normally
shapiro.test(mark43$Mark) #normally


#--------------- NOVANA ------------------

artsindex_NOVANA <- lm(NOVANA ~ markareal, data = model_data1) # Testing if there is difference among fields 
art_NOVANA_anova <- aov(artsindex_NOVANA)
summary(art_NOVANA_anova)

tukey_art_NOVANA <- TukeyHSD(art_NOVANA_anova) # Grouping fields based on these results
plot(tukey_art_NOVANA)
summary(tukey_art_NOVANA)
print(tukey_art_NOVANA)

# Species index by field based on Tukeys test
# https://www.r-graph-gallery.com/84-tukey-test.html

# Generating lables for zone groups
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  # Putting labels in same order as in boxplot
  Tukey.labels$markareal = rownames(Tukey.labels)
  Tukey.labels = Tukey.labels[order(Tukey.labels$markareal) , ]
  return(Tukey.labels)
}

# Applying the function on my dataset
LABELS <- generate_label_df(tukey_art_NOVANA, "markareal")
as.factor(LABELS[,1])

# Defining colours
my_colours <- c("white", "white", "white", "gray", "white", "white","gray")

# Drawing boxplot with my new levels
NOV <- boxplot(model_data1$NOVANA ~ as.factor(model_data1$markareal), ylim = c(min(model_data1$NOVANA),
                                                                                        1.1*max(model_data1$NOVANA)),
             col = my_colours,
             xlab = "Markareal", ylab = "Artsindeks - NOVANA", main ="")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( NOV$stats[nrow(NOV$stats),] )

#Add the labels
text( c(1:nlevels(as.factor(model_data1$markareal))) , NOV$stats[nrow(NOV$stats),]+over , LABELS[,1]  , col= "Black" )
levels(model_data1$markareal)

#------------------ NOVANA -----------------------

my_colours <- c("white", "white", "white", "gray", "white", "white","gray")

NOVANA_kruskal <- kruskal.test(NOVANA ~ markareal, data = model_data1)
print(NOVANA_kruskal) # so they are different 

NOVANA_dunn <- dunn_test(NOVANA ~ markareal, data = model_data1)

DUNN <- boxplot(model_data1$NOVANA ~ as.factor(model_data1$markareal), ylim = c(min(model_data1$NOVANA),
                                                                                 1.1*max(model_data1$NOVANA)),
               col = my_colours,
               xlab = "Markareal", ylab = "Artsindeks - NOVANA", main ="")

sig <- c("abc", "abc", "abc", "a", "bc", "c", "ab")

# Signifikants-niveau
nbGroup <- nlevels(model_data1$markareal)
text( 
  x=c(1:nbGroup), 
  y=DUNN$stats[nrow(DUNN$stats),] + 0.5, 
  paste(sig)  
)

#------------------ Naturlig ------------

artsindex_Naturlig <- lm(Naturlig ~ markareal, data = model_data1) # Testing if there is difference among fields 
art_Naturlig_anova <- aov(artsindex_Naturlig)
summary(art_Naturlig_anova)

tukey_art_Naturlig <- TukeyHSD(art_Naturlig_anova) # Grouping fields based on these results
plot(tukey_art_Naturlig)
summary(tukey_art_Naturlig)
print(tukey_art_Naturlig)

# Species index by field based on Tukeys test
# https://www.r-graph-gallery.com/84-tukey-test.html

# Generating lables for field groups
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  # Putting labels in same order as in boxplot
  Tukey.labels$markareal = rownames(Tukey.labels)
  Tukey.labels = Tukey.labels[order(Tukey.labels$markareal) , ]
  return(Tukey.labels)
}

# Applying the function on my dataset
LABELS <- generate_label_df(tukey_art_Naturlig, "markareal")
as.factor(LABELS[,1])

# Defining colours
my_colours <- c("white", "white", "white", "gray", "white", "white","gray")

# Drawing boxplot with my new levels
NAT <- boxplot(model_data1$Naturlig ~ as.factor(model_data1$markareal), ylim = c(min(model_data1$Naturlig),
                                                                               1.1*max(model_data1$Naturlig)),
               col = my_colours,
               xlab = "Markareal", ylab = "Artsindeks - Naturlig", main ="")

# Signifikants-niveau
over <- 0.1*max( NAT$stats[nrow(NAT$stats),] )

#Add the labels
text( c(1:nlevels(as.factor(model_data1$markareal))) , NAT$stats[nrow(NAT$stats),]+over , LABELS[,1]  , col= "Black" )
levels(model_data1$markareal)


#---------------- Mark ------------------------

artsindex_Mark <- lm(Mark ~ markareal, data = model_data1) # Testing if there is difference among fields 

art_Mark_anova <- aov(artsindex_Mark)

summary(art_Mark_anova)

tukey_art_Mark <- TukeyHSD(art_Mark_anova) # Grouping fields based on these results
plot(tukey_art_Mark)
summary(tukey_art_Mark)
print(tukey_art_Mark)

# Species index by field based on Tukeys test
# https://www.r-graph-gallery.com/84-tukey-test.html

# Generating lables for zone groups
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  # Putting labels in same order as in boxplot
  Tukey.labels$markareal = rownames(Tukey.labels)
  Tukey.labels = Tukey.labels[order(Tukey.labels$markareal) , ]
  return(Tukey.labels)
}

# Applying the function on my dataset
LABELS <- generate_label_df(tukey_art_Mark, "markareal")
as.factor(LABELS[,1])


# Defining colours
my_colours <- c("white", "white", "white", "gray", "white", "white","gray")

# Drawing boxplot with my new levels
M <- boxplot(model_data1$Mark ~ as.factor(model_data1$markareal), ylim = c(min(model_data1$Mark),
                                                                               1.1*max(model_data1$Mark)),
             col = my_colours, 
               xlab = "Markareal", ylab = "Artsindeks - Mark", main ="")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( M$stats[nrow(M$stats),] )

#Add the labels
text( c(1:nlevels(as.factor(model_data1$markareal))) , M$stats[nrow(M$stats),]+over , LABELS[,1]  , col= "Black" )
levels(model_data1$markareal)

#--------------- end ---------------------

