# Speciale_Ribe_Holme

# Thesis data by Cecilie Kristensen 
# For questions and more information: ceciliekristensen1009@gmail.com

This github contains the data (6 csv-files) and the scripts (5 R-scripts) 
for the analysis botanical data from Ribe Holme:

The R-scripts:
"NMDS_presence_git.R"
"NMDS_raunkiaer_git.R"
"model_selection_git.R"
"citation_packages_git.R"
"sri_calculation_JO_CE_git.R" is created by Jonathan von Oppen with modification by Cecilie Kristensen for use to the data from Ribe Holme

# csv-files:

The csv-file "total_artslist.csv" contains all the names of the plant species found at Ribe Holme

The csv-file "presence_all.csv" contains a data matrix of all the presence/absence data (1/0) 
for the floristic species found in the 33 plots (5-m circles)

The csv-file "p_raunkiaer_all.csv" contains a data matrix of alle the raunkiaer data (a frequency of 0 - 10)
for the floristic species found in the 29 plots (5-m circles with 10 raunkiaer circles)

The csv-file "NMDS_afgrode.csv" contains a weighted estimate of grass mixtures in the crop rotation system the last 5 years.
Weighted for the precence of grass mixture: 2020: 0,5, 2019: 0,2, 2018: 0,15, 2017: 0,10, 2016: 0,05.

The csv-file "scale_miljo" cointains environmental data for the 33 plot (5-circles) and 3 different calculated species index.
The species index is calculated by the NOVANA method with different weighted systems (NOVANA, Naturlig and Mark).
NOVANA is calulated with no adjustments. Naturlig, the score 0 is given to the sown species. Mark, native species with the 
score -1 is given the score 1 and sown species is given the score from NOVANA. All the data is scaled. 

The csv-file "miljo_data_indeks_sum" cointains environmental data for the 33 plot (5-circles) and 3 different calculated species index.
The species index is calculated by the NOVANA method with different weighted systems (NOVANA, Naturlig and Mark).
NOVANA is calulated with no adjustments. Naturlig, the score 0 is given to the sown species. Mark, native species with the 
score -1 is given the score 1 and sown species is given the score from NOVANA. Non of the data is scaled. 

# R-scripts: 

The script "NMDS_presence_git.R" is containing the code for runing a NMDS-ordination on the floristic presence data.
Stress-levels is also tested for choosing the correct numbers of axes. 
Environmental data is added as vectors and all is visualized by the ggplot2 package. 

The script "NMDS_raunkiaer_git.R" is containing the code for runing a NMDS-ordination on the floristic raunkiaer data.
Stress-levels is also tested for choosing the correct numbers of axes. 
Environmental data is added as vectors and all is visualized by the ggplot2 package

The script "model_selection_git.R" has the purpose of identifying the variables that explain the 3 species index 
and testing for significant difference between fields species index.
The data is tested for normality, coliniarity and correlation. 
The models is fit as a Mixed Effect Model after testing statical asumptions for other models and by comparisons. 
The best model is selected by using The Akaies Information Criterion to run a backward model selection. 

The script "sri_calculation_JO_CE_git.R" is used for calculate the SRI for the 33 plots. 
The SRI is calculated for the time of fieldwork. 

The script "citation_packages_git.R" contains all the reference of the used packages in the above scripts. 


