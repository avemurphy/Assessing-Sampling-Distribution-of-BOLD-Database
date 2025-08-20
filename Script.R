####-Importing Libraries and Data----
library(tidyverse)
library(ggplot2)
library(dplyr)

# Data from BOLD was downloaded on September 22, 2024
# First the working directory is set using "session", "set working directory" and set it with the "to source file location" option 
df_Gast <- read_tsv(file = "../data/Gastrotricha_data.tsv") 

# A a subset called df_Gast.sub is made with the key variables that could be considered for this exploration.  This includes variables such as BIN ID, location data, and taxon data (family/species) that will be used to explore sampling distributions of this phylum. 
names(df_Gast)
df_Gast.sub <- df_Gast[, c("processid", "sampleid","recordID","bin_uri","phylum_name","class_name","order_name","family_name","subfamily_name","genus_name","species_name","subspecies_name","lat","lon","country","markercode")]

####-Part.1: Looking at the Distribution of Data----
####__1. Sample Number for each BIN----
# The first step to looking at sampling bias across the different BIN's for the phylum Gastrotricha, is to identify all the  COI gene barcodes that have been identified and then the number of specimens that have been sampled for each barcode.  To do this a subset of the dataframe was taken to only include rows or specimens with a barcode that came from a COI sequence as this is the current marker most widely used.  This also removes redundant specimen that could have data for more than 1 bio marker. The data set is also filtered to remove specimens that are not assigned to a BIN. 

df_Gast.COI <- df_Gast.sub %>% 
  filter(markercode == "COI-5P") %>% 
  filter(!is.na(bin_uri))

# Now using this data set the unique barcodes were found by looking the unique BIN_uri that is used to distinguish each different COI barcode in this family.
unique(df_Gast.COI$bin_uri)
length(unique(df_Gast.COI$bin_uri))

# There is 199 unique COI sequence barcodes for this phylum, a table can be made to look at the number of specimens that have been sampled for each BIN.
sort(table(df_Gast.COI$bin_uri), decreasing = TRUE)

max(sort(table(df_Gast.COI$bin_uri), decreasing = TRUE))
# The maximum amount of samples for one BIN is 60

min(sort(table(df_Gast.COI$bin_uri), decreasing = TRUE))
# The minimum amount of samples for one BIN is 1 

median(sort(table(df_Gast.COI$bin_uri), decreasing = TRUE))
# The median sample number in a BIN among the 199 BIN's is 1, indicating that half the BIN's in the phylum only have 1 sample.

sort(table(df_Gast.COI$bin_uri), decreasing = TRUE)[1:10]
# Within just the 10 BIN's with the most samples, the sample number differs from 60 for one BIN and 4 in another BIN indicating a large distribution of samples among the BIN's. 

# Looking through the table and different values, this shows how different the amount of specimen samples are for each BIN in this phylum. The frequency of specimen values in each BIN can be looked at across the data set:

BIN_URI <- df_Gast.COI %>% 
  group_by(bin_uri) %>% 
  summarize(total = n()) 

ggplot(BIN_URI, aes(x = total)) +labs(title = "Distibustion of Sample Numbers in BINS ", x = "Number of Samples per BIN", y = "Frequnecy of Value ") + theme(plot.title = element_text(hjust = 0.5)) + geom_histogram()
# Looking at the distribution of sample numbers in BINS, it can be seen that the majority of the frequency values fall on the left side and below 20 with only a very low frequency of BINS having higher samples numbers.  At a first glance of this data without context this would look as if there was data out of range or an outlier. The composition of the BINS will be explored to determine why 2 have much higher sample numbers. 

# This bar plot shows the names of the BINS with sample numbers over 3 to identify the BINs with high sample numbers and lower frequencies on the histogram. 

BIN_URI.2 <- df_Gast.COI %>% 
  group_by(bin_uri) %>% 
  summarize(total = n()) %>% 
  filter(total >= 3)

ggplot(BIN_URI.2, aes(x = bin_uri, y = total)) + geom_bar (stat = "identity") +
  theme(text = element_text(size = 20)) + labs(title = "Number of Samples in Each BIN", x = "Barcode Index Number", y = "Number of Samples") + theme(axis.text.x = element_text(size = rel(0.5), angle = 90)) 

# The barplot shows visually that there is very different specimen numbers among the BIN's in the Gastrotricha phylum. BIN BOLD:ACB6539 and BOLD:ACB6539 specifically show a very high specimen number compared to others.  A sampling bias is occurring creating large differences in barcode sample numbers as 1 has 60, one has 34, 11 and the remaining below 10. The next section will explore what factors could be contributing to some barcode sequences (BIN's) having much more specimen data than others.

####-Part 2: Exploring Sample Bias Across BINS ----
####__1. Families in Each BIN----
# To explore the potential biases that could be contributing to the large number of sample differences among BIN's, first the family composition of BINS will be looked at.

dfGroup_BIN <- df_Gast.COI %>% 
  group_by(family_name, bin_uri) %>% 
  count(bin_uri)

# 2 charts will be made to visualize the family composition across BINS with 2.  The first will show the family composition of BINS with 2 samples and over:

dfGroup_BIN2 <- df_Gast.COI %>% 
  group_by(family_name, bin_uri) %>% 
  count(bin_uri) %>% 
  filter(n >= 2)

ggplot(dfGroup_BIN2, aes(x = bin_uri, fill = family_name )) + theme(axis.text.x= element_text(size=rel(0.5), angle = 90)) + labs(title = "Family Composotion in Each BIN", x = "Barcode Index Number", y = "Proportion of Family in Each BIN") + geom_bar(position = "fill")

# The second plot will will show the family composition of BINS with 1 sample:
dfGroup_BIN3 <- df_Gast.COI %>% 
  group_by(family_name, bin_uri) %>% 
  count(bin_uri) %>% 
  filter(n <= 1)

ggplot(dfGroup_BIN3, aes(x = bin_uri, fill = family_name )) + theme(axis.text.x= element_text(size = rel(0.5), angle = 90)) +labs(title = "Family Composotion in Each BIN", x = "Barcode Index Number", y = "Proportion of Family in Each BIN") + geom_bar(position = "fill")

# It can be seen that the BINs with higher specimen numbers come from specimens of manily 2 family members and among BINS with lower sample numbers there is more families represented. 

# The data was spread to show the composition per BIN:
dfBINs.spread <- pivot_wider(data = dfGroup_BIN, names_from  = bin_uri, values_from = n)
# After looking at the data, a subset was made from the 10 BIN's with the most samples that were identified earlier to look at the family composition of each BIN.  
family <- (dfBINs.spread[, c("family_name", "BOLD:ACB6539", "BOLD:ACB6653", "BOLD:ACB6843", "BOLD:ACB6635", "BOLD:ADX5087", "BOLD:ADW6954", "BOLD:ACB6710", "BOLD:ADW0794", "BOLD:ACB6555", "BOLD:ADW0588")])
print(family, n = Inf, width = Inf)

# The 12th and 2nd row composed of all the samples in the 10 most plentiful BIN's, this is the family Turbanellidae and Chaetonotidae.
print(family[12,], n = Inf, width = Inf)
print(family[2,], n = Inf, width = Inf)

# The family composition of the 10 bins with the least amount of samples was also explored.
sort(table(df_Gast.COI$bin_uri), decreasing = FALSE)[1:10]
print(dfBINs.spread[, c("family_name", "BOLD:ACB6475", "BOLD:ACB6478", "BOLD:ACB6496", "BOLD:ACB6517", "BOLD:ACB6634", "BOLD:ACB6705", "BOLD:ACB6707", "BOLD:ACS6640", "BOLD:ACV7545", "BOLD:ADJ4136")], n = Inf, width = Inf)
# All these BINS contained one sample, more families were represented however most of them were again from the Family Turbanellidae.

# The total amount of each family members across the BINS can be calculated
dfBINs.spread <- dfBINs.spread %>% mutate(row_sum = rowSums(across(2:199), na.rm = TRUE))
dfBINs.spread[,c(1,201)]
# Pie chart to visualize the family composition
family.number <- dfBINs.spread$row_sum
names(family.number) <- dfBINs.spread$family_name
pie(family.number, labels = NA, main="Family Distribution", col = topo.colors(length(family.number)))
legend(x = -2, y = 1, legend = names(family.number), fill = topo.colors(length(family.number)), bty = "n", cex = 0.5)

# It can be seen that the barcodes in the phlyum Gastrotricha predominantly come from specimens in the families Turbanellidae and Chaetonotidae, and that the bins BOLD:ACB6539` `BOLD:ACB6653` with the greatest samples are composed entirely of samples from the Turbanellidae family.  Other families represented such as Planodasyidaa, Xenotrichulidae and Muselliferidae come from BINS with low sample numbers. 

pie(family.number, labels = NA, main="Family Distribution", col=topo.colors(length(family.number)))
legend(x = -2, y = 1, legend = names(family.number), fill = topo.colors(length(family.number)), bty = "n", cex = 0.5)

####__2. Exploring Species in Each BIN----
# Next to continue exploring potential bias among the BINs in this phylum the species composition in each BIN can be explored. Many specimens do not have a species assignment so they were removed from this exploration.

dfGroup_BIN_spec <- df_Gast.COI %>% 
  filter(!is.na(species_name)) %>% 
  group_by(species_name, bin_uri) %>% 
  count(bin_uri)

dfBINs.spread_spec <- pivot_wider(data = dfGroup_BIN_spec, names_from  = bin_uri, values_from = n)

# After looking at the data, again a subset was made from the 10 BIN's with the most samples that were identified earlier to look at the species composition of each BIN.
species <- (dfBINs.spread_spec[, c("species_name", "BOLD:ACB6539", "BOLD:ACB6653", "BOLD:ACB6843", "BOLD:ACB6635", "BOLD:ADX5087", "BOLD:ADW6954", "BOLD:ACB6710", "BOLD:ADW0794", "BOLD:ACB6555", "BOLD:ADW0588")])
print(species, n = Inf, width = Inf)

print(dfBINs.spread_spec[, c("species_name", "BOLD:ACB6475", "BOLD:ACB6478", "BOLD:ACB6496", "BOLD:ACB6517", "BOLD:ACB6634", "BOLD:ACB6705", "BOLD:ACB6707", "BOLD:ACS6640", "BOLD:ACV7545", "BOLD:ADJ4136")], n = Inf, width = Inf)
# Majority of species came from the Turbanella genus again showing that only one family was well sampled.  

dfBINs.spread_spec <- dfBINs.spread_spec %>% mutate(row_sum = rowSums(across(2:148), na.rm = TRUE))
(dfBINs.spread_spec[,c(1,150)]) %>% print(n = Inf, width = Inf)
# Many species only have 1 sample, while those from the Turbanella have many.  The species rows with high sample numbers were looked at:
dfBINs.spread_spec[106,] %>% print(n = Inf, width = Inf)
dfBINs.spread_spec[105,] %>% print(n = Inf, width = Inf)

# There is not a equal sampling distribution among species, many of the BINS with only one or low specimen samples are from Aspidiophorus and Chaetonotus species. The majority of samples are from species in the Turbanellidae family. 
main_species <- dfGroup_BIN_spec %>% 
  filter(species_name %in% 
  c("Turbanella hyalina", 'Turbanella cornuta'))

ggplot(main_species, aes(x = bin_uri, y = n, fill = species_name)) + theme(axis.text.x= element_text(size=rel(0.5), angle = 90)) + labs(title = "Highly Sampled Species in BIN's", x = "Barcode Index Number", y = "Number of Specimens in Each BIN") + theme(plot.title = element_text(hjust = 0.5)) + geom_bar (stat = "identity")

# It is seen that there is more sampling for species in the Turbanellidae family specifically the species Turbanella hyalina and Turbanella cornuta.  The BIN BOLD:ACB6539 that has the most samples is made entirely of Turbanella hyalina samples.  This species is highly sampled as it has also been found in other BINS meaning that specimen from this same species have different sequences for their COI barcode.  The high sampling of this species has allowed for other barcodes to be found where as some BINs still only have 1 specimen and some species have only been sampled 1 time. 

####__2. Exploring Countries of Each BIN----
# Lastly after exploring the sampling bias among family and species, the country where these samples are being found can also be explored.  Now that the main species and BINS with the majority of samples have been explored, the countries the samples are coming from can be explored.  The data was organized so all the rows sharing the same species, BIN and country were grouped and those without a species and country were removed. 

country <- df_Gast.COI %>% 
  group_by(species_name, bin_uri, country) %>% 
  filter(!is.na(country)) %>% 
  filter(!is.na(species_name)) %>% 
  count(bin_uri) 

# The species that start with genus Turbanella were looked at to see where the most dominate species were located 
country %>% filter(str_detect(species_name, 'Turbanella')) %>% print(n = Inf, width = Inf)
# The countries of specimen samples in the most sampled BINs were looked at. 
country %>% filter(bin_uri == 'BOLD:ACB6539')
country %>% filter(bin_uri == 'BOLD:ACB6653')
# Saw that Brazil was not as represented in the data.
country %>% filter(country == 'Brazil')

# It appears that the species from the Turbanellidae family were sampled from countries in Europe specifically the top BINS have many specimen from Germany. Samples from Brazil are for species with less representatives and BINs with less samples.  

####-Part.3: Looking At Sampling Completeness for Gastrotricha Phylum----
# Through out this exploration it is seen that there is not equal specimen sampling across BINs and that samples are representative of 2 main families and only a few species.  All species name and the number of specimen samples for each species were made into a vector and use to model a sample completeness curve (Type = 2) based on the reading iNEXT: an R package for rarefaction and extrapolation of species diversity (Hill numbers) by Hsieh et al in 2016.  This will plot the sample coverage based on the sample size (Chao & Hu, 2024) 

library(iNEXT)
# Check the documentation of 
spec_diversity <- df_Gast.COI %>% 
  group_by(species_name) %>% 
  summarize(total= n()) 

#The sepcies and their total were made into a vector 
community_data <- spec_diversity$total
names(community_data) <- spec_diversity$species_name

# First the curve was made using q = 0 which uses a measure of the number of different species present.
iNEXT_result <- iNEXT(community_data, q = 0, datatype = "abundance")
plot(iNEXT_result, type = 2)
# Then the sample completeness curve was made using q = 1 which refers to the Shannon diversity index which takes into account the evenness across different species which is important as sample numbers were very different. 
iNEXT_result.2 <- iNEXT(community_data, q = 1, datatype = "abundance")
plot(iNEXT_result.2, type = 2)
title(main= "Sample Coverage Curve for Gastrotricha Phylum", line = 0.8)

# It is seen when taking into account the number of individuals sequenced for the COI barcode, that the sample coverage has yet to reach 1 meaning that not all the diversity has been found and as sampling increases more BINS will be created. 
 
