library(tidyverse)
library(FactoMineR)
library(factoextra)

#####
#Data pre-processing
men_CL <- merge(men_COAD, men_LUAD[1:22], by.x = "lab_id", by.y = "lab_id")
men_LP <- merge(men_LUAD, men_PRAD[1:22], by.x = "lab_id", by.y = "lab_id")
men_total <- merge(men_CL, men_LP, by.x = "lab_id", by.y = "lab_id")

men_total <- men_total |> select(-c(78:91))
men_total <- men_total |> select(lab_id,starts_with("hsa-"),22:35)

new_column_names <- gsub("\\.x", "", colnames(men_total)[85:98])
colnames(men_total)[85:98] <- new_column_names

#WOMEN
women_BC <- merge(women_BRCA[,1:28], women_COAD[,1:25], by.x = "lab_id", by.y = "lab_id")
women_total <- merge(women_LUAD[,1:25], women_BC, by.x = "lab_id", by.y = "lab_id")
#Recoding levels

men_total <- men_total %>%
  mutate(School = recode(School, "media_inf" = "middle_school", "superiori" = "high_school"))

men_total$Antidepressant <- ifelse(men_total$Antidepressant == "si","yesAntidep","noAntidep")

men_total <- men_total %>%
  mutate(Smoking = recode(Smoking, "siSmoke" = "yesSmoke"))

men_total <- men_total %>%
  mutate(Workingstatus = recode(Workingstatus, "occupato" = "Employed",
                                "disoccupato" = "Unemployed",
                                "pensionato" = "Retired"))

men_total$Passive_Smoke <- ifelse(men_total$Passive_Smoke == "si", "yesPassive_Smoke","noPassive_Smoke")

###

men_COAD <- men_COAD %>%
  mutate(School = recode(School, "media_inf" = "middle_school", "superiori" = "high_school"))

men_COAD$Antidepressant <- ifelse(men_COAD$Antidepressant == "si","yesAntidep","noAntidep")

men_COAD <- men_COAD %>%
  mutate(Smoking = recode(Smoking, "siSmoke" = "yesSmoke"))

men_COAD <- men_COAD %>%
  mutate(Workingstatus = recode(Workingstatus, "occupato" = "Employed",
                                "disoccupato" = "Unemployed",
                                "pensionato" = "Retired"))

men_COAD$Passive_Smoke <- ifelse(men_COAD$Passive_Smoke == "si", "yesPassive_Smoke","noPassive_Smoke")

data_pca_COAD <- men_COAD |> select(-lab_id) |> na.omit()
pca_COAD <- FAMD(data_pca_COAD, sup.var = 1:20, ncp = 20)
a2 <- plot(pca_COAD, "quanti",
     title = "COAD: Graph of the quantitative variables (PC1 & 2)")
b2 <- plot(pca_COAD, "quanti",
     title = "COAD: Graph of the quantitative variables (PC2 & 3)",
     axes = c(2,3))

a2 + b2

###
men_LUAD <- men_LUAD %>%
  mutate(School = recode(School, "media_inf" = "middle_school", "superiori" = "high_school"))

men_LUAD$Antidepressant <- ifelse(men_LUAD$Antidepressant == "si","yesAntidep","noAntidep")

men_LUAD <- men_LUAD %>%
  mutate(Smoking = recode(Smoking, "siSmoke" = "yesSmoke"))

men_LUAD <- men_LUAD %>%
  mutate(Workingstatus = recode(Workingstatus, "occupato" = "Employed",
                                "disoccupato" = "Unemployed",
                                "pensionato" = "Retired"))

men_LUAD$Passive_Smoke <- ifelse(men_LUAD$Passive_Smoke == "si", "yesPassive_Smoke","noPassive_Smoke")

data_pca_LUAD <- men_LUAD |> select(-lab_id) |> na.omit()
pca_LUAD <- FAMD(data_pca_LUAD, sup.var = 1:21, ncp = 20)
a2 <- plot(pca_LUAD, "quanti",
           title = "LUAD: Graph of the quantitative variables (PC1 & 2)")
b2 <- plot(pca_LUAD, "quanti",
           title = "LUAD: Graph of the quantitative variables (PC2 & 3)",
           axes = c(2,3))

a2 + b2

###

men_PRAD <- men_PRAD %>%
  mutate(School = recode(School, "media_inf" = "middle_school", "superiori" = "high_school"))

men_PRAD$Antidepressant <- ifelse(men_PRAD$Antidepressant == "si","yesAntidep","noAntidep")

men_PRAD <- men_PRAD %>%
  mutate(Smoking = recode(Smoking, "siSmoke" = "yesSmoke"))

men_PRAD <- men_PRAD %>%
  mutate(Workingstatus = recode(Workingstatus, "occupato" = "Employed",
                                "disoccupato" = "Unemployed",
                                "pensionato" = "Retired"))

men_PRAD$Passive_Smoke <- ifelse(men_PRAD$Passive_Smoke == "si", "yesPassive_Smoke","noPassive_Smoke")

data_pca_PRAD <- men_PRAD |> select(-lab_id) |> na.omit()
pca_PRAD <- FAMD(data_pca_PRAD, sup.var = 1:21, ncp = 20)
a2 <- plot(pca_PRAD, "quanti",
           title = "PRAD: Graph of the quantitative variables (PC1 & 2)")
b2 <- plot(pca_PRAD, "quanti",
           title = "PRAD: Graph of the quantitative variables (PC2 & 3)",
           axes = c(2,3))

a2 + b2

#####
#FAMD

data_pca <- men_total |> select(-lab_id) |> na.omit()
pca1 <- FAMD(data_pca, sup.var = 1:83, ncp = 20)
pca2 <- FAMD(data_pca[,84:97], ncp = 20)

fviz_screeplot(pca1, ncp = 20)

a <- plot(pca1, "quanti", title = "Graph of the quantitative variables (PC1 & 2)")
b <- plot(pca1, "quanti", axes = c(2,3), title = "Graph of the quantitative variables (PC2 & 3)")

c <- plot(pca1, "quali")
d <- plot(pca1, "quali", axes = c(2,3))

e <- plot(pca2, "var")
f <- plot(pca2, "var", axes = c(2,3))

#Figura1
(a + b) / (e + f)
#Figura2
c + d

######
#WOMEN

colnames(women_COAD)[26] <- "Age"
colnames(women_COAD)[27] <- "BMI"
colnames(women_COAD)[28] <- "PM10_Annual_Average"
colnames(women_COAD)[29] <- "School"
colnames(women_COAD)[30] <- "Antidepressant"
colnames(women_COAD)[31] <- "PCR"
colnames(women_COAD)[32] <- "Homocysteine"
colnames(women_COAD)[33] <- "Smoking"
colnames(women_COAD)[34] <- "Traffic_no"
colnames(women_COAD)[35] <- "Physical_Activity"
colnames(women_COAD)[36] <- "Workingstatus"
colnames(women_COAD)[37] <- "Passive_Smoke"
colnames(women_COAD)[38] <- "Birth_Control_Pill"
colnames(women_COAD)[39] <- "No_Pregnancies"
colnames(women_COAD)[42] <- "Menopause"


women_COAD <- women_COAD %>%
  mutate(School = recode(School, "media_inf" = "middle_school", "superiori" = "high_school"))

women_COAD$Antidepressant <- ifelse(women_COAD$Antidepressant == "antid_yes","yesAntidep","noAntidep")

women_COAD <- women_COAD %>%
  mutate(Smoking = recode(Smoking, "siSmoke" = "yesSmoke"))

women_COAD <- women_COAD %>%
  mutate(Workingstatus = recode(Workingstatus, "occupato" = "Employed",
                                "disoccupato" = "Unemployed",
                                "pensionato" = "Retired",
                                "casalinga" = "Housewife"))

women_COAD$Passive_Smoke <- ifelse(women_COAD$Passive_Smoke == "si", "yesPassive_Smoke","noPassive_Smoke")

women_COAD$Birth_Control_Pill <- ifelse(women_COAD$Birth_Control_Pill == "oral_no", "noBC_Pill","yesBC_Pill")

women_COAD$Menopause <- ifelse(women_COAD$Menopause == "menop", "yesMenopause","noMenopause")


data_pca_COAD <- women_COAD |> select(-lab_id) |> na.omit()
pca_COAD_women <- FAMD(data_pca_COAD, sup.var = 1:24, ncp = 20)
pca_COAD_women2 <- FAMD(data_pca_COAD[,-c(1:24)])
a2 <- plot(pca_COAD_women, "quanti",
           title = "Women COAD: Graph of the quantitative variables (PC1 & 2)")
b2 <- plot(pca_COAD_women, "quanti",
           title = "Women COAD: Graph of the quantitative variables (PC2 & 3)",
           axes = c(2,3))

c2 <- plot(pca_COAD_women, "quali")
d2 <- plot(pca_COAD_women, "quali", axes = c(2,3))

e2 <- plot(pca_COAD_women2, "var")
f2 <- plot(pca_COAD_women2, "var", axes = c(2,3))

(a2 + b2) / (e2 + f2)

c2 + d2

###
#LUAD

colnames(women_LUAD)[26] <- "Age"
colnames(women_LUAD)[27] <- "BMI"
colnames(women_LUAD)[28] <- "PM10_Annual_Average"
colnames(women_LUAD)[29] <- "School"
colnames(women_LUAD)[30] <- "Antidepressant"
colnames(women_LUAD)[31] <- "PCR"
colnames(women_LUAD)[32] <- "Homocysteine"
colnames(women_LUAD)[33] <- "Smoking"
colnames(women_LUAD)[34] <- "Traffic_no"
colnames(women_LUAD)[35] <- "Physical_Activity"
colnames(women_LUAD)[36] <- "Workingstatus"
colnames(women_LUAD)[37] <- "Passive_Smoke"
colnames(women_LUAD)[38] <- "Birth_Control_Pill"
colnames(women_LUAD)[39] <- "No_Pregnancies"
colnames(women_LUAD)[42] <- "Menopause"


women_LUAD <- women_LUAD %>%
  mutate(School = recode(School, "media_inf" = "middle_school", "superiori" = "high_school"))

women_LUAD$Antidepressant <- ifelse(women_LUAD$Antidepressant == "antid_yes","yesAntidep","noAntidep")

women_LUAD <- women_LUAD %>%
  mutate(Smoking = recode(Smoking, "siSmoke" = "yesSmoke"))

women_LUAD <- women_LUAD %>%
  mutate(Workingstatus = recode(Workingstatus, "occupato" = "Employed",
                                "disoccupato" = "Unemployed",
                                "pensionato" = "Retired",
                                "casalinga" = "Housewife"))

women_LUAD$Passive_Smoke <- ifelse(women_LUAD$Passive_Smoke == "si", "yesPassive_Smoke","noPassive_Smoke")

women_LUAD$Birth_Control_Pill <- ifelse(women_LUAD$Birth_Control_Pill == "oral_no", "noBC_Pill","yesBC_Pill")

women_LUAD$Menopause <- ifelse(women_LUAD$Menopause == "menop", "yesMenopause","noMenopause")


data_pca_LUAD <- women_LUAD |> select(-lab_id) |> na.omit()
pca_LUAD_women <- FAMD(data_pca_LUAD, sup.var = 1:24, ncp = 20)
pca_LUAD_women2 <- FAMD(data_pca_LUAD[,-c(1:24)])
a3 <- plot(pca_LUAD_women, "quanti",
           title = "Women LUAD: Graph of the quantitative variables (PC1 & 2)")
b3 <- plot(pca_LUAD_women, "quanti",
           title = "Women LUAD: Graph of the quantitative variables (PC2 & 3)",
           axes = c(2,3))

c3 <- plot(pca_LUAD_women, "quali")
d3 <- plot(pca_LUAD_women, "quali", axes = c(2,3))

e3 <- plot(pca_LUAD_women2, "var")
f3 <- plot(pca_LUAD_women2, "var", axes = c(2,3))

(a3 + b3) / (e3 + f3)

c3 + d3

###
#BRCA

colnames(women_BRCA)[29] <- "Age"
colnames(women_BRCA)[30] <- "BMI"
colnames(women_BRCA)[31] <- "PM10_Annual_Average"
colnames(women_BRCA)[32] <- "School"
colnames(women_BRCA)[33] <- "Antidepressant"
colnames(women_BRCA)[34] <- "PCR"
colnames(women_BRCA)[35] <- "Homocysteine"
colnames(women_BRCA)[36] <- "Smoking"
colnames(women_BRCA)[37] <- "Traffic_n"
colnames(women_BRCA)[38] <- "Physical_Activity"
colnames(women_BRCA)[39] <- "Workingstatus"
colnames(women_BRCA)[40] <- "Passive_Smoke"
colnames(women_BRCA)[41] <- "Birth_Control_Pill"
colnames(women_BRCA)[42] <- "n_Pregnancies"
colnames(women_BRCA)[45] <- "Menopause"


women_BRCA <- women_BRCA %>%
  mutate(School = recode(School, "media_inf" = "middle_school", "superiori" = "high_school"))

women_BRCA$Antidepressant <- ifelse(women_BRCA$Antidepressant == "antid_yes","yesAntidep","noAntidep")

women_BRCA <- women_BRCA %>%
  mutate(Smoking = recode(Smoking, "siSmoke" = "yesSmoke"))

women_BRCA <- women_BRCA %>%
  mutate(Workingstatus = recode(Workingstatus, "occupato" = "Employed",
                                "disoccupato" = "Unemployed",
                                "pensionato" = "Retired",
                                "casalinga" = "Housewife"))

women_BRCA$Passive_Smoke <- ifelse(women_BRCA$Passive_Smoke == "si", "yesPassive_Smoke","noPassive_Smoke")

women_BRCA$Birth_Control_Pill <- ifelse(women_BRCA$Birth_Control_Pill == "oral_no", "noBC_Pill","yesBC_Pill")

women_BRCA$Menopause <- ifelse(women_BRCA$Menopause == "menop", "yesMenopause","noMenopause")


data_pca_BRCA <- women_BRCA |> select(-lab_id) |> na.omit()
pca_BRCA_women <- FAMD(data_pca_BRCA, sup.var = 1:27, ncp = 20)
pca_BRCA_women2 <- FAMD(data_pca_BRCA[,-c(1:27)])
a4 <- plot(pca_BRCA_women, "quanti",
           title = "Women BRCA: Graph of the quantitative variables (PC1 & 2)")
b4 <- plot(pca_BRCA_women, "quanti",
           title = "Women BRCA: Graph of the quantitative variables (PC2 & 3)",
           axes = c(2,3))

c4 <- plot(pca_BRCA_women, "quali")
d4 <- plot(pca_BRCA_women, "quali", axes = c(2,3))

e4 <- plot(pca_BRCA_women2, "var")
f4 <- plot(pca_BRCA_women2, "var", axes = c(2,3))

(a4 + b4) / (e4 + f4)

c4 + d4



######
#Correlation plot

men_corr_filtered <- men_total %>%
  select(-c(85:98),-matches("\\."), -lab_id)
men_corr_filtered <- cbind(men_corr_filtered,men_LUAD$`hsa-miR-184`)
names(men_corr_filtered)[36] <- "hsa-miR-184"

women_corr_filtered <- women_total |> select(-matches("\\."), -lab_id)

#MEN
corr_men <- round(cor(men_corr_filtered, method = "spearman"), 3) # rounded to one decimal point
ggcorrplot(corr_men, title = "Correlation matrix, clinical variables")

s1 <- rank(apply(corr_men, 2, mean))
s2 <- rank(apply(corr_men, 1, mean))
corr_men <- corr_men[,order(s1)]
corr_men <- corr_men[order(s2),]

a <- ggcorrplot(corr_men,  title = "A. Correlation between miRNAs in men")


#WOMEN
corr_women <- round(cor(women_corr_filtered, method = "spearman"), 3) # rounded to one decimal point
ggcorrplot(corr_women, title = "Correlation matrix, clinical variables")

s1 <- rank(apply(corr_women, 2, mean))
s2 <- rank(apply(corr_women, 1, mean))
corr_women <- corr_women[,order(s1)]
corr_women <- corr_women[order(s2),]

b <- ggcorrplot(corr_women,  title = "B. Correlation between miRNAs in women")






