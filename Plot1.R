library(tidyverse)
library(FactoMineR)

#####
specified_mirna_COAD<-subset( Mirna2, Cancer=="COAD" & oncomirna %in% c("hsa-miR-324-3p", "hsa-miR-671-3p", "hsa-miR-139-3p",
                                                                        
                                                                        "hsa-miR-328-3p", "hsa-miR-652-3p", "hsa-miR-193a-5p", "hsa-miR-574-3p",
                                                                        
                                                                        "hsa-miR-532-3p", "hsa-miR-628-5p","hsa-miR-136-5p", "hsa-miR-29a-5p",
                                                                        
                                                                        "hsa-miR-126-5p","hsa-miR-29b-3p",
                                                                        
                                                                        "hsa-miR-142-3p", "hsa-miR-125a-5p", "hsa-miR-197-3p","hsa-miR-378a-3p",
                                                                        
                                                                        "hsa-miR-21-5p", "hsa-miR-126-3p","hsa-miR-144-5p"))
COAD1 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-125a-5p")

COAD2 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-126-3p")

COAD3 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-126-5p")

COAD4 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-136-5p")

COAD5 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-139-3p")

COAD6 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-142-3p")

COAD7 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-144-5p")

COAD8 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-193a-5p")

COAD9 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-197-3p")

COAD10 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-21-5p")

COAD11 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-29a-5p")

COAD12 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-29b-3p")

COAD13 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-324-3p")

COAD14 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-328-3p")

COAD15 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-378a-3p")

COAD16 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-532-3p")

COAD17 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-574-3p")

COAD18 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-628-5p")

COAD19 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-652-3p")

COAD20 <- subset(specified_mirna_COAD, oncomirna == "hsa-miR-671-3p")

Men2_extracted <- Men2 %>%
  distinct(lab_id, Cancer, .keep_all = TRUE) %>%
  select(lab_id, Cancer, Residence_Traffic)
Traffic_COAD <- Men2_extracted |> filter(Cancer == "COAD")

data_covariates_COAD <- data.frame(ID = COAD1$lab_id,
                                   PM_10_Annual_Average = COAD1$PM10_Annual_Average,
                                   PCR = COAD1$PCR,
                                   Age = COAD1$Age,
                                   BMI = COAD1$BMI,
                                   Alcohol = COAD1$Alcohol,
                                   Smoking = COAD1$Smoking,
                                   Workingstatus = COAD1$Workingstatus,
                                   Antidepressant = COAD1$Antidepressant,
                                   Homocysteine = COAD1$Homocysteine,
                                   Coffee = COAD1$Coffee,
                                   Residence_Traffic = Traffic_COAD$Residence_Traffic,
                                   Life_style = COAD1$Physical_Activity_New,
                                   Education = COAD1$Primary_School,
                                   Passive_Smoke_New = COAD1$Passive_Smoke_New)

dim(data_covariates_COAD)

datamirna_list_COAD <- list(COAD1,COAD2,COAD3,COAD4,COAD5,COAD6,COAD7,COAD8,COAD9,COAD10,
                            COAD11,COAD12,COAD13,COAD14,COAD15,COAD16,COAD17,COAD18, COAD19, COAD20)   
length(datamirna_list_COAD)


datafinal_COAD  <- data_covariates_COAD
for (k in 1:20) {
  tempdataset_COAD <- datamirna_list_COAD[[k]]
  response_dftemp_COAD <- data.frame(ID = tempdataset_COAD$lab_id,
                                     log2_rq = tempdataset_COAD$log2_rq )
  datafinal_COAD <- datafinal_COAD %>%
    left_join(response_dftemp_COAD, by = "ID")
  
}
colnames(datafinal_COAD)
colnames(datafinal_COAD)[16:35] <- paste(rep(c("hsa-miR-324-3p", "hsa-miR-671-3p", "hsa-miR-139-3p",
                                               
                                               "hsa-miR-328-3p", "hsa-miR-652-3p", "hsa-miR-193a-5p", "hsa-miR-574-3p",
                                               
                                               "hsa-miR-532-3p", "hsa-miR-628-5p","hsa-miR-136-5p", "hsa-miR-29a-5p",
                                               
                                               "hsa-miR-126-5p","hsa-miR-29b-3p",
                                               
                                               "hsa-miR-142-3p", "hsa-miR-125a-5p", "hsa-miR-197-3p","hsa-miR-378a-3p",
                                               
                                               "hsa-miR-21-5p", "hsa-miR-126-3p","hsa-miR-144-5p")),1:20,sep="")

colnames(datafinal_COAD)


#####
#LUAD
specified_mirna_LUAD<-subset( Mirna2, Cancer=="LUAD" & oncomirna %in% c("hsa-let-7a-5p", "hsa-miR-30d-5p", "hsa-miR-139-3p",
                                                                        
                                                                        "hsa-miR-133a-3p", "hsa-miR-30a-3p", "hsa-miR-140-3p", "hsa-let-7c-5p",
                                                                        
                                                                        "hsa-miR-184", "hsa-miR-139-5p","hsa-miR-338-5p", "hsa-miR-21-5p",
                                                                        
                                                                        "hsa-miR-708-5p","hsa-miR-9-5p","hsa-miR-212-3p",
                                                                        
                                                                        "hsa-miR-210-3p", "hsa-miR-7-1-3p", "hsa-miR-20a-3p","hsa-miR-301a-3p",
                                                                        
                                                                        "hsa-miR-128-3p", "hsa-miR-625-5p","hsa-miR-130b-3p",
                                                                        
                                                                        "hsa-miR-590-5p", "hsa-miR-151a-5p"))


#Subsetting columns of oncomirna
LUAD1 <- subset(specified_mirna_LUAD, oncomirna == "hsa-let-7a-5p",drop = FALSE)

LUAD2 <- subset(specified_mirna_LUAD, oncomirna == "hsa-let-7c-5p",drop = FALSE)

LUAD3 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-128-3p",drop = FALSE)

LUAD4 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-130b-3p",drop = FALSE)

LUAD5 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-133a-3p",drop = FALSE)

LUAD6 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-139-3p",drop = FALSE)

LUAD7 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-139-5p",drop = FALSE)

LUAD8 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-140-3p",drop = FALSE)

LUAD9 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-151a-5p",drop = FALSE)

LUAD10 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-184",drop = FALSE)

LUAD11 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-20a-3p",drop = FALSE)

LUAD12 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-210-3p",drop = FALSE)

LUAD13 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-212-3p",drop = FALSE)

LUAD14 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-21-5p",drop = FALSE)

LUAD15 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-301a-3p",drop = FALSE)

LUAD16 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-30a-3p",drop = FALSE)

LUAD17 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-30d-5p",drop = FALSE)

LUAD18 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-338-5p",drop = FALSE)

LUAD19 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-590-5p",drop = FALSE)

LUAD20 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-625-5p",drop = FALSE)

LUAD21 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-708-5p",drop = FALSE)

LUAD22 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-7-1-3p",drop = FALSE)

LUAD23 <- subset(specified_mirna_LUAD, oncomirna == "hsa-miR-9-5p",drop = FALSE)

Traffic_LUAD <- Men2_extracted |> filter(Cancer == "LUAD")


# Example data frame with IDs
data_covariates_LUAD <- data.frame(ID = LUAD1$lab_id,
                                   PM_10_Annual_Average = LUAD1$PM10_Annual_Average,
                                   PCR = LUAD1$PCR,
                                   Age = LUAD1$Age,
                                   BMI = LUAD1$BMI,
                                   Alcohol = LUAD1$Alcohol,
                                   Smoking = LUAD1$Smoking,
                                   Workingstatus = LUAD1$Workingstatus,
                                   Antidepressant = LUAD1$Antidepressant,
                                   Homocysteine = LUAD1$Homocysteine,
                                   Coffee = LUAD1$Coffee,
                                   Residence_Traffic = Traffic_LUAD$Residence_Traffic,
                                   Life_style = LUAD1$Physical_Activity_New,
                                   Education = LUAD1$Primary_School,
                                   Passive_Smoke_New = LUAD1$Passive_Smoke_New)

dim(data_covariates_LUAD)

datamirna_list_LUAD <- list(LUAD1,LUAD2,LUAD3,LUAD4,LUAD5,LUAD6,LUAD7,LUAD8,LUAD9,LUAD10,
                            LUAD11,LUAD12,LUAD13,LUAD14,LUAD15,LUAD16,LUAD17,LUAD18,LUAD19,LUAD20,LUAD21,LUAD22,LUAD23)   
length(datamirna_list_LUAD)


datafinal_LUAD  <- data_covariates_LUAD
for (k in 1:23) {
  tempdataset_LUAD <- datamirna_list_LUAD[[k]]
  response_dftemp_LUAD <- data.frame(ID = tempdataset_LUAD$lab_id,
                                     log2_rq = tempdataset_LUAD$log2_rq )
  datafinal_LUAD <- datafinal_LUAD %>%
    left_join(response_dftemp_LUAD, by = "ID")
  
}
colnames(datafinal_LUAD)
colnames(datafinal_LUAD)[16:38] <- paste(c("hsa-let-7a-5p", "hsa-miR-30d-5p", "hsa-miR-139-3p",
                                               
                                               "hsa-miR-133a-3p", "hsa-miR-30a-3p", "hsa-miR-140-3p", "hsa-let-7c-5p",
                                               
                                               "hsa-miR-184", "hsa-miR-139-5p","hsa-miR-338-5p", "hsa-miR-21-5p",
                                               
                                               "hsa-miR-708-5p","hsa-miR-9-5p","hsa-miR-212-3p",
                                               
                                               "hsa-miR-210-3p", "hsa-miR-7-1-3p", "hsa-miR-20a-3p","hsa-miR-301a-3p",
                                               
                                               "hsa-miR-128-3p", "hsa-miR-625-5p","hsa-miR-130b-3p",
                                               
                                               "hsa-miR-590-5p", "hsa-miR-151a-5p"))

colnames(datafinal_LUAD)

#####
#PRAD

specified_mirna_PRAD<-subset( Mirna2, Cancer=="PRAD" & oncomirna %in% c("hsa-miR-342-3p", "hsa-miR-143-3p", "hsa-miR-17-5p",
                                                                        "hsa-miR-708-5p", "hsa-miR-148a-3p", "hsa-miR-423-5p", "hsa-miR-7-1-3p", 
                                                                        "hsa-miR-20a-3p", "hsa-miR-769-5p","hsa-miR-142-3p", "hsa-miR-151a-5p", 
                                                                        "hsa-miR-21-5p","hsa-miR-92a-3p","hsa-miR-30d-5p",
                                                                        "hsa-miR-103a-3p", "hsa-miR-20b-5p", "hsa-miR-375","hsa-miR-126-3p", 
                                                                        "hsa-miR-93-3p", "hsa-miR-140-5p","hsa-miR-93-5p"))

PRAD1 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-103a-3p",drop = FALSE)

PRAD2 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-126-3p",drop = FALSE)

PRAD3 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-140-5p",drop = FALSE)

PRAD4 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-142-3p",drop = FALSE)

PRAD5 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-143-3p",drop = FALSE)

PRAD6 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-148a-3p",drop = FALSE)

PRAD7 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-151a-5p",drop = FALSE)

PRAD8 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-17-5p",drop = FALSE)

PRAD9 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-20a-3p",drop = FALSE)

PRAD10 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-20b-5p",drop = FALSE)

PRAD11 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-21-5p",drop = FALSE)

PRAD12 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-30d-5p",drop = FALSE)

PRAD13 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-342-3p",drop = FALSE)

PRAD14 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-375",drop = FALSE)

PRAD15 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-423-5p",drop = FALSE)

PRAD16 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-708-5p",drop = FALSE)

PRAD17 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-7-1-3p",drop = FALSE)

PRAD18 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-769-5p",drop = FALSE)

PRAD19 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-92a-3p",drop = FALSE)

PRAD20 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-93-3p",drop = FALSE)

PRAD21 <- subset(specified_mirna_PRAD, oncomirna == "hsa-miR-93-5p",drop = FALSE)

Traffic_PRAD <- Men2_extracted |> filter(Cancer == "PRAD")

data_covariates_PRAD <- data.frame(ID = PRAD1$lab_id,
                                   PM_10_Annual_Average = PRAD1$PM10_Annual_Average,
                                   PCR = PRAD1$PCR,
                                   Age = PRAD1$Age,
                                   BMI = PRAD1$BMI,
                                   Alcohol = PRAD1$Alcohol,
                                   Smoking = PRAD1$Smoking,
                                   Workingstatus = PRAD1$Workingstatus,
                                   Antidepressant = PRAD1$Antidepressant,
                                   Homocysteine = PRAD1$Homocysteine,
                                   Coffee = PRAD1$Coffee,
                                   Residence_Traffic = Traffic_PRAD$Residence_Traffic,
                                   Life_style = PRAD1$Physical_Activity_New,
                                   Education = PRAD1$Primary_School,
                                   Passive_Smoke_New = PRAD1$Passive_Smoke_New)

dim(data_covariates_PRAD)

datamirna_list_PRAD <- list(PRAD1,PRAD2,PRAD3,PRAD4,PRAD5,PRAD6,PRAD7,PRAD8,PRAD9,PRAD10,
                            PRAD11,PRAD12,PRAD13,PRAD14,PRAD15,PRAD16,PRAD17,
                            PRAD18,PRAD19, PRAD20, PRAD21)   
length(datamirna_list_PRAD)


datafinal_PRAD  <- data_covariates_PRAD
for (k in 1:21) {
  tempdataset_PRAD <- datamirna_list_PRAD[[k]]
  response_dftemp_PRAD <- data.frame(ID = tempdataset_PRAD$lab_id,
                                     log2_rq = tempdataset_PRAD$log2_rq )
  datafinal_PRAD <- datafinal_PRAD %>%
    left_join(response_dftemp_PRAD, by = "ID")
  
}
colnames(datafinal_PRAD)
colnames(datafinal_PRAD)[16:36] <- paste(c("hsa-miR-342-3p", "hsa-miR-143-3p", "hsa-miR-17-5p",
                                           "hsa-miR-708-5p", "hsa-miR-148a-3p", "hsa-miR-423-5p", "hsa-miR-7-1-3p", 
                                           "hsa-miR-20a-3p", "hsa-miR-769-5p","hsa-miR-142-3p", "hsa-miR-151a-5p", 
                                           "hsa-miR-21-5p","hsa-miR-92a-3p","hsa-miR-30d-5p",
                                           "hsa-miR-103a-3p", "hsa-miR-20b-5p", "hsa-miR-375","hsa-miR-126-3p", 
                                           "hsa-miR-93-3p", "hsa-miR-140-5p","hsa-miR-93-5p"))

colnames(datafinal_PRAD)

#####
#Women: LUAD
women_LUAD<-subset(women_data, Cancer=="LUAD" & oncomirna %in% c("hsa-miR-210-3p", "hsa-miR-151a-5p", "hsa-miR-21-5p",
                                                              
                                                              "hsa-miR-130b-3p", "hsa-miR-590-5p", "hsa-miR-625-5p", "hsa-miR-20a-3p",
                                                              
                                                              "hsa-miR-301a-3p", "hsa-miR-128-3p","hsa-miR-26b-3p", "hsa-miR-708-5p",
                                                              
                                                              "hsa-let-7a-5p","hsa-miR-133a-3p","hsa-miR-184",
                                                              
                                                              "hsa-miR-338-5p", "hsa-miR-30d-5p", "hsa-let-7c-5p","hsa-miR-139-5p",
                                                              
                                                              "hsa-miR-7-1-3p", "hsa-miR-212-3p","hsa-miR-9-5p","hsa-miR-139-3p",
                                                              
                                                              "hsa-miR-30a-3p", "hsa-miR-140-3p"))

W_LUAD1 <- subset(women_LUAD, oncomirna == "hsa-let-7a-5p",drop = FALSE)

W_LUAD2 <- subset(women_LUAD, oncomirna == "hsa-let-7c-5p",drop = FALSE)

W_LUAD3 <- subset(women_LUAD, oncomirna == "hsa-miR-128-3p",drop = FALSE)

W_LUAD4 <- subset(women_LUAD, oncomirna == "hsa-miR-130b-3p",drop = FALSE)

W_LUAD5 <- subset(women_LUAD, oncomirna == "hsa-miR-133a-3p",drop = FALSE)

W_LUAD6 <- subset(women_LUAD, oncomirna == "hsa-miR-139-3p",drop = FALSE)

W_LUAD7 <- subset(women_LUAD, oncomirna == "hsa-miR-139-5p",drop = FALSE)

W_LUAD8 <- subset(women_LUAD, oncomirna == "hsa-miR-140-3p",drop = FALSE)

W_LUAD9 <- subset(women_LUAD, oncomirna == "hsa-miR-151a-5p",drop = FALSE)

W_LUAD10 <- subset(women_LUAD, oncomirna == "hsa-miR-184",drop = FALSE)

W_LUAD11 <- subset(women_LUAD, oncomirna == "hsa-miR-20a-3p",drop = FALSE)

W_LUAD12 <- subset(women_LUAD, oncomirna == "hsa-miR-210-3p",drop = FALSE)

W_LUAD13 <- subset(women_LUAD, oncomirna == "hsa-miR-212-3p",drop = FALSE)

W_LUAD14 <- subset(women_LUAD, oncomirna == "hsa-miR-21-5p",drop = FALSE)

W_LUAD15 <- subset(women_LUAD, oncomirna == "hsa-miR-26b-3p",drop = FALSE)

W_LUAD16 <- subset(women_LUAD, oncomirna == "hsa-miR-301a-3p",drop = FALSE)

W_LUAD17 <- subset(women_LUAD, oncomirna == "hsa-miR-30a-3p",drop = FALSE)

W_LUAD18 <- subset(women_LUAD, oncomirna == "hsa-miR-30d-5p",drop = FALSE)

W_LUAD19 <- subset(women_LUAD, oncomirna == "hsa-miR-338-5p",drop = FALSE)

W_LUAD20 <- subset(women_LUAD, oncomirna == "hsa-miR-590-5p",drop = FALSE)

W_LUAD21 <- subset(women_LUAD, oncomirna == "hsa-miR-625-5p",drop = FALSE)

W_LUAD22 <- subset(women_LUAD, oncomirna == "hsa-miR-708-5p",drop = FALSE)

W_LUAD23 <- subset(women_LUAD, oncomirna == "hsa-miR-7-1-3p",drop = FALSE)

W_LUAD24 <- subset(women_LUAD, oncomirna == "hsa-miR-9-5p",drop = FALSE)

Women2_extracted <- women_data %>%
  distinct(lab_id, Cancer, .keep_all = TRUE) %>%
  select(lab_id, Cancer, Residence_Traffic)

Traffic_LUAD_women <- Women2_extracted |> filter(Cancer == "LUAD")

data_covariates_W_LUAD <- data.frame(ID = W_LUAD1$lab_id,
                                     PM_10_Annual_Average = W_LUAD1$PM10_Annual_Average,
                                     PCR = W_LUAD1$PCR,
                                     Age = W_LUAD1$age,
                                     BMI = W_LUAD1$bmi,
                                     Alcohol = W_LUAD1$Alcohol,
                                     Workingstatus = W_LUAD1$Working_Status,
                                     Antidepressant = W_LUAD1$Antidepressent,
                                     Homocysteine = W_LUAD1$Homocysteine,
                                     Menopause = W_LUAD1$Menopause,
                                     No_Pregnancy = W_LUAD1$No_Pregency,
                                     Oral_Contraceptives = W_LUAD1$Oral_Contraceptives,
                                     Smoking = W_LUAD1$Smoking_New,
                                     Coffee = W_LUAD1$Coffee,
                                     Residence_Traffic = Traffic_LUAD_women$Residence_Traffic,
                                     Life_style = W_LUAD1$Life_style,
                                     Primary_School = W_LUAD1$Prima,
                                     Passive_Smoke_New = W_LUAD1$Passive_Smoke_New)

dim(data_covariates_W_LUAD)
length(unique(data_covariates_W_LUAD$ID))

datamirna_list_W_LUAD <- list(W_LUAD1,W_LUAD2,W_LUAD3,W_LUAD4,W_LUAD5,W_LUAD6,W_LUAD7,W_LUAD8,W_LUAD9,W_LUAD10,
                              W_LUAD11,W_LUAD12,W_LUAD13,W_LUAD14,W_LUAD15,W_LUAD16,W_LUAD17,W_LUAD18,W_LUAD19,
                              W_LUAD20,W_LUAD21,W_LUAD22,W_LUAD23, W_LUAD24)   
length(datamirna_list_W_LUAD)


datafinal_W_LUAD  <- data_covariates_W_LUAD
for (k in 1:24) {
  tempdataset_W_LUAD <- datamirna_list_W_LUAD[[k]]
  response_dftemp_W_LUAD <- data.frame(ID = tempdataset_W_LUAD$lab_id,
                                       log2_rq = tempdataset_W_LUAD$log2_rq )
  datafinal_W_LUAD <- datafinal_W_LUAD %>%
    left_join(response_dftemp_W_LUAD, by = "ID")
  
}
colnames(datafinal_W_LUAD)
colnames(datafinal_W_LUAD)[21:44] <- paste(rep("log2_rq",24),1:24,sep="")

colnames(datafinal_W_LUAD)





#####
#FAMD

str(datafinal_COAD)
datafinal_COAD$Smoking <- factor(datafinal_COAD$Smoking)
datafinal_COAD$Workingstatus <- factor(datafinal_COAD$Workingstatus)
datafinal_COAD$Antidepressant <- factor(datafinal_COAD$Antidepressant)
datafinal_COAD[,13:15] <- lapply(datafinal_COAD[,13:15], factor)

tmp <- na.omit(datafinal_COAD)
tmp2 <- datafinal_COAD[,-1]
tmp <- tmp |> select(-ID)

prova1 <- FAMD(tmp, sup.var = 15:34)
FAMD(tmp2)


library(readr)
library(cluster)
library(Rtsne)
data <- read_delim("data/data_MANCOAD.csv", 
                   delim = ";", escape_double = FALSE, trim_ws = TRUE)

data$Antidepressant <- ifelse(data$Antidepressant == "si","Antidepressant_sì","Antidepressant_no")
tmp3 <- rep("",nrow(data))

tmp3[data$Smoking == "si"] <- "Fumatore"
tmp3[data$Smoking == "ex"] <- "ex_Fumatore"
tmp3[data$Smoking == "no"] <- "Non_Fumatore"

data$Smoking <- factor(tmp3, levels = c("Non_Fumatore","ex_Fumatore","Fumatore"))
data$Passive_Smoke <- ifelse(data$Passive_Smoke == "si","Fumo_passivo_sì","Fumo_passivo_no")
tmp2 <- na.omit(data)
tmp2 <- tmp2 |> select(-lab_id)
prova1 <- FAMD(tmp2, sup.var = 1:20, axes = c(2,3))

fviz_famd_var(prova1, geom = c("point"),
              alpha.var = 0.5,
              axes = c(1,2))

tmp2[,c(24,25,29,30,31,32)] <- lapply(tmp2[,c(24,25,29,30,31,32)], factor)
gower_dist <- daisy(tmp2,
                    metric = "gower",
                    type = list(logratio = c(21:23,26:27,33:34)))

gower_mat <- as.matrix(gower_dist)
tmp2[
  which(gower_mat == min(gower_mat[gower_mat != min(gower_mat)]),
        arr.ind = TRUE)[1, ], ]

sil_width <- c(NA)
for(i in 2:10){
  
  pam_fit <- pam(gower_dist,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}

plot(1:10, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:10, sil_width)



pam_fit <- pam(gower_dist, diss = TRUE, k = 3)

pam_results <- tmp2 %>%
  mutate(cluster = pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))

pam_results$the_summary

tsne_obj <- Rtsne(gower_dist, is_distance = TRUE)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering))

ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster))


tmp3 <- cbind(tmp2,pam_fit$clustering)
tmp3 <- tmp3 |> rename(Cluster = `pam_fit$clustering`)

clust1 <- tmp3 |> filter(Cluster == 1)
clust2 <- tmp3 |> filter(Cluster == 2)
clust3 <- tmp3 |> filter(Cluster == 3)

clust1_miRNA <- clust1 |> select(starts_with("hsa-"), Traffic_no)
clust1_miRNA_long <- pivot_longer(clust1_miRNA, cols = 1:20, names_to = "miRNA")
clust1_miRNA_long$Traffic_no <- factor(clust1_miRNA_long$Traffic_no, levels = c("mild","moderate","heavy"))

ggplot(clust1_miRNA_long, aes(x = miRNA, y = value)) +
  geom_boxplot() +
  facet_wrap(~Traffic_no)

clust1_miRNA_long$group <- interaction(clust1_miRNA_long$miRNA, clust1_miRNA_long$Traffic_no)

ggplot(clust1_miRNA_long, aes(x = miRNA, y = value, fill = Traffic_no)) +
  geom_boxplot(position = "dodge") +
  xlab("miRNA") +
  ylab("Expression Value") +
  ggtitle("Boxplot of miRNA Expression by Traffic_no") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

