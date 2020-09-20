#Library initiation 
library("readxl")
library("xlsx")
library("Rmpfr")

#Load the main excel files - stble1 (master table), stable 2(track reference) and stable 3 (ENIGMA+ClinVar reference variants) 
stable1 = as.data.frame(read_excel("Sup_tables_GIM_AUG_2020 REV FINAL.xlsx", sheet = "STable 1", col_names = TRUE, na="", guess_max = 15000))
stable2 = as.data.frame(read_excel("Sup_tables_GIM_AUG_2020 REV FINAL.xlsx", sheet = "STable 2", col_names = TRUE, na=""))
stable3 = as.data.frame(read_excel("Sup_tables_GIM_AUG_2020 REV FINAL.xlsx", sheet = "STable 3", col_names = TRUE, na=""))

#Table with sum of numbers from each variant in master table (row by row)
stable1_counts <- data.frame()[1:nrow(stable1),]
stable1_counts$"no_functional_impact" <- rowSums(stable1[,10:ncol(stable1)] == 0, na.rm= TRUE)  
stable1_counts$"functional_impact" <- rowSums(stable1[,10:ncol(stable1)] == 1, na.rm= TRUE)
stable1_counts$"number of tests" <- stable1_counts$"no_functional_impact"  + stable1_counts$"functional_impact"

#Preparing reference only panel table [E+C] 
reference_panel_E_C <- stable1
reference_panel_E_C$"T7"[!is.na(reference_panel_E_C$"T7")] <- NA
for (i in 1:nrow(reference_panel_E_C)){
  var_aa_st1 <- reference_panel_E_C[i,5]
  for (z in 1:nrow(stable3)-8){
    var_aa_st3 <- stable3[z,3]
    if (var_aa_st1 == toString(var_aa_st3)){
      reference_panel_E_C[i,7] <- stable3[z,8] 
    }
  }
}
reference_panel_E_C <- as.data.frame(reference_panel_E_C[!is.na(reference_panel_E_C$"T7"), ])

#Table with sum of variants from each track in reference panel table (column by column)
count_reference_panel <- function (reference_panel, master_table){
  x <- matrix(c(colSums(reference_panel[,10:length(reference_panel)] == 0, na.rm= TRUE),
                colSums(reference_panel[,10:length(reference_panel)] == 1, na.rm= TRUE),
                colSums(master_table[,10:length(master_table)] <= 1, na.rm= TRUE),
                colSums(reference_panel[,10:length(reference_panel)] <= 1, na.rm= TRUE)), ncol=4)
  rownames(x) <- c("T10",	"T11",	"T12",	"T13",	"T14",	"T15",	"T16",	"T17",	"T18",	
                                "T19",	"T20",	"T21",	"T22",	"T23",	"T24",	"T25",	"T26",	"T27",	
                                "T28",	"T29",	"T30",	"T31",	"T32",	"T33",	"T34",	"T35",	"T36",	
                                "T37",	"T38",	"T39",	"T40",	"T41",	"T42",	"T43",	"T44",	"T45",	
                                "T46",	"T47",	"T48",	"T49",	"T50",	"T51",	"T52",	"T53",	"T54",	
                                "T55",	"T56",	"T57",	"T58",	"T59",	"T60",	"T61",	"T62",	"T63",	
                                "T64",	"T65",	"T66",	"T67",	"T68",	"T69",	"T70",	"T71",	"T72",	
                                "T73",	"T74",	"T75",	"T76",	"T77",	"T78",	"T79",	"T80",	"T81",
                                "T82",	"T83",	"T84",	"T85",	"T86",	"T87",	"T88",	"T89",	"T90",	
                                "T91",	"T92",	"T93",	"T94",	"T95",	"T96",	"T97",	"T98",	"T99",	
                                "T100",	"T101",	"T102",	"T103",	"T104",	"T105",	"T106",	"T107",	"T108",	
                                "T109",	"T110",	"T111",	"T112",	"T113",	"T114",	"T115",	"T116",	"T117",	
                                "T118",	"T119",	"T120",	"T121",	"T122",	"T123",	"T124",	"T125",	"T126",	
                                "T127",	"T128",	"T129",	"T130",	"T131",	"T132",	"T133",	"T134",	"T135",	
                                "T136",	"T137",	"T138",	"T139",	"T140")
  colnames(x) <- c("Assay_no_functional_impact", "Assay_functional_impact","#_of_missense_variants_tested", "Total_classified_variants_tested")
  
  return(x)
}
                                       
#SUPPLEMENTARY TABLE 5: SUMMARY STATISTICS FOR TESTED VARIANTS (NUMBER OF ASSAYS PERFORMED FOR EACH VARIANT)
#Missense variants - Calculation for ALL FUNCTIONAL TRACKS of number of times variants were tested
count_variant_tested <- function (master_table, counts){
  x<-matrix(c("N", "Percentage ",
                            nrow(master_table),(nrow(master_table)*100)/nrow(master_table),
                            sum(counts$"number of tests" == 0),round((sum(counts$"number of tests" == 0)*100)/nrow(master_table),digits = 2),
                            sum(counts$"number of tests" >= 1),round((sum(counts$"number of tests" >= 1)*100)/nrow(master_table),digits =2),
                            sum(counts$"number of tests" >= 2),round((sum(counts$"number of tests" >= 2)*100)/nrow(master_table),digits = 2),
                            sum(counts$"number of tests" >= 10),round((sum(counts$"number of tests" >= 10)*100)/nrow(master_table),digits = 2),
                            sum(counts$"number of tests" >= 15),round((sum(counts$"number of tests" >= 15)*100)/nrow(master_table),digits = 2),
                            sum(counts$"number of tests" >= 20),round((sum(counts$"number of tests" >= 20)*100)/nrow(master_table),digits = 2),
                            sum(counts$"number of tests" >= 25),round((sum(counts$"number of tests" >= 25)*100)/nrow(master_table),digits = 2),
                            sum(counts$"number of tests" >= 30),round((sum(counts$"number of tests" >= 30)*100)/nrow(master_table),digits = 2),
                            sum(counts$"number of tests" >= 35),round((sum(counts$"number of tests" >= 35)*100)/nrow(master_table),digits = 2),
                            sum(counts$"number of tests" >= 40),round((sum(counts$"number of tests" >= 40)*100)/nrow(master_table),digits = 2),
                            sum(counts$"number of tests" >= 50),round((sum(counts$"number of tests" >= 50)*100)/nrow(master_table),digits = 2),
                            sum(counts$"number of tests" >= 60),round((sum(counts$"number of tests" >= 60)*100)/nrow(master_table),digits = 2)),
                          ncol = 2, byrow=TRUE)
  
  rownames(x)<-c("Missense Variants","Possible unique (BRCA1)","not yet tested","tested >= 1 time","tested >= 2 times","tested >= 10 times",
                               "tested >= 15 times","tested >= 20 times","tested >= 25 times","tested >= 30 times","tested >= 35 times","tested >= 40 times",
                               "tested >= 50 times","tested >= 60 times")
  colnames(x)<-c("All_Functional_Tracks", "All_Functional_Tracks(Relative")
  return(x)
}

#SUPPLEMENTARY TABLE 5: Calculation for Documented variants (BRCAExchange)
count_BRCAExchange <- function(master_table, counts, reference_panel){
  x <- matrix(c(sum(master_table$"T9" == 1),(sum(master_table$"T9" == 1)*100)/sum(master_table$"T9" == 1), 
                sum((master_table$"T9" == 1) & (counts$"number of tests" >= 1)),
           round((sum((master_table$"T9" == 1) & (counts$"number of tests" >= 1))*100)/sum(master_table$"T9" == 1),2)),
           ncol=2, byrow = TRUE)
  #Percentage  of documented VUS 
  documented_vus <- matrix(c("-", round((100*(x[2,1])/(x[1,1]-nrow(reference_panel))),2)))
  x <- cbind(x, documented_vus[c(1:2),])
  rownames(x)<- c("Total documented missense variants","Total documented missense variants tested") 
  colnames(x)<-c("N", "Percentage", "Percentage_VUS" )
  return(x)
} 

#SUPPLEMENTARY TABLE 6: ASSAY THROUGHPUT
#Number of variants tested by track
count_variant_track <- function(reference_table_counts){
  Tracks <- matrix(c(nrow(reference_table_counts),round(nrow(reference_table_counts)*100/nrow(reference_table_counts),2), 
                     sum(reference_table_counts[,3] == 1), round(sum(reference_table_counts[,3] == 1)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 5),round(sum(reference_table_counts[,3] >= 5)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 10),round(sum(reference_table_counts[,3] >= 10)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 20),round(sum(reference_table_counts[,3] >= 20)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 30),round(sum(reference_table_counts[,3] >= 30)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 40),round(sum(reference_table_counts[,3] >= 40)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 50),round(sum(reference_table_counts[,3] >= 50)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 60),round(sum(reference_table_counts[,3] >= 60)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 70),round(sum(reference_table_counts[,3] >= 70)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 80),round(sum(reference_table_counts[,3] >= 80)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 90),round(sum(reference_table_counts[,3] >= 90)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 100),round(sum(reference_table_counts[,3] >= 100)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 200),round(sum(reference_table_counts[,3] >= 200)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 300),round(sum(reference_table_counts[,3] >= 300)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 400),round(sum(reference_table_counts[,3] >= 400)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 500),round(sum(reference_table_counts[,3] >= 500)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 1000),round(sum(reference_table_counts[,3] >= 1000)*100/nrow(reference_table_counts),2),
                     sum(reference_table_counts[,3] >= 1500),round(sum(reference_table_counts[,3] >= 1500)*100/nrow(reference_table_counts),2)),
                   ncol=2,byrow = TRUE)
  rownames(Tracks) <- c("Total number of functional tracks","Tracks testing 1 variant only","Tracks testing >= 5 variants",
                        "Tracks testing >= 10 variants", "Tracks testing >= 20 variants","Tracks testing >= 30 variants",
                        "Tracks testing >= 40 variants","Tracks testing >= 50 variants","Tracks testing >= 60 variants",
                        "Tracks testing >= 70 variants","Tracks testing >= 80 variants","Tracks testing >= 90 variants",
                        "Tracks testing >= 100 variants","Tracks testing >= 200 variants","Tracks testing >= 300 variants",
                        "Tracks testing >= 400 variants","Tracks testing >= 500 variants","Tracks testing >= 1000 variants",
                        "Tracks testing >= 1500 variants")
  colnames(Tracks) <- c("N","Percentage ")
  return(Tracks)
}

#SUPPLEMENTARY TABLE 9: SUMMARY STATISTICS AND SPECIFICITY AND SENSITIVITY FOR EACH ASSAY FULLFILLING BOTH CRITERIA  ( #VARIANTS TESTED > 10 AND # VARIANTS IN REFERENCE PANEL [NON-PATHOGENIC; PATHOGENIC] [3:3])***********************************************************************************
#Calculating and creating data matrix
espec_sens_calculation <- function(reference_panel){
  Espec_Sensi <- as.data.frame(matrix(c(colSums((reference_panel[,10:ncol(reference_panel)] == 0) & (reference_panel$"T7" > 3),na.rm=TRUE),
                                        colSums((reference_panel[,10:ncol(reference_panel)] == 1) & (reference_panel$"T7" > 3),na.rm=TRUE),
                                        colSums((reference_panel[,10:ncol(reference_panel)] == 0) & (reference_panel$"T7" < 3),na.rm=TRUE),
                                        colSums((reference_panel[,10:ncol(reference_panel)] == 1) & (reference_panel$"T7" < 3),na.rm=TRUE)),
                                      ncol = 4))
  rownames(Espec_Sensi) <- c("T10",	"T11",	"T12",	"T13",	"T14",	"T15",	"T16",	"T17",	"T18",	
                             "T19",	"T20",	"T21",	"T22",	"T23",	"T24",	"T25",	"T26",	"T27",	
                             "T28",	"T29",	"T30",	"T31",	"T32",	"T33",	"T34",	"T35",	"T36",	
                             "T37",	"T38",	"T39",	"T40",	"T41",	"T42",	"T43",	"T44",	"T45",	
                             "T46",	"T47",	"T48",	"T49",	"T50",	"T51",	"T52",	"T53",	"T54",	
                             "T55",	"T56",	"T57",	"T58",	"T59",	"T60",	"T61",	"T62",	"T63",	
                             "T64",	"T65",	"T66",	"T67",	"T68",	"T69",	"T70",	"T71",	"T72",	
                             "T73",	"T74",	"T75",	"T76",	"T77",	"T78",	"T79",	"T80",	"T81",
                             "T82",	"T83",	"T84",	"T85",	"T86",	"T87",	"T88",	"T89",	"T90",	
                             "T91",	"T92",	"T93",	"T94",	"T95",	"T96",	"T97",	"T98",	"T99",	
                             "T100",	"T101",	"T102",	"T103",	"T104",	"T105",	"T106",	"T107",	"T108",	
                             "T109",	"T110",	"T111",	"T112",	"T113",	"T114",	"T115",	"T116",	"T117",	
                             "T118",	"T119",	"T120",	"T121",	"T122",	"T123",	"T124",	"T125",	"T126",	
                             "T127",	"T128",	"T129",	"T130",	"T131",	"T132",	"T133",	"T134",	"T135",	
                             "T136",	"T137",	"T138",	"T139",	"T140")
  colnames(Espec_Sensi) <- c("classified_false_no_impact", "classified_true_impact","classified_true_no_impact", "classified_false_impact")
  Espec_Sensi$"number of reference variants as bening" <- rowSums(Espec_Sensi[,1:2],na.rm=TRUE)
  Espec_Sensi$"number of reference variants as pathogenic" <- rowSums(Espec_Sensi[,3:4],na.rm=TRUE)
  Espec_Sensi$"total" <- rowSums(Espec_Sensi[,c(5,6)],na.rm=TRUE)
  
  #Sensitivity and Specificity calculation and CI bounds
  Espec_Sensi$"Sensitivity" <- round(Espec_Sensi[1:131,2] / Espec_Sensi[1:131,5], digits = 2)
  Sens_l95b <- 0
  Sens_u95b <- 0
  z <- 1*1.95996
  zsq <- z*z
  for(i in 1:nrow(Espec_Sensi)) {
    #SENSITIVITY CI
    n <- Espec_Sensi[i,5]
    p <- Espec_Sensi[i,8]
    q <- 1-p
    #lower bound
    num <- (2*n*p)+zsq-1-(z*sqrt(zsq-2-(1/n)+4*p*((n*q)+1)))
    denom <- 2*(n+zsq)
    l95b <- num/denom
    Sens_l95b[i] <- l95b
    #upper bound
    num <- (2*n*p)+zsq+1+(z*sqrt(zsq+2-(1/n)+4*p*((n*q)-1)))
    denom <- 2*(n+zsq)
    u95b <- num/denom
    Sens_u95b[i] <- u95b
  }
  Espec_Sensi <- as.data.frame(cbind(Espec_Sensi,(round(Sens_l95b, digits =2))))
  Espec_Sensi <- as.data.frame(cbind(Espec_Sensi,(round(Sens_u95b,digits=2))))
  
  Espec_Sensi$"specificity" <- round(Espec_Sensi[1:131,3] / Espec_Sensi[1:131,6], digits = 2)
  Espec_l95b <- 0
  Espec_u95b <- 0
  for(i in 1:nrow(Espec_Sensi)) {
    #SPECIFICITY CI
    n <- Espec_Sensi[i,6]
    p <- Espec_Sensi[i,11]
    q <- 1-p
    #lower bound
    num <- (2*n*p)+zsq-1-(z*sqrt(zsq-2-(1/n)+4*p*((n*q)+1)))
    denom <- 2*(n+zsq)
    l95b <- num/denom
    Espec_l95b[i] <- l95b
    #upper bound
    num <- (2*n*p)+zsq+1+(z*sqrt(zsq+2-(1/n)+4*p*((n*q)-1)))
    denom <- 2*(n+zsq)
    u95b <- num/denom
    Espec_u95b[i] <- u95b
  }
  Espec_Sensi <- as.data.frame(cbind(Espec_Sensi,(round(Espec_l95b, digits =2))))
  Espec_Sensi <- as.data.frame(cbind(Espec_Sensi,(round(Espec_u95b,digits=2))))
  
  #Changing columns names
  colnames(Espec_Sensi)[9] <- "95Percentage  CI lower"
  colnames(Espec_Sensi)[10] <- "95Percentage  CI upper"
  colnames(Espec_Sensi)[12] <- "95Percentage  CI lower"
  colnames(Espec_Sensi)[13] <- "95Percentage  CI upper"
  
  #Prior, P1 calculation
  Espec_Sensi$"Proportion pathogenic (Prior, P1)" <- round(Espec_Sensi[,5]/Espec_Sensi[,7], digits=2)
  #"+1 proportions" (Posterior, P2)
  Espec_Sensi$"+1 proportions (Posterior, P2) - Pathogenic" <- round((Espec_Sensi[,2]+Espec_Sensi[,4])/((Espec_Sensi[,2]+Espec_Sensi[,4])+1), digits=3)
  Espec_Sensi$"+1 proportions (Posterior, P2) - Benign" <- round(1/((Espec_Sensi[,1]+Espec_Sensi[,3])+1), digits=3)
  #OddsPath calculation
  Espec_Sensi$"Odds path - Pathogenic" <- round((Espec_Sensi[,15]*(1-Espec_Sensi[,14]))/((1-Espec_Sensi[,15])*Espec_Sensi[,14]), digits =3)
  Espec_Sensi$"Odds path - Benign" <- round((Espec_Sensi[,16]*(1-Espec_Sensi[,14]))/((1-Espec_Sensi[,16])*Espec_Sensi[,14]), digits =3)
  
  #Evidence strenght equivalence 
  Espec_Sensi$"Evidence_Path" <- ifelse(Espec_Sensi$"Odds path - Pathogenic" > 350 & Espec_Sensi$"Odds path - Pathogenic" < 800, "PS3_very_strong", 
                                 ifelse(Espec_Sensi$"Odds path - Pathogenic" > 18.7 & Espec_Sensi$"Odds path - Pathogenic" <= 350, "PS3",
                                 ifelse(Espec_Sensi$"Odds path - Pathogenic" > 4.3 & Espec_Sensi$"Odds path - Pathogenic" <= 18.7,"PS3_moderate",
                                 ifelse(Espec_Sensi$"Odds path - Pathogenic" > 2.1 & Espec_Sensi$"Odds path - Pathogenic" <= 4.3,"PS3_supporting", "indeterminate"))))
  
  Espec_Sensi$"Evidence_Benign" <- ifelse(Espec_Sensi$"Odds path - Benign" < 0.053 & Espec_Sensi$"Odds path - Benign" > 0, "BS3", 
                                   ifelse(Espec_Sensi$"Odds path - Benign" < 0.23 & Espec_Sensi$"Odds path - Benign" >=0.053, "BS3_moderate",
                                   ifelse(Espec_Sensi$"Odds path - Benign" < 0.48 & Espec_Sensi$"Odds path - Benign" >=0.23,"BS3_supporting","Indeterminate")))
 return(Espec_Sensi) 
}

#SUPPLEMENTARY TABLE 6: CRITERIAS
#Tracks with number of variants in reference panel non-pathogenic; pathogenic  X:Y
ratio_criteria <- function(reference_panel_counts, reference_panel){
  Ratios <- matrix(c(sum((reference_panel_counts[,5] >= 3) & (reference_panel_counts[,6] >= 3)),round(sum((reference_panel_counts[,5] >= 3) & (reference_panel_counts[,6] >= 3))*100/nrow(reference_panel),2),
                     sum((reference_panel_counts[,5] >= 4) & (reference_panel_counts[,6] >= 4)),round(sum((reference_panel_counts[,5] >= 4) & (reference_panel_counts[,6] >= 4))*100/nrow(reference_panel),2),
                     sum((reference_panel_counts[,5] >= 5) & (reference_panel_counts[,6] >= 5)),round(sum((reference_panel_counts[,5] >= 5) & (reference_panel_counts[,6] >= 5))*100/nrow(reference_panel),2),
                     sum((reference_panel_counts[,5] >= 10) & (reference_panel_counts[,6] >= 10)),round(sum((reference_panel_counts[,5] >= 10) & (reference_panel_counts[,6] >= 10))*100/nrow(reference_panel),2)),
                   ncol=2, byrow = TRUE)
  rownames(Ratios) <- c("Tracks with # variants in ENIGMA reference panel [non-pathogenic; pathogenic] >= [3:3]",
                        "Tracks with # variants in ENIGMA reference panel [non-pathogenic; pathogenic] >= [4:4]",
                        "Tracks with # variants in ENIGMA reference panel [non-pathogenic; pathogenic] >= [5:5]",
                        "Tracks with # variants in ENIGMA reference panel [non-pathogenic; pathogenic] >= [10:10]")
  colnames(Ratios) <- c("N","Percentage")
  return(Ratios)
}
#Tracks meeting the two criteria ( number of variants tested >= 10 and number of variants in reference panel equal or greater non-pathogenic; pathogenic x:y)
ratio_criteria_2 <- function(reference_panel_counts, reference_panel){
  Criteria <- matrix(c(sum((reference_panel[,3] >= 10) & (reference_panel_counts[,5] >= 3) & (reference_panel_counts[,6] >= 3)),
                       round(sum((reference_panel[,3] >= 10) & (reference_panel_counts[,5] >= 3) & (reference_panel_counts[,6] >= 3))*100/nrow(reference_panel),2),
                       sum((reference_panel[,3] >= 10) & (reference_panel_counts[,5] >= 4) & (reference_panel_counts[,6] >= 4)),
                       round(sum((reference_panel[,3] >= 10) & (reference_panel_counts[,5] >= 4) & (reference_panel_counts[,6] >= 4))*100/nrow(reference_panel),2),
                       sum((reference_panel[,3] >= 10) & (reference_panel_counts[,5] >= 5) & (reference_panel_counts[,6] >= 5)),
                       round(sum((reference_panel[,3] >= 10) & (reference_panel_counts[,5] >= 5) & (reference_panel_counts[,6] >= 5))*100/nrow(reference_panel),2),
                       sum((reference_panel[,3] >= 10) & (reference_panel_counts[,5] >= 10) & (reference_panel_counts[,6] >= 10)),
                       round(sum((reference_panel[,3] >= 10) & (reference_panel_counts[,5] >= 10) & (reference_panel_counts[,6] >= 10))*100/nrow(reference_panel),2)),
                     ncol=2, byrow=TRUE)
  rownames(Criteria) <- c("Tracks meeting the two criteria ( #variants tested >= 10 and # variants in ENIGMA reference panel >= [non-pathogenic; pathogenic] [3:3])",
                          "Tracks meeting the two criteria ( #variants tested >= 10 and # variants in ENIGMA reference panel >= [non-pathogenic; pathogenic] [4:4])",
                          "Tracks meeting the two criteria ( #variants tested >= 10 and # variants in ENIGMA reference panel >= [non-pathogenic; pathogenic] [5:5])",
                          "Tracks meeting the two criteria ( #variants tested >= 10 and # variants in ENIGMA reference panel >= [non-pathogenic; pathogenic] [10:10])")
  colnames(Criteria) <- c("N","Percentage")
  return(Criteria)
}

#Creating the HI-SET Table
hi_set_table <- function(Espec_Sensi, reference_panel){
  hi_set_table <- data.frame()[1:11009,]
  for (i in 1:nrow(Espec_Sensi)) {
    if ((Espec_Sensi[i,5])>= 4 & (Espec_Sensi[i,6])>= 4 & (reference_panel[i,3]) >= 10 & (Espec_Sensi[i,8]) >= 0.8 & (Espec_Sensi[i,11]) >= 0.8){
      hi_set_table <- as.data.frame(rbind(hi_set_table, Espec_Sensi[i,]))
    }
  }
  return(hi_set_table)
}

#Creating the HI-SET Tracks
hi_set_tracks <- function(hi_set_table, master_table){
  hi_set_tracks <- data.frame()[1:nrow(master_table),]
  tracks <- list()
  for (i in rownames(hi_set_table)){
    for (x in colnames(master_table)){
      if (i == x){
        hi_set_tracks <- as.data.frame(cbind(hi_set_tracks, master_table[,x]))
        q <- colnames(master_table[x])
        tracks <- c(tracks, q)
      }
    }
  }
  colnames(hi_set_tracks) <- tracks 
  row.names(hi_set_tracks) <- 1:11009
  hi_set_tracks <- data.frame(cbind(master_table[,c(5,7,9)],hi_set_tracks))
  return(hi_set_tracks)
}

#SUPPLEMENTARY TABLE 7: ASSAY CLASSES
assay_classes <- function(assay_classes, df_variants_reference_panel, Espec_Sensi){
  assay_classe <- matrix(c(sum(assay_classes[,5] == 2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==2)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==2))/sum(assay_classes[,5] == 2))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==2) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==2) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,5] == 2))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==2) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==2)))*100,2),
                           
                           sum(assay_classes[,5] == 3),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==3)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==3))/sum(assay_classes[,5] == 3))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==3) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==3) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,5] == 3))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==3) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==3)))*100,2),
                           
                           sum(assay_classes[,5] == 4),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==4)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==4))/sum(assay_classes[,5] == 4))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==4) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==4) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,5] == 4))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==4) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==4)))*100,2),
                           
                           sum(assay_classes[,5] == 5),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==5)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==5))/sum(assay_classes[,5] == 5))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==5) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==5) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,5] == 5))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==5) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==5)))*100,2),
                           
                           sum(assay_classes[,5] == 6),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==6)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==6))/sum(assay_classes[,5] == 6))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==6) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==6) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,5] == 6))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==6) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==6)))*100,2),
                           
                           sum(assay_classes[,5] == 7),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==7)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==7))/sum(assay_classes[,5] == 7))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==7) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==7) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,5] == 7))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==7) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==7)))*100,2),
                           
                           sum(assay_classes[,5] == 8),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==8))/sum(assay_classes[,5] == 8))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==8) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==8) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,5] == 8))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==8) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==8)))*100,2),
                           
                           sum(assay_classes[,5] == 9),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==9)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==9))/sum(assay_classes[,5] == 9))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==9) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==9) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,5] == 9))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==9) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==9)))*100,2),
                           
                           sum(assay_classes[,5] == 10),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==10)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==10))/sum(assay_classes[,5] == 10))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==10) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==10) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,5] == 10))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==10) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==10)))*100,2),
                           
                           sum(assay_classes[,5] == 11),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==11)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==11))/sum(assay_classes[,5] == 11))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==11) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==11) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,5] == 11))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==11) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==11)))*100,2),
                           
                           sum(assay_classes[,5] == 12),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==12)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==12))/sum(assay_classes[,5] == 12))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==12) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==12) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,5] == 12))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==12) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==12)))*100,2),
                           
                           sum(assay_classes[,5] == 13),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==13)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==13))/sum(assay_classes[,5] == 13))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==13) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==13) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,5] == 13))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==13) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] ==13)))*100,2),
                           
                           nrow(df_variants_reference_panel),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5]>=1)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5]>=1))/nrow(df_variants_reference_panel))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] >=1) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] >=1) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,5] >= 1))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] >=1) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,5] >=1)))*100,2)),
                         ncol = 6, byrow = TRUE)
  rownames(assay_classe) <- c("Protein expression, stability and folding",	"Binding",	"Recombination",
                              "Chromosome/Mitotic apparatus", "Proliferation/Growth",	"Sensitivity",
                              "Localization",	"Transcriptional activation",	"Focus formation",
                              "Catalytic activity",	"Cell viability", 	"Cell cycle checkpoint", "Total")
  colnames(assay_classe) <- c("N", "Met_criteria","Percentage_met_criteria", "Validated", "Percentage_validated", "Percentage_of_all_meeting_criteria_that_validated")
  return(assay_classe)  
}
#By host species
host_species <- function(assay_classes, df_variants_reference_panel, Espec_Sensi){
  host_species <- matrix(c(sum(assay_classes[,9] == 1),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==1)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==1))/sum(assay_classes[,9] == 1))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==1) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==1) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,9] == 1))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==1) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==1)))*100,2),
                           
                           sum(assay_classes[,9] == 2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==2)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==2))/sum(assay_classes[,9] == 2))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==2) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==2) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,9] == 2))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==2) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==2)))*100,2),
                           
                           sum(assay_classes[,9] == 3),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==3)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==3))/sum(assay_classes[,9] == 3))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==3) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==3) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,9] == 3))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==3) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==3)))*100,2),
                           
                           sum(assay_classes[,9] == 4),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==4)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==4))/sum(assay_classes[,9] == 4))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==4) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==4) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,9] == 4))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==4) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==4)))*100,2),
                           
                           sum(assay_classes[,9] == 5),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==5)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==5))/sum(assay_classes[,9] == 5))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==5) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==5) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,9] == 5))*100,2),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==5) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==5)))*100,2),
                           
                           sum(assay_classes[,9] == 6),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==6)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==6))/sum(assay_classes[,9] == 6))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==6) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==6) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,9] == 6))*100,2),  
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==6) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] ==6)))*100,2),
                           
                           nrow(df_variants_reference_panel),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] >=1)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] >=1))/sum(assay_classes[,9] >= 1))*100,2),
                           sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] >=1) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8)),
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] >=1) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum(assay_classes[,9] >= 1))*100,2),  
                           round((sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] >=1) & (Espec_Sensi[,8]>= 0.8) & (Espec_Sensi[,11]>=0.8))/sum((df_variants_reference_panel[,3] >= 10) & (Espec_Sensi[,5] >= 4) & (Espec_Sensi[,6] >= 4) & (assay_classes[,9] >=1)))*100,2)),
                         ncol = 6, byrow = TRUE)
  rownames(host_species) <- c("Homo sapiens",	"Saccharomyces cerevisiae",	"Mus musculus",
                              "Escherichia coli", "Oryctolagus cuniculus",	"Mesocricetus auratus", "Total")
  colnames(host_species) <- c("N", "Met_criteria","Percentage_met_criteria", "Validated", "Percentage_validated", "Percentage_of_all_meeting_criteria_that_validated")
  return(host_species)  
}

#Hi-Set Counts
hi_set_counts <- function(hi_set_tracks){
  hi_set_counts <- data.frame()[1:nrow(hi_set_tracks),]
  #Adds two columns to the end of the table with numbers of "no functional impact" and "functional impact", respectively, for each missense variant 
  hi_set_counts <- cbind(hi_set_counts, rowSums(hi_set_tracks[,4:ncol(hi_set_tracks)] == 0, na.rm= TRUE))
  hi_set_counts <- cbind(hi_set_counts, rowSums(hi_set_tracks[,4:ncol(hi_set_tracks)] == 1, na.rm= TRUE))
  #Adds one column to the end of the table after sum of the columns "no functional impact" and "functional impact"
  hi_set_counts <- cbind(hi_set_counts, rowSums(hi_set_tracks[,4:ncol(hi_set_tracks)] <= 1, na.rm= TRUE))
  colnames(hi_set_counts) <- c("no_functional_impact_Hi_set", "functional_impact_Hi_set", "number_of_tests_Hi_set")
  rownames(hi_set_counts) <- 1:nrow(hi_set_tracks)
  return(hi_set_counts)
}

#Numbers cited in the paper
hi_set_numbers <- function (hi_set_track){
  hi_set_track$"number_of_tests_Hi_set" <- rowSums(hi_set_track[,4:ncol(hi_set_track)] <= 1, na.rm= TRUE)
  y <- sum(hi_set_track$"number_of_tests_Hi_set" >= 1)
  x <- sum((hi_set_track$"T7" > 0)  & (hi_set_track$"number_of_tests_Hi_set" >= 1), na.rm= TRUE) 
  z <- y - x
  hi_set_numbers <- matrix(c(y,
                             z,
                             x,
                             sum(hi_set_track$"number_of_tests_Hi_set" == 1)))
  rownames(hi_set_numbers) <- c("Number of variants tested in Hi Set", "Number of VUS tested in Hi Set",
                                "Number of variants tested present in reference panel", "Number of variants tested 1 time only") 
  colnames(hi_set_numbers)<-c("Final numbers")
  return(hi_set_numbers)
}

#SUPPLEMENTARY TABLE 5C: VUS only (excluding reference variants)
vus_tested <- function(master_table, hi_set_tracks, reference_panel){
  master_table$"T7"[is.na(master_table$"T7")] <- "NA"
  master_table$"number of tests" <- rowSums(master_table[,10:ncol(master_table)] <= 1, na.rm= TRUE)
  hi_set_tracks$"number_of_tests_Hi_set" <- rowSums(hi_set_tracks[,4:ncol(hi_set_tracks)] <= 1, na.rm= TRUE)
  hi_set_tracks$"T7"[is.na(hi_set_tracks$"T7")] <- "NA"
  
  VUS_tested<-matrix(c(sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 1)),
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 1))*100)/nrow(master_table),2), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 1))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 1)),  
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 1))*100)/nrow(master_table),2), 
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 1))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 2)), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 2))*100)/nrow(master_table),2), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 2))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 2)),  
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 2))*100)/nrow(master_table),2), 
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 2))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 10)), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 10))*100)/nrow(master_table),2), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 10))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 10)),  
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 10))*100)/nrow(master_table),2), 
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 10))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 15)), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 15))*100)/nrow(master_table),2), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 15))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 15)),  
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 15))*100)/nrow(master_table),2), 
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 15))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 20)), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 20))*100)/nrow(master_table),2), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 20))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 20)),  
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 20))*100)/nrow(master_table),2), 
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 20))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 25)), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 25))*100)/nrow(master_table),2), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 25))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 25)),  
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 25))*100)/nrow(master_table),2), 
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 25))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 30)), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 30))*100)/nrow(master_table),2), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 30))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 30)),  
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 30))*100)/nrow(master_table),2), 
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 30))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 50)), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 50))*100)/nrow(master_table),2), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 50))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 50)),  
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 50))*100)/nrow(master_table),2), 
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 50))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 60)), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 60))*100)/nrow(master_table),2), 
                       round((sum((master_table$"T7" == "NA") & (master_table$"number of tests" >= 60))*100)/(nrow(master_table)- nrow(reference_panel)),2),
                       
                       sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 60)),  
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 60))*100)/nrow(master_table),2), 
                       round((sum((hi_set_tracks$"T7" == "NA") & (hi_set_tracks$"number_of_tests_Hi_set" >= 60))*100)/(nrow(master_table)- nrow(reference_panel)),2)),
                     ncol = 6, byrow=TRUE)
  rownames(VUS_tested) <- c("VUS Tested >= one time","VUS tested >= two times",
                            "VUS tested >= 10 times", "VUS tested >= 15 times", "VUS tested >= 20 times","VUS tested >= 25 times",
                            "VUS tested >= 30 times","VUS tested >= 50 times","VUS tested >= 60 times") 
  colnames(VUS_tested)<-c("N_All_Functional_tracks","Percentage", "Percentage_VUS", "N_Hi_Set_E", "Percentage", "Percentage_VUS")
  return(VUS_tested)  
}

#Generating the Stable9 ordered
stable8_form <- function(df_variants_reference_panel,my_data, assay_classes, Espec_Sensi, stable2){
  #dataframe with all missense variants tested (including the indeterminated/uncertain)
  whole_variants <- matrix(c(23,	4,	25,	32,	2,	35,	15,	2,	6,	6,	8,	25,	32,	35,	25,	4,	41,	117,	117,	9,	5,	5,
                             15,	16,	25,	80,	32,	5,	6,	2,	70,	13,	3,	2,	65,	9,	11,	15,	4,	5,	8,	13,	77,	36,
                             76,	1,	8,	11,	4,	3,	36,	5,	3,	2,	11,	8,	5,	3,	13,	13,	8,	2,	8,	1,	8,	8,	4,
                             9,	35,	12,	65,	77,	15,	22,	575,	5,	13,	1,	12,	12,	9,	36,	28,	8,	2,	2,	77,	4,	5,	5,
                             18,	18,	333,	3,	2,	3,	12,	28,	7,	11,	25,	3,	2,	35,	11,	117,	2,	2,	2,	10,	25,	7,
                             10,	10,	3,	5,	4,	25,	10,	25,	7,	2086,	6,	1677,	1677,	5,	3,	2,	2,	2,	2))
  #% of tested missense variants classified
  variants_classified <- matrix(c(round((df_variants_reference_panel[,3]/whole_variants)*100,2)),byrow = TRUE)
  #% of all possible missense (n = 11,009) tested
  possible_missense <- matrix(c(round((whole_variants/11009)*100,1)),byrow = TRUE)
  #Number of all documented missense variants classified
  documented_classified <- matrix(c(colSums((my_data[,10:140] <= 1) & (my_data$"T9" == 1),na.rm=TRUE)),byrow=TRUE)
  #% of all documented missense variants classified
  documented_classified <- cbind(documented_classified,matrix(c(round(documented_classified/df_variants_reference_panel[,3]*100,2)),byrow=TRUE)) 
  #gathering all data into a new table
  stable8 <- data.frame()[0:131,]
  stable8 <- cbind(stable8, stable2[,1:2],whole_variants, df_variants_reference_panel[,3],
                   variants_classified,possible_missense,documented_classified, Espec_Sensi[,c(7,5,6,8:14)],
                   df_variants_reference_panel[,c(2,1)],Espec_Sensi[,15:20])
  colnames(stable8)[c(2:9)] <- c("assay", "number of missense variants tested", "number of missense variants classified",
                                 "percentage of tested missense variants classified",
                                 "percentage of all possible missense (n = 11,009) tested",
                                 "number of all documented missense variants classified",
                                 "percentage of all documented missense variants classified",
                                 "number total of Classified")
  return(stable8)
}

#Remove variants not tested in Hi-Set
stable9_form <- function(hi_set_tracks, hi_set_table){
  n_col <- ncol(hi_set_tracks)
  set_numb <- ncol(hi_set_tracks)+5
  functional_col <- ncol(hi_set_tracks)+2
  no_functional_col <- ncol(hi_set_tracks)+1
  
  hi_set_tracks$"no_functional_impact_Hi_set" <- rowSums(hi_set_tracks[,4:n_col] == 0, na.rm = TRUE)
  hi_set_tracks$"functional_impact_Hi_set" <- rowSums(hi_set_tracks[,4:n_col] == 1, na.rm = TRUE)
  hi_set_tracks$"number_of_tests_Hi_set" <- rowSums(hi_set_tracks[,4:n_col] <= 1, na.rm = TRUE)
  clean_hi_set <- data.frame()[1:11009,]
  for (i in 1:nrow(hi_set_tracks)){
    select <- hi_set_tracks[i,ncol(hi_set_tracks)]
    if (select > 0){
      clean_hi_set <- as.data.frame(rbind(clean_hi_set, hi_set_tracks[i,]))
    }
  }
  
  clean_hi_set$"Call" <- ifelse(clean_hi_set$"no_functional_impact_Hi_set" >0 & clean_hi_set$"functional_impact_Hi_set" > 0, "Discordant", "")
  clean_hi_set[is.na(clean_hi_set)] <- "NA"
  for (x in colnames(clean_hi_set)){
    for (y in rownames(hi_set_table)){
      if (x == y){
        evidence_path <- hi_set_table[y,19]
        evidence_benign <- hi_set_table[y,20]
        for (z in 1:nrow(clean_hi_set)){
          
          if (clean_hi_set[z,x] == 1){
            if (clean_hi_set[z,ncol(clean_hi_set)-1]== 1){
              clean_hi_set[z,ncol(clean_hi_set)] <- evidence_path
            }
          }
          if (clean_hi_set[z,x] == 0){
            if (clean_hi_set[z,ncol(clean_hi_set)-1]== 1){
              clean_hi_set[z,ncol(clean_hi_set)] <- evidence_benign
            }
          }
        }
      }
    }
  }
  
  for (z in 1:nrow(clean_hi_set)){
    p <- set_numb
    if (clean_hi_set[z,no_functional_col] > 1 & clean_hi_set[z,functional_col] == 0){
      for (q in 4:n_col){
        if (clean_hi_set[z,q]==0){    
          y <- colnames(clean_hi_set)[q]  
          for (x in rownames(hi_set_table)){
            if (x==y){
              clean_hi_set[z,p] <- hi_set_table [y,20]
              p <- p+1
            }
          }
        } 
      }  
    }
    if (clean_hi_set[z,functional_col] > 1 & clean_hi_set[z,no_functional_col] == 0){
      for (q in 4:n_col){
        if (clean_hi_set[z,q]==1){    
          y <- colnames(clean_hi_set)[q]  
          for (x in rownames(hi_set_table)){
            if (x==y){
              clean_hi_set[z,p] <- hi_set_table [y,19]
              p <- p+1
            }
          }
        } 
      }  
    }  
  }
  return(clean_hi_set)
}

stable10_form <- function(stable9_E_C, hi_set_table){
  set_numb <- 33
  clean_discordant <- data.frame()[1:11009,]
  for (i in 1:nrow(stable9_E_C)){
    select <- stable9_E_C[i,29]
    if (select == "Discordant"){
      clean_discordant <- as.data.frame(rbind(clean_discordant, stable9_E_C[i,1:29]))
    }
  }
  clean_discordant$"Ratio Benign/Path" <- (-1*((clean_discordant[,27]-clean_discordant[,28])/clean_discordant[,27]))
  clean_discordant$"Log2ratio" <- log2(clean_discordant$"Ratio Benign/Path")
  
  clean_discordant$Call <- NA
  clean_discordant[is.na(clean_discordant)] <- "NA"
  for (x in colnames(clean_discordant)){
    for (y in rownames(Espec_Sens_table_E_C)){
      if (x == y){
        evidence_path <- hi_set_table[y,19]
        evidence_benign <- hi_set_table[y,20]
        
        for (z in 1:nrow(clean_discordant)){ #Now moving through the rows to find the variant tested in the track
          
          if (clean_discordant[z,x] == 1){ #Find the functional impact mark
            if (clean_discordant[z,30]== 1){ #Check if the variant was tested only once
              clean_discordant[z,32] <- "Discordant"
            }
          }
          if (clean_discordant[z,x] == 0){
            if (clean_discordant[z,30]== 2){
              clean_discordant[z,32] <- "Discordant"
            }
          }
        }
      }
    }
  }
  return (clean_discordant)
}

stable11_form <- function(stable1, hi_set_table){
  #Setting important numbers before adding columns to the table
  n_col <- ncol(stable1)
  set_numb <- 147
  #Getting a new table with only tested variants from table 1
  stable1_modified <- stable1
  stable1_modified$"no_functional_impact_Hi_set" <- rowSums(stable1_modified[,10:140] == 0, na.rm = TRUE)
  stable1_modified$"functional_impact_Hi_set" <- rowSums(stable1_modified[,10:140] == 1, na.rm = TRUE)
  stable1_modified$"number_of_tests_Hi_set" <- rowSums(stable1_modified[,10:140] <= 1, na.rm = TRUE)
  clean_stable1 <- data.frame()[1:11009,]
  for (i in 1:nrow(stable1_modified)){
    select <- stable1_modified[i,ncol(stable1_modified)]
    if (select > 0){
      clean_stable1 <- as.data.frame(rbind(clean_stable1, stable1_modified[i,]))
    }
  }
  
  clean_stable1$"Ratio Benign/Path" <- (-1*((clean_stable1[,142]-clean_stable1[,143])/clean_stable1[,142]))
  clean_stable1$"Log2ratio" <- log2(clean_stable1$"Ratio Benign/Path")
  
  #Catching the evidence path for each track (column) - one by one
  clean_stable1$Call <- NA
  clean_stable1[is.na(clean_stable1)] <- "NA"
  for (x in colnames(clean_stable1)){
    for (y in rownames(Espec_Sens_table_E_C)){
      if (x == y){
        evidence_path <- hi_set_table[y,19]
        evidence_benign <- hi_set_table[y,20]
        
        for (z in 1:nrow(clean_stable1)){ #Now moving through the rows to find the variant tested in the track
          
          if (clean_stable1[z,x] == 1){ #Find the functional impact mark
            if (clean_stable1[z,143]== 1){ #Check if the variant was tested only once
              clean_stable1[z,146] <- evidence_path
            }
          }
          if (clean_stable1[z,x] == 0){
            if (clean_stable1[z,143]== 1){
              clean_stable1[z,146] <- evidence_benign
            }
          }
        }
      }
    }
  }
  
  for (z in 1:2701){
    p <- set_numb
    if (clean_stable1[z,141] > 1 & clean_stable1[z,142] == 0){
      for (q in 10:140){
        if (clean_stable1[z,q]==0){    
          y <- colnames(clean_stable1)[q]  
          for (x in rownames(hi_set_table)){
            if (x==y){
              clean_stable1[z,p] <- hi_set_table [y,20]
              p <- p+1
            }
          }
        } 
      }  
    }
    if (clean_stable1[z,142] > 1 & clean_stable1[z,141] == 0){
      for (q in 10:140){
        if (clean_stable1[z,q]==1){    
          y <- colnames(clean_stable1)[q]  
          for (x in rownames(hi_set_table)){
            if (x==y){
              clean_stable1[z,p] <- hi_set_table [y,19]
              p <- p+1
            }
          }
        } 
      }  
    }
    if (clean_stable1[z,145] > 0){
      for (q in 10:140){
        if (clean_stable1[z,q]==0){    
          y <- colnames(clean_stable1)[q]  
          for (x in rownames(hi_set_table)){
            if (x==y){
              clean_stable1[z,p] <- hi_set_table [y,20]
              p <- p+1
            }
          }
        } 
      }  
    }
    if (clean_stable1[z,145] < 0){
      for (q in 10:140){
        if (clean_stable1[z,q]==1){    
          y <- colnames(clean_stable1)[q]  
          for (x in rownames(hi_set_table)){
            if (x==y){
              clean_stable1[z,p] <- hi_set_table [y,19]
              p <- p+1
            }
          }
        } 
      }  
    }
  }
  
  return(clean_stable1)
}


#OUTPUTS
#Espec_sens_Calculation
Espec_Sens_table_E_C <- espec_sens_calculation(reference_panel_E_C)
#Counts from Reference panel Enigma and Reference panel Enigma + ClinVar
reference_panel_E_C_counts <- count_reference_panel(reference_panel_E_C, stable1) 
#Hi-Set Table
hi_set_table_E_C <- hi_set_table(Espec_Sens_table_E_C, reference_panel_E_C_counts)
#Hi-set tracks
hi_set_tracks_E_C <- hi_set_tracks(hi_set_table_E_C, stable1)
#Hi-set Counts
hi_set_counts_E_C <- hi_set_counts(hi_set_tracks_E_C)
#Hi-set numbers cited on the paper
hi_set_numbers_E_C <- hi_set_numbers(hi_set_tracks_E_C)

#######################################################|||| TABLES ||||#######################################################

#STABL5_A_Number of variants tested per 
variants_tested <- count_variant_tested(stable1,stable1_counts)
#STABLE5_B_Documented variants tested
BRCAExchange_E_C_counts <- count_BRCAExchange(stable1,stable1_counts,reference_panel_E_C)
#STABLE5_C_VUS_ONLY(Excluding reference panel)
VUS_only_E_C <- vus_tested(stable1, hi_set_tracks_E_C, reference_panel_E_C)

#STABLE6_A_Assay throughput
tested_variants_E_C <- count_variant_track(reference_panel_E_C_counts)
#STABLE6_B_Number of tracks with number in reference panel [non-pathogenic; pathogenic]>= [x:y]
ratio_criteria_1_E_C <- ratio_criteria(Espec_Sens_table_E_C,reference_panel_E_C_counts)
#STABLE6_C_Number of tracks with number in reference panel tested variants >= 10 and [non-pathogenic; pathogenic]>= [x:y]
ratio_criteria_2_E_C <- ratio_criteria_2(Espec_Sens_table_E_C, reference_panel_E_C_counts)

#STABLE7_A_Assay_classes
assay_class_table_E_C <- assay_classes(stable2, reference_panel_E_C_counts, Espec_Sens_table_E_C)
#STABLE7_B_Host_Species_classes
host_species_table_E_C <- host_species(stable2, reference_panel_E_C_counts, Espec_Sens_table_E_C)

#STABLE8 
stable8_E_C <- stable8_form (reference_panel_E_C_counts, stable1, assay_class_table_E_C, Espec_Sens_table_E_C, stable2)

#STABLE9
stable9_E_C <- stable9_form (hi_set_tracks_E_C, hi_set_table_E_C)

#STABLE10
stable10_E_C <- stable10_form (stable9_E_C)

#STABLE11
stable11_E_C <- stable11_form (stable1, Espec_Sens_table_E_C)

#######################################################|||| OUTPUT ||||#######################################################

#Output excel file for reference panel [ENIGMA + CLIN VAR]
write.xlsx(variants_tested,file="output_E_C.xlsx", sheetName="STable6_Missense_Variants_tested")
write.xlsx(BRCAExchange_E_C_counts,file="output_E_C.xlsx", sheetName="STable6_Documented_Variants(BRCAExchange)", append=TRUE)
write.xlsx(VUS_only_E_C,file="output_E_C.xlsx", sheetName="STable6_VUS_only_(excluding reference variants)", append=TRUE)
write.xlsx(tested_variants_E_C,file="output_E_C.xlsx", sheetName="STable7_Functional_track_throughput", append=TRUE)
write.xlsx(ratio_criteria_1_E_C,file="output_E_C.xlsx", sheetName="STable7_#_variants_classified_by_track", append=TRUE)
write.xlsx(ratio_criteria_2_E_C,file="output_E_C.xlsx", sheetName="STable7_#_tracks_meeting_criteria", append=TRUE)
write.xlsx(assay_class_table_E_C,file="output_E_C.xlsx", sheetName="STable8_Assay_classes", append=TRUE)
write.xlsx(host_species_table_E_C,file="output_E_C.xlsx", sheetName="STable8_host_classes", append=TRUE)
write.xlsx(stable8_E_C,file="output_E_C.xlsx", sheetName="STable9_Sensitivity_specificity", append=TRUE, row.names=FALSE)
write.xlsx(stable9_E_C,file="output_E_C.xlsx", sheetName="STable 10_Hi_Set_Approach", append=TRUE, row.names=FALSE)
write.xlsx(stable11_E_C,file="output_E_C.xlsx", sheetName="STable 11_Majority_Voting_approach", append=TRUE)
#Extra numbers
write.xlsx(hi_set_numbers_E_C,file="output_E_C.xlsx", sheetName="Numbers_from_Hi_Set_Pick", append=TRUE)


