### Modelling pSIA + oSIA strategy ###

# the code below specifies parameters for the pSIA + oSIA model  
# to run the model, you require the csv files titled "random_importations_2_per_year_SS.csv" 


## IMPORTANT! ##
# One limitation of SimInf is that the package can only process "events" with pre-determined times
# For the pSIA + oSIA model, we need to take inventory of which nodes had at least 1 case every 90 days (SOP for polio)
# Therefore, the model needs to be stopped every 90 days in order to save a vector of all nodes requiring 
# an outbreak response. From there, we re-specify the initial conditions for the model running for the next 90 days
# Then, the model can resume, conducting an oSIA only in nodes with one case. 
# If we assume a time horizon of 5 years, this requires quite a few 90-day intervals
# with potential for an oSIA, hence the longer and complex code below. 

# Outputs from this code are the same as the other scenarios and 
# can be further processed using the accompanying script titled "post_model_processing_costs.R"
rm(list=ls())	
require(SimInf)
library(tidyverse)
library(dplyr)
library(data.table)

beta <- 0.375
gamma <- 0.125
r <- beta/gamma
caseinfection <- 1/200
I0 <- 0
sims <- 1:1000 #change for number of simulations
pop <- 8000000
reps <- 1000
ri_coverage <- 0.25
sia_coverage <- 0.25
ipv_coverage <- ri_coverage
ve <- 0.5 #vaccine efficacy for bOPV
doses <- 3 
seroconversion <- (1 -(1-ve)^doses)
birth_rate <- 4000 #or 5x10^-4 (4000 births per 8 Million) = avg across africa
death_rate = birth_rate
births_S <- round(birth_rate * (1-ri_coverage), digits=0)
births_S3 <- round(birth_rate * ri_coverage, digits=0)

########################################################################################################
### ONCE STEADY STATE IS REACHED, INTRODUCE 1 INFECTION, IMPORTATIONS AND MODEL DIFFERENT STRATEGIES ###
########################################################################################################

# ANNUAL PREVENTATIVE CAMPAIGNS
psia_annual <- data.frame(event = "intTrans", 
                          #campaigns introduced at year 30 to reach steady state
                          time = rep(c(1, 366, 731, 1096, 1461), each=length(sims)),
                          #for biannual SIAs, change all subsequent SIA code to include only the time points below
                          #time = rep(c(1, 731, 1461), each=length(sims)), #biannial SIAs - 
                          node= 1:length(sims), 
                          dest = 0,
                          n = 0,
                          proportion = sia_coverage,
                          select = 5, shift = 1)

psia_no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(c(3, 368, 733, 1098, 1463), each=length(sims)),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = 1, # move to next S tier(Sn+1)
                                select = 6, shift = 2)

psia_seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(c(2, 367, 732, 1097, 1462), each=length(sims)),
                               node= 1:length(sims), 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) , # seroconvert and move to Rv
                               select = 7, shift = 3)

psia_seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(c(2, 367, 732, 1097, 1462), each=length(sims)),
                               node= 1:length(sims), 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) , # another dose worth of protection, seroconvert and move to Rv
                               select = 8, shift = 3)

psia_seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(c(2, 367, 732, 1097, 1462), each=length(sims)),
                               node= 1:length(sims), 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , # another dose worth of protection, seroconvert and move to Rv
                               select = 9, shift = 3)

psia_seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(c(2, 367, 732, 1097, 1462), each=length(sims)),
                               node= 1:length(sims), 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

psia_seroMAX <- data.frame(event = "intTrans", 
                      time = rep(c(3, 368, 733, 1098, 1463), each=length(sims)),
                      node= 1:length(sims), 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)


#############################################
############## pSIAs + oSIAs ################
#############################################
# oSIA occurs only in the node where CASE >=1, within 90 days from the first case

# NEW f-time (5 year time horizon)
f_time <- 5*365 #time horizon for steady state

birthsS <- data.frame(event = "enter", 
                      time = c(rep(1:89, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(1:89, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(1:89, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

infection <- data.frame(event = "intTrans", 
                        time = rep((1), each=length(sims)),
                        node= 1:length(sims), 
                        dest = 0, 
                        n = 1,
                        proportion = 0, 
                        select = 10, shift = 4)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(1:89, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(1:89, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3) , # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)


random_importations <- read.csv("random_importations_2_per_year_SS.csv")
random_importations$time <- random_importations$time - (50*365)
# importations only for the first 89 days
importation_sia <- data.frame(event = "intTrans", 
                              time = random_importations, #same vector of random dates as RI only scenario
                              node= 1:length(sims), 
                              dest = 0, 
                              n = 1,
                              proportion = 0, 
                              select = 4, shift = 4)
importation_sia0 <-filter(importation_sia, time<=89)

### preventative campaigns just for this time interval
psia1 <-filter(psia_annual, time<=89)
psia1_no_seroconversion <-filter(psia_no_seroconversion, time<=89)
psia1_seroconversionS1 <-filter(psia_seroconversionS1, time<=89)
psia1_seroconversionS2 <-filter(psia_seroconversionS2, time<=89)
psia1_seroconversionS3 <-filter(psia_seroconversionS3, time<=89)
psia1_seroconversionS4 <-filter(psia_seroconversionS4, time<=89)
psia1_seroMAX <-filter(psia_seroMAX, time<=89)

# "Select" matrix
Et <- matrix(c(1, rep(0,15), #births R0
               rep(0,8),1,rep(0,7), #births S3
               rep(1,11),0,1,1,1,0, #deaths
               1, rep(0,15), #dynamic importation
               1,0,1,0,1,0,1,0,1,rep(0, 7), #vaccination
               0,1,0,1,0,1,0,1,rep(0,8), #NO seroconversion
               0,1,0,rep(0,13), #seroconversionS1
               0,0,0,1,0,rep(0,11), #seroconversionS2
               0,0,0,0,0,1,0,rep(0,9), #seroconversionS3
               1, rep(0,15), #infection SS
               rep(0,8),1,0,0,0,0,0,0,0,#seroconversionS3r 
               rep(0,7),1,rep(0,8), #seroconversionS4
               rep(0,9),1,rep(0,6), #seroMAX
               rep(0,8),1,rep(0,7)), #ipv
             nrow = 16, ncol = 14, 
             dimnames = list(c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v","S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum"), 
                             c("birthsS0", "birthsS3", "deaths", "importations", "vaccination",
                               "no_seroconversion", "seroconversionS1", "seroconversionS2", "seroconversionS3", "infection",
                               "seroconversionS3r", "seroconversionS4", "seroMAX", "ipv")))

# "Shift" matrix: column 1 = vaccination (shift to Snv vaccinated compartment), column 2 = seroconversion (shift to Rv), 
# column 3 = no seroconversion (shift to next vaccination tier)
Nt <- matrix(c(1,0,1,0,1,0,1,0,1,rep(0, 7), #vaccination
               0,1,0,1,0,1,0,2,rep(0,8), #NO seroconversion
               0,9,0,7,0,5,0,2,0,1,rep(0,6), #seroconversion
               12, rep(0,15), #importation
               rep(0,8),2,rep(0,7)), #ipv vaccination
             nrow = 16, ncol = 5,
             dimnames = list(c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v","S3r","S4", "Rv", "Ii", "I", "C", "Rn", "Ccum"), 
                             c("vaccination", "no_seroconversion", "seroconversion", "importation", "ipv")))

u0psia_osia <- read.csv(file="u0psia_ri25.csv") #initial conditions from steady state
create.StocSIRobr0 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- u0psia_osia
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v","S3r","S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  # accounting for the dynamic process that happens regardless of SIAs
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  
  tspan <- seq(from = 1, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, infection, importation_sia0,
                               psia1, #preventative campaign
                               psia1_no_seroconversion, psia1_seroconversionS1,
                               psia1_seroconversionS2, psia1_seroconversionS3,
                               psia1_seroconversionS4, psia1_seroMAX,
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs from first 89 days
SIRmodelobr0 <- create.StocSIRobr0(I0, beta, gamma, caseinfection, 89) 
out89 <- run(model = SIRmodelobr0) 
df89 <- trajectory(out89) # df for the first 89 days simulated, check and see which nodes had cases
obr_df <- subset(df89, df89$C>=1)
obr1_nodes <- sort(unique(obr_df$node)) # nodes with 1 or more cases

#############################################
################### OBR 1 ###################
#############################################
############## days 90 to 179 ###############
#############################################
importation_sia1 <-filter(importation_sia, time>=90 & time<=179)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(90:179, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(90:179, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(90:179, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(90:179, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(90:179, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr1 <- data.frame(event = "intTrans", 
                   time = rep(90, length(obr1_nodes)),
                   node= obr1_nodes, 
                   dest = 0, 
                   n = 0,
                   proportion = sia_coverage,
                   select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(92, length(obr1_nodes)),
                                node= obr1_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, # move to next S tier(Sn+1)
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(91, length(obr1_nodes)),
                               node= obr1_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) , # seroconvert and move to Rv
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(91, length(obr1_nodes)),
                               node= obr1_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) , # another dose worth of protection, seroconvert and move to Rv
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(91, length(obr1_nodes)),
                               node= obr1_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , # another dose worth of protection, seroconvert and move to Rv
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(91, length(obr1_nodes)),
                               node= obr1_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(92, length(obr1_nodes)),
                      node= obr1_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(df89, time==89) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame() 

# now, schedule the oSIA on day 90 and allow model to run for an additional 89 days
create.StocSIRobr1 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 90, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia1,
                               obr1, seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               seroconversionS4, seroMAX,
                               seroconversionS3r, ipv), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs from second 89 days
SIRmodelobr1 <- create.StocSIRobr1(I0, beta, gamma, caseinfection, 179) 
out179 <- run(model = SIRmodelobr1) 
df179 <- trajectory(out179) # df for the second 89 days simulated, check and see which nodes had cases
obr_df <- subset(df179, df179$C>=1)
obr2_nodes <- sort(unique(obr_df$node)) # nodes with 1 or more cases
end_time <- 179
end_dat <- df179

#############################################
################### OBR 2 ###################
#############################################
############# days 180 to 269 ###############
#############################################
importation_sia2 <-filter(importation_sia, time>=180 & time<=269)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(180:269, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(180:269, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(180:269, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(180:269, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(180:269, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr2 <- data.frame(event = "intTrans", 
                   time = rep(180, length(obr2_nodes)),
                   node= obr2_nodes, 
                   dest = 0, n = 0, 
                   proportion = sia_coverage,
                   select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(182, length(obr2_nodes)),
                                node= obr2_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(181, length(obr2_nodes)),
                               node= obr2_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(181, length(obr2_nodes)),
                               node= obr2_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(181, length(obr2_nodes)),
                               node= obr2_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(181, length(obr2_nodes)),
                               node= obr2_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(182, length(obr2_nodes)),
                      node= obr2_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(df179, time==179) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame() 

create.StocSIRobr2 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 180, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia2,
                               obr2, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               seroconversionS3r, ipv), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs from next interval of 89 days
SIRmodelobr2 <- create.StocSIRobr2(I0, beta, gamma, caseinfection, 269) 
out269 <- run(model = SIRmodelobr2) 
df269 <- trajectory(out269) 
obr_df <- subset(df269, df269$C>=1)
obr3_nodes <- sort(unique(obr_df$node)) 
end_time <- 269
end_dat <- df269

#############################################
################### OBR 3 ###################
#############################################
############# days 270 to 359 ###############
#############################################
importation_sia3 <-filter(importation_sia, time>=270 & time<=359)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(270:359, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(270:359, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(270:359, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(270:359, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(270:359, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr3 <- data.frame(event = "intTrans", 
                   time = rep(270, length(obr3_nodes)),
                   node= obr3_nodes, 
                   dest = 0, n = 0, 
                   proportion = sia_coverage,
                   select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(272, length(obr3_nodes)),
                                node= obr3_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(271, length(obr3_nodes)),
                               node= obr3_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(271, length(obr3_nodes)),
                               node= obr3_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(271, length(obr3_nodes)),
                               node= obr3_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(271, length(obr3_nodes)),
                               node= obr3_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(272, length(obr3_nodes)),
                      node= obr3_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame() 

create.StocSIRobr3 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 270, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia3,
                               obr3, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               ipv, seroconversionS3r),  
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr3 <- create.StocSIRobr3(I0, beta, gamma, caseinfection, 359) 
out359 <- run(model = SIRmodelobr3) 
df359 <- trajectory(out359) 
obr_df <- subset(df359, df359$C>=1)
obr4_nodes <- sort(unique(obr_df$node)) 
end_time <- 359
end_dat <- df359

#############################################
################# OBR 4 #####################
#############################################
############# days 360 to 449 ###############
#############################################
importation_sia4 <-filter(importation_sia, time>=360 & time<=449)
psia2 <-filter(psia_annual, time>=360 & time<=449)
psia2_no_seroconversion <-filter(psia_no_seroconversion,time>=360 & time<=449)
psia2_seroconversionS1 <-filter(psia_seroconversionS1, time>=360 & time<=449)
psia2_seroconversionS2 <-filter(psia_seroconversionS2, time>=360 & time<=449)
psia2_seroconversionS3 <-filter(psia_seroconversionS3, time>=360 & time<=449)
psia2_seroconversionS4 <-filter(psia_seroconversionS4, time>=360 & time<=449)
psia2_seroMAX <-filter(psia_seroMAX,time>=360 & time<=449)


birthsS <- data.frame(event = "enter", 
                      time = c(rep(360:449, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(360:449, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(360:449, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(360:449, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(360:449, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr4 <- data.frame(event = "intTrans", 
                   time = rep(360, length(obr4_nodes)),
                   node= obr4_nodes, 
                   dest = 0, n = 0, 
                   proportion = sia_coverage,
                   select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(362, length(obr4_nodes)),
                                node= obr4_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(361, length(obr4_nodes)),
                               node= obr4_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(361, length(obr4_nodes)),
                               node= obr4_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(361, length(obr4_nodes)),
                               node= obr4_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(361, length(obr4_nodes)),
                               node= obr4_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(362, length(obr4_nodes)),
                      node= obr4_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame() 

create.StocSIRobr4 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 360, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia4,
                               obr4, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               psia2,
                               psia2_no_seroconversion, psia2_seroconversionS1, 
                               psia2_seroconversionS2, psia2_seroconversionS3,
                               psia2_seroconversionS4, psia2_seroMAX,
                               ipv, seroconversionS3r),  
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr4 <- create.StocSIRobr4(I0, beta, gamma, caseinfection, 449) 
out449 <- run(model = SIRmodelobr4) 
df449 <- trajectory(out449) 
obr_df <- subset(df449, df449$C>=1)
obr5_nodes <- sort(unique(obr_df$node)) 
end_time <- 449
end_dat <- df449

#############################################
################# OBR 5 #####################
#############################################
############# days 450 to 539 ###############
#############################################
importation_sia5 <-filter(importation_sia, time>=450 & time<=539)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(450:539, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(450:539, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(450:539, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(450:539, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(450:539, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr5 <- data.frame(event = "intTrans", 
                   time = rep(450, length(obr5_nodes)),
                   node= obr5_nodes, 
                   dest = 0, n = 0, 
                   proportion = sia_coverage,
                   select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(452, length(obr5_nodes)),
                                node= obr5_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(451, length(obr5_nodes)),
                               node= obr5_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(451, length(obr5_nodes)),
                               node= obr5_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(451, length(obr5_nodes)),
                               node= obr5_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(451, length(obr5_nodes)),
                               node= obr5_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(452, length(obr5_nodes)),
                      node= obr5_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)


# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame() 

create.StocSIRobr5 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 450, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia5,
                               obr5, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               ipv, seroconversionS3r),  
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr5 <- create.StocSIRobr5(I0, beta, gamma, caseinfection, 539) 
out539 <- run(model = SIRmodelobr5) 
df539 <- trajectory(out539) 
obr_df <- subset(df539, df539$C>=1)
obr6_nodes <- sort(unique(obr_df$node)) 
end_time <- 539
end_dat <- df539

#############################################
############ Outbreak response 6 ############
#############################################
############# days 540 to 629 ###############
#############################################
importation_sia6 <-filter(importation_sia, time>=540 & time<=629)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(540:629, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(540:629, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(540:629, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(540:629, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(540:629, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr6 <- data.frame(event = "intTrans", 
                   time = rep(540, length(obr6_nodes)),
                   node= obr6_nodes, 
                   dest = 0, n = 0, 
                   proportion = sia_coverage,
                   select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(542, length(obr6_nodes)),
                                node= obr6_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(541, length(obr6_nodes)),
                               node= obr6_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(541, length(obr6_nodes)),
                               node= obr6_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(541, length(obr6_nodes)),
                               node= obr6_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(541, length(obr6_nodes)),
                               node= obr6_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(542, length(obr6_nodes)),
                      node= obr6_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame()  

create.StocSIRobr6 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 540, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia6,
                               obr6, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr6 <- create.StocSIRobr6(I0, beta, gamma, caseinfection, 629) 
out629 <- run(model = SIRmodelobr6) 
df629 <- trajectory(out629) 
obr_df <- subset(df629, df629$C>=1)
obr7_nodes <- sort(unique(obr_df$node)) 
end_time <- 629
end_dat <- df629

#############################################
############ Outbreak response 7 ############
#############################################
############# days 630 to 719 ###############
#############################################
importation_sia7 <-filter(importation_sia, time>=630 & time<=719)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(630:719, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(630:719, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(630:719, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(630:719, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(630:719, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr7 <- data.frame(event = "intTrans", 
                   time = rep(630, length(obr7_nodes)),
                   node= obr7_nodes, 
                   dest = 0, n = 0, 
                   proportion = sia_coverage,
                   select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(632, length(obr7_nodes)),
                                node= obr7_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(631, length(obr7_nodes)),
                               node= obr7_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(631, length(obr7_nodes)),
                               node= obr7_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(631, length(obr7_nodes)),
                               node= obr7_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(631, length(obr7_nodes)),
                               node= obr7_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(632, length(obr7_nodes)),
                      node= obr7_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame()   

create.StocSIRobr7 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 630, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia7,
                               obr7, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr7 <- create.StocSIRobr7(I0, beta, gamma, caseinfection, 719) 
out719 <- run(model = SIRmodelobr7) 
df719 <- trajectory(out719) 
obr_df <- subset(df719, df719$C>=1)
obr8_nodes <- sort(unique(obr_df$node)) 
end_time <- 719
end_dat <- df719

#############################################
############ Outbreak response 8 ############
#############################################
############# days 720 to 809 ###############
#############################################
importation_sia8 <-filter(importation_sia, time>=720 & time<=809)
psia3 <-filter(psia_annual, time>=720 & time<=809)
psia3_no_seroconversion <-filter(psia_no_seroconversion, time>=720 & time<=809)
psia3_seroconversionS1 <-filter(psia_seroconversionS1, time>=720 & time<=809)
psia3_seroconversionS2 <-filter(psia_seroconversionS2, time>=720 & time<=809)
psia3_seroconversionS3 <-filter(psia_seroconversionS3, time>=720 & time<=809)
psia3_seroconversionS4 <-filter(psia_seroconversionS4, time>=720 & time<=809)
psia3_seroMAX <-filter(psia_seroMAX, time>=720 & time<=809)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(720:809, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(720:809, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(720:809, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(720:809, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(720:809, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr8 <- data.frame(event = "intTrans", 
                   time = rep(720, length(obr8_nodes)),
                   node= obr8_nodes, 
                   dest = 0, n = 0, 
                   proportion = sia_coverage,
                   select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(722, length(obr8_nodes)),
                                node= obr8_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(721, length(obr8_nodes)),
                               node= obr8_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(721, length(obr8_nodes)),
                               node= obr8_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(721, length(obr8_nodes)),
                               node= obr8_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(721, length(obr8_nodes)),
                               node= obr8_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(722, length(obr8_nodes)),
                      node= obr8_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame() 

create.StocSIRobr8 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 720, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia8,
                               obr8, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               psia3,
                               psia3_no_seroconversion, psia3_seroconversionS1,
                               psia3_seroconversionS2, psia3_seroconversionS3,
                               psia3_seroconversionS4, psia3_seroMAX,
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr8 <- create.StocSIRobr8(I0, beta, gamma, caseinfection, 809) 
out809 <- run(model = SIRmodelobr8) 
df809 <- trajectory(out809) 
obr_df <- subset(df809, df809$C>=1)
obr9_nodes <- sort(unique(obr_df$node)) 
end_time <- 809
end_dat <- df809

#############################################
############ Outbreak response 9 ############
#############################################
############# days 810 to 899 ###############
#############################################
importation_sia9 <-filter(importation_sia, time>=810 & time<=899)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(810:899, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(810:899, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(810:899, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(810:899, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(810:899, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr9 <- data.frame(event = "intTrans", 
                   time = rep(810, length(obr9_nodes)),
                   node= obr9_nodes, 
                   dest = 0, n = 0, 
                   proportion = sia_coverage,
                   select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(811, length(obr9_nodes)),
                                node= obr9_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(811, length(obr9_nodes)),
                               node= obr9_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(811, length(obr9_nodes)),
                               node= obr9_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(811, length(obr9_nodes)),
                               node= obr9_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(811, length(obr9_nodes)),
                               node= obr9_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(812, length(obr9_nodes)),
                      node= obr9_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)


# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame() 

# now, schedule the oSIA on day 90 and allow model to run for an additional 89 days
create.StocSIRobr9 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 810, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia9,
                               obr9, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr9 <- create.StocSIRobr9(I0, beta, gamma, caseinfection, 899) 
out899 <- run(model = SIRmodelobr9) 
df899 <- trajectory(out899) # df for the second 89 days simulated, check and see which nodes had cases
obr_df <- subset(df899, df899$C>=1)
obr10_nodes <- sort(unique(obr_df$node)) # nodes with 1 or more cases
end_time <- 899
end_dat <- df899

#############################################
########### Outbreak response 10 ############
#############################################
############# days 900 to 989 ###############
#############################################
importation_sia10 <-filter(importation_sia, time>=900 & time<=989)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(900:989, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(900:989, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(900:989, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(900:989, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(900:989, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr10 <- data.frame(event = "intTrans", 
                    time = rep(900, length(obr10_nodes)),
                    node= obr10_nodes, 
                    dest = 0, n = 0, 
                    proportion = sia_coverage,
                    select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(902, length(obr10_nodes)),
                                node= obr10_nodes,  
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(901, length(obr10_nodes)),
                               node= obr10_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(901, length(obr10_nodes)),
                               node= obr10_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(901, length(obr10_nodes)),
                               node= obr10_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(901, length(obr10_nodes)),
                               node= obr10_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(902, length(obr10_nodes)),
                      node= obr10_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)


# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame()  

# now, schedule the oSIA on day 90 and allow model to run for an additional 89 days
create.StocSIRobr10 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 900, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia10,
                               obr10, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr10 <- create.StocSIRobr10(I0, beta, gamma, caseinfection, 989) 
out989 <- run(model = SIRmodelobr10) 
df989 <- trajectory(out989) # df for the second 89 days simulated, check and see which nodes had cases
obr_df <- subset(df989, df989$C>=1)
obr11_nodes <- sort(unique(obr_df$node)) # nodes with 1 or more cases
end_time <- 989
end_dat <- df989

#############################################
########### Outbreak response 11 ############
#############################################
############# days 990 to 1079 ##############
#############################################
importation_sia11 <-filter(importation_sia, time>=990 & time<=1079)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(990:1079, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(990:1179, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(990:1179, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0,  
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(990:1179, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(990:1179, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr11 <- data.frame(event = "intTrans", 
                    time = rep(990, length(obr11_nodes)),
                    node= obr11_nodes, 
                    dest = 0, n = 0, 
                    proportion = sia_coverage,
                    select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(992, length(obr11_nodes)),
                                node= obr11_nodes,  
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(991, length(obr11_nodes)),
                               node= obr11_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(991, length(obr11_nodes)),
                               node= obr11_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(991, length(obr11_nodes)),
                               node= obr11_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(991, length(obr11_nodes)),
                               node= obr11_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(992, length(obr11_nodes)),
                      node= obr11_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame()  

# now, schedule the oSIA on day 90 and allow model to run for an additional 89 days
create.StocSIRobr11 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 990, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia11,
                               obr11, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr11 <- create.StocSIRobr11(I0, beta, gamma, caseinfection, 1079) 
out1079 <- run(model = SIRmodelobr11) 
df1079 <- trajectory(out1079) # df for the second 89 days simulated, check and see which nodes had cases
obr_df <- subset(df1079, df1079$C>=1)
obr12_nodes <- sort(unique(obr_df$node)) # nodes with 1 or more cases
end_time <- 1079
end_dat <- df1079

#############################################
########### Outbreak response 12 ############
#############################################
############# days 1080 to 1169 #############
#############################################
importation_sia12 <-filter(importation_sia, time>=1080 & time<=1169)
psia4 <-filter(psia_annual, time>=1080 & time<=1169)
psia4_no_seroconversion <-filter(psia_no_seroconversion, time>=1080 & time<=1169)
psia4_seroconversionS1 <-filter(psia_seroconversionS1, time>=1080 & time<=1169)
psia4_seroconversionS2 <-filter(psia_seroconversionS2, time>=1080 & time<=1169)
psia4_seroconversionS3 <-filter(psia_seroconversionS3, time>=1080 & time<=1169)
psia4_seroconversionS4 <-filter(psia_seroconversionS4, time>=1080 & time<=1169)
psia4_seroMAX <-filter(psia_seroMAX, time>=1080 & time<=1169)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(1080:1169, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(1080:1169, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(1080:1169, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(1080:1169, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(1080:1169, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr12 <- data.frame(event = "intTrans", 
                    time = rep(1080, length(obr12_nodes)),
                    node= obr12_nodes, 
                    dest = 0, n = 0, 
                    proportion = sia_coverage,
                    select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(1082, length(obr12_nodes)),
                                node= obr12_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(1081, length(obr12_nodes)),
                               node= obr12_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(1081, length(obr12_nodes)),
                               node= obr12_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(1081, length(obr12_nodes)),
                               node= obr12_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(1081, length(obr12_nodes)),
                               node= obr12_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(1082, length(obr12_nodes)),
                      node= obr12_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)


# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame()  

# now, schedule the oSIA on day 90 and allow model to run for an additional 89 days
create.StocSIRobr12 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 1080, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia12,
                               obr12, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               psia4,
                               psia4_no_seroconversion, psia4_seroconversionS1, 
                               psia4_seroconversionS2, psia4_seroconversionS3, 
                               psia4_seroconversionS4, psia4_seroMAX, 
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr12 <- create.StocSIRobr12(I0, beta, gamma, caseinfection, 1169) 
out1169 <- run(model = SIRmodelobr12) 
df1169 <- trajectory(out1169) 
obr_df <- subset(df1169, df1169$C>=1)
obr13_nodes <- sort(unique(obr_df$node)) # nodes with 1 or more cases
end_time <- 1169
end_dat <- df1169

#############################################
########### Outbreak response 13 ############
#############################################
############# days 1170 to 1259 #############
#############################################
importation_sia13 <-filter(importation_sia, time>=1170 & time<=1259)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(1170:1259, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(1170:1259, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(1170:1259, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(1170:1259, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(1170:1259, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr13 <- data.frame(event = "intTrans", 
                    time = rep(1170, length(obr13_nodes)),
                    node= obr13_nodes, 
                    dest = 0, n = 0, 
                    proportion = sia_coverage,
                    select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(1172, length(obr13_nodes)),
                                node= obr13_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(1171, length(obr13_nodes)),
                               node= obr13_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(1171, length(obr13_nodes)),
                               node= obr13_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(1171, length(obr13_nodes)),
                               node= obr13_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(1171, length(obr13_nodes)),
                               node= obr13_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(1172, length(obr13_nodes)),
                      node= obr13_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)


# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame()  

# now, schedule the oSIA on day 90 and allow model to run for an additional 89 days
create.StocSIRobr13 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 1170, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia13,
                               obr13, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr13 <- create.StocSIRobr13(I0, beta, gamma, caseinfection, 1259) 
out1259 <- run(model = SIRmodelobr13) 
df1259 <- trajectory(out1259) # df for the second 89 days simulated, check and see which nodes had cases
obr_df <- subset(df1259, df1259$C>=1)
obr14_nodes <- sort(unique(obr_df$node)) # nodes with 1 or more cases
end_time <- 1259
end_dat <- df1259

#############################################
########### Outbreak response 14 ############
#############################################
############# days 1260 to 1349 #############
#############################################
importation_sia14 <-filter(importation_sia, time>=1260 & time<=1349)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(1260:1349, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(1260:1439, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(1260:1439, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(1260:1439, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(1260:1439, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr14 <- data.frame(event = "intTrans", 
                    time = rep(1260, length(obr14_nodes)),
                    node= obr14_nodes, 
                    dest = 0, n = 0, 
                    proportion = sia_coverage,
                    select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(1262, length(obr14_nodes)),
                                node= obr14_nodes,
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(1261, length(obr14_nodes)),
                               node= obr14_nodes,
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(1261, length(obr14_nodes)),
                               node= obr14_nodes,
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(1261, length(obr14_nodes)),
                               node= obr14_nodes,
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(1261, length(obr14_nodes)),
                               node= obr14_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(1262, length(obr14_nodes)),
                      node= obr14_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame() 

# now, schedule the oSIA on day 90 and allow model to run for an additional 89 days
create.StocSIRobr14 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 1260, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia14,
                               obr14, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               ipv, seroconversionS3r),  
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr14 <- create.StocSIRobr14(I0, beta, gamma, caseinfection, 1349) 
out1349 <- run(model = SIRmodelobr14) 
df1349 <- trajectory(out1349) # df for the second 89 days simulated, check and see which nodes had cases
obr_df <- subset(df1349, df1349$C>=1)
obr15_nodes <- sort(unique(obr_df$node)) # nodes with 1 or more cases
end_time <- 1349
end_dat <- df1349

#############################################
########### Outbreak response 15 ############
#############################################
############# days 1350 to 1439 #############
#############################################
importation_sia15 <-filter(importation_sia, time>=1350 & time<=1439)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(1350:1439, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(1350:1439, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(1350:1439, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(1350:1439, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(1350:1439, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr15 <- data.frame(event = "intTrans", 
                    time = rep(1350, length(obr15_nodes)),
                    node= obr15_nodes, 
                    dest = 0, n = 0, 
                    proportion = sia_coverage,
                    select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(1352, length(obr15_nodes)),
                                node= obr15_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(1351, length(obr15_nodes)),
                               node= obr15_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(1351, length(obr15_nodes)),
                               node= obr15_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(1351, length(obr15_nodes)),
                               node= obr15_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(1351, length(obr15_nodes)),
                               node= obr15_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(1352, length(obr15_nodes)),
                      node= obr15_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame() 

# now, schedule the oSIA on day 90 and allow model to run for an additional 89 days
create.StocSIRobr15 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 1350, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia15,
                               obr15, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr15 <- create.StocSIRobr15(I0, beta, gamma, caseinfection, 1439) 
out1439 <- run(model = SIRmodelobr15) 
df1439 <- trajectory(out1439) # df for the second 89 days simulated, check and see which nodes had cases
obr_df <- subset(df1439, df1439$C>=1)
obr16_nodes <- sort(unique(obr_df$node)) # nodes with 1 or more cases
end_time <- 1439
end_dat <- df1439

#############################################
########### Outbreak response 16 ############
#############################################
############# days 1440 to 1529 #############
#############################################
importation_sia16 <-filter(importation_sia, time>=1440 & time<=1529)
psia5 <-filter(psia_annual, time>=1440 & time<=1529)
psia5_no_seroconversion <-filter(psia_no_seroconversion, time>=1440 & time<=1529)
psia5_seroconversionS1 <-filter(psia_seroconversionS1, time>=1440 & time<=1529)
psia5_seroconversionS2 <-filter(psia_seroconversionS2, time>=1440 & time<=1529)
psia5_seroconversionS3 <-filter(psia_seroconversionS3, time>=1440 & time<=1529)
psia5_seroconversionS4 <-filter(psia_seroconversionS4, time>=1440 & time<=1529)
psia5_seroMAX <-filter(psia_seroMAX, time>=1440 & time<=1529)


birthsS <- data.frame(event = "enter", 
                      time = c(rep(1440:1529, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(1440:1529, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(1440:1529, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(1440:1529, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(1440:1529, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr16 <- data.frame(event = "intTrans", 
                    time = rep(1440, length(obr16_nodes)),
                    node= obr16_nodes, 
                    dest = 0, n = 0, 
                    proportion = sia_coverage,
                    select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(1442, length(obr16_nodes)),
                                node= obr16_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(1441, length(obr16_nodes)),
                               node= obr16_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(1441, length(obr16_nodes)),
                               node= obr16_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(1441, length(obr16_nodes)),
                               node= obr16_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(1441, length(obr16_nodes)),
                               node= obr16_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(1442, length(obr16_nodes)),
                      node= obr16_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame() 

# now, schedule the oSIA on day 90 and allow model to run for an additional 89 days
create.StocSIRobr16 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 1440, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia16,
                               obr16, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               psia5,
                               psia5_no_seroconversion, psia5_seroconversionS1,
                               psia5_seroconversionS2, psia5_seroconversionS3, 
                               psia5_seroconversionS4, psia5_seroMAX, 
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr16 <- create.StocSIRobr16(I0, beta, gamma, caseinfection, 1529) 
out1529 <- run(model = SIRmodelobr16) 
df1529 <- trajectory(out1529) 
obr_df <- subset(df1529, df1529$C>=1)
obr17_nodes <- sort(unique(obr_df$node)) # nodes with 1 or more cases
end_time <- 1529
end_dat <- df1529

#############################################
########### Outbreak response 17 ############
#############################################
############# days 1530 to 1619 #############
#############################################
importation_sia17 <-filter(importation_sia, time>=1530 & time<=1619)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(1530:1619, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(1530:1619, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(1530:1619, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(1530:1619, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(1530:1619, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr17 <- data.frame(event = "intTrans", 
                    time = rep(1530, length(obr17_nodes)),
                    node= obr17_nodes, 
                    dest = 0, n = 0, 
                    proportion = sia_coverage,
                    select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(1532, length(obr17_nodes)),
                                node= obr17_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(1531, length(obr17_nodes)),
                               node= obr17_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(1531, length(obr17_nodes)),
                               node= obr17_nodes,  
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(1531, length(obr17_nodes)),
                               node= obr17_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(1531, length(obr17_nodes)),
                               node= obr17_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(1532, length(obr17_nodes)),
                      node= obr17_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame() 

create.StocSIRobr17 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 1530, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia17,
                               obr17, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr17 <- create.StocSIRobr17(I0, beta, gamma, caseinfection, 1619) 
out1619 <- run(model = SIRmodelobr17) 
df1619 <- trajectory(out1619) 
obr_df <- subset(df1619, df1619$C>=1)
obr18_nodes <- sort(unique(obr_df$node)) # nodes with 1 or more cases
end_time <- 1619
end_dat <- df1619

#############################################
########### Outbreak response 18 ############
#############################################
############# days 1620 to 1709 #############
#############################################
importation_sia18 <-filter(importation_sia, time>=1620 & time<=1709)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(1620:1709, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(1620:1709, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(1620:1709, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(1620:1709, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(1620:1709, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr18 <- data.frame(event = "intTrans", 
                    time = rep(1620, length(obr18_nodes)),
                    node= obr18_nodes, 
                    dest = 0, n = 0, 
                    proportion = sia_coverage,
                    select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(1622, length(obr18_nodes)),
                                node= obr18_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(1621, length(obr18_nodes)),
                               node= obr18_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(1621, length(obr18_nodes)),
                               node= obr18_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(1621, length(obr18_nodes)),
                               node= obr18_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(1621, length(obr18_nodes)),
                               node= obr18_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(1622, length(obr18_nodes)),
                      node= obr18_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame() 

create.StocSIRobr18 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 1620, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia18,
                               obr18, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr18 <- create.StocSIRobr18(I0, beta, gamma, caseinfection, 1709) 
out1709 <- run(model = SIRmodelobr18) 
df1709 <- trajectory(out1709) 
obr_df <- subset(df1709, df1709$C>=1)
obr19_nodes <- sort(unique(obr_df$node)) # nodes with 1 or more cases
end_time <- 1709
end_dat <- df1709

#############################################
########### Outbreak response 19 ############
#############################################
############# days 1710 to 1799 #############
#############################################
importation_sia19 <-filter(importation_sia, time>=1710 & time<=1799)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(1710:1799, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(1710:1799, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(1710:1799, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

ipv <- data.frame(event = "intTrans", 
                  time = c(rep(1710:1799, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(1710:1799, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr19 <- data.frame(event = "intTrans", 
                    time = rep(1710, length(obr19_nodes)),
                    node= obr19_nodes, 
                    dest = 0, n = 0, 
                    proportion = sia_coverage,
                    select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(1712, length(obr19_nodes)),
                                node= obr19_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(1711, length(obr19_nodes)),
                               node= obr19_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(1711, length(obr19_nodes)),
                               node= obr19_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(1711, length(obr19_nodes)),
                               node= obr19_nodes,  
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(1711, length(obr19_nodes)),
                               node= obr19_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(1712, length(obr19_nodes)),
                      node= obr19_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame() 

create.StocSIRobr19 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 1710, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia19,
                               obr19, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr19 <- create.StocSIRobr19(I0, beta, gamma, caseinfection, 1799) 
out1799 <- run(model = SIRmodelobr19) 
df1799 <- trajectory(out1799) 
obr_df <- subset(df1799, df1799$C>=1)
obr20_nodes <- sort(unique(obr_df$node)) 
end_time <- 1799
end_dat <- df1799

#############################################
########### Outbreak response 20 ############
#############################################
############# days 1800 to 1825 #############
#############################################
importation_sia20 <-filter(importation_sia, time>=1800 & time<=1825)

birthsS <- data.frame(event = "enter", 
                      time = c(rep(1800:1825, each=length(sims))),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = c(rep(1800:1825, each=length(sims))),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = c(rep(1800:1825, each=length(sims))),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)


ipv <- data.frame(event = "intTrans", 
                  time = c(rep(1800:1825, each=length(sims))),
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = ipv_coverage, 
                  select = 14, shift = 5) 

seroconversionS3r <- data.frame(event = "intTrans", 
                                time = c(rep(1800:1825, each=length(sims))),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3), # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

obr20 <- data.frame(event = "intTrans", 
                    time = rep(1800, length(obr20_nodes)),
                    node= obr20_nodes, 
                    dest = 0, n = 0, 
                    proportion = sia_coverage,
                    select = 5, shift = 1)

no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(1802, length(obr20_nodes)),
                                node= obr20_nodes, 
                                dest = 0,
                                n = 0, 
                                proportion = 1, 
                                select = 6, shift = 2)

seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(1801, length(obr20_nodes)),
                               node= obr20_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1) ,
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(1801, length(obr20_nodes)),
                               node= obr20_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) ,
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(1801, length(obr20_nodes)),
                               node= obr20_nodes,  
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , 
                               select = 9, shift = 3)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(1801, length(obr20_nodes)),
                               node= obr20_nodes, 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

seroMAX <- data.frame(event = "intTrans", 
                      time = rep(1802, length(obr20_nodes)),
                      node= obr20_nodes, 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

# specify the new initial conditions
new_initial <- filter(end_dat, time==end_time) %>% select(S0, S0v, S1, S1v, S2, S2v, S3, S3v, S3r, S4, Rv, Ii, I, C, Rn, Ccum) %>% as.data.frame()  

# now, schedule the oSIA on day 90 and allow model to run for an additional 89 days
create.StocSIRobr20 <- function(I0, beta, gamma, rho, f_time){
  initial_state <- new_initial
  compartments <- c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v", "S3r", "S4", "Rv", "Ii", "I", "C", "Rn", "Ccum")
  transitions <- c("S0 -> (1-rho)*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S0 -> rho*(beta*S0*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S1 -> (1-rho)*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S1 -> rho*(beta*S1*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S2 -> (1-rho)*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S2 -> rho*(beta*S2*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3 -> (1-rho)*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3 -> rho*(beta*S3*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "S3r -> (1-rho)*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> I",
                   "S3r -> rho*(beta*S3r*(I+C))/(S0+S0v+S1+S1v+S2+S2v+S3+S3v+S3r+Rv+I+C+Rn) -> C + Ccum",
                   "I -> gamma*I -> Rn",
                   "C -> gamma*C -> Rn") 
  tspan <- seq(from = 1800, to = f_time, by = 1) 
  return(mparse(transitions = transitions, 
                compartments = compartments, 
                events = rbind(birthsS, birthsS3, deaths, importation_sia20,
                               obr20, seroconversionS4, seroMAX,
                               seroconversionS1, seroconversionS2, 
                               seroconversionS3, no_seroconversion, 
                               ipv, seroconversionS3r), 
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}

# outputs next interval of 89 days
SIRmodelobr20 <- create.StocSIRobr20(I0, beta, gamma, caseinfection, 1825) 
out1825 <- run(model = SIRmodelobr20) 
df1825 <- trajectory(out1825) 
obr_df <- subset(df1825, df1825$C>=1)
obr21_nodes <- sort(unique(obr_df$node)) 
end_time <- df1825
end_dat <- df1825

#################################################
# When outbreaks end, put all the data together #
#################################################
#################################################
out_psia <- rbind(df89, df179, df269, df359, df449, 
                  df539, df629, df719, df809, df899, 
                  df989, df1079, df1169, df1259, 
                  df1349, df1439, df1529, df1619, 
                  df1709, df1799, df1825)

#write.csv(out_psia, file="psia_osia_ri25_2imports_per_year.csv")

x <- out_psia %>% group_by(node) %>% dplyr::summarize(total_cases = max(Ccum)) 
nonzero <- subset(x, x$total_cases>0)
mean(nonzero$total_cases)
(nrow(nonzero)/1000)*100 #probability of an outbreak

## nodes needing oSIAs
ll <- list(obr1_nodes, obr2_nodes, obr3_nodes, obr4_nodes, obr5_nodes, obr6_nodes, obr7_nodes, 
           obr8_nodes, obr9_nodes, obr10_nodes, obr11_nodes, obr12_nodes, obr13_nodes, obr14_nodes, 
           obr15_nodes, obr16_nodes, obr17_nodes, obr18_nodes, obr19_nodes, obr20_nodes)
nodes <- lapply(ll, function(x) x[1: max(sapply(ll, length))]) %>% do.call(cbind, .) 

colnames(nodes) <- c("obr1_nodes", "obr2_nodes", "obr3_nodes", "obr4_nodes", "obr5_nodes",
                     "obr6_nodes", "obr7_nodes", "obr8_nodes", "obr9_nodes", "obr10_nodes",
                     "obr11_nodes", "obr12_nodes", "obr13_nodes", "obr14_nodes", "obr15_nodes",
                     "obr16_nodes", "obr17_nodes", "obr18_nodes", "obr19_nodes", "obr20_nodes")
nodes <- as.data.frame(nodes)
all_nodes <- gather(nodes, obr, sims, obr1_nodes:obr20_nodes, factor_key=TRUE)
all_nodes$obr_num <- as.numeric(all_nodes$obr)
all_nodes$year <- ""
all_nodes$year[all_nodes$obr_num<=3]<- 1
all_nodes$year[all_nodes$obr_num>=4 & all_nodes$obr_num<=7]<- 2
all_nodes$year[all_nodes$obr_num>=8 & all_nodes$obr_num<=12]<- 3
all_nodes$year[all_nodes$obr_num>=13 & all_nodes$obr_num<=16]<- 4
all_nodes$year[all_nodes$obr_num>=17]<- 5

freq_nodes <- as.data.frame(table(sims=all_nodes$sims, year=all_nodes$year))
#summary(freq_nodes$Freq) #only for counting total outbreaks over entire time horizon
#write.csv(freq_nodes, file="freq_nodes_psia_osia_ri25.csv")
#eng
