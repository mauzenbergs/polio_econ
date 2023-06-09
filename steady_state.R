# This code achieves a steady state before introducing any SIAs or the first infection. 
# After steady state is reached, model outputs from the steady state model are used as 
# initial conditions in subsequent models with SIAs
# The pSIA model continues in the code "psia_model.R"

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
ve <- 0.5 # vaccine efficacy for bOPV
doses <- 3 
seroconversion <- (1 -(1-ve)^doses)
birth_rate <- 4000 # or 5x10^-4 (4000 births per 8 Million) = avg across africa
death_rate = birth_rate
births_S <- round(birth_rate * (1-ri_coverage), digits=0)
births_S3 <- round(birth_rate * ri_coverage, digits=0)

## First, allowing the model to reach steady state
# Specifying data for birth and death compartments 
f_time <- (50*365)-1 #time horizon for steady state
birthsS <- data.frame(event = "enter", 
                      time = rep(1:f_time, each=length(sims)),
                      node= 1:length(sims), 
                      dest = 0, 
                      n = births_S,
                      proportion = 0, 
                      select = 1, shift = 0)

birthsS3 <- data.frame(event = "enter", 
                       time = rep(1:f_time, each=length(sims)),
                       node= 1:length(sims), 
                       dest = 0, 
                       n = births_S3,
                       proportion = 0, 
                       select = 2, shift = 0)

deaths <- data.frame(event = "exit", 
                     time = rep(1:f_time, each=length(sims)),
                     node= 1:length(sims), 
                     dest = 0, 
                     n = death_rate,
                     proportion = 0, 
                     select = 3, shift = 0)

vaccination <- data.frame(event = "intTrans", 
                          #campaigns introduced at year 30 to reach steady state
                          time = rep(c(10950, 11315, 11680, 12045, 12410, 
                                       12775, 13140, 13505, 13870, 14235, 
                                       14600, 14965, 15330, 15695, 16060, 
                                       16425, 16790, 17155, 17520, 17885), each=length(sims)),
                          node= 1:length(sims), 
                          dest = 0,
                          n = 0,
                          proportion = sia_coverage,
                          select = 5, shift = 1)

#can't have two internal transfer events on same day in same compartment, so on SIA + 2 days, 
#those who did not seroconvert, move up to next tier
no_seroconversion <- data.frame(event = "intTrans", 
                                time = rep(c(10952, 11317, 11682, 12047, 12412, 
                                             12777, 13142, 13507, 13872, 14237, 
                                             14602, 14967, 15332, 15697, 16062, 
                                             16427, 16792, 17157, 17522, 17887), each=length(sims)),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = 1, # move to next S tier(Sn+1)
                                select = 6, shift = 2)

#can't have two internal transfer events on same day in same compartment, so on SIA + 1 day, 
#seroconversion happens depending on VE and doses
seroconversionS1 <- data.frame(event = "intTrans", 
                               time = rep(c(10951, 11316, 11681, 12046, 12411, 
                                            12776, 13141, 13506, 13871, 14236, 
                                            14601, 14966, 15331, 15696, 16061, 
                                            16426, 16791, 17156, 17521, 17886), each=length(sims)),
                               node= 1:length(sims), 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^1), # seroconvert and move to Rv
                               select = 7, shift = 3)

seroconversionS2 <- data.frame(event = "intTrans", 
                               time = rep(c(10951, 11316, 11681, 12046, 12411, 
                                            12776, 13141, 13506, 13871, 14236, 
                                            14601, 14966, 15331, 15696, 16061, 
                                            16426, 16791, 17156, 17521, 17886), each=length(sims)),
                               node= 1:length(sims), 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^2) , # another dose worth of protection, seroconvert and move to Rv
                               select = 8, shift = 3)

seroconversionS3 <- data.frame(event = "intTrans", 
                               time = rep(c(10951, 11316, 11681, 12046, 12411, 
                                            12776, 13141, 13506, 13871, 14236, 
                                            14601, 14966, 15331, 15696, 16061, 
                                            16426, 16791, 17156, 17521, 17886), each=length(sims)),
                               node= 1:length(sims), 
                               dest = 0,
                               n = 0, 
                               proportion = (1 -(1-ve)^3) , # another dose worth of protection, seroconvert and move to Rv
                               select = 9, shift = 3)


seroconversionS3r <- data.frame(event = "intTrans", 
                                time = rep((50*365):f_time, each=length(sims)),
                                node= 1:length(sims), 
                                dest = 0,
                                n = 0, 
                                proportion = (1 -(1-ve)^3) , # another dose worth of protection, seroconvert and move to Rv
                                select = 11, shift = 5)

seroconversionS4 <- data.frame(event = "intTrans", 
                               time = rep(c(10951, 11316, 11681, 12046, 12411, 
                                            12776, 13141, 13506, 13871, 14236, 
                                            14601, 14966, 15331, 15696, 16061, 
                                            16426, 16791, 17156, 17521, 17886), each=length(sims)),
                               node= 1:length(sims), 
                               dest = 0,
                               n = 0, 
                               proportion = 1, # assume 100% move to S4 after 4th dose
                               select = 12, shift = 3)

# a bit of an unnecessary internal transfer event as we assume 100% immune after 4 doses, 
# but this just moves individuals from S4 to Rv to keep track of movement
seroMAX <- data.frame(event = "intTrans", 
                      time = rep(c(10952, 11317, 11682, 12047, 12412, 
                                   12777, 13142, 13507, 13872, 14237, 
                                   14602, 14967, 15332, 15697, 16062, 
                                   16427, 16792, 17157, 17522, 17887), each=length(sims)),
                      node= 1:length(sims), 
                      dest = 0,
                      n = 0, 
                      proportion = 1, # assume 100% seroconvert after 4 bOPV doses
                      select = 13, shift = 3)

ipv <- data.frame(event = "intTrans", 
                  time = rep((50*365):f_time, each=length(sims)), #get to steady state without IPV
                  node= 1:length(sims), 
                  dest = 0, 
                  n = 0,
                  proportion = 1, 
                  select = 14, shift = 5)

# "Select" matrix
Et <- matrix(c(1, rep(0,15), #births R0
               rep(0,8),1,rep(0,7), #births S3
               rep(1,11),0,1,1,1,0, #deaths
               1, rep(0,15), #dynamic importation
               1,0,1,0,1,0,1,0,1,rep(0,7), #vaccination
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
               0,1,0,1,0,1,0,2,0,rep(0,7), #NO seroconversion
               0,9,0,7,0,5,0,3,0,1,rep(0,6), #seroconversion
               12, rep(0,15), #importation
               rep(0,8),2,rep(0,7)), #ipv vaccination
             nrow = 16, ncol = 5,
             dimnames = list(c("S0", "S0v", "S1", "S1v", "S2", "S2v", "S3", "S3v","S3r","S4", "Rv", "Ii", "I", "C", "Rn", "Ccum"), 
                             c("vaccination", "no_seroconversion", "seroconversion", "importation", "ipv")))


# initial state assumptions
u0psia <- data.frame(S0=rep(pop*(1-ri_coverage)-I0,reps), 
                     S0v=rep(0, reps),
                     S1=rep(0, reps),
                     S1v=rep(0, reps),
                     S2=rep(0, reps),
                     S2v=rep(0, reps),
                     S3=rep(0, reps),
                     S3v=rep(0, reps),
                     S3r=rep(pop*(ri_coverage), reps),
                     S4=rep(0, reps),
                     Rv=rep(0, reps),
                     Ii=rep(0, reps),
                     I=rep(I0,reps),
                     C=rep(0, reps),
                     Rn=rep(0, reps),
                     Ccum=rep(0, reps))


create.StocSIRpsia <- function(I0, beta, gamma, rho, f_time){
  initial_state <- u0psia
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
                events = rbind(birthsS, birthsS3, deaths, 
                               vaccination, 
                               seroconversionS1, seroconversionS2,
                               seroconversionS3, no_seroconversion, 
                               seroconversionS3r, 
                               ipv, 
                               seroconversionS4, 
                               seroMAX),
                E = Et, N = Nt,
                gdata = c(beta=beta, gamma=gamma, rho=rho), 
                u0 = initial_state, 
                tspan = tspan))
}


SIRmodelpsia <- create.StocSIRpsia(I0, beta, gamma, caseinfection, f_time)
out_psia_model <- run(model = SIRmodelpsia) 
out_psia <- trajectory(out_psia_model)

# IMPORTANT!! save the last day of the steady state pSIA model for initial conditions 
u0psia_osia <- out_psia %>% filter(time==((50*365)-1)) #one day before the first infection is introduced
u0psia_osia <- u0psia_osia[ -c(1:2) ]
write.csv(u0psia_osia, file="u0psia_ri25.csv")

u0osia_osia <- out_psia %>% filter(time==((30*365)-1)) #one day before the first SIA is introduced
u0osia_osia <- u0osia_osia[ -c(1:2) ]
write.csv(u0osia_osia, file="u0osia_ri25.csv")
#end
