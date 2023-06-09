################################################################
### Processing model outputs into burden estimates and costs ###
################################################################

# The data formatted in this script will be used to estimate economic outputs
# All economic outputs are done for health system and GPEI perspective
# this script will require the files "psia_all.csv" "psia2_all.csv" and "osia_all.csv"

library(tidyverse)
library(dplyr)
library(data.table)
library(Rmisc)
library(ggallin)

# Defining cost parameters and assumptions
pop <- 8000000 #rough average of target population across AFRO countries
birth_rate <- 4000 #enter the population every day
osia_cost_per_child <- 1.0 #admin, social mob, procurement, data from GPEI
psia_cost_per_child <- 0.50 #admin, social mob, procurement, data from GPEI
discount_rate3 <- 0.03
discount_rate0 <- 0
discount_factor <- data.table(year = c(1:5)) %>% 
  dplyr::mutate(discount_factor3 = (1+discount_rate3)^(year-1),
                discount_factor0= (1+discount_rate0)^(year-1))

bopv <- 0.19 # cost per dose in 10 dose vile, UNICEF 2023
ipv <- 2.00 # UNICEF 2023 cost in 10 dose vile

## wastage assumptions can be varied 
opv_wastage_ri <- 1.15 #range 10-20%
opv_wastage_sia <- 1.10 #range 5-15%
ipv_wastage <- 1.15 #range 5-20%, only given via RI
## wastage high assumptions ##
#opv_wastage_ri <- 1.20 #range 10-20%
#opv_wastage_sia <- 1.15 #range 5-15%
#ipv_wastage <- 1.20 #range 5-20%, only given via RI
## wastage low assumptions ##
#opv_wastage_ri <- 1.10 #range 10-20%
#opv_wastage_sia <- 1.05 #range 5-15%
#ipv_wastage <- 1.05 #range 5-20%, only given via RI

doses <- 3 #bOPV doses
case_cost <- 700 #AFP = VAPP
ri_opv_admin <- 0.95 * ((1+discount_rate3)^(2023-2019)) #admin cost in USD$2019, so apply 3% discounting per year
ri_ipv_admin <- 1.78 * ((1+discount_rate3)^(2023-2019)) #admin cost in USD$2019, so apply 3% discounting per year
daly_rate <- 14
psia_num <- 5 #annual campaigns
psia2_num <- 3 #reduced frequency, t=0, t=2 years, t=4 years
sia_coverage <- 0.25

###############################################################################
############################ pSIA scenario ####################################
###############################################################################
psia_all_by_year <- read.csv(file="psia_all.csv") 

psia_all_by_year <- psia_all_by_year %>% 
  dplyr::mutate(births = birth_rate * 365,
                births_vax = births * (ri_coverage/100), #births per year actually vaccinated
                ri_bopv_doses = births_vax * doses, #bOPV doses actually administered via RI per year 
                vapp = vapp,
                dalys = (total_cases + vapp) * daly_rate)

psia_all_by_year <- merge(psia_all_by_year, discount_factor, by="year") #discounting per year

### HEALTH CARE SYSTEM PERSPECTIVE ###
psia_costs_hc <- psia_all_by_year %>% 
  dplyr::mutate(costs =
                  # treatment costs
                  (total_cases * case_cost) + # cost per AFP case
                  (vapp * case_cost) + # cost per VAPP case
                  # vaccine and admin costs depend on births
                  (births * opv_wastage_ri * bopv * doses * ri_opv_admin),
                
                # discount factor
                costs_discounted = (costs / discount_factor3),
                dalys_discounted0 = dalys / discount_factor0,
                dalys_discounted3 = round((dalys / discount_factor3),0)) 

## GPEI PERSPECTIVE ##
psia_costs_gpei <- psia_all_by_year %>% 
  dplyr::mutate(costs = (pop * opv_wastage_sia * bopv * psia_num * psia_cost_per_child) +
                  (pop * opv_wastage_sia * bopv * outbreaks * osia_cost_per_child) +
                  (births * ipv_wastage * ipv * ri_ipv_admin),
                costs_discounted = (costs / discount_factor3),
                dalys_discounted0 = dalys / discount_factor0,
                dalys_discounted3 = round((dalys / discount_factor3),0)) 


## COST EFFECTIVENESS ANALYSIS - COST PER DALY AVERTED ##
### Costs discounted ###
psia_costs_hc_by_year <-psia_costs_hc %>% 
  group_by(ri_coverage, year) %>%
  dplyr::summarise(scenario = "pSIA",
                   cases = mean(total_cases), 
                   vapp = mean(vapp),
                   dalys = mean(dalys),
                   costs = mean(costs),
                   costs_discounted = mean(costs_discounted),
                   dalys_discounted0 = mean(dalys_discounted0),
                   dalys_discounted3 = mean(dalys_discounted3)) %>%
  dplyr::mutate_if(is.numeric, round, 0)

psia_costs_gpei_by_year <-psia_costs_gpei %>% 
  group_by(ri_coverage, year) %>%
  dplyr::summarise(scenario = "pSIA",
                   cases = mean(total_cases), 
                   vapp = mean(vapp),
                   dalys = mean(dalys),
                   costs = mean(costs),
                   costs_discounted = mean(costs_discounted),
                   dalys_discounted0 = mean(dalys_discounted0),
                   dalys_discounted3 = mean(dalys_discounted3)) %>%
  dplyr::mutate_if(is.numeric, round, 0)


## BUDGET IMPACT ANALYSIS ##
# aggregating costs across 5 years, don't need discounting
psia_costs_hc_5_years <-psia_costs_hc %>% 
  group_by(ri_coverage, node) %>%
  dplyr::summarise(scenario = "pSIA", #summing across all years
                   cases = sum(total_cases), 
                   vapp = sum(vapp),
                   dalys = sum(dalys),
                   costs = sum(costs))

psia_costs_hc_5_years <-psia_costs_hc_5_years %>%
  group_by(ri_coverage) %>%
  dplyr::summarise(scenario = "pSIA",
                   cases_mean = mean(cases),
                   cases_lwr = CI(cases)[3],
                   cases_upr = CI(cases)[1],
                   vapp_mean = mean(vapp),
                   dalys_mean = mean(dalys),
                   costs_mean = mean(costs),
                   costs_lwr = CI(costs)[3],
                   costs_upr = CI(costs)[1]) %>%
  dplyr::mutate_if(is.numeric, round, 0)

psia_costs_gpei_5_years <-psia_costs_gpei %>% 
  group_by(ri_coverage, node) %>%
  dplyr::summarise(scenario = "pSIA", #summing across all years
                   cases = sum(total_cases), 
                   vapp = sum(vapp),
                   dalys = sum(dalys),
                   costs = sum(costs))

psia_costs_gpei_5_years <-psia_costs_gpei_5_years %>%
  group_by(ri_coverage) %>%
  dplyr::summarise(scenario = "pSIA",
                   cases_mean = mean(cases),
                   cases_lwr = CI(cases)[3],
                   cases_upr = CI(cases)[1],
                   vapp_mean = mean(vapp),
                   dalys_mean = mean(dalys),
                   costs_mean = mean(costs),
                   costs_lwr = CI(costs)[3],
                   costs_upr = CI(costs)[1]) %>%
  dplyr::mutate_if(is.numeric, round, 0)


###############################################################################
############################ oSIA scenario ####################################
###############################################################################
osia_all_by_year <- read.csv(file="osia_all.csv") 

osia_all_by_year <- osia_all_by_year %>% 
  dplyr::mutate(births = birth_rate * 365,
                births_vax = births * (ri_coverage/100), #births per year actually vaccinated
                ri_bopv_doses = births_vax * doses, #bOPV doses actually administered via RI per year
                vapp = vapp, 
                dalys = (total_cases + vapp) * daly_rate)

osia_all_by_year <- merge(osia_all_by_year, discount_factor, by="year") #discounting per year

### HEALTH CARE SYSTEM PERSPECTIVE ###
osia_costs_hc <- osia_all_by_year %>% 
  dplyr::mutate(costs =
                  # treatment costs
                  (total_cases * case_cost) + # cost per AFP case
                  (vapp * case_cost) + # cost per VAPP case
                  # vaccine and admin costs depend on births
                  (births * opv_wastage_ri * bopv * doses * ri_opv_admin),
                
                ## discount factor
                costs_discounted = (costs / discount_factor3),
                dalys_discounted0 = dalys / discount_factor0,
                dalys_discounted3 = round((dalys / discount_factor3),0)) 

## GPEI PERSPECTIVE ##
osia_costs_gpei <- osia_all_by_year %>% 
  dplyr::mutate(costs = (pop * opv_wastage_sia * bopv * outbreaks * osia_cost_per_child) +
                  (births * ipv_wastage * ipv * ri_ipv_admin),
                costs_discounted = (costs / discount_factor3),
                dalys_discounted0 = dalys / discount_factor0,
                dalys_discounted3 = round((dalys / discount_factor3),0)) 



## COST EFFECTIVENESS ANALYSIS - COST PER DALY AVERTED ##
### Costs discounted ###
osia_costs_hc_by_year <-osia_costs_hc %>% 
  group_by(ri_coverage, year) %>%
  dplyr::summarise(scenario = "oSIA",
                   cases = mean(total_cases), 
                   vapp = mean(vapp),
                   dalys = mean(dalys),
                   costs = mean(costs),
                   costs_discounted = mean(costs_discounted),
                   dalys_discounted0 = mean(dalys_discounted0),
                   dalys_discounted3 = mean(dalys_discounted3)) %>%
  dplyr::mutate_if(is.numeric, round, 0)

osia_costs_gpei_by_year <-osia_costs_gpei %>% 
  group_by(ri_coverage, year) %>%
  dplyr::summarise(scenario = "oSIA",
                   cases = mean(total_cases), 
                   vapp = mean(vapp),
                   dalys = mean(dalys),
                   costs = mean(costs),
                   costs_discounted = mean(costs_discounted),
                   dalys_discounted0 = mean(dalys_discounted0),
                   dalys_discounted3 = mean(dalys_discounted3)) %>%
  dplyr::mutate_if(is.numeric, round, 0)


## BUDGET IMPACT ANALYSIS ##
# aggregating costs across 5 years, don't need discounting
osia_costs_hc_5_years <-osia_costs_hc %>% 
  group_by(ri_coverage, node) %>%
  dplyr::summarise(scenario = "oSIA", #summing across all years
                   cases = sum(total_cases), 
                   vapp = sum(vapp),
                   dalys = sum(dalys),
                   outbreaks = sum(outbreaks),
                   costs = sum(costs))

osia_costs_hc_5_years <-osia_costs_hc_5_years %>%
  group_by(ri_coverage) %>%
  dplyr::summarise(scenario = "oSIA",
                   cases_mean = mean(cases),
                   cases_lwr = CI(cases)[3],
                   cases_upr = CI(cases)[1],
                   vapp_mean = mean(vapp),
                   dalys_mean = mean(dalys),
                   costs_mean = mean(costs),
                   costs_lwr = CI(costs)[3],
                   costs_upr = CI(costs)[1]) %>%
  dplyr::mutate_if(is.numeric, round, 0)

osia_costs_gpei_5_years <-osia_costs_gpei %>% 
  group_by(ri_coverage, node) %>%
  dplyr::summarise(scenario = "oSIA", #summing across all years
                   cases = sum(total_cases), 
                   vapp = sum(vapp),
                   dalys = sum(dalys),
                   costs = sum(costs))

osia_costs_gpei_5_years <-osia_costs_gpei_5_years %>%
  group_by(ri_coverage) %>%
  dplyr::summarise(scenario = "oSIA",
                   cases_mean = mean(cases),
                   cases_lwr = CI(cases)[3],
                   cases_upr = CI(cases)[1],
                   vapp_mean = mean(vapp),
                   dalys_mean = mean(dalys),
                   costs_mean = mean(costs),
                   costs_lwr = CI(costs)[3],
                   costs_upr = CI(costs)[1]) %>%
  dplyr::mutate_if(is.numeric, round, 0)


###############################################################################
############################ pSIA2 scenario ###################################
###############################################################################
psia2_all_by_year <- read.csv(file="psia2_all.csv") 

psia2_all_by_year <- psia2_all_by_year %>% 
  dplyr::mutate(births = birth_rate * 365,
                births_vax = births * (ri_coverage/100), #births per year actually vaccinated
                ri_bopv_doses = births_vax * doses, #bOPV doses actually administered via RI per year 
                vapp = vapp,
                dalys = (total_cases + vapp) * daly_rate)

psia2_all_by_year <- merge(psia2_all_by_year, discount_factor, by="year") #discounting per year

### HEALTH CARE SYSTEM PERSPECTIVE ###
psia2_costs_hc <- psia2_all_by_year %>% 
  dplyr::mutate(costs =
                  # treatment costs
                  (total_cases * case_cost) + # cost per AFP case
                  (vapp * case_cost) + # cost per VAPP case
                  # vaccine and admin costs depend on births
                  (births * opv_wastage_ri * bopv * doses * ri_opv_admin),
                
                # discount factor
                costs_discounted = (costs / discount_factor3),
                dalys_discounted0 = dalys / discount_factor0,
                dalys_discounted3 = round((dalys / discount_factor3),0)) 

## GPEI PERSPECTIVE ##
psia2_costs_gpei <- psia2_all_by_year %>% 
  dplyr::mutate(costs = (pop * opv_wastage_sia * bopv * psia2_num * psia_cost_per_child) +
                  (pop * opv_wastage_sia * bopv * outbreaks * osia_cost_per_child) +
                  (births * ipv_wastage * ipv * ri_ipv_admin),
                costs_discounted = (costs / discount_factor3),
                dalys_discounted0 = dalys / discount_factor0,
                dalys_discounted3 = round((dalys / discount_factor3),0)) 


## COST EFFECTIVENESS ANALYSIS - COST PER DALY AVERTED ##
### Costs discounted ###
psia2_costs_hc_by_year <-psia2_costs_hc %>% 
  group_by(ri_coverage, year) %>%
  dplyr::summarise(scenario = "pSIA2",
                   cases = mean(total_cases), 
                   vapp = mean(vapp),
                   dalys = mean(dalys),
                   costs = mean(costs),
                   costs_discounted = mean(costs_discounted),
                   dalys_discounted0 = mean(dalys_discounted0),
                   dalys_discounted3 = mean(dalys_discounted3)) %>%
  dplyr::mutate_if(is.numeric, round, 0)

psia2_costs_gpei_by_year <-psia2_costs_gpei %>% 
  group_by(ri_coverage, year) %>%
  dplyr::summarise(scenario = "pSIA2",
                   cases = mean(total_cases), 
                   vapp = mean(vapp),
                   dalys = mean(dalys),
                   costs = mean(costs),
                   costs_discounted = mean(costs_discounted),
                   dalys_discounted0 = mean(dalys_discounted0),
                   dalys_discounted3 = mean(dalys_discounted3)) %>%
  dplyr::mutate_if(is.numeric, round, 0)


## BUDGET IMPACT ANALYSIS ##
# aggregating costs across 5 years, don't need discounting
psia2_costs_hc_5_years <-psia2_costs_hc %>% 
  group_by(ri_coverage, node) %>%
  dplyr::summarise(scenario = "pSIA2", #summing across all years
                   cases = sum(total_cases), 
                   vapp = sum(vapp),
                   dalys = sum(dalys),
                   costs = sum(costs))

psia2_costs_hc_5_years <-psia2_costs_hc_5_years %>%
  group_by(ri_coverage) %>%
  dplyr::summarise(scenario = "pSIA2",
                   cases_mean = mean(cases),
                   cases_lwr = CI(cases)[3],
                   cases_upr = CI(cases)[1],
                   vapp_mean = mean(vapp),
                   dalys_mean = mean(dalys),
                   costs_mean = mean(costs),
                   costs_lwr = CI(costs)[3],
                   costs_upr = CI(costs)[1]) %>%
  dplyr::mutate_if(is.numeric, round, 0)

psia2_costs_gpei_5_years <-psia2_costs_gpei %>% 
  group_by(ri_coverage, node) %>%
  dplyr::summarise(scenario = "pSIA2", #summing across all years
                   cases = sum(total_cases), 
                   vapp = sum(vapp),
                   dalys = sum(dalys),
                   costs = sum(costs))

psia2_costs_gpei_5_years <-psia2_costs_gpei_5_years %>%
  group_by(ri_coverage) %>%
  dplyr::summarise(scenario = "pSIA2",
                   cases_mean = mean(cases),
                   cases_lwr = CI(cases)[3],
                   cases_upr = CI(cases)[1],
                   vapp_mean = mean(vapp),
                   dalys_mean = mean(dalys),
                   costs_mean = mean(costs),
                   costs_lwr = CI(costs)[3],
                   costs_upr = CI(costs)[1]) %>%
  dplyr::mutate_if(is.numeric, round, 0)

###############################################################################
###############################################################################
###############################################################################

## TOTAL ESTIMATED COSTS OVER 5 YEARS ##
budget_impact_hc <- rbind(psia_costs_hc_5_years, psia2_costs_hc_5_years, osia_costs_hc_5_years)
budget_impact_gpei <- rbind(psia_costs_gpei_5_years, psia2_costs_gpei_5_years, osia_costs_gpei_5_years)


###############################################################################
#########  INCREMENTAL COSTS - oSIA scenario as baseline comparator ###########
###############################################################################

names(psia_costs_hc)[names(psia_costs_hc) == "outbreaks"] <- "outbreaks_sia"
names(psia_costs_gpei)[names(psia_costs_gpei) == "outbreaks"] <- "outbreaks_sia"

names(psia2_costs_hc)[names(psia2_costs_hc) == "outbreaks"] <- "outbreaks_sia"
names(psia2_costs_gpei)[names(psia2_costs_gpei) == "outbreaks"] <- "outbreaks_sia"


### pSIA ##
osia_costs_hc$osia_dalys_discounted0 <- osia_costs_hc$dalys_discounted0
osia_costs_hc$osia_dalys_discounted3 <- osia_costs_hc$dalys_discounted3
osia_costs_hc$osia_costs_discounted <- osia_costs_hc$costs_discounted

psia_icer_hc <- merge(psia_costs_hc, 
                            y = osia_costs_hc[ , c("year", "node", "outbreaks", "ri_coverage","osia_dalys_discounted0",
                                                   "osia_dalys_discounted3","osia_costs_discounted" )],
                            by=c("year", "node", "ri_coverage"))

psia_icer_hc <- psia_icer_hc %>% 
  group_by(ri_coverage, node) %>%
  dplyr::summarise(scenario="pSIA",
                   sia_costs_discounted = sum(costs_discounted),
                   osia_costs_discounted = sum(osia_costs_discounted),
                   dalys0 = sum(dalys_discounted0),
                   dalys3 = sum(dalys_discounted3),
                   osia_dalys0 = sum(osia_dalys_discounted0),
                   osia_dalys3 = sum(osia_dalys_discounted3),
                   outbreaks_sia = sum(outbreaks_sia),
                   outbreaks_osia = sum(outbreaks))

psia_icer_hc <- psia_icer_hc %>% 
  dplyr::mutate(scenario="pSIA",
                dalys_averted_discounted0 = osia_dalys0 - dalys0,
                dalys_averted_discounted3 = osia_dalys3 - dalys3,
                outbreaks_averted = outbreaks_osia - outbreaks_sia,
                cost_diff = sia_costs_discounted - osia_costs_discounted)

# GPEI
osia_costs_gpei$osia_dalys_discounted0 <- osia_costs_gpei$dalys_discounted0
osia_costs_gpei$osia_dalys_discounted3 <- osia_costs_gpei$dalys_discounted3
osia_costs_gpei$osia_costs_discounted <- osia_costs_gpei$costs_discounted

psia_icer_gpei <- merge(psia_costs_gpei, 
                              y = osia_costs_gpei[ , c("year", "node", "outbreaks", "ri_coverage","osia_dalys_discounted0",
                                                       "osia_dalys_discounted3","osia_costs_discounted" )],
                              by=c("year", "node", "ri_coverage"))

psia_icer_gpei <- psia_icer_gpei %>% 
  group_by(ri_coverage, node) %>%
  dplyr::summarise(scenario="pSIA",
                   sia_costs_discounted = sum(costs_discounted),
                   osia_costs_discounted = sum(osia_costs_discounted),
                   dalys0 = sum(dalys_discounted0),
                   dalys3 = sum(dalys_discounted3),
                   osia_dalys0 = sum(osia_dalys_discounted0),
                   osia_dalys3 = sum(osia_dalys_discounted3),
                   outbreaks_sia = sum(outbreaks_sia),
                   outbreaks_osia = sum(outbreaks))

psia_icer_gpei <- psia_icer_gpei %>% 
  dplyr::mutate(scenario="pSIA",
                dalys_averted_discounted0 = osia_dalys0 - dalys0,
                dalys_averted_discounted3 = osia_dalys3 - dalys3,
                outbreaks_averted = outbreaks_osia - outbreaks_sia,
                cost_diff = sia_costs_discounted - osia_costs_discounted)



### pSIA2 - campaign evey 2 years ##
osia_costs_hc$outbreaks_osia <- osia_costs_hc$outbreaks
psia2_icer_hc <- merge(psia2_costs_hc, 
                             y = osia_costs_hc[ , c("year", "node", "ri_coverage","outbreaks_osia", "osia_dalys_discounted0",
                                                    "osia_dalys_discounted3","osia_costs_discounted" )],
                             by=c("year", "node", "ri_coverage"))

psia2_icer_hc <- psia2_icer_hc %>% 
  group_by(ri_coverage, node) %>%
  dplyr::summarise(scenario="pSIA2",
                   sia_costs_discounted = sum(costs_discounted),
                   osia_costs_discounted = sum(osia_costs_discounted),
                   dalys0 = sum(dalys_discounted0),
                   dalys3 = sum(dalys_discounted3),
                   osia_dalys0 = sum(osia_dalys_discounted0),
                   osia_dalys3 = sum(osia_dalys_discounted3),
                   outbreaks_sia = sum(outbreaks_sia),
                   outbreaks_osia = sum(outbreaks_osia))


psia2_icer_hc <- psia2_icer_hc %>% 
  dplyr::mutate(scenario="pSIA2",
                dalys_averted_discounted0 = osia_dalys0 - dalys0,
                dalys_averted_discounted3 = osia_dalys3 - dalys3,
                outbreaks_averted = outbreaks_osia - outbreaks_sia,
                cost_diff = sia_costs_discounted - osia_costs_discounted)

# GPEI
osia_costs_gpei$outbreaks_osia <- osia_costs_gpei$outbreaks
psia2_icer_gpei <- merge(psia2_costs_gpei, 
                               y = osia_costs_gpei[ , c("year", "node", "ri_coverage","outbreaks_osia", "osia_dalys_discounted0",
                                                        "osia_dalys_discounted3","osia_costs_discounted" )],
                               by=c("year", "node", "ri_coverage"))

psia2_icer_gpei <- psia2_icer_gpei %>% 
  group_by(ri_coverage, node) %>%
  dplyr::summarise(scenario="pSIA2",
                   sia_costs_discounted = sum(costs_discounted),
                   osia_costs_discounted = sum(osia_costs_discounted),
                   dalys0 = sum(dalys_discounted0),
                   dalys3 = sum(dalys_discounted3),
                   osia_dalys0 = sum(osia_dalys_discounted0),
                   osia_dalys3 = sum(osia_dalys_discounted3),
                   outbreaks_sia = sum(outbreaks_sia),
                   outbreaks_osia = sum(outbreaks_osia))

psia2_icer_gpei <- psia2_icer_gpei %>% 
  dplyr::mutate(scenario="pSIA2",
                dalys_averted_discounted0 = osia_dalys0 - dalys0,
                dalys_averted_discounted3 = osia_dalys3 - dalys3,
                outbreaks_averted = outbreaks_osia - outbreaks_sia,
                cost_diff = sia_costs_discounted - osia_costs_discounted)


## binding all data for INCREMENTAL results ##
icer_hc <- rbind(psia_icer_hc, psia2_icer_hc)
icer_gpei <- rbind(psia_icer_gpei, psia2_icer_gpei)

icer_hc$ri <- paste(icer_hc$ri_coverage, " RI coverage", sep="%")
icer_gpei$ri <- paste(icer_gpei$ri_coverage, " RI coverage", sep="%")
#end
