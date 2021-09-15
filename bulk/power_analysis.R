### Performed mixed-model power analysis for ST data #### 
library(simr)
simrOptions(progress=FALSE)

# Read the fuctions file 
setwd("/Users/sanjavickovic/Desktop/morphoSPOT/3dst_repo/3dst/data")
source('../functions/Read_Functions.R')

# Which samples do you want to use? 
norm_samples = c("RA1", "RA2", "RA3", "RA4", "RA5", "RA6")

# Where are you norm expression R files located? 
path_samples = "../data/"

# Load data as R objects for selected RA patient
section_labels = ""
all_ann_all = ""
id_samples = ""
sero_status = ""
treatment_status_dmards = ""
treatment_status_tnf = ""
treatment_status_il6 = ""
location_status = ""

for (norm_sample in norm_samples){
  files_norm = list.files(pattern = glob2rx(paste0(norm_sample, "_norm.exp.values.*")), path = path_samples)
  all_ann_rest_tmp = ""
  all_ann_inf_tmp = ""
  counter = 0
  for (i in 1:length(files_norm)){
    load(paste0(path_samples, files_norm[i]))
    load(paste0(path_samples, norm_sample, "_", i, "_selected_adjusted_spots_3D_manual_app"))
    assign(paste0("inf_s", i), read.csv(paste0("../data/", norm_sample, "_", i, "_all_inf.csv"), header = T, row.names = 1)) # files found on SCP
    section_labels = c(section_labels, sapply(strsplit(colnames(get(paste0("m", i))), "_"), "[[", 1))
    all_ann_inf_tmp = c(all_ann_inf_tmp, rep("Inf", length(row.names(get(paste0("inf_s", i))))))
    all_ann_rest_tmp = c(all_ann_rest_tmp, rep("Rest", length(sapply(strsplit(colnames(get(paste0("m", i))), "_"), "[[", 1))-length(row.names(get(paste0("inf_s", i))))))
    counter = counter + length(sapply(strsplit(colnames(get(paste0("m", i))), "_"), "[[", 1))
  }
  all_ann_inf_tmp = all_ann_inf_tmp[-1]
  all_ann_rest_tmp = all_ann_rest_tmp[-1]
  all_ann_all = c(all_ann_all, c(all_ann_inf_tmp, all_ann_rest_tmp))
  id_samples = c(id_samples, rep(norm_sample, counter))
  
  # collect data on sero status
  if ((norm_sample == "RA1") | (norm_sample == "RA2") | (norm_sample == "RA3")){
    sero_status = c(sero_status, rep("positive", counter))
  }
  else sero_status = c(sero_status, rep("negative", counter))
  
  # collect data on biopsy location and treatment status
  if ((norm_sample == "RA1") | (norm_sample == "RA3") | (norm_sample == "RA4")){
    location_status = c(location_status, rep("knee", counter))
    treatment_status_dmards = c(treatment_status_dmards, rep("yes", counter))
  }
  else {
    location_status = c(location_status, rep("hip", counter))
    treatment_status_dmards = c(treatment_status_dmards, rep("yes", counter))
  }
    
  if ((norm_sample == "RA1") | (norm_sample == "RA2") | (norm_sample == "RA3")| (norm_sample == "RA4")){
    treatment_status_tnf = c(treatment_status_tnf, rep("no", counter))}
  else treatment_status_tnf = c(treatment_status_tnf, rep("yes", counter))
  
  if ((norm_sample == "RA4") | (norm_sample == "RA5")){
    treatment_status_il6 = c(treatment_status_il6, rep("yes", counter))}
  else treatment_status_il6  = c(treatment_status_il6, rep("no", counter))
}
id_samples = id_samples[-1]
section_labels = section_labels[-1]
all_ann_all = all_ann_all[-1]
location_status = location_status[-1]
sero_status = sero_status[-1]
treatment_status_dmards = treatment_status_dmards[-1]
treatment_status_tnf = treatment_status_tnf[-1]
treatment_status_il6 = treatment_status_il6[-1]

# set up covariates data frame
covars = data.frame(cbind(as.numeric(as.character(str_replace(id_samples, "RA", ""))), as.numeric(as.character(str_replace(section_labels, "X", ""))), all_ann_all,location_status,sero_status,treatment_status_dmards,treatment_status_tnf,treatment_status_il6))
colnames(covars) = c("id_samples", "section_labels", "all_ann_all","location_status","sero_status",
                     "treatment_status_dmards,treatment_status_tnf","treatment_status_il6")

# perform power analysis 
## Intercept and slopes for intervention, time1, time2, intervention:time1, intervention:time2
fixed <- c(5, 0, 0.1, 0.2, 1, 0.9)

## Random intercepts for participants clustered by class
rand <- list(0.5)

## residual variance
res <- 1

# numbers of covariates
sero_status_numbers = length(unique(sero_status))
section_labels_numbers = length(unique(section_labels))
location_status_numbers = length(unique(location_status))
id_samples_numbers = length(unique(id_samples))

## Create model
# Simulating without data based on pseudo-bukl ie. do not take spot replication into account
samples_id <- c(rep(1,4), rep(2,7), rep(3,5), rep(4,4), rep(5,3), rep(6,4))
section_id <- c(c(1:4), c(1:7), c(1:5), c(1:4), c(1:3), c(1:4))
location <- c(rep("knee",4), rep("hip",7), rep( "knee",5), rep("hip",4), rep("knee",3), rep("hip",4))
group<- c(rep("positive",4), rep("positive",7), rep( "positive",5), rep("negative",4), rep("negative",3), rep("negative",4))
covars <- data.frame(samples_id, section_id, location, group)

#model
model = makeLmer(y ~ location*section_id + (1|group/samples_id), fixef=c(rep(0.5,4)), VarCorr=list(0.5, 0.5), sigma=2, data=covars)

## power sim for time
sim_location <- powerSim(model, nsim=100, test = fcompare(y~location))
sim_location

## power sim for sero status
sim_sero <- powerSim(model, nsim=100, test = fcompare(y~group))
sim_sero

# ## power sim for section labels
sim_section <- powerSim(model, nsim=100, test = fcompare(y~samples_id))
sim_section

# power sim for section labels
sim_inf <- powerSim(model, nsim=100, test = fcompare(y~section_id))

# extend model # this is pic c
model_ext_class <- extend(model, within = "group+location", n=500) # possible also within samples instrad of along
p_curve_treat <- powerCurve(model_ext_class, test=fcompare(y~section_id), along = "section_id", breaks = seq(1,7, by = 1))
plot(p_curve_treat)

p_curve_treat <- powerCurve(model_ext_class, test=fcompare(y~samples_id), along = "samples_id", breaks = seq(1,5, by = 1))
plot(p_curve_treat)

model_ext_class <- extend(model, along = "samples_id", n=250) # possible also within samples instrad of along
p_curve_treat <- powerCurve(model_ext_class, test=fcompare(y~samples_id), along = "samples_id", breaks = seq(1,260, by = 10))
plot(p_curve_treat)


## Create model
# Simulating without data based spot replication into account
samples_id <- as.numeric(sapply(strsplit(id_samples, "RA"), "[[", 2))
section_id <- as.numeric(sapply(strsplit(section_labels, "X"), "[[", 2))
location <- location_status
group<- sero_status
covars <- data.frame(samples_id, section_id, location, group)
covars

#model
model = makeLmer(y ~ location*section_id + (1|group/samples_id), fixef=c(rep(0.5, 4)), VarCorr=list(0.5, 0.5), sigma=2, data=covars)

## power sim for time
sim_location <- powerSim(model, nsim=100, test = fcompare(y~location))
sim_location

## power sim for sero status
sim_sero <- powerSim(model, nsim=100, test = fcompare(y~group))
sim_sero

# ## power sim for section labels
sim_section <- powerSim(model, nsim=100, test = fcompare(y~samples_id))
sim_section

# power sim for section labels
sim_inf <- powerSim(model, nsim=100, test = fcompare(y~section_id))
sim_inf

p_curve_treat <- powerCurve(model, test=fcompare(y~section_id), along = "section_id", breaks = seq(1,7, by = 1))
plot(p_curve_treat)

p_curve_treat <- powerCurve(model, test=fcompare(y~samples_id), along = "samples_id", breaks = seq(1,5, by = 1))
plot(p_curve_treat)
