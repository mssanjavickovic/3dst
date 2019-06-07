#Load libraries
library(scran)

# Normalize counts RA1 samples
m1 = read.csv("RA1_Raw.exp_1.csv", header = T, row.names = 1)
m2 = read.csv("RA1_Raw.exp_2.csv", header = T, row.names = 1)
m3 = read.csv("RA1_Raw.exp_3.csv", header = T, row.names = 1)
m4 = read.csv("RA1_Raw.exp_4.csv", header = T, row.names = 1)

### Combine all the inf norm data in one matrix
exp.values = mbind(m1, m2, m3, m4)

# Normalize using Scran
all.exp.values = newSCESet(countData=data.frame(exp.values))
clusters = quickCluster(all.exp.values,  min.size=20)
all.exp.values.all = computeSumFactors(all.exp.values, sizes = c(5,10,15,20,25),clusters=clusters, positive = TRUE)
all.exp.values.norm = normalize(all.exp.values.all)
fit = trendVar(all.exp.values.norm, use.spikes=FALSE)
RA1.norm = all.exp.values.norm

# Read in counts RA2 samples
m1 = read.csv("RA2_Raw.exp_1.csv", header = T, row.names = 1)
m2 = read.csv("RA2_Raw.exp_2.csv", header = T, row.names = 1)
m3 = read.csv("RA2_Raw.exp_3.csv", header = T, row.names = 1)
m4 = read.csv("RA2_Raw.exp_4.csv", header = T, row.names = 1)
m5 = read.csv("RA2_Raw.exp_5.csv", header = T, row.names = 1)
m6 = read.csv("RA2_Raw.exp_6.csv", header = T, row.names = 1)
m7 = read.csv("RA2_Raw.exp_7.csv", header = T, row.names = 1)

### Combine all the inf norm data in one matrix
exp.values = mbind(m1, m2, m3, m4, m5, m6, m7)

# Normalize counts for RA2 
all.exp.values = newSCESet(countData=data.frame(exp.values))
clusters = quickCluster(all.exp.values,  min.size=60)
all.exp.values.all = computeSumFactors(all.exp.values, clusters=clusters, sizes = c(2,5,7,10,12,15,17,20,22,25,27,30,35,40), positive = TRUE, assay="counts")
all.exp.values.norm = normalize(all.exp.values.all, log=TRUE)
fit = trendVar(all.exp.values.norm, use.spikes=FALSE)
RA2.norm = all.exp.values.norm

