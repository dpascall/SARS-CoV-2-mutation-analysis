##R version 3.6.0
##rstan version 2.21.2
library(brms) ##2.13.5
library(ape) ##5.3
library(lubridate) ##1.7.4
options(mc.cores = parallel::detectCores())
rm(list=ls())

##read in data
data <- read.csv("COVID19-DavidPascall10082020.2.csv", stringsAsFactors = F)

##subset out rows required for the analysis
data <- data[,c("ordinal_scale", "age1", "sex", "d_614_g", "n_439_k", "cog_tree_identifier", "cluster", "date")]

##give any row not associated with a cluster a unique cluster ID
k <- 1 
for  (i in 1:nrow(data)) {
  if (data$cluster[i] == "") {
    data$cluster[i] <- k
    k <- k + 1
  }
}

##format dates into days since first record
data$date <- as.POSIXct(data$date, format = "%d/%m/%Y")
data$dayssincefirstrecord <- as.numeric(round(difftime(data$date, min(data$date), units = "days")))


##exclude rows with missing data (n = 1087)
data <- data[!data$d_614_g == "na",]
data <- data[!data$n_439_k == "na",]
data <- data[complete.cases(data),]

##generate mutation variable
data$mutation <- paste(as.character(data$d_614_g),as.character(data$n_439_k),sep=".")

##read tree
tree <- read.tree("scottish_cog_alignment_wdates.fasta.timetree.nwk")
tree$tip.label <- do.call("c", lapply(tree$tip.label, function(x){return(strsplit(x, "[|]")[[1]][[1]])}))
tree <- drop.tip(tree,tree$tip.label[!tree$tip.label%in%data$cog_tree_identifier], root.edge = 1)

##exclude all rows without corresponding tips in the tree
data <- data[data$cog_tree_identifier%in%tree$tip.label,]

##set data types and centre numeric variables
data$age1 <- scale(data$age1, scale = F)
data$sex <- as.factor(data$sex)
data$mutation <- as.factor(data$mutation)
data$cluster <- as.factor(data$cluster)
data$dayssincefirstrecord <- scale(data$dayssincefirstrecord, scale = F)

##generate phylogenetic correlation matrix under brownian motion assumption
A <- vcv.phylo(tree, corr = T)

##set priors
bprior <- c(prior(student_t(3, 0, 2.5), class = Intercept), prior(normal(0, 2.5), class = b), prior(exponential(0.4), class = sd), prior(exponential(0.4), class = sds))

##rescale so that 1 is best and 4 is worst
##set death to 0
data$ordinal_scale[data$ordinal_scale == 1] <- 0
##set none of the above to 1
data$ordinal_scale[data$ordinal_scale == 4] <- 1
##set death to 4
data$ordinal_scale[data$ordinal_scale == 0] <- 4
##set vent to 0
data$ordinal_scale[data$ordinal_scale == 2] <- 0
##set O2 to 2
data$ordinal_scale[data$ordinal_scale == 3] <- 2
##set vent to 3
data$ordinal_scale[data$ordinal_scale == 0] <- 3

##model
model <- brm(ordinal_scale ~ mutation + sex + s(age1, k = 30) + s(dayssincefirstrecord, k = 30) 
    + (1|cluster) + (1|gr(cog_tree_identifier, cov = A)), family = cumulative(link = logit), 
    data = data, data2 = list(A = A), iter = 4000, warmup = 1500, control = list(adapt_delta = 0.99999), 
    chains = 4, prior = bprior, refresh = 5)
save(model, file = "Severitymodel.Rdata")
load("~/Desktop/Severitymodel.Rdata")


##posterior predictive check
pp_check(model, type = "ecdf_overlay")
##looks good

library(SPIn)
library(rstan)

##check diagnostics for all parameters
diagnostics <- monitor(model$fit)

##extract draws for parameters to report
draws <- as.matrix(model)[,c("b_mutation1.0", "b_mutation1.1", "b_sex1", "b_Intercept[1]", "b_Intercept[2]", "b_Intercept[3]", "sd_cluster__Intercept", "sd_cog_tree_identifier__Intercept", "sds_sage1_1", "sds_sdayssincefirstrecord_1")]

##get credible intervals accounting for limit of 0 for scales
results <- matrix(NA, nrow = 10, ncol = 3)
results[,1] <- colMeans(draws)

for (i in 1:10) {
  if (grepl("sd", colnames(draws)[i])) {
    results[i, c(2:3)] <- SPIn(draws[,i], conf = 0.95, lb = 0)$spin
  } else {
    results[i, c(2:3)] <- SPIn(draws[,i], conf = 0.95)$spin
  }
}

results <- round(results, digits = 2)

colnames(results) <- c("Posterior mean", "Lower bound of 95% shortest probability interval", "Upper bound of 95% shortest probability interval")
row.names(results) <- c("Presence of 614G", "Presence of 614G and 439K", "Female biological sex", "Threshold 1", "Threshold 2", "Threshold 3", "Cluster random effect standard deviation", "Phylogenetic random effect standard deviation", "Age penalised spline standard deviation", "Days since first case in the dataset penalised spline standard deviation")

write.csv(results, file = "~/Desktop/SupResultsTableSeverity.csv")
