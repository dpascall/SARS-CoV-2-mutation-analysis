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
data <- data[,c("ct", "age1", "sex", "d_614_g", "n_439_k", "cog_tree_identifier", "cluster", "date")]

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

##exclude rows with missing data
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
bprior <- c(prior(student_t(3, 20, 10), class = Intercept), prior(normal(0, 10), class = b), prior(exponential(0.1), class = sd), prior(exponential(0.1), class = sds))

##model
model <- brm(ct ~ mutation + sex + s(age1, k = 30) + s(dayssincefirstrecord, k = 30) 
    + (1|cluster) + (1|gr(cog_tree_identifier, cov = A)), family = gaussian, 
    data = data, data2 = list(A = A), iter = 4000, warmup = 1500, control = list(adapt_delta = 0.99999), 
    chains = 4, prior = bprior, refresh = 5)

save(model, file = "CTmodel.Rdata")
load("~/Desktop/CTmodel.Rdata")

##posterior predictive check
pp_check(model, type = "ecdf_overlay")
##looks good

library(SPIn)
library(rstan)

##check diagnostics for all parameters
diagnostics <- monitor(model$fit)

##extract draws for parameters to report
draws <- as.matrix(model)[,c("b_mutation1.0", "b_mutation1.1", "b_sex1", "b_Intercept", "sd_cluster__Intercept", "sd_cog_tree_identifier__Intercept", "sds_sage1_1", "sds_sdayssincefirstrecord_1", "sigma")]

##get credible intervals accounting for limit of 0 for scales
results <- matrix(NA, nrow = 9, ncol = 3)
results[,1] <- colMeans(draws)

for (i in 1:9) {
  if (grepl("sd", colnames(draws)[i])|grepl("sigma", colnames(draws)[i])) {
    results[i, c(2:3)] <- SPIn(draws[,i], conf = 0.95, lb = 0)$spin
  } else {
    results[i, c(2:3)] <- SPIn(draws[,i], conf = 0.95)$spin
  }
}

results <- round(results, digits = 2)

colnames(results) <- c("Posterior mean", "Lower bound of 95% shortest probability interval", "Upper bound of 95% shortest probability interval")
row.names(results) <- c("Presence of 614G", "Presence of 614G and 439K", "Female biological sex", "Intercept", "Cluster random effect standard deviation", "Phylogenetic random effect standard deviation", "Age penalised spline standard deviation", "Days since first case in the dataset penalised spline standard deviation", "Residual standard deviation")

SPIn(draws[,2]-draws[,1], conf = 0.95)$spin

write.csv(results, file = "~/Desktop/SupResultsTableCT.csv")
