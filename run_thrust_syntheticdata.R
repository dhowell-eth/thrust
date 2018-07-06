## Script for running thrust model
# DH - TODO 7/5/18
#   * Refine output plotting
#   * Add convergence_rate latent variable to Stan model [done]
#   * Modify Stan model to use lognormal dist. for uplift latent variable

## ------------- Load Libraries------------------------
library("rstan", lib.loc="~/R/win-library/3.4")
rstan_options(auto_write = TRUE)
library(ggplot2)
library(shinystan)
library(Matrix)
library(readr)
library(mclust)
library(factoextra)
library(ggpubr)
library(RColorBrewer)
## --------------- Run Parameters----------------------
centfile = "./inputs/SyntheticUpliftCentroids-LN.csv"      # Data for drainage basin centroids
linefile = "./inputs/SyntheticUpliftRefLine-LN.csv"        # Reference line
faultfile = "./inputs/SyntheticUpliftFaultGeometry-LN.csv" # Control fault geometry

wkdir = "~/ETH/Thesis/Main Project Tasks/Fault Modeling/thrust/"
model_src = "./uplift2dip.stan" # Stan source file

x_mft = 1e3                 # x-coordinate of where fault reaches surface [m]
convergence_rate = 21.5e-3  # convergence rate, (along the deepest portion of the decollement) [m/yr]
convergence_rate_unc = 0.25e-3 # convergence rate uncertainty [m/yr]

M = 3                      # (Optional: NA or int), forces clustering to a certain number of fault segments

nChains = 2                 # Number of Markov Chains
nCores = 1                  # Number of cores to use for processing

n_plot_draws = 100

options(mc.cores = nCores)
setwd(wkdir)
source("./bin/thrust_utils.R")
## -------------- Pre-processing --------------------
# Read data files
cent = read.csv(centfile)
refline = read.csv(linefile)
controlfault = read.csv(faultfile)
# Make sure elevations have positive axis below MSL
#refline$z = refline$z*-1
#cent$zobs = cent$zobs*-1

# Trim points past range front
if (any(cent$xobs < x_mft)) {
  cent = cent[-which(cent$xobs < x_mft),]
}

rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}
elevation.profile = approxfun(refline$x,refline$z)

# Normalize data for clustering
cluster_data = scale(array(cbind(cent$xobs,exp(cent$u_location)),dim=c(length(cent$u),2)),center=TRUE,scale=TRUE)

## ------------- Begin Modeling ----------------------

## 1. Run GMM Clustering
if (is.na(M)) {
  cluster_fit = Mclust(cluster_data)
  #plot(cluster_fit,what="BIC")
  #plot(cluster_fit,what="classification")
  p_bic = factoextra::fviz_mclust(cluster_fit,"BIC",geom="line")
  p_class = factoextra::fviz_mclust(cluster_fit,"classification",geom="point",xlab="Normalized Distance",ylab="Normalized Uplift")
} else {
  cluster_fit = Mclust(cluster_data, G = M)
  p_class = factoextra::fviz_mclust(cluster_fit,"classification",geom="point",xlab="Normalized Distance",ylab="Normalized Uplift")
}

## 2. Extract the best fitting number of populations according to BIC
n_segments = cluster_fit$G

## 3. Extract cluster means and back-transform to unnormalized space
mu = cluster_fit$parameters$mean
scaling = attr(cluster_data,"scaled:scale")
centering = attr(cluster_data,"scaled:center")
if (is.null(centering)) {
  centering = array(0)
}
if (is.null(scaling)) {
  scaling = array(1)
}
mu_t = mu
mu_t[1,] = mu_t[1,] * scaling[1] + centering[1] 
mu_t[2,] = mu_t[2,] * scaling[2] + centering[2] 
mu_t = t(mu_t)
# Sort clusters based on x position (ascending order) 
mu_LR = order(mu_t[,1])
mu_t = mu_t[mu_LR,]

## 4. Use Stan to fit a fault dip to clustered Uplift Data 
categories = cluster_fit$classification
bp = array(0,dim=c(n_segments-1))
for (i in 1:(n_segments-1)) {
   thisx = cent$xobs[which(categories==mu_LR[i])]
   nextx = cent$xobs[which(categories==mu_LR[i+1])]
   bp[i] = (min(nextx)-max(thisx))/2 + max(thisx)
 }
standata = list(
    x = refline$x,
    z = refline$z,
    x_obs = cent$xobs,
    z_obs = cent$zobs,
    u_location = cent$u_location,
    u_scale = cent$u_scale,
    groupid = categories,
    nobs=nrow(cent),
    n_segments=n_segments,
    nx=nrow(refline),
    gamma_guess = array(1.5),
    convergence=convergence_rate,
    convergence_std=convergence_rate_unc,
    x1u =  x_mft,
    z1u = elevation.profile(x_mft),
    N=nrow(cent),                     
    modelization_error_unc = 5        # Weakly informative prior on modelization uncertainty (mm/yr)
)

sm = stan_model(model_src)

fit = sampling(sm,data=standata)
print("Sampling Complete, HMC Diagnostics:")
check_hmc_diagnostics(fit)

samples = rstan::extract(fit)
dip_samples = samples$theta
dip_samples = dip_samples[,mu_LR] # Reording theta to match clusters
mean_dip = apply(dip_samples,c(2),mean)
u_samples = samples$predicted_uplift_obs
mean_predicted_u = apply(u_samples,c(2),mean)

## 5. Pass theta samples and cluster means to geometry estimator
geoms = array(0,dim=c(dim(dip_samples)[1],n_segments+1,2))
for (i in 1:dim(dip_samples)[1]) {
  geoms[i,,] = posteriordraw2geometry_bp(dip_samples[i,],bp,standata$x,elevation.profile,standata$x1u,standata$z1u)
}
mean_geom = posteriordraw2geometry_bp(mean_dip,bp,standata$x,elevation.profile,standata$x1u,standata$z1u)

## 6. Plot geometries
all_bp = array(0,dim=c(dim(geoms)[1]*dim(geoms)[2],3))
for (i in 1:dim(geoms)[1]) {
  nr = dim(geoms)[2]
  j = (i-1)*nr + 1
  all_bp[j:(j+nr-1),1:2] = geoms[i,,]
  all_bp[j:(j+nr-1),3] = 1:nr
}
bpdata = data.frame(
  x = all_bp[,1] ,
  depth = all_bp[,2],
  c = as.factor(all_bp[,3])
)
meanlinedf = data.frame(
  x = mean_geom[,1],
  depth = mean_geom[,2]
)

# Posterior Draws
plot_indices = round(seq(1,dim(geoms)[1],length.out=n_plot_draws))
p_xsection = ggplot() 
for (i in 1:length(plot_indices)) {
  thisgeom = geoms[plot_indices[i],,]
  toplot = data.frame(
    x = thisgeom[,1],
    z = thisgeom[,2]
  )
 p_xsection = p_xsection + geom_line(data=toplot,aes(x=x,y=z,colour="draw",alpha="draw"))
}
# Mean Geometry
p_xsection = p_xsection + 
              geom_point(data=meanlinedf,aes(x=x,y=depth,color="mean"),size=1) + 
              geom_line(data=meanlinedf,aes(x=x,y=depth,color="mean"),size=1) + 
              geom_point(data=controlfault,aes(x=x,y=z,color="control"),shape=8) + 
              geom_line(data=controlfault,aes(x=x,y=z,color="control")) +
              scale_y_reverse() + 
              coord_fixed() + 
              ggtitle("Fault Geometries") + 
              labs(colour="Fault Breakpoint") + 
              xlab("Along Profile Distance [m]") + 
              ylab("Depth [m]") +
              scale_colour_manual(name="",
                                  values=c(mean="red",
                                           draw="blue",
                                           control="gray40"),
                                  labels=c(mean="Mean Posterior Prediction",
                                           draw="Posterior Model Draws",
                                           control="Control Fault Geometry")) + 
              scale_alpha_manual(name="",
                                  values=c(draw=0.05),
                                  labels=c(),guide=FALSE)

print(p_xsection)


uprofiledf = data.frame(
  x = standata$x_obs,
  u_observed = exp(standata$u_location),
  unc_up = exp(standata$u_location + standata$u_scale),
  unc_down = exp(standata$u_location - standata$u_scale),
  mean_prediction = mean_predicted_u
)
p_uprofile = ggplot(data=uprofiledf,aes(x=x)) + 
              geom_point(aes(y=u_observed,color="obs")) + 
              geom_errorbar(aes(ymin=unc_down,ymax=unc_up,color="obs")) + 
              geom_point(aes(y=mean_prediction,color="pred")) + 
              scale_colour_manual(name="",
                                  values=c(
                                    obs="gray1",
                                    pred="red"
                                  ),
                                  labels=c(
                                    obs="Observed",
                                    pred="Mean Posterior Prediction")) + 
              ylab("Uplift [mm/yr]") + 
              xlab("") + 
              ggtitle("Uplift Field")
              
p_combined = ggarrange(p_uprofile, p_xsection, heights = c(4, 4), nrow = 2, align = "v")




## 7. Display and/or print plots and other outputs
print(p_bic)
print(p_class)
print(p_combined)
