


## ------------- Load Libraries------------------------
library("rstan", lib.loc="~/R/win-library/3.4")
rstan_options(auto_write = TRUE)
library(ggplot2)
library(shinystan)
library(Matrix)
library(readr)
library(mclust)

## --------------- Run Parameters----------------------
centfile = "./inputs/SyntheticUpliftCentroids.csv" # Data for drainage basin centroids
linefile = "./inputs/SyntheticUpliftRefLine.csv"   # Reference line

model_src = "./uplift2dip.stan" # Stan source file

x_mft = 1e3                 # x-coordinate of where fault reaches surface [m]
convergence_rate = 21.5e-3  # convergence rate, (along the deepest portion of the decollement) [m/yr]
convergence_rate_unc = 1e-3 # convergence rate uncertainty [m/yr]

M = NA                      # (Optional: NA or int), forces clustering to a certain number of fault segments

nChains = 2                 # Number of Markov Chains
nCores = 1                  # Number of cores to use for processing

options(mc.cores = nCores)

## -------------- Pre-processing --------------------
# Read data files
cent = read.csv(centfile)
refline = read.csv(linefile)

# Convert elevations to z (positive axis is below MSL)
refline$z = refline$z*-1
cent$zobs = cent$zobs*-1
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}
elevation.profile = approxfun(refline$x,refline$z)

# Normalize data for clustering
cluster_data = scale(array(cbind(cent$xobs,cent$u),dim=c(length(cent$u),2)),center=TRUE,scale=TRUE)

## ------------- Begin Modeling ----------------------

## 1. Run GMM Clustering
if (is.na(M)) {
  cluster_fit = Mclust(cluster_data)
  plot(cluster_fit,what="BIC")
  plot(cluster_fit,what="classification")
  
} else {
  cluster_fit = Mclust(cluster_data, G = M)
  plot(cluster_fit,what="classification")
  
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
standata = list(
    x = refline$x,
    z = refline$z,
    x_obs = cent$xobs,
    z_obs = cent$zobs,
    u_location = cent$u,
    u_scale = rep(0.05,nrow(cent)),
    groupid = categories,
    nobs=nrow(cent),
    n_segments=n_segments,
    nx=nrow(refline),
    gamma_guess = array(1.5),
    convergence=convergence_rate,
    x1u =  x_mft,
    z1u = elevation.profile(x_mft),
    N=nrow(cent),                     
    modelization_error_unc = 5        # Weakly informative prior on modelization uncertainty (mm/yr)
)

sm = stan_model(model_src)

fit = sampling(sm,data=standata)
samples = rstan::extract(fit)
dip_samples = samples$theta
dip_samples = dip_samples[,mu_LR] # Reording theta to match clusters
mean_dip = apply(dip_samples,c(2),mean)

## 5. Pass theta samples and cluster means to geometry estimator
geoms = array(0,dim=c(dim(dip_samples)[1],n_segments+1,2))
for (i in 1:dim(dip_samples)[1]) {
  geoms[i,,] = posteriordraw2geometry(dip_samples[i,],mu_t,data$x,elevation.profile,data$x1u,data$z1u)
}
mean_geom = posteriordraw2geometry(mean_dip,mu_t,data$x,elevation.profile,data$x1u,data$z1u)

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

p_xsection = ggplot(data=bpdata,aes(x=x,y=depth,color=c),alpha=0.1) + 
              geom_point() + 
              geom_line(data=meanlinedf,aes(x=x,y=depth,color="Mean Geometry")) + 
              scale_y_reverse() + 
              coord_fixed() + 
              ggtitle("Posterior Fault Geometries") + 
              labs(colour="Fault Breakpoint") + 
              xlab("Along Profile Distance [m]") + 
              ylab("Depth [m]") 
print(p_xsection)
