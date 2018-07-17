## Script for running thrust model
# DH - TODO 7/5/18
#   * Refine output plotting
#   * Add convergence_rate latent variable to Stan model [done]
#   * Modify Stan model to use lognormal dist. for uplift latent variable

## ------------- Load Libraries------------------------
library("rstan")
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
centfile = "./inputs/Line9_PreprocessedCentroids.csv"      # Data for drainage basin centroids
linefile = "./inputs/Line9_ReferenceLine.csv"        # Reference line

outtag = "Line9"
outdir = "./results/"
xcol = "xobs"
zcol = "zobs"
uloccol = "u_location"
uscalecol = "u_scale"


wkdir = "~/ETH/Thesis/Main Project Tasks/Fault Modeling/thrust/"
model_src = "./uplift2dip.stan" # Stan source file

x_mft = 0                 # x-coordinate of where fault reaches surface [m]
convergence_rate = 21.2  # convergence rate, (along the deepest portion of the decollement) [m/yr]
convergence_rate_unc = 2 # convergence rate uncertainty [m/yr]

M = NA                     # (Optional: NA or int), forces clustering to a certain number of fault segments

ONLY_CLUSTER = FALSE


nChains = 4                 # Number of Markov Chains
nCores = 1                  # Number of cores to use for processing

n_plot_draws = 100
plot_margin = 5000
vertical_exag = 2

options(mc.cores = nCores)
setwd(wkdir)
source("./bin/thrust_utils.R")
## -------------- Pre-processing --------------------
# Read data files
cent = read.csv(centfile)
refline = read.csv(linefile)

cent["xobs"] = cent[xcol]
cent["zobs"] = cent[zcol]
cent["u_location"] = cent[uloccol]
cent["u_scale"] = cent[uscalecol]

# Trim points past range front
if (any(cent$xobs < x_mft)) {
  cent = cent[-which(cent$xobs < x_mft),]
}
cent = cent[order(cent$xobs),]


rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}
elevation.profile = approxfun(refline$x,refline$z)

# Normalize data for clustering
cluster_data = scale(array(cbind(cent$xobs,exp(cent$u_location)),dim=c(length(cent$u_location),2)),center=TRUE,scale=TRUE)

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

if (ONLY_CLUSTER) {
  if (exists("p_bic")) {
    print(p_bic)
  }
  print(p_class)
  stop("Stopping After Clustering Step. Adjust CLUSTER_ONLY parameter to run full script.")
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
# Prep data structure to pass to stan
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
    modelization_error_unc = 1.5        # Weakly informative prior on modelization uncertainty (mm/yr)
)

# Compile Model Code
sm = stan_model(model_src)

## Specify Any Non-Default Parameter Initialization Points Here
# init_pt = list()
# for (i in 1:nChains){
#   this_init = list(
#     u_latent = array(rep(1,standata$nobs)))
#   init_pt[[length(init_pt)+1]] = this_init
# }

fit = sampling(sm,data=standata,iter=4000,chains=nChains,control=list(adapt_delta=0.95,max_treedepth=14))
writeLines("Sampling Complete, HMC Diagnostics:")
check_hmc_diagnostics(fit)

samples = rstan::extract(fit)
dip_samples = samples$theta
dip_samples = dip_samples[,mu_LR] # Reording theta to match clusters
mean_dip = apply(dip_samples,c(2),mean)
u_samples = samples$predicted_uplift_obs
mean_predicted_u = apply(u_samples,c(2),mean)
sigma_samples = exp(samples$sigma_modelization)
mean_sigma = mean(sigma_samples)


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

xmin = min(cent$xobs) - plot_margin
xmax = max(cent$xobs) + plot_margin

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
              geom_line(data=meanlinedf,aes(x=x,y=depth,color="mean"),size=1) + 
              geom_point(data=meanlinedf,aes(x=x,y=depth,color="meanbp"),size=1) + 
              scale_y_reverse() + 
              coord_fixed(ratio=vertical_exag,xlim=c(xmin,xmax)) + 
              ggtitle("Fault Geometries") + 
              labs(colour="Fault Breakpoint") + 
              xlab("Along Profile Distance [m]") + 
              ylab("Depth [m]") +
              scale_colour_manual(name="",
                                  values=c(mean="red",
                                           draw="blue",
                                           meanbp='gray10'),
                                  labels=c(mean="Mean Posterior Prediction",
                                           draw="Posterior Model Draws",
                                           meanbp="Mean Fault Breakpoints")) + 
              scale_alpha_manual(name="",
                                  values=c(draw=0.05),
                                  labels=c(),guide=FALSE)

uprofiledf = data.frame(
  x = standata$x_obs,
  u_observed = exp(standata$u_location),
  unc_up = exp(standata$u_location + standata$u_scale),
  unc_down = exp(standata$u_location - standata$u_scale),
  mean_prediction = mean_predicted_u
)
p_uprofile = ggplot(data=uprofiledf,aes(x=x)) + 
              geom_point(aes(y=u_observed,color="obs"),size=0.5) + 
              geom_errorbar(aes(ymin=unc_down,ymax=unc_up,color="obs"),size=0.5) + 
              geom_point(aes(y=mean_prediction,color="pred"),size=0.25) + 
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
              ggtitle("Uplift Field") + 
              coord_cartesian(xlim=c(xmin,xmax))
              
p_combined = ggarrange(p_uprofile, p_xsection, heights = c(4, 4), nrow = 2, align = "v")




## 7. Display and/or print plots and other outputs
if (exists("p_bic")) {
  print(p_bic)
}
print(p_class)
print(p_combined)

outtext = sprintf(
  "Run Complete! Fit Model for %i segments.
  Results:
  Theta (degrees):
  %s
  Modelization Error:
  %f"
,n_segments,paste(rad2deg(mean_dip),collapse=", "),mean_sigma)
writeLines(outtext)

fitfile = paste(outdir,outtag,"_fit.RDS",sep="")
summaryfile = paste(outdir,outtag,"_summary.txt",sep="")
bicfile = paste(outdir,outtag,"_BIC.pdf",sep="")
clusteringfile = paste(outdir,outtag,"_Clustering.pdf",sep="")
profilefile = paste(outdir,outtag,"_GeometryProfile.pdf",sep="")
geomfile = paste(outdir,outtag,"_MeanGeometry.csv",sep="")

fitsum = summary(fit,pars=c("theta","sigma_modelization"))
fitsum = fitsum$summary
# transform sigma from log to real space
fitsum[n_segments+1,1:(dim(fitsum)[2]-2)] = exp(fitsum[n_segments+1,1:(dim(fitsum)[2]-2)])
sumtext = paste(capture.output(print(fitsum)),collapse="\n")

meanxy = distance2xy(refline,mean_geom)
meangeomdf = data.frame(
  x_m = meanxy[,1],
  y_m = meanxy[,2],
  d_m = mean_geom[,1],
  z_m = mean_geom[,2]
)
geomtext = paste(capture.output(print(meangeomdf,row.names=FALSE)),collapse="\n")

summarytext = sprintf("
\"thrust\" Model Summary:
%s

Name: %s
N_Segments: %i

Mean Geometry:
%s

Segment Dips:
%s
",format(Sys.time(), "%a %b %d %X %Y"),outtag,n_segments,geomtext,sumtext)

# Write outputs
saveRDS(fit,file=fitfile)
write(summarytext, file = summaryfile)
write.csv(meangeomdf,file=geomfile,row.names=FALSE)
if (exists("p_bic")) {
  # Remove Auto-Generated Best Fit Line
  p_bic$layers[[3]] = NULL
  ggsave(bicfile,plot=p_bic,height=4,width=5,units="in",dpi=200,device="pdf")
}
ggsave(clusteringfile,plot=p_class,height=4,width=5,units="in",dpi=200,device="pdf")
ggsave(profilefile,plot=p_combined,height=6,width=8,units="in",dpi=200,device="pdf")


