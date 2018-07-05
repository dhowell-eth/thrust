Package: thrust
Written By: Dorran Howell - dorran.howell@gmail.com
Date: 2018

"thrust" is a set of scripts/models for converting estimated uplift rates to a fault geometry assuming an emergent thrust. This package was written for use in my Master's Thesis, which examines the geometry of the Main Himalayan Thrust using uplift rate estimates. 

The main script, "run_thrust.R", is designed to be run in R Studio and takes two input csv files:
	1.) .csv file of a reference profile line on which your data is projected, containing fields:
		- x			 : x coordinates along ref line (fault is assumed to dip towards greater coordinates)
		- z			 : elevation of topography along reference line
	2.) .csv file of uplift data with the following fields:
		- xobs		 : x coordinate of data projected to reference line 
		- zobs		 : elevation of data projected to reference line
		- u_location : estimated location parameter for uplift distribution
		- u_scale	 : estimated scale parameter for uplift distribution

run_thrust.R finds the optimal number of populations in your uplift profile using the Bayesian Information Criterion (BIC) and uses this to cluster uplift data into M populations. These populations are then fed into an MCMC sampler for a Bayesian inverse model for fault dip written in STAN.  
