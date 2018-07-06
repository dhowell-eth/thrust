functions {
  
    vector suppe_algebra_system(vector y, vector theta, real[] x_r, int[] x_i) {
      // Algebra system for solving Suppe Equation (Following Allmendinger Script)
      // x_r[1] = data values, theta, used in evaluating system (i.e. not latent)
      // y[1] = gamma
      // theta[1] = theta
      // theta[2] = phi
      vector[1] result;
      result[1] = atan((-1*sin(y[1]-theta[1])*(sin(2*y[1]-theta[1])-sin(theta[1]))) / 
                       (cos(y[1]-theta[1])*(sin(2*y[1]-theta[1])-sin(theta[1])) - sin(y[1]))) 
          - theta[2];
      return result;
    }

    vector predict_uplift_theta_groups(vector theta, real s, int K, vector gamma_guess,
                                       real[] x_r_placeholder,int[] x_i_placeholder) {
    vector[K] R;
    vector[K] phi;
    vector[2] theta_solver;
    vector[1] solver_result;
    vector[K] upred;
    vector[K] slip;
    real gamma;
      
    for (c in 1:K) {
        // For each fault segment...
        // 1. Compute phi, the difference in dip between the current segment and the next one
        if (c==K) {
          phi[c] = theta[c]; // last step is onto an assumed flat (at MFT), so phi=theta
        }
        else {
          phi[c] = theta[c] - theta[c+1];
        }
        // Provide required variables for algebra solver
        theta_solver[1] = theta[c];
        theta_solver[2] = phi[c];
        // Compute gamma from phi and theta
        solver_result = algebra_solver(suppe_algebra_system,gamma_guess,
                          theta_solver,x_r_placeholder,x_i_placeholder);
        gamma = solver_result[1];
        // Use theta, phi, and gamma to compute slip ratios
        R[c] = sin(gamma - theta[c]) / sin(phi[c] + gamma - theta[c]);
    }
    
    // Use convergence rate and slip ratios to compute slip rates along each segment
    slip[1] = s;
    for (c in 2:K) {
        slip[c] = slip[c-1]*R[c-1];
    }
    
    // Use each segments slip rate and dip to compute an uplift rate
    upred = slip .* sin(theta);
    
    return upred;
  }


// --- END FUNCTIONS SECTION --- 

}
data{
  int nobs;                  // Number of Observations
  int nx;                    // Number of grid nodes in x (horizontal) direction
  int n_segments;            // Number of fault segments (i.e. how many parameters to solve for)
  
  vector[nx] x;              // Horizontal Grid Nodes (model) [m]
  vector[nx] z;              // Elevation @ Grid Nodes (model) [m]
  vector[nobs] x_obs;        // x coordinates of observations [m]
  vector[nobs] z_obs;        // z coordinates of observations [m]
  
  vector[nobs] u_location;   // observed uplift rates [mm/yr]
  vector[nobs] u_scale; 
  
  int groupid[nobs];
  
  real convergence;          // horizontal velocity, [m/yr]
  real convergence_std;
  real x1u;                  // x coord. of fault tip at surface [m]
  real z1u;                  // z coord. of fault tip at surface [m]    
  vector[1] gamma_guess;
  
  real modelization_error_unc;

}

transformed data {
  real x_r_placeholder[1];
  int x_i_placeholder[1];
  // Placeholders for required inputs for the algebra solver 
  x_r_placeholder[1] = 0.0;
  x_i_placeholder[1] = 1;
}

parameters {
  vector<lower=-0.01745329,upper=0.5061455>[n_segments] theta; // fault segment angles [radians]
  real<lower=0.0> sigma_modelization; // modelization uncertainty
  vector<lower=0.0>[nobs] u_latent;
  real<lower=0.0> s_latent;
}
transformed parameters {
  vector[n_segments] upred; // Predicted Uplift by Zone
  vector[nobs] predicted_uplift_obs; // Predicted Uplift at Obs Points Using Cluster Classifications
  
  upred = predict_uplift_theta_groups(theta,s_latent,n_segments, gamma_guess,
                                      x_r_placeholder,x_i_placeholder);
  upred  = upred*1000; // m/yr->mm/yr
  
  predicted_uplift_obs = upred[groupid];
}
model {
  
  
  // --------------------- PRE-PROCESSING ---------------------------------
 
  // ----------------------- PRIORS ---------------------------------------
  u_latent ~ lognormal(u_location,u_scale);
  s_latent ~ normal(convergence,convergence_std);
  
  theta ~ normal(0.2443461,0.2443461);
  sigma_modelization ~ normal(0,modelization_error_unc);
  
  // ----------------------LIKELIHOOD FUNCTIONS----------------------------
  // Uplift Rates
  target += normal_lpdf(u_latent | predicted_uplift_obs, sigma_modelization);
}
