posteriordraw2geometry <- function(theta,mu,x,elevation.profile,x1u,z1u) {
  # mu should be n_segments*2 matrix
  # theta is n_segements length vector
  M = length(theta)[1]
  
  # 1. Get x coordinates of surface breakpoints
  surface_bpx = array(0,dim=c(M-1))
  for (j in 1:M-1) {
    surface_bpx[j] = (mu[j,1] + mu[j+1,1]) / 2
  }
  fault_bp = array(0,dim=c(M+1,2))
  
  # First point = MFT @ surface
  fault_bp[1,1] = x1u
  fault_bp[1,2] = z1u
  # Next points are intersection between fault dip line and kink band line
  for (j in 1:(M-1)){
    # Get the equation for the line that bisects fault segments (tau)
    bisector = ((pi + theta[j+1]) - theta[j]) / 2
    bisector_slope_angle = -1*(bisector - theta[j+1])
    bisector_slope = tan(bisector_slope_angle)
    bisector_yint = elevation.profile(surface_bpx[j]) - surface_bpx[j]*bisector_slope
    
    # Get the equation for the fault line
    fault_slope = tan(theta[j])
    fault_yint = fault_bp[j,2] - fault_bp[j,1]*fault_slope
    # The fault breakpoint is at the intersection of these lines:
    fault_bp[j+1,1] = (fault_yint-bisector_yint) / (bisector_slope-fault_slope)
    fault_bp[j+1,2] = (bisector_slope*fault_yint - fault_slope*bisector_yint) / (bisector_slope-fault_slope)
  }
  # Final point is at edge of profile (xmax)
  fault_bp[M+1,1] = max(x) 
  fault_bp[M+1,2] = fault_bp[M,2] + tan(theta[M])*(fault_bp[M+1,1]-fault_bp[M,1])
  
  return(fault_bp)
}

posteriordraw2geometry_bp <- function(theta,surface_bpx,x,elevation.profile,x1u,z1u) {
  # mu should be n_segments*2 matrix
  # theta is n_segements length vector
  M = length(theta)[1]
  fault_bp = array(0,dim=c(M+1,2))
  
  # First point = MFT @ surface
  fault_bp[1,1] = x1u
  fault_bp[1,2] = z1u
  # Next points are intersection between fault dip line and kink band line
  for (j in 1:(M-1)){
    # Get the equation for the line that bisects fault segments (tau)
    bisector = ((pi + theta[j+1]) - theta[j]) / 2
    bisector_slope_angle = -1*(bisector - theta[j+1])
    bisector_slope = tan(bisector_slope_angle)
    bisector_yint = elevation.profile(surface_bpx[j]) - surface_bpx[j]*bisector_slope
    
    # Get the equation for the fault line
    fault_slope = tan(theta[j])
    fault_yint = fault_bp[j,2] - fault_bp[j,1]*fault_slope
    # The fault breakpoint is at the intersection of these lines:
    fault_bp[j+1,1] = (fault_yint-bisector_yint) / (bisector_slope-fault_slope)
    fault_bp[j+1,2] = (bisector_slope*fault_yint - fault_slope*bisector_yint) / (bisector_slope-fault_slope)
  }
  # Final point is at edge of profile (xmax)
  fault_bp[M+1,1] = max(x) 
  fault_bp[M+1,2] = fault_bp[M,2] + tan(theta[M])*(fault_bp[M+1,1]-fault_bp[M,1])
  
  return(fault_bp)
}
