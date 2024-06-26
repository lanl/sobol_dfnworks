######################################################################
## Custom functions for loq quantile density discrepancy regression ##
######################################################################
# Load libraries
library(ggplot2)
library(dplyr)
library(KernSmooth)
library(pracma)
library(gplots)
library(tidyr)
library(refund)
library(fda)
library(fdapace)
library(laGP)
library(mvtnorm)

# Convert pdf to cdf
pdf_to_cdf <- function(pdf, x_grid, norm=TRUE){
  cdf = cumtrapz(x_grid, pdf)
  if(norm){ 
    cdf = cdf/cdf[length(cdf)] 
  }
  return(cdf)
}

# Convert cdf to pdf
cdf_to_pdf <- function(cdf, x_grid){
  pdf = gradient(as.numeric(cdf), x_grid)
  return(pdf)
}

# Convert cdf to quantile function
# Uses linear interpolation
cdf_to_quant <- function(cdf, x_grid){
  interp_fn = approxfun(cdf, x_grid)
  quant = interp_fn(x_grid)
  return(quant)
}

# Convert log quantile density to cdf
# Uses linear interpolation
lqd_to_cdf <- function(lqd, x_grid, lb=0){
  p_grid = seq(0, 1, length.out=length(lqd))
  F_inv = lb + as.numeric(cumtrapz(p_grid, exp(lqd))/trapz(p_grid, exp(lqd)))
  cdf = approx(F_inv, p_grid, xout=x_grid, yleft=0, yright=1)
  return(cdf$y)
}

# Convert pdf + quantile to LQD
# Uses linear interpolation
pdf_to_lqd <- function(pdf, quant, x_grid){
  interp_fn = approxfun(x_grid, pdf)
  lqd = -log(interp_fn(quant))
  return(lqd)
}

# alpha-mixture pdf with uniform
alpha_mix <- function(pdf, alpha){
  pdf_star = (1-alpha)*pdf + alpha*1
  return(pdf_star)
}

# Recover original pdf from alpha-mixed pdf
recover_pdf <- function(pdf_star, alpha, x_grid){
  W = trapz(x_grid, abs((pdf_star-alpha)/(1-alpha)))
  pdf = abs((pdf_star-alpha)/(1-alpha))/W
  return(pdf)
}

# Standard FPCA (fdapace)
std_fpca <- function(data_list, x, optns = list(dataType="Dense", error=FALSE,
                                                FVEthreshold=0.99)){
  n = length(data_list)
  m = length(data_list[[1]]$delta_star)
  
  Ly = list()
  Lt = list()
  for(idx in 1:n){
    Ly[[idx]] = data_list[[idx]]$delta_star
    Lt[[idx]] = x
  }
  
  # FPCA
  fpca_out = FPCA(Ly, Lt, optns=optns)
  n_pc = fpca_out$selectK
  delta_mean = fpca_out$mu
  Phi = fpca_out$phi
  
  # Append scores and FPC reconstructions of discrepancies
  for(idx in 1:n){
    data_list[[idx]]$bc = fpca_out$xiEst[idx,]
    data_list[[idx]]$delta_rec = delta_mean
    for(j in 1:n_pc){
      data_list[[idx]]$delta_rec = data_list[[idx]]$delta_rec + data_list[[idx]]$bc[j]*Phi[,j]
    }
    
    # Map to reconstructed DFN LQD # This is getting back to the LQD!!!!  NOT descrepencies.
    data_list[[idx]]$psiD_star_rec = data_list[[idx]]$psiG_star + data_list[[idx]]$delta_rec
  }
  
  fpca_obj = list(data=data_list, delta_mean=delta_mean, Phi=Phi, fpca_out=fpca_out)
  return(fpca_obj)
}


# Murph added function --> I need the same fPCA calculation, except now on the
# LQD of the graph features.
std_fpca_forfeats <- function(data_list, x, optns = list(dataType="Dense", error=FALSE,
                                                FVEthreshold=0.99), n_pc = NULL){
  n = length(data_list)
  m = length(data_list[[1]]$psife_star)
  
  Ly = list()
  Lt = list()
  for(idx in 1:n){
    Ly[[idx]] = data_list[[idx]]$psife_star
    Lt[[idx]] = x
  }
  
  # FPCA
  fpca_out = FPCA(Ly, Lt, optns=optns)
  if(is.null(n_pc)){
    n_pc = fpca_out$selectK
  }
  # 
  delta_mean = fpca_out$mu
  Phi = fpca_out$phi
  
  # Append scores and FPC reconstructions of discrepancies
  for(idx in 1:n){
    data_list[[idx]]$bc = fpca_out$xiEst[idx,]
    data_list[[idx]]$delta_rec = delta_mean
    for(j in 1:n_pc){
      data_list[[idx]]$delta_rec = data_list[[idx]]$delta_rec + data_list[[idx]]$bc[j]*Phi[,j]
    }
    
    # Murph 7/5/2023: We don't need this anymore, since we aren't working with descrepencies
    # with the graph features.
    # Map to reconstructed DFN LQD # This is getting back to the LQD!!!!  NOT descrepencies.
    # data_list[[idx]]$psiD_star_rec = data_list[[idx]]$psiG_star + data_list[[idx]]$delta_rec
  }
  
  fpca_obj = list(data=data_list, delta_mean=delta_mean, Phi=Phi, fpca_out=fpca_out)
  return(fpca_obj)
}


# Fit independent GP model given FPCA object
fit_gp <- function(fpca_obj, feat_type="quantile", n_feat=25, isotropic=TRUE, d=0.01, 
                   g_scl=0.01){
  data_list = fpca_obj$data
  n = length(data_list)
  m = length(data_list[[1]]$delta_star)
  delta_mean = fpca_obj$delta_mean
  Phi = fpca_obj$Phi
  n_pc = dim(Phi)[2]
  
  # Features (inputs)
  if(feat_type=="quantile"){
    # Equally-spaced quantiles across feature vector
    X = matrix(NA, nrow=n, ncol=n_feat)
    # input_idx = round(seq(from=1, to=length(data_list[[1]]$feat), length.out=n_feat))
    for(idx in 1:n){
      X[idx,] = quantile(log(data_list[[idx]]$feat), probs=seq(0,1,length.out=n_feat))
    }
    
    feat_fpca = NULL
  } else if(feat_type=="quantile_fpc"){
    # FPCA on quantile functions of feature vector
    Ly = list()
    Lt = list()
    for(idx in 1:n){
      Ly[[idx]] = log(data_list[[idx]]$feat)
      Lt[[idx]] = seq(0, 1, length.out=length(data_list[[1]]$feat))
    }
    
    feat_fpca = FPCA(Ly, Lt, optns = list(dataType="Dense", error=FALSE,
                                                FVEthreshold=0.999))
    nf_pc = feat_fpca$selectK
    feat_mean = feat_fpca$mu
    Phi_feat = feat_fpca$phi
    
    X = matrix(NA, nrow=n, ncol=nf_pc)
    for(idx in 1:n){
      X[idx,] = feat_fpca$xiEst[idx,]
    }
  }

  # Standardize
  cent = colMeans(X)
  scl = apply(X,2,sd)
  X_std = (X-cent)/scl
  
  # Outputs (PC basis coefficients in separate lists)
  pZ = list()
  p_mean = rep(NA,n_pc)
  p_sd = rep(NA,n_pc)
  for(k in 1:n_pc){
    pZ[[k]] = rep(NA, n)
    for(i in 1:n){
      pZ[[k]][i] = data_list[[i]]$bc[k]
    }
    p_mean[k] = mean(pZ[[k]])
    p_sd[k] = sd(pZ[[k]])
    pZ[[k]] = (pZ[[k]]-p_mean[k])/p_sd[k]
  }
  
  # Fit separate GPs
  pgp_obj = list()
  mle_pgp_obj = list()
  pred_pgp_obj = list()
  for(k in 1:n_pc){
    if(isotropic){  # isotropic
      pgp_obj[[k]] = newGP(X_std, pZ[[k]], d=d, g=g_scl*var(pZ[[k]]), dK=TRUE)
      mle_pgp_obj[[k]] = mleGP(pgp_obj[[k]], param="d")
    } else{  # anisotropic
      pgp_obj[[k]] = newGPsep(X_std, pZ[[k]], d=d, g=g_scl*var(pZ[[k]]), dK=TRUE)
      mle_pgp_obj[[k]] = mleGPsep(pgp_obj[[k]], param="d")
    }
    
    print(paste0("GP for PC ", k, " of ", n_pc, " fit"))
  }

  gp_obj = list(gp=pgp_obj, cent=cent, scl=scl, p_mean=p_mean, p_sd=p_sd, 
                X_std=X_std, isotropic=isotropic, n_feat=n_feat, feat_fpca=feat_fpca)
  return(gp_obj)
}

# Predictions using fitted GP model
pred_gp <- function(data_list, gp_obj, fpca_obj, x_grid, alpha, samp=FALSE, 
                    N_samp=200, feat_type="quantile",
                    trunc_err=matrix(0, nrow=length(data_list), ncol=length(data_list[[1]]$delta_star))){
  n = length(data_list)
  m = length(data_list[[1]]$delta_star)
  n_feat = gp_obj$n_feat
  
  pgp_obj = gp_obj$gp
  cent = gp_obj$cent
  scl = gp_obj$scl
  p_mean = gp_obj$p_mean
  p_sd = gp_obj$p_sd
  isotropic = gp_obj$isotropic
  n_pc = length(pgp_obj)
  
  # Construct features standardized in consistent way with model fitting
  if(feat_type=="quantile"){
    input_idx = round(seq(from=1, to=length(data_list[[1]]$feat), length.out=n_feat))
    X = matrix(NA, nrow=n, ncol=n_feat)
    for(idx in 1:n){
      X[idx,] = log(data_list[[idx]]$feat[input_idx])
    }
    X_std = (X-cent)/scl
  } else if(feat_type=="quantile_fpc"){
    Ly = list()
    Lt = list()
    for(idx in 1:n){
      Ly[[idx]] = log(data_list[[idx]]$feat)
      Lt[[idx]] = seq(0, 1, length.out=length(data_list[[1]]$feat))
    }
    X = predict(gp_obj$feat_fpca, Ly, Lt, xiMethod="IN")$scores
    X_std = (X-cent)/scl
  } else if(feat_type=="lqd_fpc"){
    Ly = list()
    Lt = list()
    for(idx in 1:n){
      Ly[[idx]] = log(data_list[[idx]]$feat)
      Lt[[idx]] = seq(0, 1, length.out=length(data_list[[1]]$feat))
    }
    X = predict(gp_obj$feat_fpca, Ly, Lt, xiMethod="IN")$scores
    X_std = (X-cent)/scl
  }
  
  # Predicted mean PC scores
  pred_pgp_obj = list()
  if(isotropic){
    for(k in 1:n_pc){
      pred_pgp_obj[[k]] = predGP(pgp_obj[[k]], X_std)
    }
  } else{
    for(k in 1:n_pc){
      pred_pgp_obj[[k]] = predGPsep(pgp_obj[[k]], X_std)
    }
  }
  
  # Predicted mean discrepancy, and convert to DFN LQD and cdf
  delta_mean = fpca_obj$delta_mean
  Phi = fpca_obj$Phi
  
  for(idx in 1:n){
    data_list[[idx]]$pred_bc_mean = rep(0,n_pc)
    tmp_del_mean = delta_mean
    for(k in 1:n_pc){
      bc_k = pred_pgp_obj[[k]]$mean[idx]*p_sd[k] + p_mean[k]
      data_list[[idx]]$pred_bc_mean[k] = bc_k
      tmp_del_mean = tmp_del_mean + bc_k*Phi[,k]
    }
    data_list[[idx]]$pred_delta_mean = tmp_del_mean + colMeans(trunc_err)
    
    # Map to DFN mixture LQD/CDF/PDF
    data_list[[idx]]$pred_psiD_star_mean = data_list[[idx]]$psiG_star + data_list[[idx]]$pred_delta_mean
    data_list[[idx]]$pred_FD_star_mean = lqd_to_cdf(data_list[[idx]]$pred_psiD_star_mean, x_grid)
    data_list[[idx]]$pred_fD_star_mean = cdf_to_pdf(data_list[[idx]]$pred_FD_star_mean, x_grid)
    
    # Map to DFN PDF/CDF
    data_list[[idx]]$pred_fD_mean = recover_pdf(data_list[[idx]]$pred_fD_star_mean, alpha, x_grid)
    data_list[[idx]]$pred_FD_mean = pdf_to_cdf(data_list[[idx]]$pred_fD_mean, x_grid)
  }
  
  # Predicted PC score samples, if desired
  if(samp){
    pbc_samp = array(NA, dim=c(N_samp,n,n_pc))
    for(k in 1:n_pc){
      tmp = rmvt(N_samp, pred_pgp_obj[[k]]$Sigma, pred_pgp_obj[[k]]$df)
      pbc_samp[,,k] = tmp + t(matrix(rep(pred_pgp_obj[[k]]$mean,N_samp), ncol=N_samp))
    }
    
    for(idx in 1:n){
      data_list[[idx]]$pred_bc_samp = matrix(0, nrow=N_samp, ncol=n_pc)
      data_list[[idx]]$pred_delta_samp = matrix(NA, nrow=N_samp, ncol=m)
      data_list[[idx]]$pred_psiD_star_samp = matrix(NA, nrow=N_samp, ncol=m)
      data_list[[idx]]$pred_FD_star_samp = matrix(NA, nrow=N_samp, ncol=m)
      data_list[[idx]]$pred_fD_star_samp = matrix(NA, nrow=N_samp, ncol=m)
      
      data_list[[idx]]$pred_FD_samp = matrix(NA, nrow=N_samp, ncol=m)
      data_list[[idx]]$pred_fD_samp = matrix(NA, nrow=N_samp, ncol=m)
      
      for(j in 1:N_samp){
        tmp_del_samp = delta_mean
        for(k in 1:n_pc){
          bc_k = pbc_samp[j,idx,k]*p_sd[k] + p_mean[k]
          data_list[[idx]]$pred_bc_samp[j,k] = bc_k
          tmp_del_samp = tmp_del_samp + bc_k*Phi[,k]
        }
        data_list[[idx]]$pred_delta_samp[j,] = tmp_del_samp + trunc_err[sample(1:n,1),]
        
        data_list[[idx]]$pred_psiD_star_samp[j,] = data_list[[idx]]$psiG_star + data_list[[idx]]$pred_delta_samp[j,]
        data_list[[idx]]$pred_FD_star_samp[j,] = lqd_to_cdf(data_list[[idx]]$pred_psiD_star_samp[j,], x_grid)
        data_list[[idx]]$pred_fD_star_samp[j,] = cdf_to_pdf(data_list[[idx]]$pred_FD_star_samp[j,], x_grid)
        
        data_list[[idx]]$pred_fD_samp[j,] = recover_pdf(data_list[[idx]]$pred_fD_star_samp[j,], alpha, x_grid)
        data_list[[idx]]$pred_FD_samp[j,] = pdf_to_cdf(data_list[[idx]]$pred_fD_samp[j,], x_grid)
      }
      
      if(idx%%25==0){
        print(paste0("observation ", idx, " of ", n, " completed"))
      }
    }
  }
  
  return(data_list)
}

##### Other functions not used currently #####
#####
# FPCA by smoothed covariance of discrepancies, and use to reconstruct
smooth_fpca <- function(data_list, n_bs=15, thresh=0.995){
  n = length(data_list)
  m = length(data_list[[1]]$delta)
  
  delta = matrix(NA, nrow=m, ncol=n)
  for(idx in 1:n){
    delta[,idx] = data_list[[idx]]$delta
  }
  
  # FPCA by smoothed covariance
  fpca_out = fpca.sc(Y=t(delta), nbasis=n_bs, pve=thresh)
  n_pc = fpca_out$npc
  delta_mean = as.numeric(fpca_out$mu)
  Phi = fpca_out$efunctions
  
  # Append scores and FPC reconstructions of discrepancies
  for(idx in 1:n){
    data_list[[idx]]$bc = fpca_out$scores[idx,]
    data_list[[idx]]$delta_rec = fpca_out$Yhat[idx,]
    
    # Map to reconstructed DFN LQD
    data_list[[idx]]$lqdD_rec = data_list[[idx]]$lqdG + data_list[[idx]]$delta_rec
  }
  
  fpca_obj = list(data=data_list, delta_mean=delta_mean, Phi=Phi, fpca_out=fpca_out)
  return(fpca_obj)
}


# Weighted FPCA
weighted_fpca <- function(data_list, x, alpha=0.98, optns = list(dataType="Dense", error=FALSE,
                                                                 FVEthreshold=0.99)){
  n = length(data_list)
  m = length(data_list[[1]]$delta)
  
  # Weighting
  W = alpha^(1:m)
  W = W/sum(W)
  W = diag(W)
  
  Ly = list()
  Lt = list()
  for(idx in 1:n){
    Ly[[idx]] = diag(W)*data_list[[idx]]$delta
    Lt[[idx]] = x
  }
  
  # FPCA
  fpca_out = FPCA(Ly, Lt, optns=optns)
  n_pc = fpca_out$selectK
  delta_mean = fpca_out$mu
  Phi = fpca_out$phi
  
  # Append scores and FPC reconstructions of discrepancies
  for(idx in 1:n){
    data_list[[idx]]$bc = fpca_out$xiEst[idx,]
    data_list[[idx]]$delta_rec = delta_mean
    for(j in 1:n_pc){
      data_list[[idx]]$delta_rec = data_list[[idx]]$delta_rec + data_list[[idx]]$bc[j]*Phi[,j]
    }
    data_list[[idx]]$delta_rec = (1/diag(W))*data_list[[idx]]$delta_rec
    
    # Map to reconstructed DFN LQD
    data_list[[idx]]$lqdD_rec = data_list[[idx]]$lqdG + data_list[[idx]]$delta_rec
  }
  
  fpca_obj = list(data=data_list, delta_mean=delta_mean, Phi=Phi, W=diag(W), fpca_out=fpca_out)
  return(fpca_obj)
}

# Weighted FPCA
jump_fpca <- function(data_list, x, alpha=0.98, optns = list(dataType="Dense", error=FALSE,
                                                             FVEthreshold=0.99)){
  n = length(data_list)
  m = length(data_list[[1]]$delta)
  
  # Initial mean estimate
  delta = matrix(NA, nrow=m, ncol=n)
  for(idx in 1:n){
    delta[,idx] = data_list[[idx]]$delta
  }
  init_delta_mean = rowMeans(delta)
  
  # Construct initial jump basis function
  Phi_jump = c(1,rep(0,m-1))
  
  # Basis coefficient for this initial jump
  bc_jump = rep(NA,n)
  for(idx in 1:n){
    bc_jump[idx] = trapz(x_full_std, (delta[,idx]-init_delta_mean)*Phi_jump)
  }
  
  
  Ly = list()
  Lt = list()
  for(idx in 1:n){
    Ly[[idx]] = data_list[[idx]]$delta
    Lt[[idx]] = x
  }
  
  # FPCA
  fpca_out = FPCA(Ly, Lt, optns=optns)
  n_pc = fpca_out$selectK
  delta_mean = fpca_out$mu
  Phi = fpca_out$phi
  
  # Append scores and FPC reconstructions of discrepancies
  for(idx in 1:n){
    data_list[[idx]]$bc = fpca_out$xiEst[idx,]
    data_list[[idx]]$delta_rec = delta_mean
    for(j in 1:n_pc){
      data_list[[idx]]$delta_rec = data_list[[idx]]$delta_rec + data_list[[idx]]$bc[j]*Phi[,j]
    }
    data_list[[idx]]$delta_rec = (1/diag(W))*data_list[[idx]]$delta_rec
    
    # Map to reconstructed DFN LQD
    data_list[[idx]]$lqdD_rec = data_list[[idx]]$lqdG + data_list[[idx]]$delta_rec
  }
  
  fpca_obj = list(data=data_list, delta_mean=delta_mean, Phi=Phi, W=diag(W), fpca_out=fpca_out)
  return(fpca_obj)
}