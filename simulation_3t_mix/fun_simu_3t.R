

run_simu <- function(expr, scr.mat, scr.mat.true, t.net, sm, tnpar, settings){
  
  res.gla.t1 = glasso(var(expr[[1]]),rho = tnpar$lam1,penalize.diagonal = tnpar$penalize.diag,thr = 1e-4)
  res.gla.t2 = glasso(var(expr[[2]]),rho = tnpar$lam1,penalize.diagonal = tnpar$penalize.diag,thr = 1e-4)
  res.gla.t3 = glasso(var(expr[[3]]),rho = tnpar$lam1,penalize.diagonal = tnpar$penalize.diag,thr = 1e-4)
  theta.GLA = list(t1=res.gla.t1$wi,t2=res.gla.t2$wi,t3=res.gla.t3$wi)
  
  res.JGL.F = JGL(expr,penalty="fused",lambda1 = tnpar$lam1,lambda2 = tnpar$lam2,penalize.diagonal = tnpar$penalize.diag,return.whole.theta=T,tol = 1e-4)
  theta.JGLF = res.JGL.F$theta
  
  res.JGL.S = JGLS(Y = expr, lambda1 = tnpar$lam1, lambda2 = tnpar$lam2, btc.screening = scr.mat, penalize.diag = c(T,T),tol = 1e-4)
  theta.JGLS = res.JGL.S$theta
  
  res.JGL.O = JGLS(Y = expr, lambda1 = tnpar$lam1, lambda2 = tnpar$lam2, btc.screening = scr.mat.true, penalize.diag = c(T,T),tol = 1e-4)
  theta.JGLO = res.JGL.O$theta
  
  
  ####################
  resg <- list()
  resf <- list()
  ress <- list()
  reso <- list()
  
  # S
  ##############
  S <- cal_sample_cov(expr)
  
  ##############
  # loglik
  resg$loglik <- cal_loglik(theta.GLA,S,n = sm$sampn)
  resf$loglik <- cal_loglik(theta.JGLF,S,n = sm$sampn)
  ress$loglik <- cal_loglik(theta.JGLS,S,n = sm$sampn)
  reso$loglik <- cal_loglik(theta.JGLO,S,n = sm$sampn)
  
  ##############
  # AIC
  resg$aic <- cal_AIC(theta.GLA,resg$loglik)
  resf$aic <- cal_AIC(theta.JGLF,resf$loglik)
  ress$aic <- cal_AIC(theta.JGLS,ress$loglik)
  reso$aic <- cal_AIC(theta.JGLO,ress$loglik)
  
  ##############
  # BIC
  resg$bic <- cal_BIC(theta.GLA,resg$loglik,sm$sampn)
  resf$bic <- cal_BIC(theta.JGLF,resf$loglik,sm$sampn)
  ress$bic <- cal_BIC(theta.JGLS,ress$loglik,sm$sampn)
  reso$bic <- cal_BIC(theta.JGLO,ress$loglik,sm$sampn)
  
  ##############
  # SSE
  resg$sse <- get_SSE(theta.GLA,t.net$pcsm)
  resf$sse <- get_SSE(theta.JGLF,t.net$pcsm)
  ress$sse <- get_SSE(theta.JGLS,t.net$pcsm)
  reso$sse <- get_SSE(theta.JGLO,t.net$pcsm)
  
  ##############
  # ROC
  resg$roc <- get_FPR_TPR(theta.GLA,t.net$pcsm)
  resf$roc <- get_FPR_TPR(theta.JGLF,t.net$pcsm)
  ress$roc <- get_FPR_TPR(theta.JGLS,t.net$pcsm)
  reso$roc <- get_FPR_TPR(theta.JGLO,t.net$pcsm)
  
  ##############
  # s selection
  ress$s.selected <- scrmats[[i]]$s.sl
  ress$scr.n  <- scrmats[[i]]$scr.n
  
  simu.res <- list(GLA=resg, FGL=resf, FGLS=ress, FGLO=reso)
  return(simu.res)
}





