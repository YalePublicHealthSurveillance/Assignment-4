step_func <- function(ds,
                      time_points=time_points1,
                      post_period1=c('2009-08-01', '2010-08-01'),
                      post_period2=c('2010-08-01', '2012-05-01'),
                      post_period3=c('2012-05-01', '2017-12-01'),
                      vax.vars=c('post1','post2', 'post3') ,
                      other.covars='none' ,
                      mod='pois',
                      denom, 
                      outcome_name ){
  
  ds <- ds[order(ds$date),]
  ds$time_index<-(1:nrow(ds))/nrow(ds)
  
  #Create the dummy variables for 3 post-vaccine time periods
  ds$post1<-as.numeric(ds$date >= post_period1[1] & ds$date <  post_period1[2]) 
  ds$post2<-as.numeric(ds$date >= post_period2[1] & ds$date <  post_period2[2]) 
  ds$post3<-as.numeric(ds$date >= post_period3[1] ) 
  
  ds$time_index<-1:nrow(ds)
  ds$sin12<-sin(2*pi* ds$time_index/12)
  ds$cos12<-cos(2*pi* ds$time_index/12)
  ds$sin6<-sin(2*pi* ds$time_index/6)
  ds$cos6<-cos(2*pi* ds$time_index/6)
  ds$month<-month(ds$date)
  #month.dummies<-as.data.frame(dummies::dummy(ds$month))[,1:11]
  month.dummies <- model.matrix(~ as.factor(month), data=ds)
  month.dummies <- as.data.frame(month.dummies[,-1])
  names(month.dummies) <- paste0('month', 2:(ncol(month.dummies)+1))
  
  ds$one<-1
  
  ds<-cbind.data.frame(ds, month.dummies)
  
  ds$obs <- as.factor(1:nrow(ds))
  ds$log.offset<-log(ds[,denom]+0.5)
  
  seas.vars<- c(names(month.dummies) )
  offset.vars<-'offset(log.offset)'
  if((other.covars=='none')[1]){
    mod.vars<-c(seas.vars,vax.vars)
    mod.covars.cf <- seas.vars
  }else{
    mod.vars<-c(seas.vars,other.covars,vax.vars)
    mod.covars.cf<-c(seas.vars,other.covars)
    for(i in 1:length(other.covars)){
      ds[,other.covars[i]]<-scale(log( ds[,other.covars[i]]+0.5)) 
    }
  }
  
  #Re-define time_index (it gets overwritten in loop above if other covars present)
  ds$time_index<-(1:nrow(ds))/nrow(ds)
  
  form1<-as.formula(paste0(outcome_name, '~'   ,paste0(c(mod.vars,offset.vars), collapse='+')))
  
  #Fit model
  if(mod=='pois'){
    mod1 <-
      glm(form1,
          data = ds, family='poisson'
      )
  }else if(mod=='negbin'){
    mod1 <-
      glm.nb(form1,
             data = ds
      )
  }
  aic1<-AIC(mod1)
  deviance1<-summary(mod1)$deviance
  df1<-summary(mod1)$df[2]
  overdispersion<-deviance1/df1
  
  #GENERATE PREDICTIONS
  covars3 <-
    as.matrix(cbind(ds[, mod.vars])) 
  covars3 <- cbind.data.frame(rep(1, times = nrow(covars3)), covars3)
  names(covars3)[1] <- "Intercept"
  
  pred.coefs.reg.mean <-
    mvrnorm(n = 1000,
            mu = coef(mod1),
            Sigma = vcov(mod1))
  
  preds.stage1.regmean <-
    exp(as.matrix(covars3) %*% t(pred.coefs.reg.mean) +ds$log.offset)
  
  preds.q<-t(apply(preds.stage1.regmean,1,quantile, probs=c(0.025,0.5,0.975)))
  
  #Then for counterfactual, set post-vax effects to 0.
  covars3.cf <-
    as.matrix(cbind(ds[, c(
      mod.covars.cf
    )], matrix(
      0, nrow = nrow(ds), ncol = length(vax.vars)
    )))
  covars3.cf <-
    cbind.data.frame(rep(1, times = nrow(covars3.cf)), covars3.cf)
  preds.stage1.regmean.cf <-    exp(as.matrix(covars3.cf) %*% t(pred.coefs.reg.mean)+ ds$log.offset)
  preds.cf.q<-t(apply(preds.stage1.regmean.cf,1,quantile, probs=c(0.025,0.5,0.975)))
  
  rr.t <- preds.stage1.regmean / preds.stage1.regmean.cf
  rr.q.t <- t(apply(rr.t, 1, quantile, probs = c(0.025, 0.5, 0.975)))
  
  last.t<-nrow(rr.t) #evaluate at last time point
  preds.stage1.regmean.SUM <-   preds.stage1.regmean[last.t, ]
  preds.stage1.regmean.cf.SUM <-preds.stage1.regmean.cf[last.t, ]
  rr.post <- preds.stage1.regmean.SUM / preds.stage1.regmean.cf.SUM
  rr.q.post <- quantile(rr.post, probs = c(0.025, 0.5, 0.975))
  
  rr.out <- list('rr.q.post' = rr.q.post, 'aic1'=aic1,'outcome'=ds[,outcome_name],
                 'preds.cf.q'=preds.cf.q,'preds.q'=preds.q,
                 'rr.q.t'=rr.q.t,'overdispersion'=overdispersion, 'dates'=ds$date)
  return(rr.out)
}

plot.step.func<-function(ds){
  par(mfrow=c(1,2))
  ylimits<-range(c(0,ds$outcome,ds$preds.q[,'50%'],ds$preds.cf.q[,'50%'] ))
  plot(ds$dates,ds$outcome , bty='l', ylab='N cases',xlab='date', pch=16, ylim=ylimits)
  points(ds$dates, ds$preds.q[,'50%'], type='l', col='red')
  points(ds$dates, ds$preds.cf.q[,'50%'], type='l', col='blue')
  
  plot(ds$dates,ds$rr.q.t[,'50%'] , bty='l', ylab='Rate ratio',xlab='date', type='l',ylim=c(0.2,2))
  points(ds$dates,ds$rr.q.t[,'2.5%'], lty=2, col='gray', type='l')
  points(ds$dates,ds$rr.q.t[,'97.5%'], lty=2, col='gray', type='l')
  abline(h=1, lty=2, col='gray')
}


########################################
spline_func <- function(ds,
                        time_points=time_points1,
                        post_period1=c('2009-08-01', '2010-08-01'),
                        post_period2=c('2010-08-01', '2012-05-01'),
                        post_period3=c('2012-05-01', '2017-12-01'),
                        vax.vars=c('post1','post2', 'post3') ,
                        other.covars='none' ,
                        mod='pois',
                        denom, 
                        outcome_name ){
  ds <- ds[order(ds$date),]
  
  ds$time_index<-(1:nrow(ds))/nrow(ds)
  
  #Create the spline variables for 3 post-vaccine time periods
  ds$post1<-ds$time_index - ds$time_index[ds$date == post_period1[1]] 
  ds$post1[ds$post1<0]<-0
  ds$post2<-ds$time_index - ds$time_index[ds$date == post_period2[1]] 
  ds$post2[ds$post2<0]<-0
  ds$post3<-ds$time_index - ds$time_index[ds$date == post_period3[1]] 
  ds$post3[ds$post3<0]<-0
  #level out after last time point--hold all 3 to same value as at end of period 3
  #ds$post1[time_points >= post_period3[2]]<- ds$post1[time_points == post_period3[2]] 
  #ds$post2[time_points >= post_period3[2]]<- ds$post2[time_points == post_period3[2]] 
  #ds$post3[time_points >= post_period3[2]]<- ds$post3[time_points == post_period3[2]] 
  
  
  ds$time_index<-1:nrow(ds)
  ds$sin12<-sin(2*pi* ds$time_index/12)
  ds$cos12<-cos(2*pi* ds$time_index/12)
  ds$sin6<-sin(2*pi* ds$time_index/6)
  ds$cos6<-cos(2*pi* ds$time_index/6)
  ds$month<-month(ds$date)
  # month.dummies<-as.data.frame(dummies::dummy(ds$month))[,1:11]
  
  #month.dummies<-as.data.frame(dummies::dummy(ds$month))[,1:11]
  month.dummies <- model.matrix(~ as.factor(month), data=ds)
  month.dummies <- as.data.frame(month.dummies[,-1])
  names(month.dummies) <- paste0('month', 2:(ncol(month.dummies)+1))
  
  
  ds$one<-1
  ds<-cbind.data.frame(ds, month.dummies)
  
  ds$obs <- as.factor(1:nrow(ds))
  ds$log.offset<-log(ds[,denom]+0.5)
  
  seas.vars<- c(names(month.dummies) )
  offset.vars<-'offset(log.offset)'
  if((other.covars=='none')[1]){
    mod.vars<-c(seas.vars,vax.vars)
    mod.covars.cf <- seas.vars
  }else{
    mod.vars<-c(seas.vars,other.covars,vax.vars)
    mod.covars.cf<-c(seas.vars,other.covars)
    for(i in 1:length(other.covars)){
      ds[,other.covars[i]]<-scale(log( ds[,other.covars[i]]+0.5)) 
    }
  }
  
  #Re-define time_index (it gets overwritten in loop above if other covars present)
  ds$time_index<-(1:nrow(ds))/nrow(ds)
  
  form1<-as.formula(paste0(outcome_name, '~'   ,paste0(c(mod.vars,offset.vars), collapse='+')))
  
  #Fit model
  if(mod=='pois'){
    mod1 <-
      glm(form1,
          data = ds, family='poisson'
      )
  }else if(mod=='negbin'){
    mod1 <-
      glm.nb(form1,
             data = ds
      )
  }
  aic1<-AIC(mod1)
  deviance1<-summary(mod1)$deviance
  df1<-summary(mod1)$df[2]
  overdispersion<-deviance1/df1
  #GENERATE PREDICTIONS
  covars3 <-
    as.matrix(cbind(ds[, mod.vars])) 
  covars3 <- cbind.data.frame(rep(1, times = nrow(covars3)), covars3)
  names(covars3)[1] <- "Intercept"
  pred.coefs.reg.mean <-
    mvrnorm(n = 1000,
            mu = coef(mod1),
            Sigma = vcov(mod1))
  preds.stage1.regmean <-
    exp(as.matrix(covars3) %*% t(pred.coefs.reg.mean) +ds$log.offset)
  preds.q<-t(apply(preds.stage1.regmean,1,quantile, probs=c(0.025,0.5,0.975)))
  
  #Then for counterfactual, set post-vax effects to 0.
  covars3.cf <-
    as.matrix(cbind(ds[, c(
      mod.covars.cf
    )], matrix(
      0, nrow = nrow(ds), ncol = length(vax.vars)
    )))
  covars3.cf <-
    cbind.data.frame(rep(1, times = nrow(covars3.cf)), covars3.cf)
  preds.stage1.regmean.cf <-    exp(as.matrix(covars3.cf) %*% t(pred.coefs.reg.mean)+ ds$log.offset)
  preds.cf.q<-t(apply(preds.stage1.regmean.cf,1,quantile, probs=c(0.025,0.5,0.975)))
  
  rr.t <- preds.stage1.regmean / preds.stage1.regmean.cf
  rr.q.t <- t(apply(rr.t, 1, quantile, probs = c(0.025, 0.5, 0.975)))
  
  last.t<-nrow(rr.t) #evaluate at last time point
  preds.stage1.regmean.SUM <-   preds.stage1.regmean[last.t, ]
  preds.stage1.regmean.cf.SUM <-preds.stage1.regmean.cf[last.t, ]
  rr.post <- preds.stage1.regmean.SUM / preds.stage1.regmean.cf.SUM
  rr.q.post <- quantile(rr.post, probs = c(0.025, 0.5, 0.975))
  
  rr.out <- list('rr.q.post' = rr.q.post, 'aic1'=aic1,'outcome'=ds[,outcome_name],
                 'preds.cf.q'=preds.cf.q,'preds.q'=preds.q,
                 'rr.q.t'=rr.q.t,'overdispersion'=overdispersion, 'dates'=ds$date)
  return(rr.out)
}

# plot.step.func<-function(ds){
#   par(mfrow=c(1,2))
#   ylimits<-range(c(ds$outcome,ds$preds.q[,'50%'],ds$preds.cf.q[,'50%'] ))
#   plot(time_points1,ds$outcome , bty='l', ylab='N cases',xlab='date', pch=16, ylim=ylimits)
#   points(time_points1, ds$preds.q[,'50%'], type='l', col='red')
#   points(time_points1, ds$preds.cf.q[,'50%'], type='l', col='blue')
#   
#   plot(time_points1,ds$rr.q.t[,'50%'] , bty='l', ylab='Rate ratio',xlab='date', type='l',ylim=c(0.2,2))
#   points(time_points1,ds$rr.q.t[,'2.5%'], lty=2, col='gray', type='l')
#   points(time_points1,ds$rr.q.t[,'97.5%'], lty=2, col='gray', type='l')
#   abline(h=1, lty=2, col='gray')
# }
