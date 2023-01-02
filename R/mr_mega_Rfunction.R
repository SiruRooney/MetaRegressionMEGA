#'MR_mega
#'@description Meta regression
#'@param data_name: an input file which contains the list of all study files
#'@param data_loc:  a character vector of working directory
#'@param qt: an indicator of whether the trait is quatitative or not
#'@param pcCount: the number of principal components
#'@param nCores: the number of cores which are used for running parallel
#'@param parallel: an indicator of whether proceeed parallel
#'@import doParallel
#'@import parallel
#'@return a csv file which contains statistics and p values under different hyporthesis testings
#'@export
MR_mega<-function(data_name,data_loc=NULL,qt=TRUE,pcCount=1,nCores,parallel=TRUE){
  
  #sum_stat$Pop_i is data.frame
  sum_stat_pop<-readInputfile(data_name,data_loc=NULL)
  sum_stat_pop=check_eaf(sum_stat_pop)
  
  marker_name_pop=lapply(sum_stat_pop,function(ssp){
    return(unique(ssp$MARKERNAME))
  })
  
  num_snp=length(unique(unlist(marker_name_pop)))
  cat("There are ",num_snp,"different SNPs among",length(sum_stat_pop),"populations\n")
  snp_pop=unique(unlist(marker_name_pop))
  
  #Creat a distance matrix
  sel_stat_snp<-sel_MafChrPos(sum_stat_pop,snp_pop)
  sel_markernum=foreach(i_chr=1:23,.combine=rbind)%do%{
    sel_markernum=rep(0,300)
    sel_snp_pos=(sel_stat_snp[(sel_stat_snp$CHR==i_chr),]$POS)%/%(10^6)
    dup_last=(!duplicated(sel_snp_pos,fromLast=TRUE))
    sel_markernum[unique(sel_snp_pos,fromLast=TRUE)]=sel_stat_snp[(sel_stat_snp$CHR==i_chr),]$MARKERNUM[dup_last]
    return(sel_markernum)
  }
  
  cat("There are ",length(sel_markernum[sel_markernum!=0]),"independent variants for EAF correlation calculation.\n")
  
  #forDistance<-array(dim=c(length(sel_markernum[sel_markernum!=0]),length(sum_stat_pop)),dimnames=list(sel_stat_snp$MARKERNAME[sort(sel_markernum[sel_markernum!=0])],1:length(sum_stat_pop)))
  exp3=paste("forDistance=cbind(",paste("sel_stat_snp$EAF",1:4,"[(sel_markernum[sel_markernum!=0])]",sep="",collapse=","),")",sep="")
  #exp3=paste("forDistance=cbind(",paste("sel_stat_snp$EAF",1:4,"[sel_mrmega+1]",sep="",collapse=","),")",sep="")
  eval(parse(text=exp3))
  
  dim(forDistance)#25*4
  dist<-dist_pop(t(forDistance))
  dist<-double_center(dist)
  eigen_d<-eigen(dist)
  eigen_d$sqrt.value=sqrt(ifelse(eigen_d$values>0,eigen_d$values,0))
  pcs<-eigen_d$vectors%*%diag(eigen_d$sqrt.value)[,1:pcCount]

  
  marker_cohort_count=array(0,dim=c(num_snp,1),dimnames=list(snp_pop,"Count_cohorts"))
  eaf_filt=list()
  
  for (i_pop in 1:length(sum_stat_pop)){
    sub_snp=snp_pop[(snp_pop%in%sum_stat_pop[[i_pop]]$MARKERNAME)]
    eaf_filt[[i_pop]]=(sum_stat_pop[[i_pop]]$EAF!=-1)
    marker_cohort_count[(snp_pop%in%(sum_stat_pop[[i_pop]]$MARKERNAME[eaf_filt[[i_pop]]]))]=marker_cohort_count[(snp_pop%in%(sum_stat_pop[[i_pop]]$MARKERNAME[eaf_filt[[i_pop]]]))]+1
  }
  
  #filtering the SNPs
  
  snp_pop_filt=rownames(marker_cohort_count)[marker_cohort_count-2>pcCount]
  for (i_pop in 1:length(sum_stat_pop)){
    filt_threshold=(sum_stat_pop[[i_pop]]$MARKERNAME%in%snp_pop_filt)
    sum_stat_pop[[i_pop]]=sum_stat_pop[[i_pop]][filt_threshold,]
  }
  
  n_avgeaf<-getN_AverageEAF(sum_stat_pop,snp_pop_filt)
 
  beta_pop=data.frame(MARKERNAME=snp_pop_filt)
  invse2_pop=data.frame(MARKERNAME=snp_pop_filt)
  marker_inf_pop=data.frame(MARKERNAME=snp_pop_filt)
 
  for(i_pop in 1:length(sum_stat_pop)){
    
    beta_pop=merge(beta_pop,sum_stat_pop[[i_pop]][,c("MARKERNAME","BETA")],by="MARKERNAME",all.x=TRUE)
    invse2_pop=merge(invse2_pop,sum_stat_pop[[i_pop]][,c("MARKERNAME","SE")],by="MARKERNAME",all.x=TRUE)
    
    if(i_pop==1){
      marker_inf_pop= merge(marker_inf_pop,sum_stat_pop[[i_pop]][,c("MARKERNAME","CHROMOSOME","POSITION","EA","NEA")],by="MARKERNAME",all.x=TRUE)  
    }else{
      #ind_sam=(marker_inf_pop$MARKERNAME%in%sum_stat_pop[[i_pop]]$MARKERNAME)
      ind_marker=!(marker_inf_pop$MARKERNAME%in%sum_stat_pop[[i_pop-1]]$MARKERNAME)
      if(any(isTRUE(ind_marker))){
        marker_inf_pop[ind_marker,]= merge(marker_inf_pop[ind_marker,],sum_stat_pop[[i_pop]][,c("MARKERNAME","CHROMOSOME","POSITION","EA","NEA")],by="MARKERNAME",all.x=TRUE)  
      }
    }
  }
  
  names(beta_pop)=c("MARKERNAME",sprintf("Pop%s",1:length(sum_stat_pop)))
  names(invse2_pop)=c("MARKERNAME",sprintf("Pop%s",1:length(sum_stat_pop)))
  
  invse2_pop[-1]=round(invse2_pop[-1]^{-2})
  
  #Run linear regression in parallel and return results of hyporthese testings
  cl<- makeCluster(nCores,type="FORK")
  registerDoParallel(cl)
  #################################################################
  lr_out<-foreach(i_snp=1:length(snp_pop_filt),.combine=rbind)%dopar%{
    cat("i_snp=",i_snp,"\n")
    lr_out<-lr_w(unlist(beta_pop[i_snp,-1]),pcs,unlist(invse2_pop[i_snp,-1]))
    return(lr_out)
  }
  
  TSS_RSS=abs(lr_out$TSS-lr_out$RSS)
  RSS=abs(lr_out$RSS)
  TSS0_RSS=abs(lr_out$TSS0-lr_out$RSS)
  
  pModelHet=pchisq(TSS_RSS,pcCount,lower.tail = FALSE)
  #df_pModelHet=pcCount
  pResidHet=pchisq(RSS,lr_out$df_RSS,lower.tail = FALSE)
  #df_pResidHet=lr_out$df_RSS
  pModelTest=pchisq(TSS0_RSS,pcCount+1,lower.tail = FALSE)
  #df_pModelTest=pcCount+1
  
  logBF=0.5*(TSS0_RSS-(pcCount+1)*log(marker_cohort_count[marker_cohort_count-2>pcCount]))
  
  ##################################################################
  logout=data.frame(marker_inf_pop,n_avgeaf[,-1],
                    lr_out[,1:((pcCount+1)*2)],
                    chisq_association=TSS0_RSS,
                    ndf_association=pcCount+1,
                    pvalue_association=pModelTest,
                    chisq_anceheter=TSS_RSS,
                    ndf_chisq_anceheter=pcCount,
                    pvalue_anceheter=pModelHet,
                    chisq_residual=RSS,
                    ndf_residual=lr_out$df_RSS,
                    pvalue_residual=pResidHet,
                    logBF=logBF)
  
  write.csv(logout,file="mr_mega_out.txt")
  
}


#' getN_AverageEAF
#' @description summarize the sample size and eaf
#' @param ssp: the revised list of all study files
#' @param ssf: the filtered variants names
#' @return a list which contains the sum of sample size among all studies and the average of eaf
#' @export
getN_AverageEAF<-function(ssp,ssf){
  cl<- detectCores()
  registerDoParallel(cl)
  
  n_eaf_pop=foreach(i_pop=1:length(ssp))%do%{
    ind_marker=(ssp[[i_pop]]$MARKERNAME%in%ssf)
    dtf=data.frame(MARKERNAME=ssp[[i_pop]][ind_marker,]$MARKERNAME,N=ssp[[i_pop]][ind_marker,]$N,sum_eaf=ssp[[i_pop]][ind_marker,]$N*ssp[[i_pop]][ind_marker,]$EAF)
    return(dtf)  
  }
  n_eaf<-data.frame(MARKERNAME=ssf,N=0,EAF_avg=0)
  for(i_pop in 1:length(ssp)){
    loc=match(n_eaf$MARKERNAME,n_eaf_pop[[i_pop]]$MARKERNAME)
    n_eaf$N=n_eaf$N+n_eaf_pop[[i_pop]]$N[loc]
    n_eaf$EAF_avg=n_eaf$EAF_avg+n_eaf_pop[[i_pop]]$sum_eaf[loc]
  }
  
  n_eaf$EAF_avg=n_eaf$EAF_avg/n_eaf$N
  return(n_eaf)
}



#' linear regression
#'@description linear regression for each variants
#'@param Y: BETA
#'@param X: pcs principal components
#'@param W: SE standard errors
#'@return a list of regression result for one variant
#'@export

lr_w<-function(Y,X,W){
  #Y referes to BETA
  #X refers to pcs
  #W refers to SE
  m0<-lm(Y~1,weights=W)
  m1<-lm(Y~X,weights=W)
  
  RSS=sum(W*(m1$residuals^2))
  TSS0=sum(W*Y^2)
  TSS=sum(W*m0$residuals^2)
  
  sum_m1<-summary(m1)
  coef_se<-c(t(sum_m1$coefficients[,1:2]))
  names(coef_se)=sprintf(rep(c("beta%s","se%s"),NROW(sum_m1$coefficients)),rep(0:(NROW(sum_m1$coefficients)-1),each=2))
  coef.se=as.data.frame(as.list(coef_se),col.names=names(coef_se))
  
  return(data.frame(coef.se,RSS=RSS,df_RSS=m1$df.residual,TSS0=TSS0,TSS=TSS))  
}
