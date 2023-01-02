#' sel_MafChrPos
#' @description filt out some variants of which summary statistics are not included
#' @param ssp: a list of all study files which contain summary statistics data in different ancestries
#' @param sp: a vector which contains all variants names
#' @return a filtered list of all study files
#' @export

sel_MafChrPos<-function(ssp,sp){
  
  exp=paste("sel_stat_snp=data.frame(MARKERNAME=sp,CHR=NA,POS=NA,MARKERNUM=1:length(sp),",paste("EAF",1:4,"=NA",sep ="",collapse=","),")",sep="")
  eval(parse(text=exp))
  
  for(i_pop in 1:length(ssp)){
    ind_maf=as.matrix((ssp[[i_pop]]$EAF<=0.99)&(ssp[[i_pop]]$EAF>=0.01))#rownames(ind_maf)=ssp[[i_pop]]$MARKNERNAME
    ind_chr=as.matrix(ssp[[i_pop]]$CHROMOSOME>0)
    ind_pos=as.matrix(ssp[[i_pop]]$POSITION>0)
    sub_snp=sp[(sp%in%ssp[[i_pop]]$MARKERNAME)]
    loc=match(sub_snp,ssp[[i_pop]]$MARKERNAME)
    exp2=paste("sel_stat_snp$EAF",i_pop,"[(sp%in%ssp[[i_pop]]$MARKERNAME)]=(ssp[[i_pop]]$EAF[ind_maf])[loc]",sep="")
    eval(parse(text=exp2))
    sel_stat_snp$CHR[(sp%in%ssp[[i_pop]]$MARKERNAME)]=(ssp[[i_pop]]$CHROMOSOME[ind_chr])[loc]   
    sel_stat_snp$POS[(sp%in%ssp[[i_pop]]$MARKERNAME)]=(ssp[[i_pop]]$POSITION[ind_pos])[loc]      
  }
  
  sel=complete.cases(sel_stat_snp)
  return(sel_stat_snp[sel,])
}

#' dist_pop
#' @description Calculating the multi dimensional scaling
#' @param X: a matrix (N * the number of study files) of which entry is EAF
#' @return a matrix named double centering for multi-dimension scaling
#' @export
dist_pop<-function (X){
  
  ind_X=(!is.na(X))
  denomi<-ind_X%*%t(ind_X)
  nomi <- X %*% t(X)
  vec <- apply(X, 1, function(x) sum(x * x))
  d <- -2 * nomi + vec[col(nomi)] + vec[row(nomi)]
  diag(d) <- 0
  d <- d/denomi
  return(d)
}


#'double_center
#'@description double centering for multi-dimension scaling
#'@param D: matrix of multi-dimension scaling
#'@return a matrix through double centering
#'@export
double_center<-function(D){
  row_meanD=rowMeans(D)
  mean_dist=mean(row_meanD)
  D_center=(-0.5)*(D-row_meanD[col(D)]-row_meanD[row(D)]+mean_dist)
  return(D_center)
}
