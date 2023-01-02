#' readInputfile
#' @description read file which contains names of different ancestry groups summary statistics
#' @param data_name: an input file which contains the list of all study files
#' @param qt: an indicator of whether the trait is quatitative or dichotomous
#' @param data_loc: a character vector of working directory
#' @import R.utils
#' @import data.table
#' @return a list of all study files which contain summary statistics data in different ancestries
#' @export
readInputfile<-function(data_name,qt=TRUE,data_loc=NULL){
  sum_dat<-list()
  
  if(is.null(data_loc)){
    data_loc=getwd()
  }
  
  if(isTRUE(qt)){
    sum_stat_nd=c("MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","BETA","SE","N")
  }else{
    sum_stat_nd=c("MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","OR","OR_95L","OR_95U","N")
  }
  
  if (!is.null(data_name)){
    for (i in 1:length(data_name)){
      sum_dat[[i]]=data.table::fread(paste(data_loc,"/",data_name[i],sep=""),select=sum_stat_nd)#[,pop_sum_stat%in%sum_stat_nd]
    }
  }else{
    print("data_name should be a file containing names of several summary statistics files")
  }
  names(sum_dat)=paste("Pop",1:length(data_name),sep="")
  return(sum_dat)
}




#' check_eaf

#'@description Check whether all efficient groups share the same efficient allele
#'@param ssp: a list of all study files which contain summary statistics data in different ancestries

#'@return the revised list of all study files
#'@export

check_eaf<-function(ssp){
  
  for(i_pop in 1:(length(ssp)-1)){
    
    cover=ssp[[i_pop+1]]$MARKERNAME%in%ssp[[i_pop]]$MARKERNAME
    loc=match(ssp[[i_pop]]$MARKERNAME,ssp[[i_pop+1]]$MARKERNAME[cover])
    same_ea=(ssp[[i_pop+1]]$EA[cover][loc]!=ssp[[i_pop]]$EA)
    if(any(same_ea)){
      ref_ea=ssp[[i_pop]]$EA[same_ea]
      ssp[[i_pop+1]]$NEA[cover][loc][same_ea]=ssp[[i_pop+1]]$EA[cover][loc][same_ea]
      ssp[[i_pop+1]]$EA[cover][loc][same_ea]=ref_ea
      ssp[[i_pop+1]]$EAF[cover][loc][same_ea]=1-ssp[[i_pop+1]]$EAF[cover][loc][same_ea]
      ssp[[i_pop+1]]$BETA[cover][loc][same_ea]=-ssp[[i_pop+1]]$BETA[cover][loc][same_ea]
    }
  }
  return(ssp)
}




