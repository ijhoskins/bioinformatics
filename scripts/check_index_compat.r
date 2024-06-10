#!/usr/bin/env Rscript

# Checks index sequences for compatibility in sequencing

library(argparser, quietly=TRUE)
library(DescTools, quietly=TRUE, warn.conflicts=FALSE)
library(data.table, quietly=TRUE)

parse_cline_args<- function(){

  p<- arg_parser("Check index compatibility")
  p<- add_argument(p, "--indices", help="Text file with columns {Library, i5, i7}")
  p<- add_argument(p, "--levenshtein", help="Use edit distance instead of Hamming distance", flag=TRUE)
  p<- add_argument(p, "--dist", help="Distance at/below which indices are flagged", type="integer", default=1)
  argv<- parse_args(p)
  return(argv)
}

workflow<- function(){
  
  message("Starting index compatibility check")
  
  argv<- parse_cline_args()
  
  index_dt<- fread(argv$indices, header=TRUE)
  index_dt[,i5_i7:=paste(i5, i7, sep=",")]
  
  dist_method<- "hamming"
  if(argv$levenshtein){
    dist_method<- "levenshtein"
  }
  
  # Compute all pairwise distances
  nlibs<- nrow(index_dt)
  dist_mat<- matrix(nrow=nlibs, ncol=nlibs)
  for(i in 1:nlibs){
    for(j in 1:nlibs){
      
      if(i==j){
        next
      }
      
      dist_mat[i,j]<- StrDist(
        x=index_dt[i, i5_i7], y=index_dt[j, i5_i7], method=dist_method, ignore.case=TRUE)
    }
  }
  
  rownames(dist_mat)<- index_dt[,Library]
  colnames(dist_mat)<- index_dt[,Library]
  
  is_incompat_index<- function(x){
    res<- which(x <= argv$dist)
    return(res)
  }
  
  # Determine if indices are conflicting, and if so, which
  dist_mat_incompat_list<- apply(dist_mat, 1, is_incompat_index)
  dist_mat_incompat<- sapply(dist_mat_incompat_list, function(x){length(x) > 0})
  
  if(all(!dist_mat_incompat)){
    message("All index sequences are compatible.")
    return(NULL)
  }
    
  
  dist_mat_incompat_list_filt<- dist_mat_incompat_list[dist_mat_incompat]
  
  # Combine names of conflicting libraries
  get_lib_name<- function(x){
    lib_name<- index_dt[x, Library]
    return(lib_name)
  }
  
  dist_mat_incompat_list_filt_pair<- sapply(dist_mat_incompat_list_filt, get_lib_name)
  dist_mat_incompat_list_filt_pair_names<- names(dist_mat_incompat_list_filt_pair)
  
  # Finally uniquify the pairs
  # One might think our set should always be even due to symmetry, so one could take the first half of values
  # However, this is not the case for indices that may match more than one other index, which may be nonsymmetric
  cat_ragged_libs<- function(key_index){
    
    key_lib<- dist_mat_incompat_list_filt_pair_names[key_index]
    val_libs<- dist_mat_incompat_list_filt_pair[[key_index]]
    
    libs<- vector(mode="character", length=length(val_libs))
    for(i in 1:length(val_libs)){
      val_lib<- val_libs[i]
      libs[i]<- paste(key_lib, val_lib, sep=",")
    }
    return(libs)
  }

  incompat_libs<- unlist(lapply(seq_along(dist_mat_incompat_list_filt_pair), cat_ragged_libs))
  
  # Finally filter out redundant pairs
  incompat_libs_filt<- NULL
  for(i in 1:length(incompat_libs)){
    
    incompat_lib<- incompat_libs[i]
    incompat_lib_split<- unlist(strsplit(incompat_lib, ","))
    incompat_lib_rev<- paste(incompat_lib_split[2], incompat_lib_split[1], sep=",")
    
    if(!incompat_lib_rev%in%incompat_libs_filt){
      incompat_libs_filt<- c(incompat_libs_filt, incompat_lib)
    }
  }
  
  # Print the results
  msg<- paste("Libraries with index incompatibility", paste(
    incompat_libs_filt, collapse=" | "), sep=": ")
  
  message(msg)
  message("Completed index compatibility check")
  
  return(dist_mat)
}


if(getOption('run.workflow', default=TRUE)){
  dist_mat<- workflow()
}

