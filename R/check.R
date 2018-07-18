#################################################################################
##
## Author:  Nat Goodman
## Created: 18-05-04
##          from repwr.R restarted 18-02-15, created 17-10-05
##          from swfdr.R created 17-10-02
##          from swfdr_base.R created 17-04-10
##          from fdr.R created 17-04-03
##          from pval_error.R created 16-11-01
## Restart: 18-02-15
##          from scope.R created 17-12-04
##
## Copyright (C) 2018 Nat Goodman.
## 
## Error checking code for repwr.R
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################

## --- Error Checking Functions ---
## make sure smry.type initialized and provided types are valid
##   limit type to mesr
##   fix single valued one
check_type=function(type,mesr=parent(mesr,NULL),multiok=F) {
  ## initialize smry.type unless already done
  init_smry();
  ## make sure all types in smry.type
  bad=!(type %in% smry.type);
  if (any(bad)) {
    bad=paste(collapse=', ',type[bad]);
    stop(paste('Invalid smry types:',bad));
  }
  ## make sure from.type and relto.type singleton or vectors keyed by mesr
  if (length(unique(type))>1&&is.null(names(type)))
    stop('Cannot have multiple from or relto types unless keyed by mesr');
  ## extend singleton to vector of correct size keyed by mesr 
  if (is.null(names(type))) type=setNames(rep(type[1],length=length(mesr)),mesr)
  else {
    type=type[mesr];
    if (length(unique(type))>1&&!multiok) {
      bad=paste(collapse=', ',unique(type));
      stop(paste('Cannot have mutiple from or relto types unless multiok is TRUE:',bad));
    }}
  type;
}  
## make sure mesr initialized and provided mesrs are valid
check_mesr=function(mesr=parent(mesr,NULL)) {
  ## initialize mesr unless already done
  init_mesr();
  ## make sure all provided measures are good
  bad=!(mesr %in% mesr.all);
  if (any(bad)) {
    bad=paste(collapse=', ',mesr[bad]);
    stop(paste('Invalid measures:',bad));
  }
  T;
}
## make sure 'm' and 'datadir' consistent
## to defend against changing m on-the-fly without rerunning init
check_m=function(m=parent(m),datadir=parent(datadir)) {
  mpat=paste(sep='',paste_nv(m,m_pretty(m)),'$');
  if (!grepl(mpat,datadir))
    stop(paste(sep='','Inconsistent m and datadir: ',
               paste(sep=', ',paste_nv(m,m_pretty(m)),paste_nv(datadir)),
               '. Did you forget to rerun init after changing m?'));
  T;
}
