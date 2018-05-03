#################################################################################
##
## Author:  Nat Goodman
## Created: 18-05-03
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
## Utility functions for repwr.R
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################

## ---- Utility Functions ----
## generate name=value
paste_nv=function(name,value,sep='=') {
  name=as.character(pryr::subs(name));
  if (missing(value))
    if (exists(name,envir=parent.frame(n=1))) value=get(name,envir=parent.frame(n=1))
    else stop(paste('no value for',name,'in function call or parent environment'));
  paste(sep=sep,name,value); 
}
## generate list of name=value using values from parent environment. code adapted from base::rm
nvq=function(...,sep=' ') {
  dots=match.call(expand.dots=FALSE)$...
   if (length(dots) &&
     !all(vapply(dots,function(x) is.symbol(x) || is.character(x),NA,USE.NAMES=FALSE))) 
     stop("... must contain names or character strings");
  ## CAUTION: for some reason, doesn't work to use 'parent.frame(n=1)' inside sapply
  env=parent.frame(n=1);
  names=vapply(dots,as.character,"");
  values=sapply(names,function(name) {
    if (exists(name,envir=env)) get(name,envir=env)
    else stop(paste('no value for',name,'in parent environment'));
  });
  ## values=sapply(names,function(name)
  ##   if (exists(name,envir=parent.frame(n=2))) get(name,envir=parent.frame(n=2))
  ##   else stop(paste('no value for',name,'in parent environment')));
  paste(collapse=sep,mapply(function(name,value) paste(sep='=',name,value),names,values));
}
## tack id onto filebase if not NULL or NA
paste_id=function(base,id=NULL,sep='.') {
  ## test id this way to avoid running is.na when id=NULL 
  if (is.null(id)) return(base);
  if (is.na(id)) return(base);
  paste(sep=sep,base,id);
}  
## pretty print typical values of n, d & m
n_pretty=function(n) as.character(n);
d_pretty=function(d) sprintf('%3.2f',d);
m_pretty=function(m) {
  if (round(log10(m))==log10(m)) sub('e\\+0{0,1}','e',sprintf("%0.0e",m),perl=TRUE)
  else as.character(m);
}
## get value of variable from parent or set to default
## call with quoted or unquoted variable name
## if default missing, throws error if variable not found
parent=function(what,default) {
  what=as.character(pryr::subs(what));
  if (exists(what,envir=parent.frame(n=2))) return(get(what,envir=parent.frame(n=2)));
  if (!missing(default)) return(default);
  stop(paste(sep='',"object '",what,"' not found"));
}
## copy local variables to global - to simplify init
assign_global=function() {
  env=parent.frame(n=1);
  sapply(ls(envir=env),function(what) assign(what,get(what,envir=env),envir=.GlobalEnv));
}
## quote names in paramter list. code adapted from base::rm
cq=function(...) {
 dots=match.call(expand.dots=FALSE)$...
 if (length(dots) &&
     !all(vapply(dots,function(x) is.symbol(x) || is.character(x),NA,USE.NAMES=FALSE))) 
   stop("... must contain names or character strings");
return(vapply(dots,as.character,""));
}
## extend akima::aspline for matrix
asplinem=function(x,y,xout,...) {
  if (is.vector(y)) return(akima::aspline(x,y,xout,...));
  if (length(dim(y))!=2) stop('y must be vector or 2-dimensional matrix-like object');
  ## yout=apply(y,2,function(y) akima::aspline(x,y,xout,...)$y);
  yout=apply(y,2,function(y) akima::aspline(x,y,xout,...)$y);
  ## if yout has single row (ie, xout has one element), R turns it into a vector...
  if (length(xout)==1) yout=t(yout);
  yout;
}
## extend loess.smoth for matrix - probably only useful for plotting
loessm=function(x,y,xout,...) {
  if (is.vector(y)) y=data.frame(y=y);
  if (length(dim(y))!=2) stop('y must be vector or 2-dimensional matrix-like object');
  data=data.frame(x=x,y);
  yout=do.call(data.frame,lapply(colnames(y),function(name) {
    fmla=as.formula(paste(name,'~ x'));
    loess.obj=loess(fmla,data=data,...);
    yout=predict(loess(fmla,data=data),xout)
  }));
  ## if yout has single row (ie, xout has one element), R turns it into a vector...
  ## if (length(xout)==1) yout=t(yout);
  colnames(yout)=colnames(y);
  yout;
}
## not in - based on example in RefMan - more intutive than !%in%
"%notin%"=function(x,table) match(x,table,nomatch=0)==0
## between, near - to subset sim results
between=function(x,lo,hi) x>=lo&x<hi
near=function(x,target,tol=.01) between(x,target-tol,target+tol)

