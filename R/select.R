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
## Code to select data for analysis and plotting for repwr.R
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################

## --- Data Selection and Interpolation Functions ---
## filter and interpolate smry to get the data needed for plotting and sort based on x
smry_select=
  function(smry=parent(smry),n1=parent(n1),n2=parent(n2),d1=parent(d1),d2=parent(d2),
           mesr=parent(mesr),x=parent(x),
           from.type=parent(from.type),relto.type=parent(relto.type)) {
    smry=smry_prune();
    smry=smry_interp();
    ## filter smry and sort based on x
    ## in database parlance, the select part is a natural semijoin
    ## CAUTION: need cbind in case params have incompatible lengths. sigh...
    xdata=suppressWarnings(data.frame(cbind(n1,n2,d1,d2)));
    xdata$i=seq_len(nrow(xdata));
    smry=merge(smry,xdata);
    ## sort based on x
   if (x=='auto') {
      ## set x to first var that can drive loop, else 'asis'
      x=do.call(
        c,sapply(cq(n1,n2,d1,d2),simplify=F,
                 function(x) if (length(unique(xdata[[x]]))==nrow(xdata)) x));
      if (is.null(x)) x='asis' else x=x[1];
    }
    smry=switch(x,
                asis=smry[with(smry,order(i,type)),],
                n1=smry[with(smry,order(n1,n2,d1,d2,type)),],
                n2=smry[with(smry,order(n2,n1,d1,d2,type)),],
                d1=smry[with(smry,order(d1,d2,n1,n2,type)),],
                d2=smry[with(smry,order(d2,d1,n1,n2,type)),]);
    list(smry=subset(smry,select=-i),x=x);
  }
## prune smry so interp will be faster
smry_prune=
  function(smry=parent(smry),n1=parent(n1),n2=parent(n2),d1=parent(d1),d2=parent(d2),
           mesr=parent(mesr),
           from.type=parent(from.type),relto.type=parent(relto.type)) {
    n1.lo=max(smry$n1[smry$n1<=min(n1)]);
    n1.hi=min(smry$n1[smry$n1>=max(n1)]);
    n2.lo=max(smry$n2[smry$n2<=min(n2)]);
    n2.hi=min(smry$n2[smry$n2>=max(n2)])
    d1.lo=max(smry$d1[smry$d1<=min(d1)]);
    d1.hi=min(smry$d1[smry$d1>=max(d1)]);
    d2.lo=max(smry$d2[smry$d2<=min(d2)]);
    d2.hi=min(smry$d2[smry$d2>=max(d2)]);
    types=unique(c(from.type,relto.type));
    smry=subset(smry,subset=(n1>=n1.lo&n1<=n1.hi&n2>=n2.lo&n2<=n2.hi&
                             d1>=d1.lo&d1<=d1.hi&d2>=d2.lo&d2<=d2.hi&
                             type%in%types),
                select=c(cq(n1,n2,d1,d2,type,m),mesr));
  }
## interpolate smry one parameter at a time
smry_interp=
  function(smry=parent(smry),n1=parent(n1),n2=parent(n2),d1=parent(d1),d2=parent(d2),
           mesr=parent(mesr)) {
    smry=smry_interp1(n1);
    smry=smry_interp1(n2);
    smry=smry_interp1(d1);
    smry=smry_interp1(d2);
    ## remove rownames contructed by R and put in order we want
    rownames(smry)=NULL;
    smry[with(smry,order(n1,n2,d1,d2,type)),c(cq(n1,n2,d1,d2,type,m),mesr)];
  }
smry_interp1=function(x,xout,smry=parent(smry)) {
  x=as.character(pryr::subs(x));
  if (missing(xout))
    if (exists(x,envir=parent.frame(n=1))) xout=get(x,envir=parent.frame(n=1))
    else stop(paste('no value for',x,'in function call or parent environment'));
  xout=unique(xout);
  if (all(xout%in%smry[,x])) {
    ## no need to interpolate. subset will do.
    smry=smry[smry[,x]%in%xout,];
  } else {
    ## get selectors for parameter and data columns
    param=cq(n1,n2,d1,d2,type);
    fix=param[param %notin% x];
    data=colnames(smry) %notin% param;
    ## split smry by fixed params
    smry.fix=smry[,fix];
    smry.byfix=split(smry,apply(smry.fix,1,function(row) paste(collapse=' ',row)));
    ## interpolate each group over x
    smry=do.call(rbind,lapply(smry.byfix,function(smry) {
      yout=asplinem(smry[,x],smry[,data],y,xout=xout);
      smry=suppressWarnings(data.frame(xout,smry[1,fix],yout)); # tack params onto output
      colnames(smry)[1]=x;                                        # and fix names
      smry;
    }));
  }
  smry;
}
## make sure from.type & relto.type are valid
##   limit type to mesr
##   fix single valued one
check_type=function(mesr=parent(mesr,NULL),type,multiok=F) {
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
## make sure smry.type and mesr initialized and all provided measures are legal
check_mesr=
  function(mesr=parent(mesr,NULL),from.type=parent(from.type),relto.type=parent(relto.type),
           relto.multiok=parent(relto.multiok,F)) {
    ## make sure all provided measures are good
    bad=!(mesr %in% mesr.all);
    if (any(bad)) {
      bad=paste(collapse=', ',mesr[bad]);
      stop(paste('Invalid measures:',bad));
    }
    ## make sure all provided smry types are good
    type=unique(c(from.type,relto.type));
    bad=!(type %in% smry.type);
    if (any(bad)) {
      bad=paste(collapse=', ',type[bad]);
      stop(paste('Invalid smry types:',bad));
    }
    ## make sure from.type and relto.type singleton or vectors keyed by mesr
    if (length(unique(from.type))>1&&
        (is.null(names(from.type))||any(names(from.type) %notin% mesr.all)))
      stop('Cannot have multiple from.types unless keyed by mesr');
    if (length(unique(relto.type))>1&&
        (is.null(names(relto.type))||any(names(relto.type) %notin% mesr.all)))
      stop('Cannot have multiple relto.types unless keyed by mesr');
    ## make sure relto.type singleton unless relto.multiok is T
    if (length(unique(relto.type))>1) {
      relto.type=unique(relto.type[mesr]);
      if (length(relto.type)>1&!relto.multiok) {
        bad=paste(collapse=', ',relto.type);
        stop(paste('Cannot have mutiple relto.types unless relto.multiok is TRUE:',bad));
      }}
    T;
  }
## compute boolean vector of which rows are supposed to say 'yes' and which 'no'
## note various flavors of significance criteria - not sure which is best
##   nonzro: default non-zero rate. currently nonz1
##   nonz1: d1 non-zero; d2 doesn't matter
##   nonz1or2: d1 or d2 non-zero
##   nonz1and2: d1 and d2 non-zero
##   sameff: d1, d2 equal
##   uni1: d1, d2 equal and non-zero
x_yes=function(pos.rate,rate) {
  with(pos.rate$x,switch(rate,
                       nonzro=d1!=0,
                       nonz1=d1!=0,
                       nonz1or2=(d1!=0|d2!=0),
                       nonz1and2=(d1!=0&d2!=0),
                       sameff=(d1==d2),
                       uni=(d1!=0&d1==d2)));
}
