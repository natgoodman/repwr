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
## post-process posr to get the data needed for plotting or analysis
##   filter and interpolate, determine rate type, convert rates per type
## specify query (aka filter) by n1,n2,d1,d2 or xdata
##   n1,n2,d1,d2 - query is row-by-row combination
##   xdata - data.frame given desired combinations
## specify posr by explicit parameter or id or from, relto smry types
## rate.rule is keyword (eg, nonzro) or function that maps posr,rate.tol to logical vector
##   specify by keyword or function that maps posr,rate.tol to logical vector
## rate.type tells whether we want raw pos or neg rates, or error or correct rates
##   pos is no-op - uses posr as is; neg negates (really complements) everything
##   error negates  'true' cases to get false negative rate
##                  'false' cases already are false positive rate
##   correct negates 'false' cases to get true negative rate
##                   'true' cases already are true positive rate 
## truedd.multi tells what to do if rate implies mix of 'true' and 'false' cases
##   'error' means not allowed, 'asis' means ok but don't reorder,
##   'false.first','true.first' means move false (resp. true) cases first if x='auto'
## x tells which x variable drives the plot, or use x vars as given, or select automatically
data_rate=
  function(posr=NULL,posr.id=parent(posr.id,'std'),
           rate.rule=parent(rate.rule,'raw'),
           rate.type=parent(rate.type,'pos'),rate.tol=parent(rate.tol,0),
           truedd.multi=parent(truedd.multi,'false.first'),
           x=parent(x,'auto'),
           n1=parent(n1),n2=parent(n2),d1=parent(d1),d2=parent(d2),
           xdata=parent(xdata,NULL),mesr=parent(mesr)) {
    ## init(must.exist=T);            # make sure environment initialized
    check_mesr();                  # make sure measures legal
    if (is.null(posr)) posr=get_posr(posr.id);
    posr=posr_select();
    if (is.function(rate.rule)) true.dd=rate.rule(posr,rate.tol)
    else true.dd=true_dd(posr,rate.rule,rate.tol);
    drat=data.frame(posr[,cq(i,n1,n2,d1,d2)],true.dd=true.dd,posr[,mesr,drop=F]);
  }
drat_order=
  function(drat,
           rate.rule=parent(rate.rule,'raw'),
           rate.type=parent(rate.type,'pos'),rate.tol=parent(rate.tol,0),
           truedd.multi=parent(truedd.multi,'false.first'),x=parent(x,'auto'),
           mesr=parent(mesr)) {
    ## sort based on x
    if (x=='auto') {
      ## set x to first var that can drive loop, else 'asis'
      x.order=do.call(
        c,sapply(cq(n1,n2,d1,d2),simplify=F,
                 function(x) if (length(unique(drat[,x]))==nrow(drat)) x));
      if (is.null(x.order)) x.order='asis' else x.order=x.order[1];
    } else x.order=x;
    drat=switch(x.order,
                asis=drat[with(drat,order(i)),],
                n1=drat[with(drat,order(n1,n2,d1,d2)),],
                n2=drat[with(drat,order(n2,n1,d1,d2)),],
                d1=drat[with(drat,order(d1,d2,n1,n2)),],
                d2=drat[with(drat,order(d2,d1,n1,n2)),]);
    ## deal with mix of 'true' and 'false' cases
    true.dd=drat$true.dd;               # set true.dd anew since order may have changed
    if (any(true.dd)&&any(!true.dd)) {
      ## there's a mixture.
      truedd.multi=match.arg(truedd.multi,cq(false.first,true.first,asis,error));
      if (truedd.multi!='asis') {
        ## Note: if truedd.multi is 'asis', no need to worry about order
        if (truedd.multi=='error')
          stop("Cannot mix 'true' and 'false' cases (eg, false positives and false negatives) when truedd.multi is 'error'")
        else {
          ## re-order if necessary and allowed. see if already in correct order
          truedd.first=if (truedd.multi=='false.first') F else T;
          nrun=length(rle(true.dd)$lengths);
          if (nrun!=2||true.dd[1]!=truedd.first)
            ## not in correct order. re-order if allowed (x='auto')
            if (x=='auto') {
              i=order(true.dd,decreasing=truedd.first);
              drat=drat[i,];
              x.order='asis'
            } else 
              stop("Cannot re-order 'true' and 'false' cases (eg, false positives and false negatives) unless x='auto'")
          }}}
    ## convert rates based on rate.type. also start setting up output variables
    rate.type=match.arg(rate.type,cq(pos,neg,error,correct));
    true.dd=drat$true.dd;
    xdata=drat[,cq(n1,n2,d1,d2)];
    ydata=apply(drat[,mesr,drop=F],2,function(rate)
      switch(rate.type,
             pos=rate,                           # nothing to do
             neg=1-rate,                         # negate everything
             error=ifelse(true.dd,1-rate,rate),  # negate 'true' cases
             correct=ifelse(true.dd,rate,1-rate) # negate 'false' cases
      ));
    ## refine rate type if possible
    rate.type=rate_type();
    ## return everything that might be useful...
    invisible(list(xdata=xdata,ydata=ydata,
                   drat=data.frame(xdata,true.dd=true.dd,ydata),
                   true.dd=true.dd,x=x.order,rate.type=rate.type));
  }
## refine rate type if possible, eg, 'fpr' or 'fnr'
rate_type=
  function(rate.rule=parent(rate.rule),rate.type=parent(rate.type),true.dd=parent(true.dd)) {
    if (rate.type %notin% cq(pos,neg,error,correct)) stop(paste('Invalid rate type:',rate.type));
    if (rate.rule=='raw') {
      ## for raw, 'error' and 'correct' are meaningless. convert to 'neg' and 'pos'
      if (rate.type %in% cq(pos,correct)) 'pos' else 'neg';
    } else 
      switch(rate.type,
             pos=rate.type,
             neg=rate.type,
             error=if (all(true.dd)) 'fnr' else if (all(!true.dd)) 'fpr' else 'error',
             correct=if (all(true.dd)) 'tpr' else if (all(!true.dd)) 'tnr' else 'correct')
  }

## do the actual selection
posr_select=
  function(posr=parent(posr),n1=parent(n1),n2=parent(n2),d1=parent(d1),d2=parent(d2),
           xdata=parent(xdata,NULL),mesr=parent(mesr,mesr.all)) {
    ## CAUTION: need cbind in case params have incompatible lengths. sigh...
    if (is.null(xdata)) xdata=suppressWarnings(data.frame(cbind(n1,n2,d1,d2)));
    ## round d1,d2 to avoid imprecise decimals
    xdata$d1=round(xdata$d1,digits=5); xdata$d2=round(xdata$d2,digits=5); 
    ## set n1,n2,d1,d2 for sub-functions
    n1=unique(xdata$n1); n2=unique(xdata$n2); d1=unique(xdata$d1); d2=unique(xdata$d2);
    ## interpolate posr if needed. prune first to improve performance
    posr=posr_prune();
    posr=posr_interp();
    ## filter posr and sort based on x. in database parlance, the merge is natural semijoin
    xdata$i=seq_len(nrow(xdata));
    posr=merge(posr,xdata);
    ## clamp pos.rate to [0,1]. interpolation can under- or over-shoot ;
    posr[,mesr]=apply(posr[,mesr,drop=F],c(1,2),function(y) max(min(y,1),0));
    invisible(posr);
  }
## prune posr so interp will be faster
posr_prune=
  function(posr=parent(posr),n1=parent(n1),n2=parent(n2),d1=parent(d1),d2=parent(d2),
           mesr=parent(mesr)) {
    n1.lo=max(posr$n1[posr$n1<=min(n1)]);
    n1.hi=min(posr$n1[posr$n1>=max(n1)]);
    n2.lo=max(posr$n2[posr$n2<=min(n2)]);
    n2.hi=min(posr$n2[posr$n2>=max(n2)])
    d1.lo=max(posr$d1[posr$d1<=min(d1)]);
    d1.hi=min(posr$d1[posr$d1>=max(d1)]);
    d2.lo=max(posr$d2[posr$d2<=min(d2)]);
    d2.hi=min(posr$d2[posr$d2>=max(d2)]);
    posr=subset(posr,subset=(n1>=n1.lo&n1<=n1.hi&n2>=n2.lo&n2<=n2.hi&
                             d1>=d1.lo&d1<=d1.hi&d2>=d2.lo&d2<=d2.hi),
                select=c(cq(n1,n2,d1,d2),mesr));
  }
## interpolate posr one parameter at a time
posr_interp=
  function(posr=parent(posr),n1=parent(n1),n2=parent(n2),d1=parent(d1),d2=parent(d2),
           mesr=parent(mesr)) {
    posr=posr_interp1(n1);
    posr=posr_interp1(n2);
    posr=posr_interp1(d1);
    posr=posr_interp1(d2);
    ## remove rownames contructed by R and put in order we want
    rownames(posr)=NULL;
    invisible(posr);
  }
posr_interp1=function(x,xout,posr=parent(posr),mesr=parent(mesr)) {
  x=as.character(pryr::subs(x));
  if (missing(xout))
    if (exists(x,envir=parent.frame(n=1))) xout=get(x,envir=parent.frame(n=1))
    else stop(paste('no value for',x,'in function call or parent environment'));
  xout=unique(xout);
  if (all(xout%in%posr[,x])) {
    ## no need to interpolate. subset will do.
    posr=posr[posr[,x]%in%xout,];
  } else {
    ## get selectors for parameter and data columns
    param=cq(n1,n2,d1,d2);
    fix=param[param %notin% x];
    ## split posr by fixed params
    posr.fix=posr[,fix];
    posr.byfix=split(posr,apply(posr.fix,1,function(row) paste(collapse=' ',row)));
    ## interpolate each group over x
    posr=do.call(rbind,lapply(posr.byfix,function(posr) {
      yout=asplinem(posr[,x],posr[,mesr],xout=xout);
      posr=suppressWarnings(data.frame(xout,posr[1,fix],yout)); # tack params onto output
      colnames(posr)=c(x,fix,mesr);                             # and fix names
      posr;
    }));
  }
  posr;
}
## compute boolean vector of which rows are supposed to say 'yes' and which 'no'
## for predefined correctness criteria
## note various flavors of significance criteria - not sure which is best
##   nonzro: default non-zero rate. currently nonz1
##   nonz1: d1 non-zero; d2 doesn't matter
##   nonz1or2: d1 or d2 non-zero
##   nonz1and2: d1 and d2 non-zero
##   nonz2: d2 non-zero; d1 doesn't matter
##   sameff: d1, d2 equal
##   uni1: d1, d2 equal and non-zero
true_dd=function(posr,rate.rule,rate.tol=0) {
  with(posr,switch(rate.rule,
                   raw=rep(T,nrow(posr)),
                   ## criteria parameterized by rate.tol
                   ## nonzro=(d1!=0),
                   nonzro=(d1>rate.tol),
                   nonz1=(d1>rate.tol),
                   nonz1or2=(d1>rate.tol|d2>rate.tol),
                   nonz1and2=(d1>rate.tol&d2>rate.tol),
                   nonz2=(d2>rate.tol),
                   farzro=(d1>=rate.tol),                   # for backwards compatibility
                   sameff=(abs(d1-d2)<=rate.tol),
                   nearff=(abs(d1-d2)<=rate.tol),           # for backwards compatibility
                   uni=(d1>=rate.tol&abs(d1-d2)<=rate.tol), # not used
                   stop(paste('Invalid rate rule:',rate.rule))
                   ));
}
