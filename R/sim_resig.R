#################################################################################
##
## Author:  Nat Goodman
## Created: 18-09-08
##          from sim.R restarted created 18-05-03
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
## Generate simulated data for repwr.R.
## Includes simulation itself and pre- and post- simulation  pipeline
## Specializes sim.R for the very simple resig analysis which only looks at sig2
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################

## --- Generate Data ---
dodata=function(need.init=T,doc=parent(doc,'readme'),...) {
  if (need.init)
    ## for sandbox runs, use doc-specific init
    if (doc=='xperiment') init_xperiment(doc=doc,...) else init(doc=doc,...);
  ## dopre();              # precalculate or load global data. not need for resig
  dosim();                 # load saved simulations or do new ones
  dopost();                # post-simulation pipeline
}

## ## ---- Precalculate or Load Global Data ----
## dopre=function() {
##   doscope();                     # load or generate small telescopes data
##   doscopd();                     # load or generate precalculated 'scopd' data 
##   doconfvl();                    # load or generate precalculated confidence interval data
##   dopredvl();                    # load or generate precalculated prediction interval data
## }
## ---- Simulation Functions ----
## do the simulation.  do it case by case for stylistic consistency with dores
dosim=function() {
  ## do simulations! each value of n x d
  cases=expand.grid(n=n,d=d);
  apply(cases,1,function(case) {
    n=case['n']; d=case['d'];
    sim=dosim1(n,d);                    # do base simulation
    simr=dosimr(sim,n,d);               # add/substract columns to/from sim for results
  })
  invisible();
}
## simulate one case. adapted from swfdr
dosim1=function(n,d) {
  ## use saved sim if exists and args permit
  sim=get_sim(n,d,id,must.exist=F);
  if (!is.null(sim)) return(invisible(sim));
  ## no saved simulation or args say not to use it. run simulation
  if (verbose) print(paste(sep=' ','>>> dosim1:',nvq(n,d)));
  sim=data.frame(do.call(rbind,lapply(seq_len(m),function(i) {
    ## draw pair of samples. group0 is control. group1 has effect=d
    group0=rnorm(n,mean=0);
    group1=rnorm(n,mean=d);
    mean0=mean(group0);
    mean1=mean(group1);
    d.raw=mean1-mean0;
    sd0=sd(group0);
    sd1=sd(group1);
    ## pval=t.test(group0,group1,var.equal=T)$p.value;
    ## store results
    ## c(pval=pval,d.raw=d.raw,mean0=mean0,mean1=mean1,sd0=sd0,sd1=sd1);
    c(d.raw=d.raw,mean0=mean0,mean1=mean1,sd0=sd0,sd1=sd1);
  })));
  sd=with(sim,pooled_sd(sd0,sd1));
  d.sdz=sim$d.raw/sd;
  pval=d2pval(n,d.sdz);
  ## NG 18-02-07. no longer need symbolic row ids. integer ids work fine when needed
  ## sim=with(sim,data.frame(id=paste(sep=',',case,seq_len(nrow(sim))),n=n,
  ##                         d=d,d.sdz=d.sdz,sd,pval,sim,row.names=NULL,stringsAsFactors=F));
  sim=with(sim,data.frame(n=n,d=d,d.sdz=d.sdz,sd,pval,sim,row.names=NULL,stringsAsFactors=F));
  ## optionally save sim and add to in-memory list. usually don't keep - no point
  save_sim(sim,n,d,id);
  invisible(sim);
}
## add/substract columns to/from sim as needed for analysis
dosimr=function(sim=NULL,n,d) {
  ## use saved simr if exists and args permit
  simr=get_simr(n,d,id,must.exist=F);
  if (!is.null(simr)) return(invisible(simr));
  ## no saved simr or args say not to use it. construct it
  if (verbose) print(paste(sep=' ','+++ dosimr:',nvq(n,d)));
  if (is.null(sim)) sim=get_sim(n,d,id,load=T);
  simr=sim;
  ## add columns for analysis
  ## simr.ci=ci(sim,confvl[confvl$n==n,]);
  ## simr[,names(simr.ci)]=simr.ci;
  simr$d.sign=sign(simr$d.sdz);
  simr$d.abs=abs(simr$d.sdz);
  ## ## for meta-analysis
  ## simr$w.meta=with(simr,1/sd_d2t(n,d.sdz)^2);
  ## simr$d.meta=with(simr,d.sdz*w.meta)
  ## drop ones we no longer need
  simr=subset(simr,select=c(-d.raw,-mean0,-mean1,-sd0,-sd1));
  ## optionally save simr and add to in-memory list. usually do - not too big
  save_simr(simr,n,d,id);
  invisible(simr);
}
## ---- Post-simulation Pipeline ----
## dores. compute the results. do it case by case to avoid rereading large detail files
dopost=function() {
  ## construct cases
  if (dopost.allcases) cases=expand.grid(n1=n,n2=n,d1=d,d2=d)
  else {
    ## NG 18-02-22: this branch is obsolete
    ##   I originally thought it would be too slow to do 'em all, but it's okay
    stop('Trying to run obsolete branch of dopost');
  }
  ## do it!
  smry=do.call(rbind,apply(cases,1,function(case) {
    n1=case['n1']; n2=case['n2']; d1=case['d1']; d2=case['d2'];
    if (is.na(d2)) d2=d1;
    ## compute detailed results or use saved detl
    detl=dodetl(n1=n1,n2=n2,d1=d1,d2=d2);
    ## compute summary or used saved smry
    smry=dosmry(detl=detl,n1=n1,n2=n2,d1=d1,d2=d2)
  }));
  ## optionally get mesrs from smry and save. usually do - needed for analysis
  mesr=setdiff(colnames(smry),cq(n1,n2,d1,d2,type,m));
  save_data(mesr);
  ## optionally save smry and add to in-memory list. usually do - needed for analysis
  save_data(smry);
  ## generate default posr (positive rate) cases. also saves results if desired
  posr=doposr(smry=smry);
  invisible(posr);
}
## contruct one detl case
dodetl=function(s1=NULL,s2=NULL,n1,n2,d1,d2) {
  ## use saved detl if exists and args permit
  detl=get_detl(n1,n2,d1,d2,id,must.exist=F);
  if (!is.null(detl)) return(invisible(detl));
  ## no saved detl or args say not to use it. construct detl
  if (verbose) print(paste(sep=' ','+++ dodetl',nvq(n1,n2,d1,d2)));
  if (is.null(s1)) s1=get_simr(n1,d1,id,load=T);
  if (is.null(s2)) s2=get_simr(n2,d2,id,load=T);
  ## NG 18-02-23: permute s1 else self-comparisons meaningless
  ## NG 18-06-26: permute s2 also. not needed, but makes narrative easier
  ## use saved si if exists and args permit
  si=get_si(n1,n2,d1,d2,id,must.exist=F);
  if (is.null(si)) {
    m=nrow(s1);
    if (dopost.permute) {
      i1=sample.int(m);
      i2=sample.int(m);
    } else i1=i2=seq_len(m);
    si=data.frame(i1=i1,i2=i2);
    ## optionally add i1 to in-memory list and save to file (mostly for testing)
    save_si(si,n1,n2,d1,d2,id);
  }
  s1=s1[si$i1,];
  s2=s2[si$i2,];
  ## ## define variables that depend on both studies
  ## ## prediction intervals
  ## p1=pi(s1,predvl[predvl$n1==n1&predvl$n2==n2,]);
  ## p2=pi(s2,predvl[predvl$n1==n2&predvl$n2==n1,]);
  ## ## small telescope thresholds
  ## scp1=scope[as.character(n1),as.character(n2)];
  ## scp2=scope[as.character(n2),as.character(n1)];
  ## ## my 'misinterpretation' of small telescopes that considers d.sdz
  ## scpd1=scpd(s1,scopd[scopd$n1==n1&scopd$n2==n2,]);
  ## scpd2=scpd(s2,scopd[scopd$n1==n2&scopd$n2==n1,]);
  ## ## meta d, ci, pval
  ## w.sum=s1$w.meta+s2$w.meta;
  ## meta.d=(s1$d.meta+s2$d.meta)/w.sum;
  ## meta.abs=abs(meta.d);
  ## meta.sd=sqrt(1/(w.sum));
  ## meta.pval=2*pnorm(-abs(meta.d/meta.sd));
  ## meta.ci=qnorm(p=0.5+conf.level/2,sd=meta.sd)
  ## meta.ci.lo=meta.d-meta.ci;
  ## meta.ci.hi=meta.d+meta.ci;
  ## apply the rules. do 'em in-line for efficiency
  ## use concise terms so detl won't be too wide to view
  ##   d{12m}=effect size for s1, s2, meta
  ##   c{12m}=confidence interval for s1, s2, meta
  ##   scp{12}=small telescope threshold for s1, s2
  ##   term.term means compare the two terms, eg, d1.c2 means s1 d.sdz in s2 confidence interval
  ## significant: s1,s2,meta 
  sig1=s1$pval<=sig.level;
  sig2=s2$pval<=sig.level;
  ## sigm=meta.pval<=sig.level;
  ## same direction (ie, sign): s1,s2
  sdir=s1$d.sign==s2$d.sign;
  ## ## d in confidence interval: each vs other two
  ## d1.c2=between(s1$d.sdz,s2$ci.lo,s2$ci.hi);
  ## d2.c1=between(s2$d.sdz,s1$ci.lo,s1$ci.hi);
  ## d1.cm=between(s1$d.sdz,meta.ci.lo,meta.ci.hi);
  ## d2.cm=between(s2$d.sdz,meta.ci.lo,meta.ci.hi);
  ## dm.c1=between(meta.d,s1$ci.lo,s1$ci.hi);
  ## dm.c2=between(meta.d,s2$ci.lo,s2$ci.hi);
  ## ## confidence intervals overlap: all 3 unique pairs
  ## ##   thanks to https://fgiesen.wordpress.com/2011/10/16/checking-for-interval-overlap/
  ## ##   for simple calculation I should have known already :)
  ## c1.c2=s1$ci.lo<=s2$ci.hi & s2$ci.lo<=s1$ci.hi;
  ## c1.cm=s1$ci.lo<=meta.ci.hi & meta.ci.lo<=s1$ci.hi;
  ## c2.cm=s2$ci.lo<=meta.ci.hi & meta.ci.lo<=s2$ci.hi;
  ## ## d in prediction interval: each d vs two intervals we have
  ## d1.p2=between(s1$d.sdz,p2$pi.lo,p2$pi.hi);
  ## d2.p1=between(s2$d.sdz,p1$pi.lo,p1$pi.hi);
  ## dm.p2=between(meta.d,p2$pi.lo,p2$pi.hi);
  ## dm.p1=between(meta.d,p1$pi.lo,p1$pi.hi);
  ## ## prediction intervals overlap
  ## p1.p2=p1$pi.lo<=p2$pi.hi & p2$pi.lo<=p1$pi.hi;
  ## ## d >= small telescope boundary
  ## ## NG 18-02-14: to handle negative values, use d.abs instead of d.sdz
  ## d1.scp2=(s1$d.abs>=scp2)&sdir;
  ## d2.scp1=(s2$d.abs>=scp1)&sdir;
  ## dm.scp2=(meta.abs>=scp2)&sdir;
  ## dm.scp1=(meta.abs>=scp1)&sdir;
  ## ## d >= my 'misinterpretation' of small telescope boundary
  ## ##   need as.vector else R treats it as matrix. screws up naming in data.frame
  ## d1.scpd2=as.vector((s1$d.abs>=scpd2)&sdir);
  ## d2.scpd1=as.vector((s2$d.abs>=scpd1)&sdir);
  ## dm.scpd2=as.vector((meta.abs>=scpd2)&sdir);
  ## dm.scpd1=as.vector((meta.abs>=scpd1)&sdir);
  ## ## d2 bigger (actually more extreme) than d1
  ## big1=(s1$d.abs>=s2$d.abs)&sdir;
  ## big2=(s1$d.abs<=s2$d.abs)&sdir;
  ## return as data.frame so dosmry can use as environment
  ## detl=data.frame(sig1,sig2,sigm,sdir,
  ##   d1.c2,d2.c1,d1.cm,d2.cm,dm.c1,dm.c2,c1.c2,c1.cm,c2.cm,d1.p2,d2.p1,dm.p1,dm.p2,p1.p2,
  ##   d1.scp2,d2.scp1,dm.scp2,dm.scp1,d1.scpd2,d2.scpd1,dm.scpd2,dm.scpd1,
  ##   big1,big2);
  ## return as data.frame so dosmry can use as environment
  detl=data.frame(sig1,sig2,sdir);
  ## optionally save detl and add to in-memory list. usually don't keep -- too big
  save_detl(detl,n1,n2,d1,d2,id);
  invisible(detl);
}
## contruct one smry case
dosmry=function(detl=NULL,n1,n2,d1,d2) {
  ## use saved smry if exists and args permit
  smry=get_smry(n1,n2,d1,d2,id,must.exist=F);
  if (!is.null(smry)) return(invisible(smry));
  ## no saved smry or args say not to use it. construct smry
  if (verbose) print(paste(sep=' ','+++ dosmry',nvq(n1,n2,d1,d2)));
  if (is.null(detl)) detl=get_detl(n1,n2,d1,d2,id,load=T);
  ## construct summaries for several cases
  type=cq(raw,sig1,bsln,sig2,big2,s2or); # order must be same as rows added to smry!
  smry.na=rep(NA,ncol(detl));            # used when smry type is empty
  ## 1) raw detail
  smry=data.frame(type='raw',m=nrow(detl),t(colSums(detl)),stringsAsFactors=F);
  raw=c(nrow(detl),colSums(detl));
  ## 2) sig1.
  detl.sig1=with(detl,detl[sig1,,drop=F]);
  sig1=c(nrow(detl.sig1),if(nrow(detl.sig1)!=0) colSums(detl.sig1) else smry.na);
  ## 3) baseline = sig1 & sdir. all require this
  detl.bsln=with(detl.sig1,detl.sig1[sdir,,drop=F]);
  bsln=c(nrow(detl.bsln),if(nrow(detl.bsln)!=0) colSums(detl.bsln) else smry.na);
  ## ## 4) baseline & sig2. small telescopes requires this
  ## detl.sig2=with(detl.bsln,detl.bsln[sig2,,drop=F]);
  ## sig2=c(nrow(detl.sig2),if(nrow(detl.sig2)!=0) colSums(detl.sig2) else smry.na);
  ## ## 5) bsln & (sig2 | each other mesr)
  ## ## CAUTION: if detl.bsln is empty, OR'ing the terms is error
  ## detl.s2or=if(nrow(detl.bsln)!=0) detl.bsln|detl.bsln$sig2 else detl.bsln;
  ## s2or=c(nrow(detl.s2or),if(nrow(detl.s2or)!=0) colSums(detl.s2or) else smry.na);
  ## ## 6) baseline | big2. many authors also accept big2
  ## ## CAUTION: if detl.bsln is empty, OR'ing the terms is error
  ## detl.big2=if(nrow(detl.bsln)!=0) detl.bsln|detl.bsln$big2 else detl.bsln;
  ## big2=c(nrow(detl.big2),if(nrow(detl.big2)!=0) colSums(detl.big2) else smry.na);
  ## turn it all into data frame with n1,n2,d1,d2,type at front
  ## smry=data.frame(n1,n2,d1,d2,type,rbind(t(raw),t(sig1),t(bsln),t(sig2),t(big2),t(s2or)),
  ##                 row.names=NULL,stringsAsFactors=F);
  smry=data.frame(n1,n2,d1,d2,type,rbind(t(raw),t(sig1),t(bsln)),
                  row.names=NULL,stringsAsFactors=F);
  colnames(smry)=c(cq(n1,n2,d1,d2,type,m),colnames(detl));
  ## optionally save smry and add to in-memory list. usually don't keep -- no point
  save_smry(smry,n1,n2,d1,d2,id);
  invisible(smry);
}
## contruct default positive rate objects
doposr=function(smry=NULL) {
  if (is.null(smry)) smry=get_data(smry);
  ## init measures and summary types unless already done
  init_mesr();
  ## construct positive rates for cases of interest
  ## 1) standard
  ## 2) from sig1 relto sig1 - to show effect of sdir in supplement
  ## 3) bias1 - to show that sig1 bias in fnr_exact is not due to sdir
  ## ## 3) from sig2 relto sig1
  ## ## 4) from sig2 relto sig2
  do_posr(smry,from.type=mesr.fromtype,relto.type=mesr.reltotype,posr.id='std')
  do_posr(smry,from.type='sig1',relto.type='sig1');
  do_posr(smry,from.type=setNames(cq(raw,sig1,raw),cq(sig1,sig2,sdir)),
          relto.type=mesr.reltotype,posr.id='bias1');
  ## do_posr(smry,from.type='bsln',relto.type='sig1');
  ## do_posr(smry,from.type='sig2',relto.type='sig1');
  ## do_posr(smry,from.type='sig2',relto.type='sig2');
  ## do_posr(smry,from.type='s2or',relto.type='sig1');
  invisible(T);
}

## get or contruct one positive rate objects and optionally save & keep
## from, relto are smry types. should be single valued unless id set
## posr.id used to create file and saved object. if missing, set from from.type, relto.type 
do_posr=
  function(smry,from.type=parent(from.type,'bsln'),relto.type=parent(relto.type,'sig1'),
           posr.id=NULL,mesr=mesr.all) {
    if (is.null(posr.id)) posr.id=paste(sep='_',from.type[1],relto.type[1]);
    ## use saved posr if exists and args permit
    posr=get_posr(posr.id,must.exist=F);
    if (!is.null(posr)) return(invisible(posr));
    ## no saved posr or args say not to use it. construct posr
    if (verbose) print(paste(sep=' ','+++ doposr',casename_posr(posr.id)));
    ## make sure measures & types legal and limit to type we need
    from.type=check_type(from.type,multiok=T);
    relto.type=check_type(relto.type,multiok=T);
    check_mesr();
    smry.bytype=split(smry,smry$type);
    posr=do.call(cbind,lapply(mesr,function(mesr) {
      from.type=from.type[mesr];
      relto.type=relto.type[mesr];
      count=smry.bytype[[from.type]][,mesr];
      m=smry.bytype[[relto.type]]$m;
      ifelse(m==0,0,count/m);
    }));
    if (any(posr[!is.na(posr)]>1))
      stop(paste(sep='',"Some posr values are greater than 1. This usually means from.type is 'before' relto.type in the workflow: from.type=",from.type,", relto.type=",relto.type));
   if (any(posr[!is.na(posr)]<0))
      stop(paste(sep='',"Some posr values are negative. I have no idea how this could happen..."));
    colnames(posr)=mesr;
    ## tack on x-params
    xdata=smry.bytype[[1]][,cq(n1,n2,d1,d2)];
    posr=data.frame(xdata,posr);
    ## optionally save posr and add to in-memory list. usually do
    save_posr(posr,posr.id);
    invisible(posr);
  }

## ## ---- Small Telescopes Functions ----
## ## generate & optionally save or load data for small telescopes
## ## presently just scope - matrix of "close enough" small telescopes d values for n1, n2
## doscope=function() {
##   ## use saved scope if exists and args permit
##   scope=get_data(scope,must.exist=F);
##   if (!is.null(scope)) return(invisible(scope));
##   ## no saved scope or args say not to use it. calculate it
##   if (verbose) print(paste(sep='','>>> doscope'));
##   scope=do_scope();
##   ## optionally save scope and add to in-memory list. usually do - needed for dodetl
##   save_data(scope);
##   invisible(scope);  
## }
## ## get matrix of "close enough" small telescopes d values for n1, n2
## do_scope=function(n1=parent(n,20),n2=parent(n,40)) {
##   cases=expand.grid(n1=n1,n2=n2);
##   scope=apply(cases,1,function(case) {
##     n1=case['n1']; n2=case['n2'];
##     d1=d_sig(n1,sig.level);
##     d.scope=uniroot(
##       function(d0) p_d2t(n1,d=d1,d0,lower.tail=F)-scope.power,interval=c(0,10))$root;
##     ## suppress annoying warning from pt:
##     ##   full precision may not have been achieved in 'pnt{final}'
##     d.close=suppressWarnings(
##       uniroot(function(d) p_d2t(n2,d,d0=d.scope)-scope.close,c(-10,10))$root);
##     d.close;})
##   scope=matrix(data=scope,nrow=length(n),ncol=length(n))
##   rownames(scope)=colnames(scope)=n;
##   scope;
## }
## ## generate & optionally save or load precalculated scopd data
## ## for my 'misinterpretation' of small telescopes that considers d.sdz
## doscopd=function() {
##   ## use saved scopd if exists and args permit
##   scopd=get_data(scopd,must.exist=F);
##   if (!is.null(scopd)) return(invisible(scopd));
##   ## no saved scopd or args say not to use it. calculate it
##   if (verbose) print(paste(sep='','>>> doscopd'));
##   scopd=do_scopd();
##   ## optionally save scopd and add to in-memory list. usually do - needed for dodetl
##   save_data(scopd);
##   invisible(scopd);  
## }
## ## generate precalculated scopd data
## do_scopd=function(n=parent(n),d=parent(dsdz.grid)) {
##   d=d[d>=0];                            # scope works on absolute values
##   cases=expand.grid(n1=n,d=d);
##   ## scpd1 is min d.pop with acceptable power to detect d1.sdz - depends on s1 only
##   ## scpd2 is min d2.sdz close enough to scpd1 - depends on s1, s2
##   scpd1=do_scopd1(n1=cases$n1,d=cases$d);
##   scpd1=data.frame(cases,scpd1=scpd1);
##   colnames(cases)=cq(n2,d);             # rename 'n' column for stylistic consistency
##   cases=merge(scpd1,expand.grid(n2=n,d=d),by='d',suffixes=c(1,2));
##   scpd2=do_scopd2(n2=cases$n2,scpd1=cases$scpd1);
##   scpd2=data.frame(cases,scpd=scpd2);
##   scpd2[,cq(n1,n2,d,scpd)];           # select & reorder columns for stylistic consistency
## }
## ## scopd1 is min d.pop with acceptable power to detect d1.sdz - depends on s1 only
## ## scopd2 is min d2.sdz close enough to scopd1 - depends on s1, s2
## do_scopd1=
##   Vectorize(function(n1,d)
##     suppressWarnings(
##       uniroot(function(d0) p_d2t(n1,d,d0,lower.tail=F)-scope.power,interval=c(-10,10))$root));
## do_scopd2=
##   Vectorize(function(n2,scpd1)
##     suppressWarnings(
##       uniroot(function(d) p_d2t(n2,d,d0=scpd1)-scope.close,c(-10,10))$root));

## ## interpolate precalculated scopd data at d.sdz
## ## simr & scopd should already be filtered to current n1, n2
## scpd=function(simr,scopd) {
##   scopd=with(scopd,akima::aspline(x=d,y=scpd,xout=simr$d.sdz)$y)
##   data.frame(scpd=scopd);
## }
## ## ---- Confidence and Prediction Interval Functions ----
## ## generate & optionally save or load precalculated confidence intervals
## doconfvl=function() {
##   ## use saved confvl if exists and args permit
##   confvl=get_data(confvl,must.exist=F);
##   if (!is.null(confvl)) return(invisible(confvl));
##   ## no saved confvl or args say not to use it. calculate it
##   if (verbose) print(paste(sep='','>>> doconfvl'));
##   confvl=do_confvl();
##   ## optionally save confvl and add to in-memory list. usually do - needed for dodetl
##   save_data(confvl);
##   invisible(confvl);  
## }
## ## generate precalculated confidence intervals
## do_confvl=function(n=parent(n),d=parent(dsdz.grid)) {
##   cases=expand.grid(n=n,d=d);
##   confvl=t(Vectorize(ci_d2t)(n=cases$n,d=cases$d));
##   confvl=data.frame(cases,confvl);
##   colnames(confvl)=cq(n,d,lo,hi);
##   confvl;
## }
## ## interpolate precalculated confidence interval data at d.sdz
## ## sim & confvl should already be filtered to current n
## ci=function(sim,confvl) {
##   lo=with(confvl,akima::aspline(x=d,y=lo,xout=sim$d.sdz)$y);
##   hi=with(confvl,akima::aspline(x=d,y=hi,xout=sim$d.sdz)$y);
##   data.frame(ci.lo=lo,ci.hi=hi)
## }
## ## generate & optionally save or load precalculated prediction intervals
## dopredvl=function() {
##   ## use saved predvl if exists and args permit
##   predvl=get_data(predvl,must.exist=F);
##   if (!is.null(predvl)) return(invisible(predvl));
##   ## no saved predvl or args say not to use it. calculate it
##   if (verbose) print(paste(sep='','>>> dopredvl'));
##   predvl=do_predvl();
##   ## optionally save predvl and add to in-memory list. usually do - needed for dodetl
##   save_data(predvl);
##   invisible(predvl);  
## }
## ## generate precalculated prediction intervals
## ## assumes doconfvl already done
## do_predvl=function(confvl=parent(confvl)) {
##   cases=merge(confvl,confvl,by='d',suffixes=c(1,2));
##   ## since we already have conf intervals, no need to do full pi_d2t rigamarole
##   lo=with(cases,d-sqrt((d-lo1)^2+(hi2-d)^2));
##   hi=with(cases,d+sqrt((d-lo2)^2+(hi1-d)^2));
##   data.frame(n1=cases$n1,n2=cases$n2,d=cases$d,lo=lo,hi=hi);
## }
## ## interpolate precalculated prediction interval data at d.sdz
## ## simr & predvl should already be filtered to current n1, n2
## pi=function(simr,predvl) {
##   lo=with(predvl,akima::aspline(d,y=lo,xout=simr$d.sdz)$y);
##   hi=with(predvl,akima::aspline(d,y=hi,xout=simr$d.sdz)$y);
##   data.frame(pi.lo=lo,pi.hi=hi)
## }
## convert pos.rate matrix to correct and error rate matrices
##   yes tells which rows are supposed to say 'yes'
correct_rate=function(pos.rate,yes) {
  ydata=pos.rate$y;
  correct.rate=apply(ydata,2,function(y) ifelse(yes,y,1-y));
  list(x=pos.rate$x,y=correct.rate,yes=yes);
}
error_rate=function(pos.rate,yes) {
  ydata=pos.rate$y;
  error.rate=apply(ydata,2,function(y) ifelse(yes,1-y,y));
  list(x=pos.rate$x,y=error.rate,yes=yes);
}
## negate pos rate matrix
neg_rate=function(pos.rate) {
  ydata=pos.rate$y;
  neg.rate=1-ydata;
  list(x=pos.rate$x,y=neg.rate);
}

