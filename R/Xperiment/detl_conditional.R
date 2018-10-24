#################################################################################
##
## Author:  Nat Goodman
## Created: 18-08-08
##          from sim.R created 18-05-03
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
## Experimental sandbox code related to replication power blog post
## Tests having dodetl calculate just the measures used in a given doc
## Uses readme data for everything up to detl
##   DATA must be HAND-COPIED BEFORE RUNNING SCRIPT
## This version does readme mesrs under conditional
##   performance almost as good as hand-crafted and reasonably general
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## source 'real code'. repwr.R get's 'em all
source('R/repwr.R');
## --- Eperimental sandbox functions ---
## run for this sandbox
run=function(need.init=T,doc='xperiment',...) {
  if (need.init) wrap_fun(init_xperiment);
  need.init=F;
  wrap_fun(dodata);                   # generate data - ie, run simulation
  wrap_fun(dodoc,init_doc_xperiment); # generate figures, tables for doc
  cmp_detl();                                     # run test
}
## init for this sandbox. same parameters as readme
init_xperiment=
  function(doc='xperiment',subdoc='detl_conditional',
           n=20*2^(0:4),d=c(0,0.2,0.5,0.8,1),m=1e3,
           mdir=paste_nv(m,m_pretty(m)),            # m subdirectory
           datadir=filename('data','xperiment',subdoc,mdir),
           clean=F,clean.memlist=T,clean.sim=F,clean.simr=F,clean.si=F,clean.toplevel=F,
           clean.detl=T,clean.smry=T,clean.posr=T,
           ...) {
    ## call init with our arguments
    wrap_fun(init);
  }
## init_doc for this sandbox
init_doc_xperiment=
  function(subdoc='detl_conditional',docfun=doc_readme,
           figdir=filename('figure',doc,subdoc,mdir),tbldir=filename('table',doc,subdoc,mdir),
           clean.out=T,figscreen=T,...) {
    wrap_fun(init_doc);
 }
## contruct one detl case
##### do readme mesrs under conditional
##    doc needs:  sig2,d1.c2,sigm,d2.c1,c1.c2,d1.p2,d2.p1,p1.p2,d2.scp1
##    code needs: sig1,sig2,sigm,sdir,d1.c2,d2.c1,c1.c2,d1.p2,d2.p1,p1.p2,big2
dodetl=function(s1=NULL,s2=NULL,n1,n2,d1,d2) {
  mesr.need=cq(sig1,sig2,sigm,sdir,d1.c2,d2.c1,c1.c2,d1.p2,d2.p1,p1.p2,d2.scp1,big2);
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
  ## define variables that depend on both studies
  ## prediction intervals
  if (any(grepl('\\.p(1|2)$',mesr.need))) {
    p1=pi(s1,predvl[predvl$n1==n1&predvl$n2==n2,]);
    p2=pi(s2,predvl[predvl$n1==n2&predvl$n2==n1,]);
  }
  ## small telescope thresholds
  if (any(grepl('d(1|2)\\.scp(d{0,1})(1|2)',mesr.need))) {
    scp1=scope[as.character(n1),as.character(n2)];
    scp2=scope[as.character(n2),as.character(n1)];
    ## my 'misinterpretation' of small telescopes that considers d.sdz
    scpd1=scpd(s1,scopd[scopd$n1==n1&scopd$n2==n2,]);
    scpd2=scpd(s2,scopd[scopd$n1==n2&scopd$n2==n1,]);
  }
  ## meta d, ci, pval
  if (any(grepl('sigm|dm|cm',mesr.need))) {
    w.sum=s1$w.meta+s2$w.meta;
    meta.d=(s1$d.meta+s2$d.meta)/w.sum;
    meta.abs=abs(meta.d);
    meta.sd=sqrt(1/(w.sum));
    meta.pval=2*pnorm(-abs(meta.d/meta.sd));
    if (any(grepl('cm',mesr.need))) {
      meta.ci=qnorm(p=0.5+conf.level/2,sd=meta.sd)
      meta.ci.lo=meta.d-meta.ci;
      meta.ci.hi=meta.d+meta.ci;
    }}
  ## apply the rules. do 'em in-line for efficiency
  ## use concise terms so detl won't be too wide to view
  ##   d{12m}=effect size for s1, s2, meta
  ##   c{12m}=confidence interval for s1, s2, meta
  ##   scp{12}=small telescope threshold for s1, s2
  ##   term.term means compare the two terms, eg, d1.c2 means s1 d.sdz in s2 confidence interval
  ## significant: s1,s2,meta
  detl=list();
  detl$sig1=s1$pval<=sig.level;
  detl$sig2=s2$pval<=sig.level;
  detl$sdir=s1$d.sign==s2$d.sign;
  if ('sigm' %in% mesr.need) detl$sigm=meta.pval<=sig.level;
  ## same direction (ie, sign): s1,s2
  ## d in confidence interval: each vs other two
  if ('d1.c2' %in% mesr.need) detl$d1.c2=between(s1$d.sdz,s2$ci.lo,s2$ci.hi);
  if ('d2.c1' %in% mesr.need) detl$d2.c1=between(s2$d.sdz,s1$ci.lo,s1$ci.hi);
  if ('d1.cm' %in% mesr.need) detl$d1.cm=between(s1$d.sdz,meta.ci.lo,meta.ci.hi);
  if ('d2.cm' %in% mesr.need) detl$d2.cm=between(s2$d.sdz,meta.ci.lo,meta.ci.hi);
  if ('dm.c1' %in% mesr.need) detl$dm.c1=between(meta.d,s1$ci.lo,s1$ci.hi);
  if ('dm.c2' %in% mesr.need) detl$dm.c2=between(meta.d,s2$ci.lo,s2$ci.hi);
  ## confidence intervals overlap: all 3 unique pairs
  ##   thanks to https://fgiesen.wordpress.com/2011/10/16/checking-for-interval-overlap/
  ##   for simple calculation I should have known already :)
  if ('c1.c2' %in% mesr.need) detl$c1.c2=s1$ci.lo<=s2$ci.hi & s2$ci.lo<=s1$ci.hi;
  if ('c1.cm' %in% mesr.need) detl$c1.cm=s1$ci.lo<=meta.ci.hi & meta.ci.lo<=s1$ci.hi;
  if ('c2.cm' %in% mesr.need) detl$c2.cm=s2$ci.lo<=meta.ci.hi & meta.ci.lo<=s2$ci.hi;
  ## d in prediction interval: each d vs two intervals we have
  if ('d1.p2' %in% mesr.need) detl$d1.p2=between(s1$d.sdz,p2$pi.lo,p2$pi.hi);
  if ('d2.p1' %in% mesr.need) detl$d2.p1=between(s2$d.sdz,p1$pi.lo,p1$pi.hi);
  if ('dm.p2' %in% mesr.need) detl$dm.p2=between(meta.d,p2$pi.lo,p2$pi.hi);
  if ('dm.p1' %in% mesr.need) detl$dm.p1=between(meta.d,p1$pi.lo,p1$pi.hi);
  ## prediction intervals overlap
  if ('p1.p2' %in% mesr.need) detl$p1.p2=p1$pi.lo<=p2$pi.hi & p2$pi.lo<=p1$pi.hi;
  ## d >= small telescope boundary
  ## NG 18-02-14: to handle negative values, use d.abs instead of d.sdz
  if (any(grepl('d(1|2)\\.scp(d{0,1})(1|2)',mesr.need))) {
    if ('d1.scp2' %in% mesr.need) detl$d1.scp2=(s1$d.abs>=scp2)&detl$sdir;
    if ('d2.scp1' %in% mesr.need) detl$d2.scp1=(s2$d.abs>=scp1)&detl$sdir;
    if ('dm.scp2' %in% mesr.need) detl$dm.scp2=(meta.abs>=scp2)&detl$sdir;
    if ('dm.scp1' %in% mesr.need) detl$dm.scp1=(meta.abs>=scp1)&detl$sdir;
    ## d >= my 'misinterpretation' of small telescope boundary
    ##   need as.vector else R treats it as matrix. screws up naming in data.frame
    if ('d1.scpd2' %in% mesr.need) detl$d1.scpd2=as.vector((s1$d.abs>=scpd2)&detl$sdir);
    if ('d2.scpd1' %in% mesr.need) detl$d2.scpd1=as.vector((s2$d.abs>=scpd1)&detl$sdir);
    if ('dm.scpd2' %in% mesr.need) detl$dm.scpd2=as.vector((meta.abs>=scpd2)&detl$sdir);
    if ('dm.scpd1' %in% mesr.need) detl$dm.scpd1=as.vector((meta.abs>=scpd1)&detl$sdir);
  }
  ## d2 bigger (actually more extreme) than d1
  if ('big1' %in% mesr.need) detl$big1=(s1$d.abs>=s2$d.abs)&detl$sdir;
  if ('big2' %in% mesr.need) detl$big2=(s1$d.abs<=s2$d.abs)&detl$sdir;
  ## return as data.frame so dosmry can use as environment
  ## TODO: handle missing mesrs
  ## print('>>> before data.frame'); BREAKPOINT();
  ## detl=as.data.frame(do.call(cbind,sapply(mesr.need,get,simplify=F)))
  ## detl=data.frame(sig1,sig2,sigm,sdir,
  ##   d1.c2,d2.c1,d1.cm,d2.cm,dm.c1,dm.c2,c1.c2,c1.cm,c2.cm,d1.p2,d2.p1,dm.p1,dm.p2,p1.p2,
  ##   d1.scp2,d2.scp1,dm.scp2,dm.scp1,d1.scpd2,d2.scpd1,dm.scpd2,dm.scpd1,
  ##   big1,big2);
  ## optionally save detl and add to in-memory list. usually don't keep -- too big
  detl=as.data.frame(detl);
  save_detl(detl,n1,n2,d1,d2,id);
  invisible(detl);
}
## compare sandbox and real detls
## cmp_detl same for detl_conditional and detl_handcrafted
cmp_detl=function(verbose=F) {
  mesr.need=cq(sig1,sig2,sigm,sdir,d1.c2,d2.c1,c1.c2,d1.p2,d2.p1,p1.p2,d2.scp1,big2);
  mdir=paste_nv(m,m_pretty(m));
  xdir=detldir;
  rdir=file.path('data/readme',mdir,'detl');
  cases=expand.grid(n1=n,n2=n,d1=d,d2=d);
  bad=NULL;
  apply(cases,1,function(case) {
    n1=case['n1']; n2=case['n2']; d1=case['d1']; d2=case['d2'];
    if (is.na(d2)) d2=d1;
    xfile=filename(xdir,base='detl',tail=casename_nndd(n1,n2,d1,d2),suffix='RData');
    rfile=filename(rdir,base='detl',tail=casename_nndd(n1,n2,d1,d2),suffix='RData');
    if (verbose) print(paste(sep=', ',xfile,rfile));
    xdetl=load_nndd(xfile);
    rdetl=load_nndd(rfile);
    if (verbose) print(paste(sep='=','dim(xdetl)',paste(collapse=',',dim(xdetl))));
    if (verbose) print(paste(sep='=','dim(rdetl)',paste(collapse=',',dim(rdetl))));
    ## print(paste(collapse=','.mesr.need));
    xdetl=xdetl[,mesr.need];
    rdetl=rdetl[,mesr.need];
    casename=casename_nndd(n1,n2,d1,d2)
    if (any(xdetl!=rdetl)) {
      print(paste(sep=' ','bad: ',casename));
      bad=c(bad,casename);
      BREAKPOINT();
    } else if (verbose) print(paste(sep=' ','good:',casename));
  });
  if (verbose) print(paste(sep='=','length(bad)',length(bad)));
  if (!is.null(bad)) return(bad);
  'all good';
}
    
