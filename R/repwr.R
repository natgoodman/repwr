#################################################################################
##
## Author:  Nat Goodman
## Created: 17-10-05
##          from swfdr.R created 17-10-02
##          from swfdr_base.R created 17-04-10
##          from fdr.R created 17-04-03
##          from pval_error.R created 16-11-01
## Restart: 18-02-15
##          from scope.R created 17-12-04
##
## Copyright (C) 2018 Nat Goodman.
## 
## Explores relication power. More TBD
##
## This is a basic implementation meant to be run interactively.
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################

## ---- run ----
## run the program
## parameters defined in init
run=function(...) {
  init(...);                     # process parameters & other initialization
  dopre();                       # precalculate or load global data
  dosim();                       # load saved simulations or do new ones
  dores();                       # generate results data
}
## ---- init ----
## initialization.
## process parameters and store in global variables.
## create output directory if necessary.
init=function(
  ## simulation parameters 
  n=10*2^(0:9),                     # sample sizes
  d=round(seq(0,1,by=0.1),digits=5),# population effect sizes (means)
                                    # round to avoid imprecise decimals
  m=1e4,                            # number of instances per pop
  ## analysis parameters
  sig.level=0.05,                   # for conventional significance
  conf.level=0.95,                  # for confidence intervals
  pred.level=0.95,                  # for prediction intervals
                                    # grid for various precacluated data
  dsdz.grid=round(seq(min(d)-3,max(d)+3,by=.05),digits=5),
                                    # round to avoid imprecise decimals
  scope.power=0.33,                 # power for small telescope
  scope.close=0.05,                 # pval for big telescope "close enough" calculation
  dores.allcases=T,                 # do all cases (2500 with defaults) or 1500 select ones
  dores.permute=T,                  # permute s1 data. else self-comparisons meaningless
  ## pval.plot=c(.001,.01,.03,.05,.1), # pvalues for which we plot results
  ## 
  ## program parameters, eg, for output files, error messages, etc.
  scriptname='repwr',            #  
  ## datadir=file.path('data',scriptname), # directory for data files
  datadir='data',                       # directory for data files
  simdir=file.path(datadir,'sim'),      # directory for sim files
  simrdir=file.path(datadir,'simr'),    # directory for simr files
  i1dir=file.path(datadir,'i1'),        # directory for s1 permutation index files
  detldir=file.path(datadir,'detl'),    # directory for detailed results files
  smrydir=file.path(datadir,'smry'),    # directory for summary results files
  figdir=file.path('figure',scriptname), # directory for plots
  id=NULL,                       # identifying info tacked onto filenames
  verbose=F,                     # print progress messages
  ## program control
  load=NA,                       # shorthand for other load params
                                 #   NA means load if file exists
                                 #   T, F mean always or never load
  load.sim=load,                 # load saved simulations
  load.simr=load,                # load saved simr files - sim customized for analysis
  load.i1=load,                  # load saved s1 permutation index files
  load.detl=load,                # load saved analysis detail files
  load.smry=load,                # load saved analysis summary files
  load.data=load,                # load saved top level data files
  save=NA,                       # shorthand for other save params 
                                 #   NA means save unless file exists
                                 #   T, F mean always or never save
  save.sim=save,                 # save simulations (RData format)
  save.simr=save,                # save simr files - sim customized for analysis (RData format)
  save.i1=save,                  # save s1 permutation index files (RData format)
  save.detl=save,                # save analysis detail files (RData format)
  save.smry=save,                # save analysis summary files (RData)
  save.data=save,                # save top level results (RData & txt formats)
  save.plot=save,                # save plots
  save.txt=NA,                   # save results in txt format as well as RData
                                 #   NA means use default rule for type:
                                 #   F for all but top level data
  save.txt.sim=!is.na(save.txt)&save.txt,  # save txt simulations. default F
  save.txt.simr=!is.na(save.txt)&save.txt, # save txt simr data. default F
  save.txt.i1=!is.na(save.txt)&save.txt,   # save txt s1 permutation indexes. default F
  save.txt.detl=!is.na(save.txt)&save.txt, # save txt analysis details. default F
  save.txt.smry=!is.na(save.txt)&save.txt, # save txt case-by-case summaries. default F
  save.txt.data=is.na(save.txt)|save.txt,  # save txt top level results. default T
  keep=NA,                       # shorthand for other keep params 
                                 #   NA means use default keep rule for type:
                                 #   T for all but detl
                                 #   T, F mean always or never keep
  keep.sim=!is.na(keep)&keep,    # keep simulations. default F
  keep.simr=is.na(keep)|keep,    # keep simr data. default T
  keep.i1=is.na(keep)|keep,      # keep s1 permutation indexes. default F
  keep.detl=!is.na(keep)&keep,   # keep analysis details. default F
  keep.smry=!is.na(keep)&keep,   # keep case-by-case summaries. default F
  keep.data=is.na(keep)|keep,    # keep top-level data. default T
                                 #    
  clean=F,                       # remove contents of datadir and figdir and start fresh
  clean.all=clean,               # remove everything in datadir and start fresh
  clean.sim=clean.all,           # clean simulations. default F
  clean.simr=clean.all,          # clean simr data. default F
  clean.i1=clean.all,            # clean s1 permutation indexes. default F
  clean.detl=clean.all,          # clean analysis details. default F
  clean.smry=clean.all,          # clean case-by-case summaries. default F
  clean.data=clean.all,          # clean top-level data. default F
  clean.fig=clean,               # remove everything in figdir and start fresh
  end=NULL                       # placeholder for last parameter
  ) {
  ## assign parameters to global variables
  ## do it before calling any functions that rely on globals
  assign_global();
  ## clean and create output directories and internal memory values as needed
  outdir=c(datadir,simdir,simrdir,i1dir,detldir,smrydir);
  memlist=cq(sim.list,simr.list,i1.list,detl.list,smry.list,data.list);
  if (clean.all) {
    unlink(datadir,recursive=T);
    suppressWarnings(rm(list=memlist,envir=.GlobalEnv));
  } else {
    ## may want to clean specific types
    if (clean.sim) cleanq(sim);
    if (clean.simr) cleanq(simr);
    if (clean.i1) cleanq(i1);
    if (clean.detl) cleanq(detl);
    if (clean.smry) cleanq(smry);
    if (clean.data) cleanq(data,cleandir=F);
 }
  ## create subdirectories. nop if already exist
  sapply(outdir,function(dir) dir.create(dir,recursive=TRUE,showWarnings=FALSE));
  if (clean.fig) unlink(figdir,recursive=T);
  dir.create(figdir,recursive=TRUE,showWarnings=FALSE);
  ## setup in-memory lists to hold simulations, etc.. do carefully in case already setup
  sapply(memlist,function(what) if (!exists(what)) assign(what,list(),envir=.GlobalEnv));
  ## initialize summary type and measures if possible
  init_smry(must.exist=F);
  init_mesr(must.exist=F);
  invisible();
}
## initialize summary types. NOTE: needs smry object!
init_smry=function(
  must.exist=T,
  smry=get_data(smry,must.exist=F),
  smry.type=unique(smry$type),
  end=NULL                       # placeholder for last parameter
  ) {
  if (is.null(smry))
    if (must.exist) 
      stop('Cannot initialize summary type: smry object not in-memory and smry file does not exist')
    else {
      init.smry<<-F;    # so analysis functions will know init_smry not done 
      return(F);
    }
  ## at end, assign smry.type to global variable
  smry.type<<-smry.type;
  init.smry<<-T;         # so analysis functions will know init_smry done
  invisible(T);
}
## initialize measures. NOTE: needs detl object! detl header object will do
init_mesr=function(
  must.exist=T,
  detl=get_data(detl,must.exist=F),
  mesr.all=colnames(detl),
  ## measures grouped by source row in smry
  mesr.fromraw=cq(sig1,sign),
  mesr.frombsln=setdiff(grep('scp',mesr.all,invert=T,value=T),mesr.fromraw),
  mesr.fromsig2=setdiff(mesr.all,c(mesr.frombsln,mesr.fromraw)),
  ## measures grouped by relative row (denominator in rate calculation) in standard interpretation
  mesr.relraw=mesr.fromraw,
  mesr.relsig1=mesr.frombsln,
  mesr.relsig2=mesr.fromsig2,
  ## measures grouped by functional category
  mesr.sig=cq(sig2,sigm),
  ## mesr.dcc=cq(d1.c2,d2.c1,c1.c2,d1.p2,d2.p1,p1.p2),
  mesr.dcc=grep('(d(1|2)\\.(c|p)(1|2))|((c|p)1\\.(c|p)2)',mesr.all,value=T),
  ## mesr.scp=cq(d1.scp2,d2.scp1,d1.scpd2,d2.scpd1),
  mesr.scp=grep('d(1|2)\\.scp(d{0,1})(1|2)',mesr.all,value=T),
  mesr.meta=grep('dm|cm',mesr.all,value=T),
  mesr.other=setdiff(mesr.all,c(mesr.sig,mesr.dcc,mesr.scp,mesr.meta)),
  ## defaults for plotrate and heatrate
  mesr.plotdflt=cq(sig2,sigm,d1.c2,d2.c1,d1.p2,d2.p1),
  mesr.heatdflt=c(mesr.sig,mesr.dcc,mesr.scp,mesr.meta),
  ## measures ordered for legends and such
  mesr.order=c(mesr.sig,mesr.dcc,mesr.scp,mesr.other),
  ## RColorBrewer pallete names for plotrate
  pal.sig='Reds',
  pal.dcc='Blues',
  pal.meta='Greens',
  pal.scp='Purples',
  pal.other='Greys',
  end=NULL                       # placeholder for last parameter
  ) {
  if (is.null(detl))
    if (must.exist) 
      stop('Cannot initialize measures: detl object not in-memory and detl file does not exist')
    else {
      init.mesr<<-F;    # so dosmry will know init_mesr not done 
      return(F);
    }
  ## convert names into palettes of correct size
  col.mesr=do.call(c,sapply(cq(sig,dcc,meta,scp,other),USE.NAMES=F,function(what) {
    pal=get(paste(sep='.','pal',what));
    mesr=get(paste(sep='.','mesr',what));
    ## skip 1st two colors - usually too light - then reverse it so darker colors will be first
    col=colorRampPalette(rev(RColorBrewer::brewer.pal(9,pal)[3:8]))(length(mesr));
    col=setNames(col,mesr);
  }));
  ## at end, assign mesr parameters to global variables
  sapply(grep('mesr',ls(),value=T),function(what) assign(what,get(what),envir=.GlobalEnv));
  init.mesr<<-T;         # so dosmry will know init_mesr done
  invisible(T);
}
## ---- Precalculate or Load Global Data ----
dopre=function() {
  doscope();                     # load or generate small telescopes data
  doscopd();                     # load or generate precalculated 'scopd' data 
  doconfvl();                    # load or generate precalculated confidence interval data
  dopredvl();                    # load or generate precalculated prediction interval data
}
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
  simr.ci=ci(sim,confvl[confvl$n==n,]);
  simr[,names(simr.ci)]=simr.ci;
  simr$d.sign=sign(simr$d.sdz);
  simr$d.abs=abs(simr$d.sdz);
  ## for meta-analysis
  simr$w.meta=with(simr,1/sd_d2t(n,d.sdz)^2);
  simr$d.meta=with(simr,d.sdz*w.meta)
  ## drop ones we no longer need
  simr=subset(simr,select=c(-d.raw,-mean0,-mean1,-sd0,-sd1));
  ## optionally save simr and add to in-memory list. usually do - not too big
  save_simr(simr,n,d,id);
  invisible(simr);
}
## ---- Analysis Functions ----
## dores. compute the results. do it case by case to avoid rereading large detail files
dores=function() {
  ## construct cases
  if (dores.allcases) cases=expand.grid(n1=n,n2=n,d1=d,d2=d)
  else {
    ## NG 18-02-22: this branch is obsolete
    ##   I originally thought it would be too slow to do 'em all, but it's okay
    d=sort(d);                     # just in case...
    ## easy ones first
    exact=expand.grid(n1=n,n2=n,d1=d,d2=NA);         # exact replicas - each value of n x n x d
    s2null=expand.grid(n1=n,n2=n,d1=c(0.5,1),d2=0);  # study1>0, study2=0
    s1null=expand.grid(n1=n,n2=n,d1=0,d2=c(0.5,1));  # study1=0, study2>0
    ## study1 one step bigger than study2 & vice versa
    i1=3:length(d);
    s1bigger=expand.grid(n1=n,n2=n,d1=i1);
    s1bigger$d2=s1bigger$d1-1;
    s1bigger$d1=d[s1bigger$d1];
    s1bigger$d2=d[s1bigger$d2];
    s2bigger=s1bigger;
    s2bigger$d1=s1bigger$d2;
    s2bigger$d2=s1bigger$d1;
    ## combine 'em
    cases=rbind(exact,s2null,s1null,s1bigger,s2bigger);
  }
  ## do it!
  smry=do.call(rbind,apply(cases,1,function(case) {
    n1=case['n1']; n2=case['n2']; d1=case['d1']; d2=case['d2'];
    if (is.na(d2)) d2=d1;
    ## compute detailed results or use saved detl
    detl=dodetl(n1=n1,n2=n2,d1=d1,d2=d2);
    ## compute summaries unless they already exist
    smry=dosmry(detl=detl,n1=n1,n2=n2,d1=d1,d2=d2)
  }));
  ## optionally save detl head and add to in-memory list. usually do - needed for analysis
  detl=get_detl(n1=cases[1,'n1'],n2=cases[1,'n2'],d1=cases[1,'d1'],d2=cases[1,'d2']);
  detl=detl[0,];
  save_data(detl);
  ## optionally save smryand add to in-memory list. usually do - needed for analysis
  save_data(smry);
  invisible(smry);
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
  ## use saved i1 if exists and args permit
  i1=get_i1(n1,n2,d1,d2,id,must.exist=F);
  if (is.null(i1)) {
    m=nrow(s1);
    if (dores.permute) {
      i1=sample.int(m)
      s1=s1[i1,];
    } else i1=seq_len(m);
    ## optionally add i1 to in-memory list and save to file (mostly for testing)
    save_i1(i1,n1,n2,d1,d2,id);
  }
  ## define variables that depend on both studies
  ## prediction intervals
  p1=pi(s1,predvl[predvl$n1==n1&predvl$n2==n2,]);
  p2=pi(s2,predvl[predvl$n1==n2&predvl$n2==n1,]);
  ## small telescope thresholds
  scp1=scope[as.character(n1),as.character(n2)];
  scp2=scope[as.character(n2),as.character(n1)];
  ## my 'misinterpretation' of small telescopes that considers d.sdz
  scpd1=scpd(s1,scopd[scopd$n1==n1&scopd$n2==n2,]);
  scpd2=scpd(s2,scopd[scopd$n1==n2&scopd$n2==n1,]);
  ## meta d, ci, pval
  w.sum=s1$w.meta+s2$w.meta;
  meta.d=(s1$d.meta+s2$d.meta)/w.sum;
  meta.abs=abs(meta.d);
  meta.sd=sqrt(1/(w.sum));
  meta.pval=2*pnorm(-abs(meta.d/meta.sd));
  meta.ci=qnorm(p=0.5+conf.level/2,sd=meta.sd)
  meta.ci.lo=meta.d-meta.ci;
  meta.ci.hi=meta.d+meta.ci;
  ## apply the rules. do 'em in-line for efficiency
  ## use concise terms so detl won't be too wide to view
  ##   d{12m}=effect size for s1, s2, meta
  ##   c{12m}=confidence interval for s1, s2, meta
  ##   scp{12}=small telescope threshold for s1, s2
  ##   term.term means compare the two terms, eg, d1.c2 means s1 d.sdz in s2 confidence interval
  ## significant: s1,s2,meta 
  sig1=s1$pval<=sig.level;
  sig2=s2$pval<=sig.level;
  sigm=meta.pval<=sig.level;
  ## same sign: s1,s2
  sign=s1$d.sign==s2$d.sign;
  ## d in confidence interval: each vs other two
  d1.c2=between(s1$d.sdz,s2$ci.lo,s2$ci.hi);
  d2.c1=between(s2$d.sdz,s1$ci.lo,s1$ci.hi);
  d1.cm=between(s1$d.sdz,meta.ci.lo,meta.ci.hi);
  d2.cm=between(s2$d.sdz,meta.ci.lo,meta.ci.hi);
  dm.c1=between(meta.d,s1$ci.lo,s1$ci.hi);
  dm.c2=between(meta.d,s2$ci.lo,s2$ci.hi);
  ## confidence intervals overlap: all 3 unique pairs
  ##   thanks to https://fgiesen.wordpress.com/2011/10/16/checking-for-interval-overlap/
  ##   for simple calculation I should have known already :)
  c1.c2=s1$ci.lo<=s2$ci.hi & s2$ci.lo<=s1$ci.hi;
  c1.cm=s1$ci.lo<=meta.ci.hi & meta.ci.lo<=s1$ci.hi;
  c2.cm=s2$ci.lo<=meta.ci.hi & meta.ci.lo<=s2$ci.hi;
  ## d in prediction interval: each d vs two intervals we have
  d1.p2=between(s1$d.sdz,p2$pi.lo,p2$pi.hi);
  d2.p1=between(s2$d.sdz,p1$pi.lo,p1$pi.hi);
  dm.p2=between(meta.d,p2$pi.lo,p2$pi.hi);
  dm.p1=between(meta.d,p1$pi.lo,p1$pi.hi);
  ## prediction intervals overlap
  p1.p2=p1$pi.lo<=p2$pi.hi & p2$pi.lo<=p1$pi.hi;
  ## d >= small telescope boundary
  ## NG 18-02-14: to handle negative values, use d.abs instead of d.sdz
  d1.scp2=(s1$d.abs>=scp2)&sign;
  d2.scp1=(s2$d.abs>=scp1)&sign;
  dm.scp2=(meta.abs>=scp2)&sign;
  dm.scp1=(meta.abs>=scp1)&sign;
  ## d >= my 'misinterpretation' of small telescope boundary
  ##   need as.vector else R treats it as matrix. screws up naming in data.frame
  d1.scpd2=as.vector((s1$d.abs>=scpd2)&sign);
  d2.scpd1=as.vector((s2$d.abs>=scpd1)&sign);
  dm.scpd2=as.vector((meta.abs>=scpd2)&sign);
  dm.scpd1=as.vector((meta.abs>=scpd1)&sign);
  ## d2 bigger (actually more extreme) than d1
  big2=(s1$d.abs<=s2$d.abs)&sign;
  ## return as data.frame so dosmry can use as environment
  detl=data.frame(sig1,sig2,sigm,sign,
    d1.c2,d2.c1,d1.cm,d2.cm,dm.c1,dm.c2,c1.c2,c1.cm,c2.cm,d1.p2,d2.p1,dm.p1,dm.p2,p1.p2,
    d1.scp2,d2.scp1,dm.scp2,dm.scp1,d1.scpd2,d2.scpd1,dm.scpd2,dm.scpd1,
    big2);
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
  ## init measures unless already done
  if (!init.mesr) init_mesr(detl=detl);
  ## construct summaries for several cases
  ## 1) raw detail
  smry=data.frame(type='raw',m=nrow(detl),t(colSums(detl)),stringsAsFactors=F);
  ## 2) sig1.
  detl.sig1=with(detl,detl[sig1,]);
  smry[2,]=data.frame('sig1',nrow(detl.sig1),t(colSums(detl.sig1)),stringsAsFactors=F);
  ## 3) baseline = sig1 & sign. all require this
  detl.bsln=with(detl.sig1,detl.sig1[sign,]);
  smry[3,]=data.frame('bsln',nrow(detl.bsln),t(colSums(detl.bsln)),stringsAsFactors=F);
  ## 4) baseline & sig2. small telescopes requires this
  detl.sig2=with(detl.bsln,detl.bsln[sig2,]);
  smry[4,]=data.frame('sig2',nrow(detl.sig2),t(colSums(detl.sig2)),stringsAsFactors=F);
  ## NG 18-03-14: std ill-conceived. causes rates to be incomparable
  ## ## 5) standard measures
  ## 5) baseline | big2. many authors also accept big2
  detl.big2=detl.bsln|detl.bsln$big2;
  smry[5,]=data.frame('big2',nrow(detl.big2),t(colSums(detl.big2)),stringsAsFactors=F);
  ## 6) baseline & (sig2 | big2). just because..
  detl.s2bg=with(detl.bsln,detl.bsln[sig2|big2,]);
  smry[6,]=data.frame('s2bg',nrow(detl.s2bg),t(colSums(detl.s2bg)),stringsAsFactors=F);
  ## tack on n1,n2,d1,d2
  smry=data.frame(n1,n2,d1,d2,smry,row.names=NULL);
  ## optionally save smry and add to in-memory list. usually don't keep -- no point
  save_smry(smry,n1,n2,d1,d2,id);
  invisible(smry);
}
## ---- Small Telescopes Functions ----
## generate & optionally save or load data for small telescopes
## presently just scope - matrix of "close enough" small telescopes d values for n1, n2
doscope=function() {
  ## use saved scope if exists and args permit
  scope=get_data(scope,must.exist=F);
  if (!is.null(scope)) return(invisible(scope));
  ## no saved scope or args say not to use it. calculate it
  if (verbose) print(paste(sep='','>>> doscope'));
  scope=do_scope();
  ## optionally save scope and add to in-memory list. usually do - needed for dodetl
  save_data(scope);
  invisible(scope);  
}
## get matrix of "close enough" small telescopes d values for n1, n2
do_scope=function(n1=parent(n,20),n2=parent(n,40)) {
  cases=expand.grid(n1=n1,n2=n2);
  scope=apply(cases,1,function(case) {
    n1=case['n1']; n2=case['n2'];
    d1=d_sig(n1,sig.level);
    d.scope=uniroot(
      function(d0) p_d2t(n1,d=d1,d0,lower.tail=F)-scope.power,interval=c(0,10))$root;
    ## suppress annoying warning from pt:
    ##   full precision may not have been achieved in 'pnt{final}'
    d.close=suppressWarnings(
      uniroot(function(d) p_d2t(n2,d,d0=d.scope)-scope.close,c(-10,10))$root);
    d.close;})
  scope=matrix(data=scope,nrow=length(n),ncol=length(n))
  rownames(scope)=colnames(scope)=n;
  scope;
}
## generate & optionally save or load precalculated scopd data
## for my 'misinterpretation' of small telescopes that considers d.sdz
doscopd=function() {
  ## use saved scopd if exists and args permit
  scopd=get_data(scopd,must.exist=F);
  if (!is.null(scopd)) return(invisible(scopd));
  ## no saved scopd or args say not to use it. calculate it
  if (verbose) print(paste(sep='','>>> doscopd'));
  scopd=do_scopd();
  ## optionally save scopd and add to in-memory list. usually do - needed for dodetl
  save_data(scopd);
  invisible(scopd);  
}
## generate precalculated scopd data
do_scopd=function(n=parent(n),d=parent(dsdz.grid)) {
  d=d[d>=0];                            # scope works on absolute values
  cases=expand.grid(n1=n,d=d);
  ## scpd1 is min d.pop with acceptable power to detect d1.sdz - depends on s1 only
  ## scpd2 is min d2.sdz close enough to scpd1 - depends on s1, s2
  scpd1=do_scopd1(n1=cases$n1,d=cases$d);
  scpd1=data.frame(cases,scpd1=scpd1);
  colnames(cases)=cq(n2,d);             # rename 'n' column for stylistic consistency
  cases=merge(scpd1,expand.grid(n2=n,d=d),by='d',suffixes=c(1,2));
  scpd2=do_scopd2(n2=cases$n2,scpd1=cases$scpd1);
  scpd2=data.frame(cases,scpd=scpd2);
  scpd2[,cq(n1,n2,d,scpd)];           # select & reorder columns for stylistic consistency
}
## scopd1 is min d.pop with acceptable power to detect d1.sdz - depends on s1 only
## scopd2 is min d2.sdz close enough to scopd1 - depends on s1, s2
do_scopd1=
  Vectorize(function(n1,d)
    suppressWarnings(
      uniroot(function(d0) p_d2t(n1,d,d0,lower.tail=F)-scope.power,interval=c(-10,10))$root));
do_scopd2=
  Vectorize(function(n2,scpd1)
    suppressWarnings(
      uniroot(function(d) p_d2t(n2,d,d0=scpd1)-scope.close,c(-10,10))$root));

## interpolate precalculated scopd data at d.sdz
## simr & scopd should already be filtered to current n1, n2
scpd=function(simr,scopd) {
  scopd=with(scopd,akima::aspline(x=d,y=scpd,xout=simr$d.sdz)$y)
  data.frame(scpd=scopd);
}
## ---- Confidence and Prediction Interval Functions ----
## generate & optionally save or load precalculated confidence intervals
doconfvl=function() {
  ## use saved confvl if exists and args permit
  confvl=get_data(confvl,must.exist=F);
  if (!is.null(confvl)) return(invisible(confvl));
  ## no saved confvl or args say not to use it. calculate it
  if (verbose) print(paste(sep='','>>> doconfvl'));
  confvl=do_confvl();
  ## optionally save confvl and add to in-memory list. usually do - needed for dodetl
  save_data(confvl);
  invisible(confvl);  
}
## generate precalculated confidence intervals
do_confvl=function(n=parent(n),d=parent(dsdz.grid)) {
  cases=expand.grid(n=n,d=d);
  confvl=t(Vectorize(ci_d2t)(n=cases$n,d=cases$d));
  confvl=data.frame(cases,confvl);
  colnames(confvl)=cq(n,d,lo,hi);
  confvl;
}
## interpolate precalculated confidence interval data at d.sdz
## sim & confvl should already be filtered to current n
ci=function(sim,confvl) {
  lo=with(confvl,akima::aspline(x=d,y=lo,xout=sim$d.sdz)$y);
  hi=with(confvl,akima::aspline(x=d,y=hi,xout=sim$d.sdz)$y);
  data.frame(ci.lo=lo,ci.hi=hi)
}
## generate & optionally save or load precalculated prediction intervals
dopredvl=function() {
  ## use saved predvl if exists and args permit
  predvl=get_data(predvl,must.exist=F);
  if (!is.null(predvl)) return(invisible(predvl));
  ## no saved predvl or args say not to use it. calculate it
  if (verbose) print(paste(sep='','>>> dopredvl'));
  predvl=do_predvl();
  ## optionally save predvl and add to in-memory list. usually do - needed for dodetl
  save_data(predvl);
  invisible(predvl);  
}
## generate precalculated prediction intervals
## assumes doconfvl already done
do_predvl=function(confvl=parent(confvl)) {
  cases=merge(confvl,confvl,by='d',suffixes=c(1,2));
  ## since we already have conf intervals, no need to do full pi_d2t rigamarole
  lo=with(cases,d-sqrt((d-lo1)^2+(hi2-d)^2));
  hi=with(cases,d+sqrt((d-lo2)^2+(hi1-d)^2));
  data.frame(n1=cases$n1,n2=cases$n2,d=cases$d,lo=lo,hi=hi);
}
## interpolate precalculated prediction interval data at d.sdz
## simr & predvl should already be filtered to current n1, n2
pi=function(simr,predvl) {
  lo=with(predvl,akima::aspline(d,y=lo,xout=simr$d.sdz)$y);
  hi=with(predvl,akima::aspline(d,y=hi,xout=simr$d.sdz)$y);
  data.frame(pi.lo=lo,pi.hi=hi)
}

## ---- Plot Functions ----
## run plot function, save if required, label result with function name
doplot=function(what,id=parent(id,NULL)) {
  what=as.character(pryr::subs(what));
  f=get(what,envir=parent.frame(n=2),mode='function');
  dev=f();
  ## function may return multiple plots
  file=
    if (length(dev)==1) filename_plot(what,id)
      else sapply(seq_along(dev), function(i) filename_plot(what,id,i=i))
  for (i in seq_along(dev)) {
    if ((is.na(save.plot)&!file.exists(file[i]))|(!is.na(save.plot)&save.plot))
      savePlot(file[i],device=dev[i]);
  }
  setNames(dev,what)
}
## plot rate vs. any of n1, n2, d1, d2
## operates on full smry unless smryx set
##   default for n1 10, 20, 40, 80, 160
##   nx is n multiplier: n2=nx*n1
##   x tells which x variable drives the plot, or to use x vars as given
##     x='asis' also sets xtitle='none' and xaxt='n' so xaxis display will make sense
##   cutoff tells where to draw dashed line indicating 'good' rate, if plot.cutoff=T
##     code sets cutoff to 1-cutoff for 'correct' rates
##   plot.cutoff tells whether to plot cutoff at all
##     default, set in code, is T for 'error' or 'correct' rates, F for 'pos' or 'neg'
##   xtitle tells whether to move single-valued vars from xdata and put in title
##   xaxt tells whether to let R compute x-axis and x-grid:
##     n means we do it; r,s,R are synonyms and mean R does it
##     CAUTION: auto works for reasonable params but not in general...
##   smooth tells whether to smooth data to make plots prettier
##   d is synonym for d1, usually used when d1==d2
##   std means standard interpretaion of rules. CAUTION: rates not comparable in this case!
plotrate_nndd=
  function(smryx,n1=n[1:5],n2=nx*n1,nx=2,d1=d,d2=d1,d=0.5,mesr=mesr.plotdflt,
           rate=cq(pos,neg,sigcorrect,sigerror,effcorrect,efferror,unicorrect,unierror),
           x=cq(auto,n1,n2,d1,d2,asis),xtitle=cq(auto,none,n,d,n1,n2,d1,d2),
           cutoff=sig.level,plot.cutoff=T,
           xaxt=cq(auto,n,r,s,R),xaxt.max=11,smooth=T,
           std=F,from.type=unique(c('bsln',smry.type)),relto.type=unique(c('sig1',smry.type)),
           title=NULL,ylab=NULL,xlab=NULL,xlim=NULL,ylim=c(0,1)) {
    rate=match.arg(rate);
    x=match.arg(x);
    xtitle=match.arg(xtitle);
    xaxt=match.arg(xaxt);
    from.type=match.arg(from.type);
    relto.type=match.arg(relto.type);
    if (missing(plot.cutoff)) {
      if (grepl('correct|error',rate)) plot.cutoff=T else plot.cutoff=F;
    }
    if (grepl('correct',rate)) cutoff=1-cutoff;
    ## get the smry subset for the parameters unless smryx passed as paramter
    if (missing(smryx)) {
      smry=get_data(smry);
      smryx=smry_select(smry,n1,n2,d1,d2,x=x);
    }
    smry=smryx$smry; x=smryx$x;
    ## if (missing(xlim)) xlim=range(x);
    data=data_rate(smry,mesr,rate,std,from.type,relto.type);
    if (is.null(ylab)) ylab=paste(sep='',rate,' rate (relative to ',relto.type,')');
    ## if (is.null(xlab)) xlab=x;
    xdata=subset(data,select=c(n1,n2,d1,d2));
    ydata=subset(data,select=-c(n1,n2,d1,d2));
   if (x=='asis') {
      x=seq_len(nrow(data));
      xtitle='none';
      xaxt='n';
      } else x=xdata[,x];
    if (is.null(title)) {
      xdt=xdata_xtitle(xdata,xtitle);
      xdata=xdt$xdata; xtitle=xdt$xtitle;
      if (!std) title=paste(sep='',rate,' rate (relative to ',relto.type,')')
      else title=paste(sep='',rate,' rate (standard interpretation)');
      title=paste(collapse=' ',c(title,xtitle));
    }
    if (xaxt=='auto') if (nrow(xdata)<=xaxt.max) xaxt='n' else xaxt='s';
    if (xaxt=='n') {
      matplot(x,ydata,type='n',xlab=NA,ylab=ylab,main=title,xlim=xlim,ylim=ylim,xaxt='n');
      xaxis(at=x,labels=xdata,xlab=xlab);
      grid(nx=NA,ny=NULL);
      abline(v=x,lty='dotted',col='lightgray');
    } else {
      matplot(x,ydata,type='n',xlab=xlab,ylab=ylab,main=title,xlim=xlim,ylim=ylim,xaxt='s');
      grid();
    }
    col=mesr2col(ydata);
    if (!smooth) {
      matlines(x,ydata,lty='solid',col=col);
    } else {
      ## smooth ydata so the plot will look nicer
      x.smooth=seq(min(x),max(x),len=100);
      y.smooth=asplinem(x,ydata,xout=x.smooth);
      col.smooth=mesr2col(y.smooth);
      matlines(x.smooth,y.smooth,lty='solid',col=col.smooth);
    }
    matpoints(x,ydata,col=col,pch=16,cex=.75);
    mesr_legend(ydata);
    if (plot.cutoff) abline(h=cutoff,lty='dashed',lwd=0.5)
    dev.cur();
  }
## wrapper for plotrate_nndd which gets params from smry
plotrate_smry=function(smry,...) plotrate_nndd(n1=smry$n1,n2=smry$n2,d1=smry$d1,d2=smry$d2,...)

## draw heat map of rate vs. any of n1, n2, d1, d2
## operates on full smry unless smryx set
##   defaults set n1,n2,d1 to fixed values and d2 to d
##   nx is n multiplier: n2=nx*n1
##   cutoff tells where to switch colors from 'lo' to 'hi'
##     code sets cutoff to 1-cutoff for 'correct' or 'pos' rates
##   xtitle tells whether to move single-valued vars from xdata and put in title
##   d is synonym for d1, usually used when d1==d2
##   std means standard interpretaion of rules. CAUTION: rates not comparable in this case!
heatrate_nndd=
  function(smryx,n1=n[1:5],n2=nx*n1,nx=2,d1=d,d2=d1,d=0.5,mesr=mesr.heatdflt,
           rate=cq(pos,neg,sigcorrect,sigerror,effcorrect,efferror,unicorrect,unierror),
           x=cq(asis,auto,n1,n2,d1,d2),xtitle=cq(none,auto,n,d,n1,n2,d1,d2),
           cutoff=sig.level,
           ## xaxt=cq(n,auto,r,s,R),xaxt.max=11,
           std=F,from.type=unique(c('bsln',smry.type)),relto.type=unique(c('sig1',smry.type)),
           title=NULL,xlab=NULL) {
    rate=match.arg(rate);
    x=match.arg(x);
    xtitle=match.arg(xtitle);
    ## xaxt=match.arg(xaxt);
    from.type=match.arg(from.type);
    relto.type=match.arg(relto.type);
    ## setup colors and breakpoints for heatmap
    if (grepl('correct|pos',rate)) cutoff=1-cutoff;
    hi.brk=seq(0,-log10(cutoff),length.out=101);
    lo.brk=seq(-log10(cutoff),4,length.out=101);
    blues=colorRampPalette(RColorBrewer::brewer.pal(5,'Blues'))(100);
    reds=colorRampPalette(RColorBrewer::brewer.pal(4,'Reds'))(100);
    brk=unique(c(hi.brk,lo.brk));
    col=c(blues,reds);
    ## get the smry subset for the parameters unless smryx passed as paramter
    if (missing(smryx)) {
      smry=get_data(smry);
      smryx=smry_select(smry,n1,n2,d1,d2,x=x);
    }
    smry=smryx$smry; x=smryx$x;
    ## if (missing(xlim)) xlim=range(x);
    data=data_rate(smry,mesr,rate,std,from.type,relto.type);
    ## if (is.null(xlab)) xlab=x;
    xdata=subset(data,select=c(n1,n2,d1,d2));
    ydata=as.matrix(subset(data,select=-c(n1,n2,d1,d2)));
    if (is.null(title)) {
      xdt=xdata_xtitle(xdata,xtitle);
      xdata=xdt$xdata; xtitle=xdt$xtitle;
      if (!std) title=paste(sep='',rate,' rate (relative to ',relto.type,')')
      else title=paste(sep='',rate,' rate (standard interpretation)');
      title=paste(collapse=' ',c(title,xtitle));
    }
    ## if (x=='asis') x=seq_len(nrow(data))
    ## else x=xdata[,x];
    x=seq_len(nrow(data));
    y=seq_len(ncol(ydata));
 
    ## if (xaxt=='auto') if (nrow(xdata)<=xaxt.max) xaxt='n' else xaxt='s';
    ## if (xaxt=='n') {
    ## sort measures by distance from 0 to make nicer picture
    yzero=rbind(zero=0,t(ydata));
    ydist=as.matrix(dist(yzero));
    ycol=rownames(yzero)[order(ydist['zero',])][2:ncol(ydist)];
    ydata=ydata[,ycol,drop=F];
    ## rownames(ydata)=x;
    ## log transform ydata so colors will work better
    ylog=ifelse(ydata<1e-2,2,-log10(ydata));
    image(x,y,ylog,axes=F,breaks=brk,col=col,main=title,xlab=NA,ylab=NA);
    xaxis(at=x,labels=xdata,xlab=xlab,lwd.axis=0);
    axis(2,at=y,labels=ycol,cex.axis=0.75,las=1,lwd=0,line=-0.5);
    box();
    abline(h=y+0.5,lty='dotted',col='lightgray',lwd=0.75);
    ## } else {
    ##   ## TODO: this arm not yet ported! probably not useful anyway...
    ##   matplot(x,ydata,type='n',xlab=xlab,ylab=ylab,main=title,xlim=xlim,ylim=ylim,xaxt='s');
    ##   grid();
    ## }
    ## col=mesr2col(ydata);
    ## matlines(x,ydata,lty='solid',col=col);
    ## matpoints(x,ydata,col=col,pch=16,cex=.75);
    ## mesr_legend(ydata);
    dev.cur();
  }

## filter and interpolate smry to get the data needed for plotting and sort based on x
smry_select=
  function(smry=parent(smry),n1=parent(n1),n2=parent(n2),d1=parent(d1),d2=parent(d2),
           mesr=parent(mesr),x=parent(x),
           std=parent(std),from.type=parent(from.type),relto.type=parent(relto.type)) {
    smry=smry_prune();
    smry=smry_interp();
    ## filter smry and sort based on x
    ## in database parlance, the select part is a natural semijoin
    ## CAUTION: need cbind in case params have incompatible lengths. sigh...
    xdata=suppressWarnings(data.frame(cbind(n1,n2,d1,d2)));
    xdata$i=seq_len(nrow(xdata));
    smry=merge(smry,xdata);
    ## sort based on x
     BREAKPOINT();
   if (x=='auto') {
      ## set x to first var that can drive loop, else 'asis'
      x=do.call(
        c,sapply(cq(n1,n2,d1,d2),simplify=F,
                 function(x) if (length(unique(xdata[[x]]))==nrow(xdata)) x));
      if (is.null(x)) x='asis' else x=x[1];
    }
    BREAKPOINT();
    smry=switch(x,
                asis=smry[with(smry,order(i,type)),],
                n1=smry[with(smry,order(n1,n2,d1,d2,type)),],
                n2=smry[with(smry,order(n2,n1,d1,d2,type)),],
                d1=smry[with(smry,order(d1,d2,n1,n2,type)),],
                d2=smry[with(smry,order(d2,d1,n1,n2,type)),]);
    BREAKPOINT();
    list(smry=subset(smry,select=-i),x=x);
  }
## prune smry so interp will be faster
smry_prune=
  function(smry=parent(smry),n1=parent(n1),n2=parent(n2),d1=parent(d1),d2=parent(d2),
           mesr=parent(mesr),
           std=parent(std),from.type=parent(from.type),relto.type=parent(relto.type)) {
    n1.lo=max(smry$n1[smry$n1<=min(n1)]);
    n1.hi=min(smry$n1[smry$n1>=max(n1)]);
    n2.lo=max(smry$n2[smry$n2<=min(n2)]);
    n2.hi=min(smry$n2[smry$n2>=max(n2)])
    d1.lo=max(smry$d1[smry$d1<=min(d1)]);
    d1.hi=min(smry$d1[smry$d1>=max(d1)]);
    d2.lo=max(smry$d2[smry$d2<=min(d2)]);
    d2.hi=min(smry$d2[smry$d2>=max(d2)]);
    if (!std) types=c(from.type,relto.type) else types=cq(raw,sig1,bsln,sig2);
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
    if (exists(x,envir=parent.frame(n=2))) xout=get(x,envir=parent.frame(n=2))
    else stop(paste('no value for',x,'in function call or parent environment'));
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

## remove single-valued xvars from xdata and put in title
xdata_xtitle=function(xdata,xtitle) {
  if (xtitle=='none') xtitle=NULL
  else {
    ## find single valued xvars we want
    xvars=colnames(xdata);          # start with all of 'em
    if (xtitle=='n') {xvars=grep('^n',xvars)}
    else if (xtitle=='d') {xvars=grep('^d',xvars)}
    else if (xtitle!='auto') {xvars=xtitle};
    x.single=do.call(
      c,sapply(xvars,simplify=F,
               function(x) {xval=unique(xdata[[x]]); if (length(xval)==1) xval}));
    ## remove single-valued vars from xdata and put in title
    if (!is.null(x.single)) {
      ## remove vars from xdata without converting one column data frame into vector (sigh...)
      ## from https://stackoverflow.com/questions/4605206/drop-data-frame-columns-by-name
      xdata=within(xdata,rm(list=names(x.single)));
      ## coalesce pairs with same values
      x.n=x.single[grep('^n',names(x.single),value=T)];
      x.d=x.single[grep('^d',names(x.single),value=T)];
      if (length(x.n)==2 && length(unique(x.n))==1) x.n=c(n=unique(x.n));
      if (length(x.d)==2 && length(unique(x.d))==1) x.d=c(d=unique(x.d));
      ## format values as we want them
      x.n=sapply(x.n,n_pretty);
      x.d=sapply(x.d,d_pretty);
      x.single=c(x.n,x.d);
      ## finally turn it all into a string
      xtitle=paste(collapse=', ',paste(sep='=',names(x.single),x.single));
    } else {
      xtitle=NULL;
    }}
  list(xdata=xdata,xtitle=xtitle);
}

## get values for measures from smry and convert to negative, sigerr, efferr, or unierror
## from, relto are smry types
## if std is TRUE, computes standard interpretaion of rules.
##   CAUTION: rates not comparable in this case!
data_rate=
  function(smry,mesr=mesr.dflt,
           rate=cq(pos,neg,sigcorrect,sigerror,effcorrect,efferror,unicorrect,unierror),
           std=F,from.type=unique(c('bsln',smry.type)),relto.type=unique(c('sig1',smry.type))) {
    rate=match.arg(rate);
    pos.rate=pos_rate(smry,mesr,std,from.type,relto.type);
    if (rate=='pos') return(pos.rate);
    ## else construct function call and evaluate
    rate_fun=paste(sep='_',rate,'rate');
    eval(parse(text=rate_fun))(pos.rate);
  }
## get values for measures from smry
## from, relto are smry types
## std means standard interpretaion of rules. CAUTION: rates not comparable in this case!
pos_rate=
  function(smry,mesr=mesr.dflt,std=F,
           from.type=unique(c('bsln',smry.type)),relto.type=unique(c('sig1',smry.type))) {
    ## init smry.type and mesr unless already done
    if (!init.smry) init_smry();
    if (!init.mesr) init_mesr();
    ## make sure all provided measures are legal
    bad=!(mesr %in% mesr.all);
    if (any(bad)) {
      bad=paste(collapse=', ',mesr[bad]);
      stop(paste('Illegal measure names:',bad));
    }
    smry.bytype=split(smry,smry$type);
    if (!std) {
      from.type=match.arg(from.type);
      relto.type=match.arg(relto.type);
      smry=smry.bytype[[from.type]];
      m.relto=smry.bytype[[relto.type]]$m;
      ## drop=FALSE causes result to remain data.frame even with single mesr
      ## thanks StackOverflow! https://stackoverflow.com/questions/7352254
      pos.rate=smry[,mesr,drop=FALSE]/m.relto;
      param=smry.bytype[[from.type]][,cq(n1,n2,d1,d2)];
    } else {
      ## compute standard measures. CAUTION: rates not comparable in this case!
      m=nrow(smry.bytype$raw);
      count=matrix(nrow=m,ncol=length(mesr),dimnames=list(NULL,mesr))
      pos.rate=matrix(nrow=m,ncol=length(mesr),dimnames=list(NULL,mesr));
      mesr.fromraw=mesr[mesr %in% mesr.fromraw];
      mesr.frombsln=mesr[mesr %in% mesr.frombsln];
      mesr.fromsig2=mesr[mesr %in% mesr.fromsig2];
      if (length(mesr.fromraw))
        count[,mesr.fromraw]=data.matrix(smry.bytype$raw[,mesr.fromraw]);
      if (length(mesr.frombsln))
        count[,mesr.frombsln]=data.matrix(smry.bytype$bsln[,mesr.frombsln]);
      if (length(mesr.fromsig2))
        count[,mesr.fromsig2]=data.matrix(smry.bytype$sig2[,mesr.fromsig2]);
      mesr.relraw=mesr[mesr %in% mesr.relraw];
      mesr.relsig1=mesr[mesr %in% mesr.relsig1];
      mesr.relsig2=mesr[mesr %in% mesr.relsig2];
      if (length(mesr.relraw))
        pos.rate[,mesr.relraw]=data.matrix(count[,mesr.relraw]/smry.bytype$raw$m);
      if (length(mesr.relsig1))
        pos.rate[,mesr.relsig1]=data.matrix(count[,mesr.relsig1]/smry.bytype$sig1$m);
      if (length(mesr.relsig2))
        pos.rate[,mesr.relsig2]=data.matrix(count[,mesr.relsig2]/smry.bytype$sig2$m);
      param=smry.bytype$raw[,cq(n1,n2,d1,d2)];
    }
    ## use 1.01 threshold (not 1) 'cuz interpolation can overshoot
    if (any(pos.rate>1.01))
      stop(paste(sep='',"Some pos_rate values are great than 1. This usually means from.type is 'before' relto.type in the workflow: from.type=",from.type,", relto.type=",relto.type));
    ## clamp pos.rate to [0,1] 'cuz interpolation can under- or over-shoot
    pos.rate=apply(pos.rate,c(1,2),function(y) max(min(y,1),0));
    cbind(param,pos.rate);
  }
## negate pos rate matrix
neg_rate=function(pos.rate) {
  neg.rate=1-subset(pos.rate,select=-c(n1,n2,d1,d2));
  cbind(subset(pos.rate,select=c(n1,n2,d1,d2)),neg.rate);
}
## convert pos.rate matrix to correct and error rate matrices
##   d.correct tells which values of d are supposed to say 'yes'
correct_rate=function(pos.rate,d.correct) {
  data=subset(pos.rate,select=-c(n1,n2,d1,d2));
  correct.rate=apply(data,2,function(data) ifelse(d.correct,data,1-data));
  cbind(subset(pos.rate,select=c(n1,n2,d1,d2)),correct.rate);
}
error_rate=function(pos.rate,d.correct) {
  data=subset(pos.rate,select=-c(n1,n2,d1,d2));
  correct.rate=apply(data,2,function(data) ifelse(d.correct,1-data,data));
  cbind(subset(pos.rate,select=c(n1,n2,d1,d2)),correct.rate);
}
## convert pos rate matrix to significance correct and error rate matrices
## correct rate is pos.rate when d1 and d2 are both non-zero 
## error rate is 1-correct = pos.rate for d1=0 or d2=0, 1-pos.rate otherwise
sigcorrect_rate=function(pos.rate) {
  d.correct=with(pos.rate,(d1!=0&d2!=0));
  correct_rate(pos.rate,d.correct);
}
## sigerror_rate=function(pos.rate) neg_rate(sigcorrect_rate(pos.rate));
sigerror_rate=function(pos.rate) {
  d.correct=with(pos.rate,(d1!=0&d2!=0));
  error_rate(pos.rate,d.correct);
}
## convert pos rate matrix to same-effect correct and error rate matrices
## correct rate is pos.rate for d1=d2, 1-pos.rate otherwise
## error rate is 1-correct = pos.rate for d1!=d2, 1-pos.rate otherwise
effcorrect_rate=function(pos.rate) {
  d.correct=with(pos.rate,(d1==d2));
  correct_rate(pos.rate,d.correct);
}
## efferror_rate=function(pos.rate) neg_rate(effcorrect_rate(pos.rate));
efferror_rate=function(pos.rate) {
  d.correct=with(pos.rate,(d1==d2));
  error_rate(pos.rate,d.correct);
}
## convert pos rate matrix to unified correct and error rate matrices
## correct rate is pos.rate for d1!=0 & d2!=0 & d1=d2, 1-pos.rate otherwise
## error rate is 1-correct = pos.rate for d1=0 | d2=0 | d1!=d2, 1-pos.rate otherwise
unicorrect_rate=function(pos.rate) {
  d.correct=with(pos.rate,(d1!=0&d1==d2));
  correct_rate(pos.rate,d.correct);
}
## unierror_rate=function(pos.rate) neg_rate(unicorrect_rate(pos.rate));
unierror_rate=function(pos.rate) {
  d.correct=with(pos.rate,(d1!=0&d1==d2));
  error_rate(pos.rate,d.correct);
}
## get color for measures in positive or error matrix
mesr2col=function(matrix) col.mesr[colnames(matrix)];
## construct legend for measures in pos or error matrix
mesr_legend=
  function(matrix,where='bottomright',x=NULL,y=NULL,title='measure',
           bty='n',pch=19,cex=par('cex'),title.col='black') {
  mesr=colnames(matrix);
  legend.text=mesr[order(match(mesr,mesr.order))];
  legend.col=col.mesr[legend.text];
  if (is.null(x)) x=where;
  legend(x,y,bty=bty,legend=legend.text,col=legend.col,pch=pch,cex=cex,
         title=title,title.col=title.col)
}

## annotate x axis
##  at is vector of x positions for ticks and labels. if NULL, no ticks or label
##  labels is vector or matrix. if matrix, columns are lines of annotion
##  xlab is overall label for axis
##  names is vector of names for annotation lines. if T, use colnames of labels
## 
xaxis=function(at,labels=NULL,xlab=NULL,tick=T,names=T,
           tcl=-0.4,
           line.labels=seq(.4,by=0.6,len=NCOL(labels)),
           line.xlab=min(max(line.labels)+1.2,par('mar')[1]-1.2),
           cex=par('cex'),col=par('fg'),lty=par('lty'),lwd=par('lwd'),
           cex.labels=0.75*cex,cex.xlab=cex,
           col.axis=col,col.labels=col,col.xlab=col,col.ticks=col,
           lwd.axis=lwd,lwd.ticks=lwd) {
    labels.axis=if (is.logical(labels)) labels else NA;
    axis(side=1,at=at,labels=labels.axis,col=col.axis,lty=lty,lwd=lwd.axis,
         tick=tick,tcl=tcl,col.ticks=col.ticks,lwd.ticks=lwd.ticks);
    if (is.na(labels.axis)) {
      ## emit each label line 
      if (is.vector(labels)) {
        mtext(labels,at=at,line=line.labels[1],cex=cex.labels,col=col.labels,side=1)
      } else {
        if (names&!is.null(colnames(labels))) {
          ## place annotation labels in left margin
          ## formula computes (min.plot-0.05*range of plot. empirically looks good
          at=c((1.05*par('usr')[1]-0.05*par('usr')[2]),at);
          ## at=c(par('usr')[1]-0.4,at); 
          labels=rbind(colnames(labels),labels);
        }
        for (j in seq_len(ncol(labels))) 
          mtext(labels[,j],at=at,line=line.labels[j],cex=cex.labels,col=col.labels,side=1);
      }
    }
    ## add overall label (hopefully) at bottom
    if (!is.null(xlab)) mtext(xlab,line=line.xlab,cex=cex.xlab,col=col.xlab,side=1);
  }

## ---- Save and Load ----
## call with file or software attempts to construct file name
## CAUTION: when called in a fresh workspace, directory variables not yet defined
##### for data keyed by n, d
## save data in RData and optionally txt formats
save_nd=function(data,n,d,id=NULL,file=NULL,what,save,save.txt=F,keep) {
  if (is.null(file)) base=basename_nd(n,d,id,what)
  else base=desuffix(file);
  file=filename(base=base,suffix='RData');
  if ((is.na(save)&!file.exists(file))|(!is.na(save)&save)) {
    save(data,file=file);
    if (save.txt)
      write.table(data,file=filename(base=base,suffix='txt'),sep='\t',quote=F,row.names=F);
  }
  if (keep) keep_nd(data,n,d,id,what=what);
}
## load data from file
load_nd=function(file=NULL,n,d,id=NULL,what) {
  if (is.null(file)) file=filename_nd(n,d,id,what);
  what=load(file=file);               # what is name of saved data
  get(what);                          # return it
}
## get data already in memory (eg, in sim.list) or read from file
##   fail if data does not exist unless must.exist is FALSE
get_nd=function(n,d,id=NULL,what,load,keep,must.exist=T) {
  if (is.na(load)|load) {
    case=casename_nd(n,d,id,short=T);
    what.list=get(paste(sep='.',what,'list'))
    data=what.list[[case]];
    if (is.null(data)) {
      file=filename_nd(n,d,id,what);
      if (file.exists(file)) data=load_nd(file=file)
      else {
        if (must.exist)
          stop(paste(sep=' ',what,'for case',casename_nd(n,d,id),
                     'not in memory and file',file,'does not exist'));
        data=NULL;
      }
    if (keep) keep_nd(data,n,d,id,what=what);
    }}
    else {
      if (must.exist)
        stop(paste(sep=' ','must.exist is TRUE but load is FALSE for',what,'case',
                   casename_nd(n,d,id)));
      data=NULL;
    }
  invisible(data);
}
keep_nd=function(data,n,d,id=NULL,what) {
  case=casename_nd(n,d,id,short=T);
  what.list=get(paste(sep='.',what,'list'))
  what.list[[case]]=data;
  assign(paste(sep='.',what,'list'),what.list,envir=.GlobalEnv);
}
##### for data keyed by  n1, n2, d1, d2
save_nndd=function(data,n1,n2,d1,d2,id=NULL,file=NULL,what,save,save.txt=F,keep) {
  if (is.null(file)) base=basename_nndd(n1,n2,d1,d2,id,what)
  else base=desuffix(file);
  file=filename(base=base,suffix='RData');
  if ((is.na(save)&!file.exists(file))|(!is.na(save)&save)) {
    save(data,file=file);
    if (save.txt)
      write.table(data,file=filename(base=base,suffix='txt'),sep='\t',quote=F,row.names=F);
  }
  if (keep) keep_nndd(data,n1,n2,d1,d2,what=what);
}
## load data from file
load_nndd=function(file=NULL,n1,n2,d1,d2,id=NULL,what) {
  if (is.null(file)) file=filename_nndd(n1,n2,d1,d2,id,what);
  what=load(file=file);               # what is name of saved detl
  get(what);                          # return it
}
## get data already in memory (eg, in detl.list) or read from file
##   fail if data does not exist unless must.exist is FALSE
get_nndd=function(n1,n2,d1,d2,id=NULL,what,load,keep,must.exist=T) {
  if (is.na(load)|load) {
    case=casename_nndd(n1,n2,d1,d2,id,short=T);
    what.list=get(paste(sep='.',what,'list'))
    data=what.list[[case]];
    if (is.null(data)) {
      file=filename_nndd(n1,n2,d1,d2,id,what);
      if (file.exists(file)) data=load_nndd(file=file)
      else {
        if (must.exist)
          stop(paste(sep=' ',what,'for case',casename_nndd(n1,n2,d1,d2,id),
                     'not in memory and file',file,'does not exist'))
        else data=NULL;
      }
      if (keep) keep_nndd(data,n1,n2,d1,d2,what=what);
    }}
    else {
      if (must.exist)
        stop(paste(sep=' ','must.exist is TRUE but load is FALSE for',what,'case',
                   casename_nnd(n1,n2,d1,d2,id)));
      data=NULL;
    }
  invisible(data);
}
keep_nndd=function(data,n1,n2,d1,d2,id=NULL,what) {
  case=casename_nndd(n1,n2,d1,d2,id,short=T);
  what.list=get(paste(sep='.',what,'list'))
  what.list[[case]]=data;
  assign(paste(sep='.',what,'list'),what.list,envir=.GlobalEnv);
}
##### sim
save_sim=function(sim,n,d,id=NULL,file=NULL,save=save.sim,save.txt=save.txt.sim,keep=keep.sim)
  save_nd(sim,n,d,id,file,what='sim',save,save.txt,keep);
load_sim=function(file=NULL,n,d,id=NULL) load_nd(file,n,d,id,what='sim');
get_sim=function(n,d,id=NULL,load=load.sim,keep=keep.sim,must.exist=T)
  get_nd(n,d,id,what='sim',load,keep,must.exist);
##### simr
save_simr=function(simr,n,d,id=NULL,file=NULL,save=save.simr,save.txt=save.txt.simr,keep=keep.simr)
  save_nd(simr,n,d,id,file,what='simr',save,save.txt,keep);
load_simr=function(file=NULL,n,d,id=NULL) load_nd(file,n,d,id,what='simr');
get_simr=function(n,d,id=NULL,load=load.simr,keep=keep.simr,must.exist=T)
  get_nd(n,d,id,what='simr',load,keep,must.exist);

##### i1 - s1 permutation indexes (mostly for testing)
save_i1=function(i1,n1,n2,d1,d2,id=NULL,file=NULL,save=save.i1,save.txt=save.txt.i1,keep=keep.i1)
  save_nndd(i1,n1,n2,d1,d2,id,file,what='i1',save,save.txt,keep);
load_i1=function(file=NULL,n1,n2,d1,d2,id=NULL)
  load_nndd(file,n1,n2,d1,d2,id,what='i1');
get_i1=function(n1,n2,d1,d2,id=NULL,load=load.i1,keep=keep.i1,must.exist=T)
  get_nndd(n1,n2,d1,d2,id,what='i1',load,keep,must.exist);
##### detl - detailed rule results
save_detl=function(detl,n1,n2,d1,d2,id=NULL,file=NULL,
                   save=save.detl,save.txt=save.txt.detl,keep=keep.detl)
  save_nndd(detl,n1,n2,d1,d2,id,file,what='detl',save,save.txt,keep);
load_detl=function(file=NULL,n1,n2,d1,d2,id=NULL)
  load_nndd(file,n1,n2,d1,d2,id,what='detl');
get_detl=function(n1,n2,d1,d2,id=NULL,load=load.detl,keep=keep.detl,must.exist=T)
  get_nndd(n1,n2,d1,d2,id,what='detl',load,keep,must.exist);
##### smry - summary rule results
save_smry=function(smry,n1,n2,d1,d2,id=NULL,file=NULL,
                   save=save.smry,save.txt=save.txt.smry,keep=keep.smry)
  save_nndd(smry,n1,n2,d1,d2,id,file,what='smry',save,save.txt,keep);
load_smry=function(file=NULL,n1,n2,d1,d2,id=NULL)
  load_nndd(file,n1,n2,d1,d2,id,what='smry');
get_smry=function(n1,n2,d1,d2,id=NULL,load=load.smry,keep=keep.smry,must.exist=T)
  get_nndd(n1,n2,d1,d2,id,what='smry',load,keep,must.exist);

##### data - top-level data saved in datadir
## save data in RData and txt formats
save_data=function(what,file=NULL,data=NULL,id=NULL,
                   save=save.data,save.txt=save.txt.data,keep=keep.data) {
  what=as.character(pryr::subs(what));
  if (missing(data) && exists(what,envir=parent.frame(n=1)))
    data=get(what,envir=parent.frame(n=1));
  if (is.null(data)) stop('Trying to save NULL object. Is "what" set correctly?');
  if (is.null(file)) base=basename_data(what,id)
  else base=desuffix(file);
  file=filename(base=base,suffix='RData');
  if ((is.na(save)&!file.exists(file))|(!is.na(save)&save)) {
    save(data,file=file);
    if (save.txt)
      write.table(data,file=filename(base=base,suffix='txt'),sep='\t',quote=F,row.names=F);
  }
  if (keep) keep_data(name=what,data=data);
  invisible(data);
}
## load data from file
load_data=function(file=NULL,what=NULL,id=NULL) {
  if (is.null(file)&is.null(what))
    stop('Cannot load data unless file or what is set');
  if (is.null(file)) file=filename_data(what,id);
  what=load(file=file);               # what is name of saved data
  get(what);                          # return it
}
## get top-level data already in memory or read from file
##   fail if data does not exist unless must.exist is FALSE
get_data=function(what,id=NULL,load=load.data,keep=keep.data,must.exist=T,name=NULL) {
  if (is.na(load)|load) {
    if (!is.null(name)) what=name else what=as.character(pryr::subs(what));
    data=data.list[[what]];
    if (is.null(data)) {
      file=filename_data(what,id);
      if (file.exists(file)) {
        data=load_data(file=file);          # what is name of saved data
      } else {
        if (must.exist)
          stop(paste(sep=' ',what,'not in memory and file',file,'does not exist'));
        data=NULL;
      }
      if (keep&&!is.null(data)) keep_data(name=what,data=data);
    }}
    else {
      if (must.exist)
        stop(paste(sep=' ','must.exist is TRUE but load is FALSE for',what));
      data=NULL;
    }
  invisible(data);
}
## keep top-level data. assign globally and keep in global list
keep_data=function(what,data=NULL,name=NULL) {
 if (!is.null(name)) what=name else what=as.character(pryr::subs(what));
  if (missing(data) && exists(what,envir=parent.frame(n=1)))
    data=get(what,envir=parent.frame(n=1));
  if (is.null(data)) stop("Trying to keep NULL object. Is 'what' set correctly?");
  assign(what,data,envir=.GlobalEnv); # assign globally
  data.list[[what]]<<-data;  # assign to global list
}

##### plot - save one or more plots
save_plot=function(dev,file=NULL,what=NULL,id=NULL) {
  if (is.null(file)&is.null(what)) stop('Cannot save plot unless file or what is set');
  if (is.null(file)) {
    file=
      if (length(dev)==1) filename_plot(what,id)
        else sapply(seq_along(dev), function(i) filename_plot(what,id=id,i=i));
  } else file=filename(desuffix(file),suffix='png');
  for (i in seq_along(dev)) savePlot(file[i],device=dev[i]);
}

## ---- Statistical Functions ----
## t.test related functions for one and two samples
## two sample functions all assume equal samples size and variance

## My formula for pooled_sd, though independent of n, is correct. It's a simplification
##   of the standard formula for n0=n1
##   standard formula: sqrt(((n-1)*sd0^2+(n-1)*sd1^2)/(n+n-2));
pooled_sd=function(sd0,sd1) sqrt((sd0^2+sd1^2)/2);

## t statistic (not used)
tstat=function(n,mean0,mean1,sd0,sd1) (mean0-mean1)/sqrt((sd0^2+sd1^2)/n);
## t statistic to Cohen's d & pval
t2d=function(n,t) t*sqrt((2*n)/n^2)
t2pval=function(n,t) 2*pt(-abs(t),df=2*(n-1))
## Cohen's d to t statistic & pval
d2t=function(n,d) d*sqrt(n^2/(2*n))
d2pval=function(n,d) t2pval(n,d2t(n,d))
## pval to t statistic & Cohen's d
pval2t=function(n,pval) qt(pval/2,df=2*(n-1),lower.tail=F)
pval2d=function(n,pval) q_d2t(n,q=pval/2,lower.tail=F)
## confidence interval of d.raw
ci_draw=function(n,d.raw,sd,conf.level=0.95) {
  p0=(1-conf.level)/2; p1=1-p0;
  tstar=sd/sqrt(n/2)*qt(p1,df=2*n-2);
  setNames(c(d.raw-tstar,d.raw+tstar),paste(sep='',100*c(p0,p1),'%'));
}
## significance boundary for Cohen's d
d_sig=function(n,sig.level=parent(sig.level,0.05)) pval2d(n,pval=sig.level)

## probability functions for t-distribution of d
ncp=function(n,d) sqrt(n/2)*d
d_d2t=function(n,d,d0=NULL) {
  df=2*(n-1);
  t=d2t(n,d);
  if (!is.null(d0)) suppressWarnings(dt(t,df=df,ncp=ncp(n,d0)))
  else dt(t,df=df)
}
p_d2t=function(n,d,d0=NULL,lower.tail=TRUE) {
  df=2*(n-1);
  t=d2t(n,d);
  if (!is.null(d0)) suppressWarnings(pt(t,df=df,ncp=ncp(n,d0),lower.tail=lower.tail))
    else pt(t,df=df,lower.tail=lower.tail)
}
q_d2t=function(n,q,d0=NULL,lower.tail=TRUE) {
  df=2*(n-1);
  if (!is.null(d0)) t=suppressWarnings(qt(q,df=df,ncp=ncp(n,d0),lower.tail=lower.tail))
    else t=qt(q,df=df,lower.tail=lower.tail)
  t2d(n,t);
}
r_d2t=function(m,n,d0=NULL) {
  df=2*(n-1);
  if (!is.null(d0)) t=suppressWarnings(rt(m,df=df,ncp=ncp(n,d0))) else t=rt(m,df=df);
  t2d(n,t)
}
## mean and sd for t-distribution of d 
mean_d2t=function(n,d0=NULL) {
  df=2*(n-1);
  if (!is.null(d0)) {
    ncp=ncp(n,d0);
    ## NG 18-02-07. gamma blows up when n>100 or so. use lgamma instead
    ## theo.mean=sqrt(df/2)*ncp*gamma((df-1)/2)/gamma(df/2)
    theo.mean=sqrt(df/2)*ncp*exp(lgamma((df-1)/2)-lgamma(df/2))
    t2d(n,theo.mean)
  } else 0;
}
sd_d2t=function(n,d0=NULL) {
  df=2*(n-1);
  if (!is.null(d0)) {
    ncp=ncp(n,d0);
    ## NG 18-02-07. gamma blows up when n>100 or so. use lgamma instead
    ## theo.mean=sqrt(df/2)*ncp*gamma((df-1)/2)/gamma(df/2)
    theo.mean=sqrt(df/2)*ncp*exp(lgamma((df-1)/2)-lgamma(df/2))
    theo.var=(1+ncp^2)*df/(df-2)-(theo.mean^2)
    theo.sd=sqrt(theo.var)
    t2d(n,theo.sd)
  } else
    (sqrt(2*n)/n)*sdt(2*(n-1));
}
## sd of (central) t distribution
sdt=function(df) sqrt(df/(df-2))

## confidence and prediction intervals  for t-distribution of d
## my adaptation of confidence interval function from
## http://urisohn.com/sohn_files/BlogAppendix/Colada20.ConfidenceIntervalsForD.R
## can possibly make it a bit faster by unwrapping p_d2t
ci_d2t=function(n,d,conf.level=0.95) {
  p0=(1-conf.level)/2; p1=1-p0;
  ci.lo=suppressWarnings(
    uniroot(function(d0) p_d2t(n,d,d0,lower.tail=F)-p0,interval=c(-10,10))$root);
  ci.hi=suppressWarnings(
    uniroot(function(d0) p_d2t(n,d,d0,lower.tail=F)-p1,interval=c(-10,10))$root);
  c(ci.lo,ci.hi);
}
## my adaptation of prediction interval function from
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5028066/ and predictionInterval pacakge
pi_d2t=function(n1,n2,d,ci1=NULL,ci2=NULL,pred.level=0.95) {
  if (is.null(ci1)) ci1=ci_d2t(n1,d);
  if (is.null(ci2)) ci2=ci_d2t(n2,d);
  l1=ci1[1]; u1=ci1[2];
  l2=ci2[1]; u2=ci2[2];
  c(d-sqrt((d-l1)^2+(u2-d)^2),d+sqrt((d-l2)^2+(u1-d)^2));
}
## ---- File Functions ----
##### names for data keyed by n, d
filename_nd=function(n,d,id=NULL,what,suffix='RData') 
  filename(basename_nd(n,d,id,what),suffix=suffix);
basename_nd=function(n,d,id=NULL,what)
  filename(get(paste(sep='',what,'dir')),base=what,tail=casename_nd(n,d,id))
casename_nd=function(n,d,id=NULL,short=F) {
  casename=
    if (short) paste(sep=',',n,d_pretty(d))
      else paste(sep=',',paste_nv(n),paste_nv(d,d_pretty(d)));
  paste_id(casename,id);
}
##### names for data keyed by n1, n2, d1, d2
filename_nndd=function(n1,n2,d1,d2,id=NULL,what,suffix='RData') 
  filename(basename_nndd(n1,n2,d1,d2,id,what),suffix=suffix);
basename_nndd=function(n1,n2,d1,d2,id=NULL,what)
  filename(get(paste(sep='',what,'dir')),base=what,tail=casename_nndd(n1,n2,d1,d2,id))
casename_nndd=function(n1,n2,d1,d2,id=NULL,short=F) {
  if (d1==d2) {
      casename=
        if (short) paste(sep=',',n1,n2,d_pretty(d1))
          else paste(sep=',',paste_nv(n1),paste_nv(n2),paste_nv('d1=d2',d_pretty(d1)));
    } else {
      casename=
        if (short) paste(sep=',',n1,n2,d_pretty(d1),d_pretty(d2))
        else paste(sep=',',paste_nv(n1),paste_nv(n2),
                   paste_nv(d1,d_pretty(d1)),paste_nv(d2,d_pretty(d2)));
    }
  paste_id(casename,id);
}
##### sim
filename_sim=function(n,d,id=NULL,suffix='RData') filename_nd(n,d,id,what='sim',suffix);
basename_sim=function(n,d,id=NULL) basename_nd(n,d,id,what='sim');
casename_sim=casename_nd;
##### simr
filename_simr=function(n,d,id=NULL,suffix='RData') filename_nd(n,d,id,what='simr',suffix);
basename_simr=function(n,d,id=NULL) basename_nd(n,d,id,what='simr');
casename_simr=casename_nd;
##### i1 (mostly for testing)
filename_i1=function(n1,n2,d1,d2,id=NULL,suffix='RData')
  filename_nndd(n1,n2,d1,d2,id,what='i1',suffix);
basename_i1=function(n1,n2,d1,d2,id=NULL) basename_nndd(n1,n2,d1,d2,id,what='i1');
casename_i1=casename_nndd;
##### detl
filename_detl=function(n1,n2,d1,d2,id=NULL,suffix='RData')
  filename_nndd(n1,n2,d1,d2,id,what='detl',suffix);
basename_detl=function(n1,n2,d1,d2,id=NULL) basename_nndd(n1,n2,d1,d2,id,what='detl');
casename_detl=casename_nndd;
##### smry
filename_smry=function(n1,n2,d1,d2,id=NULL,suffix='RData')
  filename(basename_smry(n1,n2,d1,d2,id),suffix=suffix);
basename_smry=function(n1,n2,d1,d2,id=NULL) basename_nndd(n1,n2,d1,d2,id,what='smry');
casename_smry=casename_nndd;
##### data - arbitrary objects saved in datadir
filename_data=function(what,id=NULL,suffix='RData')
  filename(basename_data(what,id),suffix=suffix);
basename_data=function(what,id=NULL) filename(datadir,base=paste_id(what,id));
##### plot - saved in figdir. may have numeric tail
filename_plot=function(what,id=NULL,i=NULL,suffix='png')
  filename(basename_plot(what,id,i),suffix=suffix);
basename_plot=function(what,id=NULL,i=NULL) {
  if (!is.null(i)) i=sprintf("%02i",i);
  basename=filename(figdir,base=what,tail=i);
  paste_id(basename,id);
}
## generate name=value
paste_nv=function(name,value,sep='=') {
  name=as.character(pryr::subs(name));
  if (missing(value))
    if (exists(name,envir=parent.frame(n=2))) value=get(name,envir=parent.frame(n=2))
    else stop(paste('no value for',name,'in function call or parent environment'));
  paste(sep=sep,name,value); 
}
## generate list of name=value using values from parent environment. code adapted from base::rm
nvq=function(...,sep=' ') {
  dots=match.call(expand.dots=FALSE)$...
   if (length(dots) &&
     !all(vapply(dots,function(x) is.symbol(x) || is.character(x),NA,USE.NAMES=FALSE))) 
     stop("... must contain names or character strings");
  names=vapply(dots,as.character,"");
  values=sapply(names,function(name)
    if (exists(name,envir=parent.frame(n=2))) get(name,envir=parent.frame(n=2))
    else stop(paste('no value for',name,'in function call or parent environment')));
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
m_pretty=function(m) sub('e\\+0{0,1}','e',sprintf("%0.0e",m),perl=TRUE)

## construct file or directory pathname from components
## wrapper for file.path with base, tail and suffix pasted on
##  base appended with '.'
##  tail components combined with '.'
##  suffix added unless already there
filename=function(...,base=NULL,tail=NULL,suffix=NULL) {
  if (!is.null(base)||!is.null(tail)) base=paste(collapse='.',c(base,tail));
  if (is.null(base)) file=file.path(...) else file=file.path(...,base);
  if (!is.null(suffix)) {
    ## remove leading '.' if present
    suffix=sub('^\\.','',suffix,perl=T);
    suffix.pattern=paste(collapse='|',paste(sep='','\\.',suffix,'$'));
    file=ifelse(grepl(suffix.pattern,file),file,paste(sep='.',file,suffix[1]));
  }
  file;
}
## remove suffix from filename
desuffix=function(file,suffix=c('RData','txt')) {
  if (!is.null(suffix)) {
    ## remove leading '.' if present
    suffix=sub('^\\.','',suffix,perl=T);
    suffix.pattern=paste(collapse='|',paste(sep='','\\.',suffix,'$'));
    file=sub(suffix.pattern,'',file);
  }
  file;
}
## filebasename same as filename but w/o suffix
filebasename=function(...) filename(...,suffix=NULL)
## construct directory pathname. synonym for filebasename
dirname=filebasename;

####################
## ---- Utility Functions ----

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
## clean specific data type. deletes directory, any top level files and in-memory list
cleanq=function(what,cleandir=T) {
  what=as.character(pryr::subs(what));
  ## delete top level files if exist
  unlink(sapply(cq(RData,txt), function(suffix) filename_data(what,suffix=suffix)));
  ## delete in-memory list
  whatlist=paste(sep='.',what,'list');
  if (exists(whatlist,envir=.GlobalEnv)) rm(list=whatlist,envir=.GlobalEnv);
  if (cleandir) {
    whatdir=paste(sep='',what,'dir');
    ## delete directory if exists
    if (exists(whatdir,envir=.GlobalEnv)) unlink(get(whatdir,envir=.GlobalEnv),recursive=T);
  }
}
## extend akima::aspline for matrix
asplinem=function(x,y,xout,...) {
  if (is.vector(y)) return(akima::aspline(x,y,xout,...));
  if (length(dim(y))!=2) stop('y must be vector or 2-dimensional matrix-like object');
  yout=apply(y,2,function(y) akima::aspline(x,y,xout,...)$y);
  ## if yout has single row (ie, xout has one element), R turns it into a vector...
  if (length(xout)==1) yout=t(yout);
  yout;
}
## not in - based on example in RefMan - more intutive than !%in%
"%notin%"=function(x,table) match(x,table,nomatch=0)==0
## between, near - to subset sim results
between=function(x,lo,hi) x>=lo&x<hi
near=function(x,target,tol=.01) between(x,target-tol,target+tol)

## like lines, but sets color (and other graphical params) for each segment
## from https://stackoverflow.com/questions/7744379/elegant-way-to-select-the-color-for-a-particular-segment-of-a-line-plot
## clever.
##   uses head(x,-1) to get all but last element of vector, x[-1] to get all but first
##   draws segments from x[i],y[i] to x[i+1],y[i+1]
seglines=function(x,y,...) segments(head(x,-1),head(y,-1),x[-1],y[-1],...)

########################################
## color ramp for pvalues.
## adapted from Stack Overflow https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
try(library(RColorBrewer), silent=TRUE);
if ('RColorBrewer' %in% .packages()) {
  ## generate color ramps
  blues=colorRampPalette(brewer.pal(5,'Blues'));
  reds=colorRampPalette(rev(brewer.pal(4,'Reds')));
  ## greens=colorRampPalette(brewer.pal(4,'Greens'));
} else {
  ## use hardcoded values previoulsy generated
  blues=colorRampPalette(c("#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C"));
  reds=colorRampPalette(rev(c("#FEE5D9","#FCAE91","#FB6A4A","#CB181D")));
}
nonsig.breaks=seq(0,-log10(0.05),length.out=11)
sig.breaks=seq(-log10(0.05),4,length.out=11)
pval2col=function(pval,scale=F) {
  if (scale) {
    pval=ifelse(pval<1e-4,1e-4,pval);
    sig.pval=pval[pval<=0.05];
    nonsig.pval=pval[pval>0.05];
    if (length(sig.pval>0))
      sig.breaks=seq(-log10(min(sig.pval)),-log10(max(sig.pval)),length.out=11);
    if (length(nonsig.pval>0))
      nonsig.breaks=seq(-log10(min(nonsig.pval)),-log10(max(nonsig.pval)),length.out=11);
    log.pval=-log10(pval);
    ifelse(pval<=0.05,reds(11)[cut(log.pval,sig.breaks,include.lowest=T)],blues(11)[cut(log.pval,nonsig.breaks,include.lowest=T)]);
  } else {
    ## this version does not scale colors to pvals we have
    log.pval=ifelse(pval<1e-4,4,-log10(pval));
    ifelse(pval<=0.05,reds(11)[cut(log.pval,sig.breaks,include.lowest=T)],blues(11)[cut(log.pval,nonsig.breaks,include.lowest=T)]);
  }
}
## construct pval legend
pval_legend=function(x0,y0,steps=100,width=0.2,height=1) {
  sig.pval=seq(0,0.05,length=steps);
  nonsig.pval=seq(0.05,1,length=steps);
  y=seq(y0,by=height/(2*steps),length=2*steps);
  rect(x0,head(y,-1),x0+width,y[-1],col=pval2col(c(sig.pval,nonsig.pval)),border=NA);
  text(x0,y0+(1.05*height),"pvalues",adj=0,cex=.8);
  text(x0+width,y0,'0',pos=4,cex=.8);
  text(x0+width,y0+height/2,'0.05',pos=4,cex=.8);
  text(x0+width,y0+height,'1',pos=4,cex=.8);
}

## run data-making function, save if required, store result in global workspace
dodata=function(what) {
  what=as.character(pryr::subs(what));
  f=get(what,envir=parent.frame(n=2),mode='function');
  data=f();
  what=sub('^do_','',what);
  if (save.data) save_data(data,what);
  ## fiddle with name to get the form we want
  assign(what,data,envir=.GlobalEnv);
  invisible(data);
}
