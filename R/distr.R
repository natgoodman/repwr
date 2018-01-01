#################################################################################
##
## Author:  Nat Goodman
## Created: 17-12-20
##          from scope.R created 17-12-04 (actually scope.with.automode.R)
##          from repwr.R created 17-10-05
##          from swfdr.R created 17-10-02
##          from swfdr_base.R created 17-04-10
##          from fdr.R created 17-04-03
##          from pval_error.R created 16-11-01
##
## Copyright (C) 2017 Nat Goodman.
## 
## Explores sampling and populating distributions
##
## This is a basic implementation meant to be run interactively.  It
## is a minimal implementation that supports an email conversation
## circa Dec 2017
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
  dosim();                       # load saved simulations or do new ones
  dores();                       # generate results, mostly plots
}
## ---- init ----
## initialization.
## process parameters and store in global variables.
## create output directory if necessary.
init=function(
  ## simulation parameters 
  m=1e6,                         # total number of instances
  n=20,                          # sample size 
  d.pop=runif(m,-1.5,3),         # 'manual' mode: population effect sizes (means)
  d.test=1.5,                    # 'auto' mode: max or range of d.pop we intend to test
                                 #   can also be list keyed by n
  init.mode='manual',            # choices: 'manual','automatic'
                                 #   'manual': software uses d.pop as given
                                 #   'auto': software sets d.pop
  ## analysis parameters (more TBD)
  ## conf.level=0.95,            # for confidence intervals
  m.plot=1e6,                    # params to select simulations for plots
  n.plot=20,                     #
  id.plot=NULL,                  #
  d.qq=0.5,                      # d for qq plots
  tol.qq=0.01,                   # range around d.qq for qq plots
  dlim.plot=c(0,d.test),         # d limits for sd vs d plots
  sig.level=0.05,
  ## pval.plot=c(.001,.01,.03,.05,.1), # pvalues for which we plot results
  ## 
  ## program parameters, eg, for output files, error messages, etc.
  scriptname='distr',            #  
  datadir=file.path('data',scriptname), # directory for data files
  figdir=file.path('figure',scriptname), # directory for plots
  id=NA,                         # identifying info tacked onto filenames
  verbose=FALSE,                 # print progress messages
  ## program control
  load=NA,                       # shorthand to override defaults for load.sim & load.res
  load.sim=load,                 # load saved simulations in automatic mode
                                 #   NA means load if file exists
                                 #   T, F mean always or never load
  save=FALSE,                    # shorthand for save.sim & save.res 
  save.sim=save,                 # save simulations (RData format)
  save.res=save,                 # save results (RData & txt formats); plots (png)
  save.plot=save.res,            # save plots. used to override save.res
  save.simtxt=FALSE,             # save simulations in txt format - big & slow!
  clean=FALSE,                   # remove contents of datadir and figdir and start fresh
  clean.data=clean,              # remove contents of datadir and start fresh
  clean.fig=clean,               # remove contents of figdir and start fresh
  end=NULL                       # placeholder for last parameter
  ) {
  init.mode=match.arg(init.mode,c('manual','automatic'),several.ok=T);
  cases=expand.grid(m=m,n=n,id=id,load.sim=load.sim,init.mode=init.mode,stringsAsFactors=F);
  ## clean and create output directories as needed
  if (clean.data) unlink(datadir,recursive=TRUE);
  dir.create(datadir,recursive=TRUE,showWarnings=FALSE);
  if (clean.fig) unlink(figdir,recursive=TRUE);
  dir.create(figdir,recursive=TRUE,showWarnings=FALSE);
  ## setup 'sim' list to hold simulations. do it carefully in case already setup
  if (!exists('sim')) sim=list() else if (!is.list(sim)) sim=list();
  ## at end, assign global parameters to global variables
  ##   don't do it earlier because it's too easy to confuse local and
  ##   global variables with the same names
  m<<-m;
  n<<-n;
  d.pop<<-d.pop;
  d.test<<-d.test;
  init.mode<<-init.mode;
  sim<<-sim;
  fig<<-1;                              # initial figure number
  m.plot<<-m.plot;
  n.plot<<-n.plot;
  id.plot<<-id.plot;
  d.qq<<-d.qq;
  tol.qq<<-tol.qq;
  dlim.plot<<-dlim.plot;
  sig.level<<-sig.level;
  ## conf.level<<-conf.level;
  ## pval.plot<<-pval.plot;
  scriptname<<-scriptname;
  datadir<<-datadir;
  figdir<<-figdir;
  id<<-id;
  verbose<<-verbose;
  load.sim<<-load.sim;
  save.sim<<-save.sim;
  save.res<<-save.res;
  save.simtxt<<-save.simtxt;
  save.plot<<-save.plot;
  clean.data<<-clean.data;         # not really needed, since only used in init
  clean.fig<<-clean.fig;           # not really needed, since only used in init
  cases<<-cases;
  invisible(cases);
}
## ---- Simulation Functions ----
## do the simulation
## parameters in case agument and global variables set by init
dosim=function(save.sim=parent(save.sim)) {
  ## do simulations! case by case
  for (i in seq_len(nrow(cases))) {
    ## extract args this way to preserve type (numeric vs string) and avoid conversion to list
    m=cases[i,'m']; n=cases[i,'n']; id=cases[i,'id'];
    load.sim=cases[i,'load.sim']; init.mode=cases[i,'init.mode'];
    if (is.na(id)) id=NULL;
    case=casename_sim(m,n,id);
    file=filename_sim(m,n,id);
    if (init.mode=='automatic') {
     ## use saved simulation if possible
      if (file.exists(file)) {
        ## saved simulation exists. use if params permit
        if (is.na(load.sim)|load.sim) {
          sim[[case]]=load_sim(file);
          next;
        } else {
          if (load.sim) stop(paste(sep=' ','load.sim is TRUE but file',file,'does not exist'));
        }}
      ## no saved simulation or params prevent its use. run simulation
      ## sd=sd_d2t(n,d), ie, sd of Cohen's d from t-distribution
      ## data range +/- 4*sd will get (essentially) all of it
     if (is.list(d.test)) d.test=d.test[[as.character(n)]];
      if (length(d.test)==1) d.test=c(0,d.test);
      lo=min(d.test); hi=max(d.test); sd.lo=sd_d2t(n,lo); sd.hi=sd_d2t(n,hi);
      d.pop=runif(m,lo-4*sd.lo,hi+4*sd.hi);
    }
    sim2=sim[[case]]=dosim2(m,n,d.pop);
    if (save.sim) save_sim(sim2,file);
  }
  sim<<-sim;
  invisible(sim);
}       
## simulate one two-sample case. adapted from swfdr
dosim2=function(m,n,d.pop) {
  if (verbose) print(paste(sep='','>>> dosim2:',' m=',m,' n=',n));
  sim2=do.call(rbind,lapply(seq_along(d.pop),function(i) {
    if (verbose&(i%%1000==0)) print(paste(sep='','dosim2:',' i=',i));
    d=d.pop[i];
    ## draw pair of samples. group0 is control. group1 has effect=d
    group0=rnorm(n,mean=0);
    group1=rnorm(n,mean=d);
    mean0=mean(group0);
    mean1=mean(group1);
    d.raw=mean1-mean0;
    sd0=sd(group0);
    sd1=sd(group1);
    pval=t.test(group0,group1,var.equal=TRUE)$p.value;
    ## store results
    c(pval=pval,d.raw=d.raw,mean0=mean0,mean1=mean1,sd0=sd0,sd1=sd1)
  }));
  ## convert to data frame and add in input and post-computed outputs
  sim2=as.data.frame(sim2);
  sd=with(sim2,pooled_sd(sd0,sd1));
  sim2=cbind(data.frame(m=m,n=n,d.pop=d.pop,d.sdz=sim2$d.raw/sd),sim2,sd=sd);
  if (verbose) print(paste(sep='','<<< dosim2:',' m=',m,' n=',n));
  invisible(sim2);
}
## ---- Analysis and Plot Functions ----
## do the analysis
dores=function(save.res=parent(save.res)) {
  fig<<-1;                    # figure number. bumped in each plot function
  plots=c(
    doplot(plotsd_sdz_vs_pop),
    doplot(plotsd_raw_vs_pop),
    doplot(plotsd_pop_vs_sdz),
    doplot(plotsd_pop_vs_raw),
    doplot(plotqq_sdz_vs_pop),
    doplot(plotqq_raw_vs_pop),
    doplot(plotqq_pop_vs_sdz),
    doplot(plotqq_pop_vs_raw));
  if (save.res) {
    sapply(names(plots),function(name) save_plot(plots[name],name));
  }
  invisible(plots);
}
## run plot function and  label result with function name
doplot=function(what) {
  what=as.character(pryr::subs(what));
  f=get(what,envir=parent.frame(n=2));
  setNames(f(),what)
}
##### these functions plot sd of one d measure vs. another d measure
## plot sd(d.sdz) vs. d.pop
plotsd_sdz_vs_pop=
  function(fig=parent(fig,1),
           m=parent(m.plot,1e6),n=parent(n.plot,10),id=parent(id.plot,NULL),
           dlim=parent(dlim.plot,c(0,1.5))) {
    dev.new();
    title=paste(sep='','Figure ',fig,". sd of d.sdz vs. d.pop for m=",m.pretty(m),", n=",n);
    sim2=get_sim(m,n,id);
    sim2=subset(sim2,subset=(d.pop>=dlim[1]&d.pop<=dlim[2]));
    sim2=sim2[order(sim2$d.pop),];
    group=cut(sim2$d.pop,breaks=2000,labels=F,ordered_result=TRUE);
    d.pop=sapply(split(sim2$d.pop,group),mean);
    sdz.by.pop=split(sim2$d.sdz,group);
    empi=sapply(sdz.by.pop,sd);
    empi.smooth=loess.smooth(d.pop,empi);
    theo=sapply(d.pop,function(d) sd_d2t(n,d));
    theosim=sapply(d.pop,function(d) sd(r_d2t(1e3,n,d0=d)));
    theosim.smooth=loess.smooth(d.pop,theosim);
    plot(d.pop,empi,pch=20,cex=.5,col='grey',ylab='standard deviation',main=title);
    lines(d.pop,theo,col='gold',lwd=6);
    lines(empi.smooth,col='blue');
    lines(theosim.smooth,col='green');
    grid();
    ## add legend
    legend.text=
      c('empirical','smoothed empirical','theoretical (s_d2t)','simulated theo (r_d2t)');
    legend.col=c('grey','blue','gold','green');
    legend('top',bty='n',legend=legend.text,col=legend.col,pch=19,cex=.8)
    fig<<-fig+1;                # bump for next guy
    dev.cur();
  }
## plot sd(d.raw) vs. d.pop
plotsd_raw_vs_pop=
  function(fig=parent(fig,2),
           m=parent(m.plot,1e6),n=parent(n.plot,10),id=parent(id.plot,NULL),
           dlim=parent(dlim.plot,c(0,1.5))) {
    dev.new();
    title=paste(sep='','Figure ',fig,". sd of d.raw vs. d.pop for m=",m.pretty(m),", n=",n);
    sim2=get_sim(m,n,id);
    sim2=subset(sim2,subset=(d.pop>=dlim[1]&d.pop<=dlim[2]));
    sim2=sim2[order(sim2$d.pop),];
    group=cut(sim2$d.pop,breaks=2000,labels=F,ordered_result=TRUE);
    d.pop=sapply(split(sim2$d.pop,group),mean);
    raw.by.pop=split(sim2$d.raw,group);
    empi=sapply(raw.by.pop,sd);
    empi.smooth=loess.smooth(d.pop,empi);
    theo=rep(sqrt(2/n),length(d.pop));
    theosim=sapply(d.pop,function(d) sd(rnorm(1e3,mean=d,sd=sqrt(2/n))));
    theosim.smooth=loess.smooth(d.pop,theosim);
    plot(d.pop,empi,pch=20,cex=.5,col='grey',ylab='standard deviation',main=title);
    lines(d.pop,theo,col='gold',lwd=6);
    lines(empi.smooth,col='blue');
    lines(theosim.smooth,col='green');
    grid();
    ## add legend
    legend.text=
      c('empirical','smoothed empirical','theoretical (sqrt(2/n))','simulated theo (rnorm)');
    legend.col=c('grey','blue','gold','green');
    legend('top',bty='n',legend=legend.text,col=legend.col,pch=19,cex=.8)
    fig<<-fig+1;                # bump for next guy
    dev.cur();
  }
## plot sd(d.pop) vs. d.sdz
plotsd_pop_vs_sdz=
  function(fig=parent(fig,3),
           m=parent(m.plot,1e6),n=parent(n.plot,10),id=parent(id.plot,NULL),
           dlim=parent(dlim.plot,c(0,1.5))) {
    dev.new();
    title=paste(sep='','Figure ',fig,". sd of d.pop vs. d.sdz for m=",m.pretty(m),", n=",n);
    sim2=get_sim(m,n,id);
    sim2=subset(sim2,subset=(d.sdz>=dlim[1]&d.sdz<=dlim[2]));
    sim2=sim2[order(sim2$d.sdz),];
    group=cut(sim2$d.sdz,breaks=2000,labels=F,ordered_result=TRUE);
    d.sdz=sapply(split(sim2$d.sdz,group),mean);
    pop.by.sdz=split(sim2$d.pop,group);
    empi=sapply(pop.by.sdz,sd);
    empi.smooth=loess.smooth(d.sdz,empi);
    theo=sapply(d.sdz,function(d) sd_d2t(n,d));
    theosim=sapply(d.sdz,function(d) sd(r_d2t(1e3,n,d0=d)));
    theosim.smooth=loess.smooth(d.sdz,theosim);
    plot(d.sdz,empi,pch=20,cex=.5,col='grey',ylab='standard deviation',main=title);
    lines(d.sdz,theo,col='gold',lwd=6);
    lines(empi.smooth,col='blue');
    lines(theosim.smooth,col='green');
    lines(d.sdz,sdd(n,d.sdz),col='red');
    grid();
    ## add legend
    legend.text=
      c('empirical','smoothed empirical','theoretical (s_d2t)','simulated theo (r_d2t)',
        'approximation from compute.es');
    legend.col=c('grey','blue','gold','green','red');
    legend('top',bty='n',legend=legend.text,col=legend.col,pch=19,cex=.8)
    fig<<-fig+1;                # bump for next guy
    dev.cur();
  }
## plot sd(d.pop) vs. d.raw
plotsd_pop_vs_raw=
  function(fig=parent(fig,4),
           m=parent(m.plot,1e6),n=parent(n.plot,10),id=parent(id.plot,NULL),
           dlim=parent(dlim.plot,c(0,1.5))) {
    dev.new();
    title=paste(sep='','Figure ',fig,". sd of d.pop vs. d.raw for m=",m.pretty(m),", n=",n);
    sim2=get_sim(m,n,id);
    sim2=subset(sim2,subset=(d.raw>=dlim[1]&d.raw<=dlim[2]));
    sim2=sim2[order(sim2$d.raw),];
    group=cut(sim2$d.raw,breaks=2000,labels=F,ordered_result=TRUE);
    d.raw=sapply(split(sim2$d.raw,group),mean);
    pop.by.raw=split(sim2$d.pop,group);
    empi=sapply(pop.by.raw,sd);
    empi.smooth=loess.smooth(d.raw,empi);
    theo=rep(sqrt(2/n),length(d.raw));
    theosim=sapply(d.raw,function(d) sd(rnorm(1e3,mean=d,sd=sqrt(2/n))));
    theosim.smooth=loess.smooth(d.raw,theosim);
    plot(d.raw,empi,pch=20,cex=.5,col='grey',ylab='standard deviation',main=title);
    lines(d.raw,theo,col='gold',lwd=6);
    lines(empi.smooth,col='blue');
    lines(theosim.smooth,col='green');
    grid();
    ## add legend
    legend.text=
      c('empirical','smoothed empirical','theoretical (sqrt(2/n))','simulated theo (rnorm)');
    legend.col=c('grey','blue','gold','green');
    legend('top',bty='n',legend=legend.text,col=legend.col,pch=19,cex=.8)
    fig<<-fig+1;                # bump for next guy
    dev.cur();
  }
##### these functions do qq plots of specific slice of d measure
## qq plot of d.sdz vs d.pop
plotqq_sdz_vs_pop=
  function(fig=parent(fig,5),
           m=parent(m.plot,1e6),n=parent(n.plot,10),id=parent(id.plot,NULL),
           d.qq=parent(d.qq,NULL),tol.qq=parent(tol.qq,.01)) {
    dev.new();
    title=paste(sep='','Figure ',fig,". QQ TDist Plot of d.sdz vs. d.pop for m=",m.pretty(m),", n=",n,", d=",d.qq);
    sim2=get_sim(m,n,id);
    sim2=subset(sim2,subset=near(d.pop,d.qq,tol.qq));
    d=sim2$d.sdz;
    qq_d2t(n,d,d.qq,main=title,cex.main=.95);
    qqline_d2t(n,d,d.qq,col='red');
    grid();
    fig<<-fig+1;                # bump for next guy
    dev.cur();
  }
## qq plot of d.raw vs d.pop
plotqq_raw_vs_pop=
  function(fig=parent(fig,5),
           m=parent(m.plot,1e6),n=parent(n.plot,10),id=parent(id.plot,NULL),
           d.qq=parent(d.qq,NULL),tol.qq=parent(tol.qq,.01)) {
    dev.new();
    title=paste(sep='','Figure ',fig,". QQ Normal Plot of d.raw vs. d.pop for m=",m.pretty(m),", n=",n,", d=",d.qq);
    sim2=get_sim(m,n,id);
    sim2=subset(sim2,subset=near(d.pop,d.qq,tol.qq));
    d=sim2$d.raw;
    qqnorm(d,main=title,cex.main=.95);
    qqline(d,col='red');
    grid();
    fig<<-fig+1;                # bump for next guy
    dev.cur();
  }
## qq plot of d.pop vs d.sdz
plotqq_pop_vs_sdz=
  function(fig=parent(fig,5),
           m=parent(m.plot,1e6),n=parent(n.plot,10),id=parent(id.plot,NULL),
           d.qq=parent(d.qq,NULL),tol.qq=parent(tol.qq,.01)) {
    dev.new();
    title=paste(sep='','Figure ',fig,". QQ TDist Plot of d.pop vs. d.sdz for m=",m.pretty(m),", n=",n,", d=",d.qq);
    sim2=get_sim(m,n,id);
    sim2=subset(sim2,subset=near(d.sdz,d.qq,tol.qq));
    d=sim2$d.pop;
    qq_d2t(n,d,d.qq,main=title,cex.main=.95);
    qqline_d2t(n,d,d.qq,col='red');
    grid();
    fig<<-fig+1;                # bump for next guy
    dev.cur();
  }
## qq plot of d.pop vs d.raw
plotqq_pop_vs_raw=
  function(fig=parent(fig,5),
           m=parent(m.plot,1e6),n=parent(n.plot,10),id=parent(id.plot,NULL),
           d.qq=parent(d.qq,NULL),tol.qq=parent(tol.qq,.01)) {
    dev.new();
    title=paste(sep='','Figure ',fig,". QQ Normal Plot of d.pop vs. d.raw for m=",m.pretty(m),", n=",n,", d=",d.qq);
    sim2=get_sim(m,n,id);
    sim2=subset(sim2,subset=near(d.raw,d.qq,tol.qq));
    d=sim2$d.pop;
    qqnorm(d,main=title,cex.main=.95);
    qqline(d,col='red');
    grid();
    fig<<-fig+1;                # bump for next guy
    dev.cur();
  }
## ---- Save and Load ----
## save simulation in RData and optionally txt formats
## CAUTION: txt big & slow!
save_sim=function(sim,file=NULL,m,n,id=NULL,save.simtxt=parent(save.simtxt),F) {
  if (missing(file)) base=basename_sim(m,n,id)
  else base=desuffix(file);
  save(sim,file=filename(base=base,suffix='RData'));
  if (save.simtxt)
    write.table(sim,file=filename(base=base,suffix='txt'),sep='\t',quote=F,row.names=F);
}
## load saved simulation
## call with file or software attempts to construct file name
## CAUTION: when called in a fresh workspace, datadir variable not yet defined
load_sim=function(file=NULL,m,n,id=NULL) {
  if (missing(file)) file=filename_sim(m,n,id);
  what=load(file=file);               # what is name of saved sim
  get(what);                          # return it
}
## get simulation already in memory (sim list) or read from file or fail
get_sim=function(m,n,id=NULL,sim=parent(sim,list())) {
  case=casename_sim(m,n,id);
  sim2=sim[[case]];
  if (is.null(sim2)) {
    file=filename_sim(m,n,id);
    if (file.exists(file)) {
      sim2=load_sim(file=file);
      sim[[case]]<<-sim2;
    } else stop(paste(sep='','simulation results for case: ',case,
                      ' not in memory and file ',file,' does not exist'));
  }
  invisible(sim2);
}
    
## save data frame
save_data=function(data,filebase) {
  base=filebasename(datadir,filebase,tail=id);
  save(data,file=filename(base=base,suffix='RData'));
  write.table(data,file=filename(base=base,suffix='txt'),sep='\t',quote=F,row.names=F);
}
## save one or more plots
## dev is vector of devices
save_plot=function(dev,filebase) {
  base=filebasename(figdir,filebase);
  if (length(dev)==1) savePlot(filename(figdir,filebase,suffix='png'),device=dev)
    else for (i in seq_along(dev)) {
    file=filename(figdir,paste(sep='.',filebase,sprintf("%02i",i)),suffix='png');
    savePlot(file,device=dev[i]);
  }
}

## ---- Statistical Functions ----
## t.test related functions for one and two samples
## two sample functions all assume equal samples size and variance

## My formula for pooled_sd, though independent of n, is correct. It's a simplification
##   of the standard formula for n0=n1
##   standard formula: sqrt(((n-1)*sd0^2+(n-1)*sd1^2)/(n+n-2));
pooled_sd=function(sd0,sd1) sqrt((sd0^2+sd1^2)/2);
## Uri's formula for ncp
ncp=function(n,d) sqrt(n/2)*d

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

## confidence intervals for normal and  t-distribution of d
ci_norm=function(n,d,conf.level=.95) {
  p0=(1-conf.level)/2; p1=1-p0;
  zstar=qnorm(p0)*(2/sqrt(n));
  setNames(c(d-zstar,d+zstar),paste(sep='',100*c(p0,p1),'%'));
}
ci_d2t=function(n,d,d0=NULL,conf.level=.95) {
  p0=(1-conf.level)/2; p1=1-p0;
  lower=q_d2t(n,p0,d0,lower.tail=TRUE);
  upper=q_d2t(n,p0,d0,lower.tail=FALSE);
 setNames(c(d+lower,d+upper),paste(sep='',100*c(p0,p1),'%'));
}

## probability functions for t-distribution of d
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
## sd of Cohen's d for two samples
## sdd (probably wrong in general but seems to work for sd(d.pop) vs d.sdz ...)
##   from https://stats.stackexchange.com/questions/144084/variance-of-cohens-d-statistic
##   also in compute.es::des
## sd_d2t   sd of t-distribution of d
## mean_d2t mean of t-distribution of d
##   adapted from UnivRNG::draw.noncentral.t
sdd=function(n,d) sqrt(((2*n)/n^2)+(d^2/(4*n)))
sd_d2t=function(n,d0=NULL) {
  df=2*(n-1);
  if (!is.null(d0)) {
    ncp=ncp(n,d0);
    theo.mean=sqrt(df/2)*ncp*gamma((df-1)/2)/gamma(df/2)
    theo.var=(1+ncp^2)*df/(df-2)-(theo.mean^2)
    theo.sd=sqrt(theo.var)
    t2d(n,theo.sd)
  } else
    (sqrt(2*n)/n)*sdt(2*(n-1));
}
mean_d2t=function(n,d0=NULL) {
  df=2*(n-1);
  if (!is.null(d0)) {
    ncp=ncp(n,d0);
    theo.mean=sqrt(df/2)*ncp*gamma((df-1)/2)/gamma(df/2)
    t2d(n,theo.mean)
  } else 0;
}

## sd of (central) t distribution
sdt=function(df) sqrt(df/(df-2))

## like qqnorm & qqline but for t-distribution of d. as usual, assumes two sample case, etc
##   why is plot not straight with large d0?? perhaps a problem in ncp approx
qq_d2t=
  function(n,y,d0=NULL,
           main="Q-Q Plot of T Distribution of Cohen's d",
           xlab="Theoretical Quantiles",ylab="Sample Quantiles",
           qnum=1000,...) {
    probs=seq(1/qnum,1-(1/qnum),by=1/qnum);
    theo=q_d2t(n,probs,d0);
    empi=quantile(y,probs=probs);
    plot(theo,empi,main=main,xlab=xlab,ylab=ylab,...);
}
qqline_d2t=function(n,y,d0=NULL,...) {
  probs=c(0.25,0.75);
  theo=q_d2t(n,probs,d0);
  empi=quantile(y,probs=probs);
  slope=diff(empi)/diff(theo);
  int=empi[1]-slope*theo[1];
  abline(int,slope,...);
}


## ---- Utility Functions ----
## file related functions
## filename for sim files
filename_sim=function(m,n,id=NULL,suffix='RData') 
  filename(basename_sim(m,n,id=id),suffix=suffix);
basename_sim=function(m,n,id=NULL) filename(datadir,base='sim2',tail=casename_sim(m,n,id))
casename_sim=function(m,n,id=NULL) {
  ## test id this way to avoid running is.na when id=NULL 
  if (!is.null(id)) if (is.na(id)) id=NULL;
  paste(collapse='.',c(paste(sep='','m=',m.pretty(m)),paste(sep='','n=',n),id));
}
## pretty print typical values of m
m.pretty=function(m) sub('e\\+0{0,1}','e',sprintf("%0.0e",m),perl=TRUE)

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
    suffix=sub('^\\.','',suffix,perl=TRUE);
    suffix.pattern=paste(collapse='|',paste(sep='','\\.',suffix,'$'));
    file=ifelse(grepl(suffix.pattern,file),file,paste(sep='.',file,suffix[1]));
  }
  file;
}
## remove suffix from filename
desuffix=function(file,suffix=c('RData','txt')) {
  if (!is.null(suffix)) {
    ## remove leading '.' if present
    suffix=sub('^\\.','',suffix,perl=TRUE);
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
## misc functions
##
## get value of variable from parent or set to default
## call with quoted or unquoted variable name
## if default missing, throws error if variable not found
parent=function(what,default) {
  what=as.character(pryr::subs(what));
  if (exists(what,envir=parent.frame(n=2))) return(get(what,envir=parent.frame(n=2)));
  if (!missing(default)) return(default);
  stop(paste(sep='',"object '",what,"' not found"));
}
## not in - based on example in RefMan - more intutive than !%in%
"%notin%"=function(x,table) match(x,table,nomatch=0)==0
## between, near - to subset sim results
between=function(x,lo,hi) x>=lo&x<hi
near=function(x,target,tol=.01) between(x,target-tol,target+tol)

