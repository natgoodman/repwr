#################################################################################
##
## Author:  Nat Goodman
## Created: 17-12-04
##          from repwr.R created 17-10-05
##          from swfdr.R created 17-10-02
##          from swfdr_base.R created 17-04-10
##          from fdr.R created 17-04-03
##          from pval_error.R created 16-11-01
##
## Copyright (C) 2017 Nat Goodman.
## 
## Explores Uri Simonsohn's Small Telescopes replication scheme.
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
  m=1e5,                         # total number of instances
  n=20,                          # sample size 
  d.pop=runif(m,0,1.5),          # population effect sizes (means)
  ## analysis parameters
  scope.power=0.33,              # min scope telescope power
  close.level=0.05,              # pval for big telescope "close enough" calculation
  n1.res=c(5,10,20*(1:10)),      # n1 for main results
  n2.res=NULL,                   # n2 for main results. by default, computed from dcross$n2
  n2max.res=1e4,                 # max n2 for uniroot intervals in main results
  pval.res=seq(0.01,0.05,by=.01),# pval for main results
  m.dvsd=m[1],                   # m for d.sdz vs. d.pop plot (plot_dvsd)
  n.dvsd=n[1],                   # n for d.sdz vs. d.pop plot (plot_dvsd)
  id.dvsd=id[1],                 # id for d.sdz vs. d.pop plot (plot_dvsd)
  n1.bound=n[1],                 # n1 for plot_bounds
  n2.bound=2.5*n1.bound,         # n2 for plot_bounds
  n1.dcrit1=n[1],                # n1 for plot_dcrit1
  n2.dcrit1=2.5*n1.dcrit1,       # n2 for plot_dcrit1
  pval.dcrit1=0.05,              # pval for plot_dcrit1
  n1.pclose1=n[1],               # n1 for plot_pclose1
  pval.pclose1=0.05,             # pval for plot_pclose1
  dsmall.dcross=0.2,             # small effect size for plot_dcross
  pval.dcross=0.05,              # pval for plot_dcross
  ## conf.level=0.95,            # for confidence intervals
  sig.level=0.05,                # for conventional significance
  ## pval.plot=c(.001,.01,.03,.05,.1), # pvalues for which we plot results
  ## 
  ## program parameters, eg, for output files, error messages, etc.
  scriptname='scope',            #  
  datadir=file.path('data',scriptname), # directory for data files
  figdir=file.path('figure',scriptname), # directory for plots
  id=NULL,                       # identifying info tacked onto filenames
  verbose=F,                     # print progress messages
  ## program control
  load=NA,                       # shorthand to override defaults for load.sim & load.res
  load.sim=load,                 # load saved simulations
                                 #   NA means load if file exists
                                 #   T, F mean always or never load
  save=F,                        # shorthand for save.sim & save.data 
  save.sim=save,                 # save simulations (RData format)
  save.data=save,                # save data results (RData & txt formats)
  save.plot=save,                # save plots
  save.simtxt=F,                 # save simulations in txt format - big & slow!
  clean=F,                       # remove contents of datadir and figdir and start fresh
  clean.data=clean,              # remove contents of datadir and start fresh
  clean.fig=clean,               # remove contents of figdir and start fresh
  end=NULL                       # placeholder for last parameter
  ) {
  cases=expand.grid(m=m,n=n,load.sim=load.sim,stringsAsFactors=F);
  ## clean and create output directories as needed
  if (clean.data) unlink(datadir,recursive=T);
  dir.create(datadir,recursive=TRUE,showWarnings=FALSE);
  if (clean.fig) unlink(figdir,recursive=T);
  dir.create(figdir,recursive=TRUE,showWarnings=FALSE);
  ## setup 'sim' list to hold simulations. do it carefully in case already setup
  if (!exists('sim')) sim=list() else if (!is.list(sim)) sim=list();
  ## at end, assign global parameters to global variables
  ##   don't do it earlier because it's too easy to confuse local and
  ##   global variables with the same names
  m<<-m;
  n<<-n;
  sim<<-sim;
  d.pop<<-d.pop;
  scope.power<<-scope.power;
  close.level<<-close.level;
  n1.res<<-n1.res;
  n2.res<<-n2.res;
  pval.res<<-pval.res;
  m.dvsd<<-m.dvsd;
  n.dvsd<<-n.dvsd;
  id.dvsd<<-id.dvsd;
  n1.bound<<-n1.bound;
  n2.bound<<-n2.bound;
  n1.dcrit1<<-n1.dcrit1;
  n2.dcrit1<<-n2.dcrit1;
  pval.dcrit1<<-pval.dcrit1;
  dsmall.dcross<<-dsmall.dcross;
  pval.dcross<<-pval.dcross;
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
  save.data<<-save.data;
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
  sim=list();
  ## do simulations! case by case
  for (i in seq_len(nrow(cases))) {
    ## extract args this way to preserve type (numeric vs string) and avoid conversion to list
    m=cases[i,'m']; n=cases[i,'n']; load.sim=cases[i,'load.sim']; 
    case=paste(m,n);
    file=filename_sim(m,n,id);
     ## use saved simulation if possible
    if (load.sim&!file.exists(file))
      stop(paste(sep=' ','load.sim is TRUE but file',file,'does not exist'));
    if ((is.na(load.sim)|load.sim)&file.exists(file)) {
      sim[[case]]=load_sim(file);
      next;
    }
    ## no saved simulation or params say not to use it. run simulation
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
    pval=t.test(group0,group1,var.equal=T)$p.value;
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
dores=function(save.data=parent(save.data)) {
  fig<<-1;                  # figure number. bumped in each plot function
  dodata(do_d1);            # get d1 for each pval
  dodata(do_d33);           # get d33 for each d1
  dodata(do_dcross);        # get n2,d2 at crossover for each d33. also sets n2.res
  dodata(do_dsig);          # get conventional significance boundaries for each n2
  dodata(do_dclose);        # get "close enough" boundaries for each n2,d33
  dodata(do_dcrit);         # combine dsig,dclose in single table
  plots=c(
    doplot(plot_dvsd),      # plot d.pop vs. d.sdz
    doplot(plot_bounds),    # plot bounds for one case of n1, n2
    doplot(plot_dcrit1),    # plot critical d values for fixed n1
    doplot(plot_pclose1),   # plot 'close enough' pvals d values for fixed n1
    doplot(plot_dcross));   # plot min rep d, aka crossover points
  }
  ## run plot function, save if required, label result with function name
doplot=function(what) {
  what=as.character(pryr::subs(what));
  f=get(what,envir=parent.frame(n=2),mode='function');
  dev=f();
  if (save.plot) save_plot(dev,what);
  setNames(dev,what)
}
## run data-making function, save if required, store result in global workspace
dodata=function(what) {
  what=as.character(pryr::subs(what));
  f=get(what,envir=parent.frame(n=2),mode='function');
  data=f();
  if (save.data) save_data(data,what);
  ## fiddle with name to get the form we want
  what=sub('^do_','',what);
  assign(what,data,envir=.GlobalEnv);
  invisible(data);
}
## get table of d1 for each pval
do_d1=function(n1=parent(n1.res,20),pval=parent(pval.res,seq(0.01,0.05,by=.01))) {
  d1=do.call(rbind,lapply(n1,function(n1) {
    d1=pval2d(n1,pval);
    data.frame(n1,pval,d1)}));
    d1;
  }
## get table of d33.
do_d33=function(d1=parent(d1,NULL)) {
  if (is.null(d1))
    stop('d1 must be set (typically by running d1=do_d1()) before calling do_d33');
  d33=do.call(rbind,apply(d1,1,function(row) {
    n1=row['n1']; pval=row['pval']; d1=row['d1'];
    d33=uniroot(function(d0) p_d2t(n1,d=d1,d0,lower.tail=F)-0.33,interval=c(0,10))$root;
    data.frame(n1,pval,d33)}));
  d33;
}
## get table of n2,d2 at crossover point
do_dcross=function(d33=parent(d33,NULL),n2.max=parent(n2max.res,1e4)) {
  if (is.null(d33))
    stop('d33 must be set (typically by running d33=do_d33()) before calling do_dcross');
  dcross=do.call(rbind,apply(d33,1,function(row) {
    n1=row['n1']; pval=row['pval']; d33=row['d33'];
    n2=round(uniroot(function(n2) d_sig(n2)-d_close(n2,d33=d33),interval=c(5,n2.max))$root);
    d2=mean(d_sig(n2),d_close(n2,d33=d33));
    ## data.frame(n1,pval,n2.cross=round(n2.cross))}));
    data.frame(n1,n2,pval,dcross=d2)}));
  rownames(dcross)=1:nrow(dcross);
  n2.res<<-seq(5,max(dcross$n2)*1.5,by=5); # add 50% so plots will look nice
  dcross;
}
## get table of conventionally significant d values
do_dsig=
  function(n2=parent(n2.res,seq(5,1680,5)),sig.level=parent(sig.level,0.05)) 
    data.frame(n2,dsig=pval2d(n2,sig.level))  

## get table of 'close enough' d values
do_dclose=function(n2=parent(n2.res,seq(5,1680,5)),d33=parent(d33,NULL)) {
  if (is.null(d33))
    stop('d33 must be set (typically by running d33=do_d33()) before calling do_dclose');
  dclose=do.call(rbind,apply(d33,1,function(row) {
    n1=row['n1']; pval=row['pval']; d33=row['d33'];
    dclose=d_close(n2,d33);
    data.frame(n1,n2,pval,dclose,row.names=NULL)}));
  rownames(dclose)=1:nrow(dclose);
  dclose;
}
## get table of critical d values - significance & close enough - by merging dsig, dclose
do_dcrit=function(dsig=parent(dsig,NULL),dclose=parent(dclose,NULL)) {
  if (is.null(dsig))
    stop('dsig must be set (typically by running dsig=do_dsig()) before calling do_dcrit');
  if (is.null(dclose))
    stop('dclose must be set (typically by running dclose=do_dclose()) before calling do_dcrit');
  dcrit=merge(dsig,dclose,by='n2');
  ## get rows and columns in prefered order
  dcrit=subset(dcrit,select=c(n1,n2,pval,dsig,dclose));
  dcrit=dcrit[with(dcrit,order(n1,n2,pval)),];
  dcrit;
}

## plot d.pop vs. d.sdz - 2 zoom levels
plot_dvsd=
  function(fig=parent(fig,1),
           m=parent(m.dvsd,1e5),n=parent(n.dvsd,20),id=parent(id.dvsd,NULL)) {
    ## entire range. limit to (-2,2) for symmetry. sample down to 10k
    dev.new();
    title=paste(sep='','Figure ',fig,". Cohen's d vs. true population d for n=",n);
    sim2=get_sim(m,n,id);
    ## sim2=subset(sim2,subset=d.sdz>=0);
    sim2=subset(sim2,subset=between(d.sdz,-1,2));
    sim2=sim2[sample.int(nrow(sim2),1e4),];
    with(sim2,plot(d.sdz,d.pop,pch=16,cex=.6,col=pval2col(pval),xlim=c(-1,2),
                   xlab="Cohen's d",ylab='true population d',main=title));
    grid();
    pval_legend(-1,0.5,width=.25,height=.95);
    dev1=dev.cur();
    ## zoom to sig boundary
    dev.new();
    title=paste(sep='','Figure ',fig+1,". Cohen's d vs. true d for n=",n,' (zoomed to significance region)');
    sim2=get_sim(m,n,id);
    d.sig=d_sig(n);
    sim2=subset(sim2,subset=near(d.sdz,d.sig,tol=.005));
    with(sim2,plot(d.sdz,d.pop,pch=16,cex=.6,col=pval2col(pval),
                   xlab="Cohen's d",ylab='true population d',main=title,cex.main=.9));
    grid();
    ## because of roundoff error or something, d.sig is not the exact boundary
    ## compute boundary from data so line will fall between the reds and greens
    max.neg=max(sim2$d.sdz[sim2$pval>0.05]);
    min.pos=min(sim2$d.sdz[sim2$pval<=0.05]);
    abline(v=mean(max.neg,min.pos),col='red');
    fig<<-fig+2;                # bump for next guy
    dev2=dev.cur();
    c(dev1,dev2);
  }
## plot bounds for one case of n1, n2
plot_bounds=function(fig=parent(fig,3),n1=parent(n1.bound,20),n2=parent(n2.bound,50)) {
  dev.new();
  title=paste(sep='','Figure ',fig,'. dsig and dclose bounds for n1=',n1,', n2=',n2);
  d33=d_33(n1);
  x=seq(-1,2,length=1e4);
  ymax=d_d2t(n2,d33,d0=d33);
  plot(x=NULL,y=NULL,type='n',xlim=range(x),ylim=c(0,ymax),xlab='d',ylab='density',main=title);
  grid();
  seglines(x,d_d2t(n1,x,d0=d33),col=pval2col(d2pval(n1,x)),lwd=4);
  seglines(x,d_d2t(n2,x,d0=d33),col=pval2col(d2pval(n2,x)),lwd=4);
  ## abline(v=d_sig(20),col='red',lty='dashed');
  dsig=d_sig(n2);
  dclose=d_close(n2,d33)
  abline(v=dsig,col='red',lty='dashed');
  abline(v=dclose,col='grey',lty='dashed',lwd=2);
  ## thanks to https://www.datacamp.com/community/tutorials/15-questions-about-r-plots#q3
  ## for axis call below!
  at=c(dsig,dclose);
  axis(side=1,at=at,labels=round(at,2),tck=-.02,cex.axis=.75,mgp=c(3,.4,0));
  y.text=0.15;
  x.text=suppressWarnings(uniroot(function(d) d_d2t(n1,d,d0=d33)-y.text,c(d33,2))$root)
  text(x.text+.02,y.text,"sampling distribution for n=20\ncentered on d33",adj=c(0,0.5),cex=.75);
  y.text=0.05;
  x.text=suppressWarnings(uniroot(function(d) d_d2t(n2,d,d0=d33)-y.text,c(d33,2))$root)
  text(x.text+.02,y.text,"sampling distribution for n=50\ncentered on d33",adj=c(0,0.5),cex=.75);
  text(d_close(n2,d33)-.01,0.2,"replication\n'close enough'\nto d33%",adj=c(0.5,0.5),cex=.75);
  text(d_sig(n2)-.01,0.1,"replication\nsignificant",adj=c(0.5,0.5),cex=.75);
  pval_legend(1.25,0.2,height=0.15);
  fig<<-fig+1;                # bump for next guy
  dev.cur();
}
## plot critical d values for fixed n1
plot_dcrit1=
  function(fig=parent(fig,4),dcrit=parent(dcrit,NULL),
           n1=parent(n1.dcrit1,20),n2=parent(n2.dcrit1,50),pval=parent(pval.dcrit1,0.05)) {
    ## make dcrit if needed (should have been made in calling function)
    if (is.null(dcrit)) dcrit=do_dcrit();
    dev.new();
    title=paste(sep='','Figure ',fig,'. dsig and dclose for n1=',n1,', varying n2, pval1=',pval);
    dcrit=dcrit[dcrit$n1==n1&dcrit$pval==pval,];
    x=subset(dcrit,select=n2);
    y=subset(dcrit,select=c(dsig,dclose))
    ## interpolate dcrit to calculate crossover point
    n2_dsig=with(dcrit,approxfun(n2,dsig));
    n2_dclose=with(dcrit,approxfun(n2,dclose));
    n2=round(uniroot(function(n2) n2_dsig(n2)-n2_dclose(n2),interval=range(x))$root);
    ## set xlim to reasonable value > n2
    xlim=c(0,max(n2*2,round(n2*2,-2)));
    col=c('red','blue');
    matplot(x=x,y=y,type='l',xlim=xlim,xlab='replication n',ylab='replication d',main=title,
            lty='solid',col=col);
    grid();
    abline(v=n2,col='grey',lty='dashed');
    ## thanks to https://www.datacamp.com/community/tutorials/15-questions-about-r-plots#q3
    ## for axis call below!
    axis(side=1,at=n2,labels=n2,tck=-.02,cex.axis=.75,mgp=c(3,.4,0));
    ## add legend
    legend.text=c('dsig','dclose');
    legend.col=col;
    legend('topright',bty='n',legend=legend.text,col=legend.col,pch=19,cex=1);
    fig<<-fig+1;                # bump for next guy
    dev.cur();
  }
## plot pvals for 'close enough' d values for fixed n1
plot_pclose1=
  function(fig=parent(fig,5),dcrit=parent(dcrit,NULL),
           n1=parent(n1.pclose1,20),pval=parent(pval.pclose1,0.05)) {
    ## make dcrit if needed (should have been made in calling function)
    if (is.null(dcrit)) dcrit=do_dcrit();
    dev.new();
    title=paste(
      sep='','Figure ',fig,'. pval for "close enough" d for n1=',n1,', varying n2, pval1=',pval);
    dcrit=dcrit[dcrit$n1==n1&dcrit$pval==pval,];
    ## interpolate dcrit to calculate crossover point
    n2_dsig=with(dcrit,approxfun(n2,dsig));
    n2_dclose=with(dcrit,approxfun(n2,dclose));
    ## set n2.max to reasonable value > crossover point
    n2.cross=round(uniroot(function(n2) n2_dsig(n2)-n2_dclose(n2),interval=range(dcrit$n2))$root);
    n2.max=round(n2.cross*1.5,-1);
    n2=5:n2.max;
    ## calculate pval for each n2 and interpolate for nicer plot
    dcrit$pval2=with(dcrit,d2pval(n2,dclose));
    n2_pval2=with(dcrit,splinefun(n2,pval2));
    plot(n2,n2_pval2(n2),type='l',xlab='replication n',ylab='replication pval',
         main=title,cex.main=.9,col='blue');
    grid();
    abline(h=.05,col='darkgreen',lty='dashed');
    axis(side=2,at=.05,labels=.05,tck=-.02,
         col.axis='darkgreen',col.tick='darkgreen',cex.axis=.75,mgp=c(3,.4,0));
    ## space out points so text won't be too crowded
    n2=c(seq(10,round(n2.cross,-1),by=10),seq(round(n2.cross,-1)+20,n2.max,by=20));
    points(n2,n2_pval2(n2),col='red',pch=16);
    text(n2+1,n2_pval2(n2)+.015,round(n2_dclose(n2),2),col='red',cex=.75,adj=c(0,.5));
    ## add legend
    legend.text=c('"close enough" d')
    legend.col=c('red');
    legend('topright',bty='n',legend=legend.text,col=legend.col,pch=19,cex=1);
    fig<<-fig+1;                # bump for next guy
    dev.cur();
  }
## plot min replication d for varying n1, aka crossover points
plot_dcross=
  function(fig=parent(fig,6),dcross=parent(dcross,NULL),
           dsmall=parent(dsmall.dcross,0.2),pval=parent(pval.dcross,0.05)) {
    ## make dcross if needed (should have been made in calling function)
    if (is.null(dcross)) dcross=do_dcross();
    dev.new();
    title=paste(sep='','Figure ',fig,'. Minimum replication d for varying n1, pval1=',pval);
    dcross=dcross[dcross$pval==pval,];
    n1_d=splinefun(dcross$n1,dcross$dcross);
    n1_n2=splinefun(dcross$n1,dcross$n2);
    n2_d=splinefun(dcross$n2,dcross$dcross);
    n1=min(dcross$n1):max(dcross$n1);
    plot(n1,n1_d(n1),type='l',ylim=c(0,max(dcross$dcross)),
         xlab='n1',ylab='d',main=title,col='blue');
   ## space out points so text won't be too crowded
    n1=sort(unique(dcross$n1));
    points(n1,n1_d(n1),col='red',pch=16);
    text(n1+3,n1_d(n1)+.015,n1_n2(n1),col='red',cex=.75,adj=c(.5,.5));
    grid();
    ## do legend here 'cuz I reuse dsmall in later code
    legend.text=c('d2 (min replication d)','n2 to achieve d2',
                  paste(sep='','min n1,n2 to achieve d2=',dsmall),
                  paste(sep='','min n1,n2 to achieve 80% power to detect d2=',dsmall));
    legend.col=c('blue','red','darkgreen','darkolivegreen4');
    legend('topright',bty='n',legend=legend.text,col=legend.col,pch=19,cex=.9)
    
    ## n2 needed to achieve d=0.2 (193)
    abline(h=dsmall,col='darkgreen',lty='dashed')
    n2.dsmall=round(uniroot(function(n2) n2_d(n2)-dsmall,range(dcrit$n2))$root)
    text(100,dsmall+.005,n2.dsmall,col='darkgreen',cex=.75,adj=c(.5,0))
    ## n1 where dcross crosses 0.2
    n1.dsmall=round(uniroot(function(n1) n1_d(n1)-dsmall,range(dcrit$n1))$root)
    abline(v=n1.dsmall,col='darkgreen',lty='dashed')
    at=n1.dsmall;
    axis(side=1,at=at,labels=at,tck=-.02,col.axis='darkgreen',col.tick='darkgreen',
         cex.axis=.75,mgp=c(3,.4,0));
    ## n2 needed for 80% power to detect 0.2 (393) and corresponding dsig (.14)
    n2.dsmall=round(power.t.test(delta=dsmall,power=0.8)$n);
    dsmall=d_sig(n2.dsmall);
    abline(h=dsmall,col='darkolivegreen4',lty='dashed')
    text(100,dsmall+.005,n2.dsmall,col='darkolivegreen4',cex=.75,adj=c(.5,0))
    at=dsmall;
    axis(side=2,at=at,labels=round(at,2),tck=-.02,
         col.axis='darkolivegreen4',col.tick='darkolivegreen4',cex.axis=.75,mgp=c(3,.4,0));
    ## n1 where dcross crosses dsmall
    n1.dsmall=round(uniroot(function(n1) n1_d(n1)-dsmall,range(dcrit$n1))$root)
    abline(v=n1.dsmall,col='darkolivegreen4',lty='dashed')
    at=n1.dsmall;
    axis(side=1,at=at,labels=at,tck=-.02,col.axis='darkolivegreen4',col.tick='darkolivegreen4',
         cex.axis=.75,mgp=c(3,.4,0));
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
  base=filebasename(figdir,filebase,tail=id);
  if (length(dev)==1) savePlot(filename(figdir,filebase,tail=id,suffix='png'),device=dev)
    else for (i in seq_along(dev)) {
    file=filename(figdir,paste(sep='.',filebase,sprintf("%02i",i)),tail=id,suffix='png');
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

## significance boundary for Cohen's d
d_sig=function(n,sig.level=parent(sig.level,0.05)) pval2d(n,pval=sig.level)
## Uri's d.33% - value of d giving power of 33%
##  specify either pval or d1
## d_33=function(n=20) power.t.test(n,power=0.33)$delta
d_33=function(n1=20,pval=.05,d1=d_sig(n1,pval)) 
  uniroot(function(d0) p_d2t(n1,d=d1,d0,lower.tail=F)-0.33,interval=c(0,10))$root
## generalization of d_33 for arbitrary 'small scope' power
##  specify either pval or d1
d_scope=
  function(n1=20,pval=.05,d1=d_sig(n1,pval),scope.power=parent(scope.power,0.33)) 
    uniroot(function(d0) p_d2t(n1,d=d1,d0,lower.tail=F)-scope.power,interval=c(0,10))$root

## "close enough to d33" boundary.
d_close1=function(n2,d33,close.level=parent(close.level,0.05))
  ## suppress annoying warning from pt: full precision may not have been achieved in 'pnt{final}'
  suppressWarnings(uniroot(function(d) p_d2t(n2,d,d0=d33)-close.level,c(-10,10))$root)
d_close=Vectorize(d_close1,vectorize.args='n2')

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

