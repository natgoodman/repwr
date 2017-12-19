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
  small.power=0.33,              # min small telescope power
  close.level=0.05,              # pval for big telescope "close enough" calculation
  n1.dcrit=n,                    # n1 for d.crit data frame (do_dcrit, plot_dcrit)
  n2.dcrit=10:1000,              # n2 for d.crit data frame (do_dcrit, plot_dcrit)
  n.dvsd=n[1],                   # n for d.sdz vs. d.pop plot (plot_dvsd)
  n1.bound=n[1],                 # n1 for "bounds" plot (plot_bounds)
  n2.bound=2.5*n1.bound,         # n2 for "bounds" plot (plot_bounds)
  ## conf.level=0.95,            # for confidence intervals
  sig.level=0.05,
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
  save=F,                        # shorthand for save.sim & save.res 
  save.sim=save,                 # save simulations (RData format)
  save.res=save,                 # save results (RData & txt formats); plots (png)
  save.plot=save.res,            # save plots. used to override save.res
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
  ## at end, assign global parameters to global variables
  ##   don't do it earlier because it's too easy to confuse local and
  ##   global variables with the same names
  m<<-m;
  n<<-n;
  d.pop<<-d.pop;
  small.power<<-small.power;
  close.level<<-close.level;
  n1.dcrit<<-n1.dcrit;
  n2.dcrit<<-n2.dcrit;
  n.dvsd<<-n.dvsd;
  n1.bound<<-n1.bound;
  n1.bound<<-n1.bound;
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
dores=function(save.res=parent(save.res)) {
  fig<<-1;                    # figure number. bumped in each plot function
  d.crit=do_dcrit();          # get table of critical d values - significance & close enough
  plot.dvsd=plot_dvsd();      # plot d.pop vs. d.sdz
  plot.bounds=plot_bounds();  # plot bounds for one case of n1, n2
  plot.dcrit=plot_dcrit();    # plot critical d values for fixed n1
  if (save.res) {
    save_data(d.crit,'dcrit');
    save_plot(plot.dvsd,'plot_dvsd');
    save_plot(plot.bounds,'plot_bounds');
    save_plot(plot.dcrit,'plot_dcrit');
  }
}
## get table of critical d values - significance & close enough
do_dcrit=function(n1=parent(n1.dcrit,20),n2=parent(n2.dcrit,10:1000)) {
  d.crit=do.call(rbind,
    lapply(n1,function(n1) {
      d.crit=do.call(rbind,
        lapply(n2,function(n2) data.frame(n2,d.sig=d_sig(n2),d.close=d_close(n2,n1))));
      d.crit=data.frame(n1,d.crit);
    }))
}
## plot d.pop vs. d.sdz - 2 zoom levels
plot_dvsd=function(fig=parent(fig,1),n=parent(n.dvsd,20)) {
  ## entire range. limit to (-2,2) for symmetry. sample down to 10k
  dev.new();
  title=paste(sep='','Figure ',fig,". Cohen's d vs. true population d for n=",n);
  case=paste(m,n);
  sim2=sim[[case]];
  if (is.null(sim2)) stop(paste(sep='','simulation results for n=',n,' m=',m,' not in memory'));
  ## sim2=subset(sim2,subset=d.sdz>=0);
  sim2=subset(sim2,subset=between(d.sdz,-1,2));
  sim2=sim2[sample.int(nrow(sim2),1e4),];
  with(sim2,plot(d.sdz,d.pop,pch=20,col=pval2col(pval),xlim=c(-1,2),
                 xlab="Cohen's d",ylab='true population d',main=title));
  grid();
  pval_legend(-1,0.5,width=.25,height=.95);
  dev1=dev.cur();
  ## zoom to sig boundary
  dev.new();
  title=paste(sep='','Figure ',fig+1,". Cohen's d vs. true d for n=",n,' (zoomed to significance region)');
  sim2=sim[[case]];
  d.sig=d_sig(n);
  sim2=subset(sim2,subset=near(d.sdz,d.sig,tol=.01));
  with(sim2,plot(d.sdz,d.pop,pch=20,col=pval2col(pval),
                 xlab="Cohen's d",ylab='true population d',main=title,cex.main=.9));
  grid();
  abline(v=d.sig,col='red');
  fig<<-fig+2;                # bump for next guy
  dev2=dev.cur();
  c(dev1,dev2);
}
## plot bounds for one case of n1, n2
plot_bounds=function(fig=parent(fig,3),n1=parent(n1.bound,20),n2=parent(n2.bound,50)) {
  dev.new();
  title=paste(sep='','Figure ',fig,'. d.sig and d.close bounds for n1=',n1,', n2=',n2);
  d.33=d_33(n1);
  x=seq(-1,2,length=1e4);
  ymax=dd2t(n2,d.33,d0=d.33);
  plot(x=NULL,y=NULL,type='n',xlim=range(x),ylim=c(0,ymax),xlab='d',ylab='density',main=title);
  grid();
  seglines(x,dd2t(n1,x,d0=d.33),col=pval2col(d2pval(n1,x)),lwd=4);
  seglines(x,dd2t(n2,x,d0=d.33),col=pval2col(d2pval(n2,x)),lwd=4);
  ## abline(v=d_sig(20),col='red',lty='dashed');
  dsig=d_sig(n2);
  dclose=d_close(n2,n1)
  abline(v=dsig,col='red',lty='dashed');
  abline(v=dclose,col='grey',lty='dashed',lwd=2);
  ## thanks to https://www.datacamp.com/community/tutorials/15-questions-about-r-plots#q3
  ## for axis call below!
  at=c(dsig,dclose);
  axis(side=1,at=at,labels=round(at,2),tck=-.02,cex.axis=.75,mgp=c(3,.4,0));
  y.text=0.15;
  x.text=suppressWarnings(uniroot(function(d) dd2t(n1,d,d0=d.33)-y.text,c(d.33,2))$root)
  text(x.text+.02,y.text,"sampling distribution for n=20\ncentered on d.33",adj=c(0,0.5),cex=.75);
  y.text=0.05;
  x.text=suppressWarnings(uniroot(function(d) dd2t(n2,d,d0=d.33)-y.text,c(d.33,2))$root)
  text(x.text+.02,y.text,"sampling distribution for n=50\ncentered on d.33",adj=c(0,0.5),cex=.75);
  text(d_close(n2,n1)-.01,0.2,"replication\n'close enough'\nto d.33%",adj=c(0.5,0.5),cex=.75);
  text(d_sig(n2)-.01,0.1,"replication\nsignificant",adj=c(0.5,0.5),cex=.75);
  pval_legend(1.25,0.2,height=0.15);
  fig<<-fig+1;                # bump for next guy
  dev.cur();
}
## plot critical d values for fixed n1
plot_dcrit=
  function(fig=parent(fig,4),d.crit=parent(d.crit,NULL),n1=parent(n1.dcrit,20),
           n2=parent(n2.bound,50)) {
    ## make d.crit if needed (should have been made in calling function)
    d.crit=do_dcrit();        # get table of critical d values - significance & close enough
    dev.new();
    title=paste(sep='','Figure ',fig,'. d.sig and d.close bounds for n1=',n1,', varying n2');
    dcrit=d.crit[d.crit$n1==n1,];
    x=subset(dcrit,select=n2);
    y=subset(dcrit,select=c(d.sig,d.close))
    col=c('red','blue');
    matplot(x=x,y=y,type='l',xlab='replication n',ylab='replication d',main=title,
            lty='solid',col=col);
    grid();
    ## calculate crossover point
    n2.cross=
      dcrit$n2[uniroot(function(i) dcrit$d.sig[i]-dcrit$d.close[i],
                       interval=c(1,nrow(dcrit)))$root];
    abline(v=n2.cross,col='grey',lty='dashed');
    ## thanks to https://www.datacamp.com/community/tutorials/15-questions-about-r-plots#q3
    ## for axis call below!
    axis(side=1,at=n2.cross,labels=n2.cross,tck=-.02);
    ## add legend
    legend.title='';
    legend.text=c('d.sig','d.close');
    legend.col=col;
    legend('topright',bty='n',legend=legend.text,col=legend.col,pch=19,cex=1,
           title=legend.title,title.col='black')
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
## t statistic to pval
t2pval=function(n,t) 2*pt(-abs(t),df=2*(n-1))
## Cohen's d to t statistic & pval
d2t=function(n,d) d*sqrt(n^2/(2*n))
d2pval=function(n,d) t2pval(n,d2t(n,d))

## probability functions for t-distribution of d
ncp=function(n,d) sqrt(n/2)*d
pd2t=function(n,d,d0=NULL) {
  df=2*(n-1);
  t=d2t(n,d);
  if (!is.null(d0)) pt(t,df=df,ncp=ncp(n,d0)) else pt(t,df=df)
}
dd2t=function(n,d,d0=NULL) {
  df=2*(n-1);
  t=d2t(n,d);
  if (!is.null(d0)) dt(t,df=df,ncp=ncp(n,d0)) else dt(t,df=df)
}

## sd of Cohen's d for two samples
## from https://stats.stackexchange.com/questions/144084/variance-of-cohens-d-statistic
## also in compute.es::des
sdd=function(n,d) sqrt(((2*n)/n^2)+(d^2/(4*n)))

## significance boundary for Cohen's d
d_sig=function(n,sig.level=parent(sig.level,0.05))
  uniroot(function(d) d2pval(n,d)-sig.level,c(0,10))$root
## Uri's d.33% - value of d giving power of 33%
d_33=function(n=20) power.t.test(n,power=0.33)$delta
## "close enough to d.33" boundary
d_close=function(n2,n1=20,
  small.power=parent(small.power,0.33),close.level=parent(close.level,0.05),
  d.33=power.t.test(n1,power=small.power)$delta)
  ## suppress annoying warning from pt: full precision may not have been achieved in 'pnt{final}'
  suppressWarnings(uniroot(function(d) pd2t(n2,d,d0=d.33)-close.level,c(-10,10))$root)

## ---- Utility Functions ----
## file related functions
## filename for sim files
filename_sim=function(m,n,id=NULL,suffix='RData') 
  filename(basename_sim(m,n,id=id),suffix=suffix);
basename_sim=function(m,n,id=NULL) {
  filename(datadir,base='sim2',tail=
           c(paste(sep='','m=',sub('e\\+0{0,1}','e',sprintf("%0.0e",m),perl=T)),
             paste(sep='','n=',n),id));
}
  
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

