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
## Initialization code for repwr.R
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################

## ---- init ----
## initialization.
## process parameters and store in global variables.
## create output directory if necessary.
init=function(
  ## simulation parameters 
  n=10*2^(0:9),                     # sample sizes
  d=round(seq(0,1,by=0.1),digits=5),# population effect sizes (means)
                                    # round to avoid imprecise decimals
  m=1e3,                            # number of instances per pop
  ## analysis parameters
  sig.level=0.05,                   # for conventional significance
  conf.level=0.95,                  # for confidence intervals
  pred.level=0.95,                  # for prediction intervals
  pwr.level=0.80,                   # default power
  fpr.cutoff=sig.level,             # false positive rate cutoff for plots
  fnr.cutoff=1-pwr.level,           # false negative rate cutoff for plots
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
  scriptname='repwr',                       #
  data.subdir=paste_nv(m,m_pretty(m)),      # vector of subdir components
  ## datadir='data',                        # directory for data files
  ## datadir=file.path('data',scriptname),  # directory for data files
  datadir=filename('data',tail=data.subdir),# directory for data files. default eg, data/m=1e4
  simdir=file.path(datadir,'sim'),          # directory for sim files
  simrdir=file.path(datadir,'simr'),        # directory for simr files
  i1dir=file.path(datadir,'i1'),            # directory for s1 permutation index files
  detldir=file.path(datadir,'detl'),        # directory for detailed results files
  smrydir=file.path(datadir,'smry'),        # directory for summary results files
  figdir=file.path('figure',scriptname),    # directory for plots
  id=NULL,                                  # info tacked onto filenames. not used much
  verbose=F,                                # print progress messages
  ## program control
  must.exist=F,                  # must all sub-inits succeed?
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
  init_smry(must.exist=must.exist);
  init_mesr(must.exist=must.exist);
  invisible();
}
## initialize summary types. NOTE: needs smry object!
init_smry=function(
  must.exist=T,
  smry=get_data(smry,must.exist=F),
  smry.type=unique(smry$type),
  end=NULL                       # placeholder for last parameter
  ) {
  ## return immediately if already initialized
  if (exists('init.smry',envir=.GlobalEnv)&&init.smry) invisible(T);
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
  mesr.relsig1=c(mesr.frombsln,mesr.fromsig2),
  mesr.relsig2=NULL,
  ## mesr.relsig2=mesr.fromsig2,
  ## measures grouped by functional category
  mesr.sig=cq(sig2,sigm),
  ## mesr.dcc=cq(d1.c2,d2.c1,c1.c2,d1.p2,d2.p1,p1.p2),
  mesr.dcc=grep('(d(1|2)\\.(c|p)(1|2))|((c|p)1\\.(c|p)2)',mesr.all,value=T),
  ## mesr.scp=cq(d1.scp2,d2.scp1,d1.scpd2,d2.scpd1),
  mesr.scp=grep('d(1|2)\\.scp(d{0,1})(1|2)',mesr.all,value=T),
  mesr.meta=grep('dm|cm',mesr.all,value=T),
  mesr.other=setdiff(mesr.all,c(mesr.sig,mesr.dcc,mesr.scp,mesr.meta)),
  ## defaults for plotrate and heatrate
  ## mesr.plotdflt=cq(sig2,sigm,d1.c2,d2.c1,d1.p2,d2.p1),
  mesr.plotdflt=c(mesr.sig,mesr.dcc,'d1.scp2','d2.scp1'),
  mesr.heatdflt=c(mesr.sig,mesr.dcc,mesr.scp,mesr.meta),
  ## measures ordered for legends and such
  mesr.order=c(mesr.sig,mesr.dcc,mesr.scp,mesr.other),
  ## RColorBrewer pallete names for plotrate
  pal.sig='Reds',
  pal.dcc='Blues',
  pal.meta='Greens',
  pal.scp='Purples',
  pal.other='Greys',
  ## line properties
  lty.max=8,                          # max 'on' value for dashed lines
  end=NULL                       # placeholder for last parameter
  ) {
 ## return immediately if already initialized
  if (exists('init.mesr',envir=.GlobalEnv)&&init.mesr) invisible(T);
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
    ## col=colorRampPalette(rev(RColorBrewer::brewer.pal(9,pal)[3:8]))(length(mesr));
    ## col=colorRampPalette(rev(RColorBrewer::brewer.pal(9,pal))[2:7])(length(mesr));
    ## RColorBrewer sequential palettes can have 3-9 colors. use directly unless too many mesrs
    ## skip 1st two colors - usually too light - then reverse so darker colors will be first
    ##   if too many mesrs, colorRampPalette will make more
    n.mesr=length(mesr);
    if (n.mesr<=7) col=rev(RColorBrewer::brewer.pal(n.mesr+2,pal))[1:n.mesr]
    else col=colorRampPalette(rev(RColorBrewer::brewer.pal(9,pal)[3:9]))(n.mesr);
    col=setNames(col,mesr);
  }));
  ## compute line widths to further discriminate measures
  lwd.mesr=do.call(c,sapply(cq(sig,dcc,meta,scp,other),USE.NAMES=F,function(what) {
    mesr=get(paste(sep='.','mesr',what));
    n.mesr=length(mesr);
    lwd=seq(1,3,len=n.mesr);
    lwd=setNames(lwd,mesr);
  }));
  ## compute line types to further discriminate measures
  lty.mesr=do.call(c,sapply(cq(sig,dcc,meta,scp,other),USE.NAMES=F,function(what) {
    mesr=get(paste(sep='.','mesr',what));
    n.mesr=length(mesr);
    ## crude effort to create divergent line types, up to lty.max
    gap=ceiling(n.mesr/2);
    on=(do.call(c,lapply(seq_len(floor(n.mesr/2)),function(i) c(i,i+gap))));
    on=1+on%%(lty.max-1);
    off=on+1;
    lty=c('solid',paste(sep='',as.hexmode(on),as.hexmode(off)))[1:n.mesr];
    lty=setNames(lty,mesr);
  }));
  ## setup mesr.from & relto variable for default interpretation
  mesr.fromtype=sapply(mesr.all,function(mesr) {
    if (mesr %in% mesr.fromraw) 'raw'
    else if (mesr %in% mesr.frombsln) 'bsln'
    else if (mesr %in% mesr.fromsig2) 'sig2'
    else stop(paste('No from.type for mesr:',mesr))});
  mesr.reltotype=sapply(mesr.all,function(mesr) {
    if (mesr %in% mesr.relraw) 'raw'
    else if (mesr %in% mesr.relsig1) 'sig1'
    else if (mesr %in% mesr.relsig2) 'sig2'
    else stop(paste('No relto.type for mesr:',mesr))});
                                 
  ## at end, assign mesr parameters to global variables
  sapply(grep('mesr',ls(),value=T),function(what) assign(what,get(what),envir=.GlobalEnv));
  init.mesr<<-T;         # so dosmry will know init_mesr done
  invisible(T);
}
