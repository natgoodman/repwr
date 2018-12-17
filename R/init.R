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
  ## doc parameters 
  doc='readme',                     # controls sim defaults, data, figure subdirs
  docx=match.arg(doc,cq(readme,resig,repwr,tech,xperiment)), 
  ## simulation parameters 
  n=switch(docx,                            # sample sizes
           readme=20*2^(0:4),               # 20,40,80,160,320 (5 values)
           resig=c(20,seq(50,500,by=50)),   # 20,50,100,...,1000 (11 values)
           repwr=c(10,20,50*1.454215^(0:8)),# 10,20,50,...,1000 (11 values)
           tech=10*2^(0:9),                 # 10,20,40,...,5120 (10 values)
           xperiment=NA),                   # xperiment must supply values
  d=switch(docx,                            # population effect sizes
           readme=c(0,0.2,0.5,0.8,1),
           resig=seq(0,1,by=0.1),repwr=seq(0,1,by=0.1),tech=seq(0,1,by=0.1),
           xperiment=NA),                   # xperiment must supply values
  m=switch(docx,                            # instances per study (ie, population)
           readme=1e3,resig=1e4,repwr=1e4,tech=1e4,
           xperiment=NA),                   # xperiment must supply values
  ## analysis parameters
  sig.level=0.05,                   # for conventional significance
  conf.level=0.95,                  # for confidence intervals
  pred.level=0.95,                  # for prediction intervals
                                    # grid for various precacluated data
  dsdz.grid=seq(min(d)-3,max(d)+3,by=.05),
  scope.power=0.33,                 # power for small telescope
  scope.close=0.05,                 # pval for big telescope "close enough" calculation
  dopost.allcases=T,                # do all cases (2500 with defaults) or 1500 select ones
  dopost.permute=T,                 # permute s1,s2 data. else self-comparisons meaningless
  ## program parameters, eg, for output files, error messages, etc.
  scriptname='repwr',                      #
  mdir=paste_nv(m,m_pretty(m)),            # m subdirectory
  datadir=file.path('data',docx,mdir),     # directory for data files. default eg, data/repwr/m=1e4
  simdir=file.path(datadir,'sim'),         # directory for sim files
  simrdir=file.path(datadir,'simr'),       # directory for simr files
  sidir=file.path(datadir,'si'),           # directory for study permutation index files
  detldir=file.path(datadir,'detl'),       # directory for detailed results files
  smrydir=file.path(datadir,'smry'),       # directory for summary results files
  posrdir=file.path(datadir,'posr'),       # directory for positive rate files
  ## NG 18-10-18: figdir, tbldir moved to init_doc
  ## figdir=file.path('figure',docx,mdir), # directory for figures. default eg, figure/repwr/m=1e4
  ## tbldir=file.path('table',docx,mdir), # directory for tables. default eg, table/repwr/m=1e4
  id=NULL,                                 # info tacked onto filenames. not used much
  verbose=F,                               # print progress messages
  ## program control
  must.exist=F,                  # must all sub-inits succeed?
  load=NA,                       # shorthand for other load params
                                 #   NA means load if file exists
                                 #   T, F mean always or never load
  load.sim=load,                 # load saved simulations
  load.simr=load,                # load saved simr files - sim customized for analysis
  load.si=load,                  # load saved study permutation index files
  load.detl=load,                # load saved analysis detail files
  load.smry=load,                # load saved analysis summary files
  load.posr=load,                # load saved positive rate files
  load.data=load,                # load saved top level data files
  save=NA,                       # shorthand for other save params 
                                 #   NA means save unless file exists
                                 #   T, F mean always or never save
  save.sim=save,                 # save simulations (RData format)
  save.simr=save,                # save simr files - sim customized for analysis (RData format)
  save.si=save,                  # save study permutation index files (RData format)
  save.detl=save,                # save analysis detail files (RData format)
  save.smry=save,                # save analysis summary files (RData)
  save.posr=save,                # save positive rate files (RData)
  save.data=save,                # save top level results (RData & txt formats)
  ## NG 18-10-18: save.fig, save.tbl moved to init_doc
  ## save.fig=T,                 # save figures (when called via dofig)
  ## save.tbl=T,                 # save tables
  save.txt=NA,                   # save results in txt format as well as RData
                                 #   NA means use default rule for type:
                                 #   F for all but top level data
  save.txt.sim=!is.na(save.txt)&save.txt,  # save txt simulations. default F
  save.txt.simr=!is.na(save.txt)&save.txt, # save txt simr data. default F
  save.txt.si=!is.na(save.txt)&save.txt,   # save txt study permutation indexes. default F
  save.txt.detl=!is.na(save.txt)&save.txt, # save txt analysis details. default F
  save.txt.smry=!is.na(save.txt)&save.txt, # save txt case-by-case summaries. default F
  save.txt.posr=is.na(save.txt)|save.txt,  # save txt positive rate files. default T
  save.txt.data=is.na(save.txt)|save.txt,  # save txt top level results. default T
  ## NG 18-10-18: save.txt.tbl moved to init_doc
  ## save.txt.tbl=T,                # save txt tables. default T
  keep=NA,                       # shorthand for other keep params 
                                 #   NA means use default keep rule for type:
                                 #   T for all but detl
                                 #   T, F mean always or never keep
  keep.sim=!is.na(keep)&keep,    # keep simulations. default F
  keep.simr=is.na(keep)|keep,    # keep simr data. default T
  keep.si=is.na(keep)|keep,      # keep study permutation indexes. default F
  keep.detl=!is.na(keep)&keep,   # keep analysis details. default F
  keep.smry=!is.na(keep)&keep,   # keep case-by-case summaries. default F
  keep.posr=is.na(keep)|keep,    # keep positive rate data. default T
  keep.data=is.na(keep)|keep,    # keep top-level data. default T
                                 #    
  clean=F,                       # remove everything and start fresh
  clean.data=clean,              # remove datadir & memlist
  ## NG 18-10-18: clean.fig, clean.tbl moved to init_doc
  ## clean.fig=clean,               # remove figdir
  ## clean.tbl=clean,               # remove tbldir
  clean.memlist=T,               # clean memlist cache - always safe
  clean.sim=F,                   # clean simulations. default F
  clean.simr=F,                  # clean simr data. default F
  clean.si=F,                    # clean study permutation indexes. default F
  clean.detl=F,                  # clean analysis details. default F
  clean.smry=F,                  # clean case-by-case summaries. default F
  clean.posr=F,                  # clean positive rate files. default F
  clean.toplevel=F,              # clean top-level data. default F
  ## NG 18-10-18: plot control moved to init_doc
  ## plot control
  ## figscreen=if (docx=='readme') T else F, # plot figures on screen
  ## fignew=figscreen,              # plot each figure in new window
  ##
  end=NULL                       # placeholder for last parameter
  ) {
  doc=docx;                      # to avoid confusion later
  if (doc=='xperiment'&any(is.na(c(n,d,m))))
    stop('doc=xperiment but no value provided for n, d, or m');
  ## round d, dsdz.grid to avoid imprecise decimals
  ## round n since 1.454215^(0:8) generates fractions
  d=round(d,digits=5); dsdz.grid=round(dsdz.grid,digits=5); n=round(n);
  ## assign parameters to global variables
  ## do it before calling any functions that rely on globals
  assign_global();
  ## clean and create output directories and internal memory values as needed
  ## NG 18-10-18:  figdir, tbldir moved to init_doc
  ## outdir=c(datadir,simdir,simrdir,sidir,detldir,smrydir,posrdir,figdir,tbldir);
  outdir=c(datadir,simdir,simrdir,sidir,detldir,smrydir,posrdir);
  memlist=cq(sim.list,simr.list,si.list,detl.list,smry.list,posr.list,data.list);
  if (clean.data) unlink(datadir,recursive=T);
  ## NG 18-10-18: clean.fig, clean.tbl moved to init_doc
  ## if (clean.fig) unlink(figdir,recursive=T);
  ## if (clean.tbl) unlink(tbldir,recursive=T);
  if (clean.memlist) suppressWarnings(rm(list=memlist,envir=.GlobalEnv));
  ## clear init_smry & init_mesr flags so these will be rerun
  suppressWarnings(rm(list=cq(init.smry,init.mesr),envir=.GlobalEnv));
  ## clean specific types if desired
  if (clean.sim) cleanq(sim);
  if (clean.simr) cleanq(simr);
  if (clean.si) cleanq(si);
  if (clean.detl) cleanq(detl);
  if (clean.smry) cleanq(smry);
  if (clean.posr) cleanq(posr);
  if (clean.toplevel) cleanq(data,cleandir=F);
  ## create data subdirectories. nop if already exist
  sapply(outdir,function(dir) dir.create(dir,recursive=TRUE,showWarnings=FALSE));
  ## setup in-memory lists to hold simulations, etc.. do carefully in case already setup
  sapply(memlist,function(what) if (!exists(what)) assign(what,list(),envir=.GlobalEnv));
  invisible();
}
## initialize measures.
## init_mesr called by doposr to init parameters used for constructing posrs
##   can't run until smry object exists!
## also called by init_docmesr for the common case of generating doc separately from data  
## init_docmesr called by init_doc to init parameters used for doc generation
init_mesr=function(must.exist=T) {
  ## return immediately if already initialized
  if (exists('init.mesr',envir=.GlobalEnv)&&init.mesr) invisible(T);
  mesr.all=get_data(mesr,must.exist=must.exist);
  if (is.null(mesr.all)) {
    init.mesr<<-F;    # so doposr will know init_mesr not done 
    invisible(F);
  }
  ## measures grouped by source row in smry
  mesr.fromraw=cq(sig1,sdir);
  mesr.frombsln=setdiff(grep('scp',mesr.all,invert=T,value=T),mesr.fromraw);
  mesr.fromsig2=setdiff(mesr.all,c(mesr.frombsln,mesr.fromraw));
  ## measures grouped by relative row (denominator in rate calculation) in std interpretation
  mesr.relraw=mesr.fromraw;
  mesr.relsig1=c(mesr.frombsln,mesr.fromsig2);
  mesr.relsig2=NULL;
  ## mesr.relsig2=mesr.fromsig2;
  ## measures grouped by functional category
  mesr.sig=cq(sig2,sigm);
  ## mesr.dcc=cq(d1.c2,d2.c1,c1.c2,d1.p2,d2.p1,p1.p2);
  mesr.dcc=grep('(d(1|2)\\.(c|p)(1|2))|((c|p)1\\.(c|p)2)',mesr.all,value=T);
  ## mesr.scp=cq(d1.scp2,d2.scp1,d1.scpd2,d2.scpd1);
  mesr.scp=grep('d(1|2)\\.scp(d{0,1})(1|2)',mesr.all,value=T);
  mesr.meta=grep('dm|cm',mesr.all,value=T);
  mesr.other=setdiff(mesr.all,c(mesr.sig,mesr.dcc,mesr.scp,mesr.meta));
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
init_docmesr=function(must.exist=T) {
  init_mesr(must.exist=must.exist);
  ## defaults for plot functions. depends on doc
  if (doc=='readme') {
    mesr.dflt=cq(sig2,d1.c2,sigm,d2.c1,c1.c2,d1.p2,d2.p1,p1.p2,d2.scp1);
    mesr.plotdflt=mesr.ragdflt=cq(sig2,d1.c2,sigm,d2.c1);
    mesr.heatdflt=mesr.rocdflt=grep('scp',mesr.dflt,invert=T,value=T);
    mesr.order=mesr.dflt;
    n.mesr=length(mesr.dflt);
    col.mesr=c(colorRampPalette(RColorBrewer::brewer.pal(min(8,n.mesr-1),'Set1'))(n.mesr-1),
               'blue');
    ## manually fix 6th color (d1.p2) - make it darker
    ## col.mesr[6]='#FFcc00';
    ## use line widths, point cex to further discriminate measures
    ## sig2 is biggest. others gradually diminish. d2.scp1 is special - shouldn't be too small
    lwd.mesr=c(2,seq(1.5,0.75,len=n.mesr-2),1);
    cex.mesr=c(1,seq(0.9,0.5,len=n.mesr-2),0.75);
    ## some docs use line types to further discriminate measures. readme doesn't
    lty.mesr=rep('solid',n.mesr);
    ## set names in all these lists
    ## CAUTION: have to use loop (not sapply) for scoping to work
    for (name in cq(col.mesr,lwd.mesr,cex.mesr,lty.mesr))
      assign(name,setNames(get(name),mesr.dflt));
  } else if (doc=='resig') {
    mesr.dflt='sig2';
    mesr.plotdflt=mesr.heatdflt=mesr.rocdflt=mesr.ragdflt=mesr.dflt;
    mesr.order=mesr.dflt;
    n.mesr=length(mesr.dflt);
    col.mesr=colorRampPalette(RColorBrewer::brewer.pal(max(3,min(8,n.mesr)),'Set1'))(n.mesr);
    lwd.mesr=2;
    cex.mesr=1;
    ## some docs use line types to further discriminate measures. repwr doesn't
    lty.mesr=rep('solid',n.mesr);
    ## set names in all these lists
    ## CAUTION: have to use loop (not sapply) for scoping to work
    for (name in cq(col.mesr,lwd.mesr,cex.mesr,lty.mesr))
      assign(name,setNames(get(name),mesr.dflt));
  } else if (doc=='repwr'&is.null(subdoc)) {
    mesr.dflt=cq(sig2,d1.c2,sigm,d2.c1,c1.c2,d1.p2,d2.p1,p1.p2,d2.scp1);
    mesr.plotdflt=mesr.heatdflt=mesr.rocdflt=mesr.ragdflt=grep('scp',mesr.dflt,invert=T,value=T);
    mesr.order=mesr.dflt;
    n.mesr=length(mesr.dflt);
    col.mesr=c(colorRampPalette(RColorBrewer::brewer.pal(min(8,n.mesr-1),'Set1'))(n.mesr-1),
               'blue');
    ## manually fix 6th color (d1.p2) - make it darker
    col.mesr[6]='#FFcc00';
    ## use line widths, point cex to further discriminate measures
    ## sig2 is biggest. others gradually diminish. d2.scp1 is special - shouldn't be too small
    lwd.mesr=c(2,seq(1.5,0.75,len=n.mesr-2),1);
    cex.mesr=c(1,seq(0.9,0.5,len=n.mesr-2),0.75);
    ## some docs use line types to further discriminate measures. repwr doesn't
    lty.mesr=rep('solid',n.mesr);
    ## set names in all these lists
    ## CAUTION: have to use loop (not sapply) for scoping to work
    for (name in cq(col.mesr,lwd.mesr,cex.mesr,lty.mesr))
      assign(name,setNames(get(name),mesr.dflt));
  } else if (doc=='repwr'&subdoc=='supp') {
    ## TODO: these are old. refine based on experience
    ## mesr.plotdflt=cq(sig2,sigm,d1.c2,d2.c1,d1.p2,d2.p1);
    mesr.plotdflt=c(mesr.sig,mesr.dcc,'d1.scp2','d2.scp1');
    mesr.heatdflt=c(mesr.sig,mesr.dcc,mesr.scp,mesr.meta);
    mesr.rocdflt=c(mesr.sig,mesr.dcc);
    mesr.ragdflt=cq(sig2,d1.c2);
    ## measures ordered for legends and such
    mesr.order=c(mesr.sig,mesr.dcc,mesr.scp,mesr.other);
    ## RColorBrewer pallete names for plotrate
    pal.sig='Reds';
    pal.dcc='Blues';
    pal.meta='Greens';
    pal.scp='Purples';
    pal.other='Greys';
    ## line properties
    lty.max=8;                          # max 'on' value for dashed lines
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
      lwd=seq(3,1,len=n.mesr);
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
    ## compute cex for points in plotroc
    cex.mesr=do.call(c,sapply(cq(sig,dcc,meta,scp,other),USE.NAMES=F,function(what) {
      mesr=get(paste(sep='.','mesr',what));
      n.mesr=length(mesr);
      cex=seq(1,0.5,len=n.mesr);
      cex=setNames(cex,mesr);
    }));
  } else if (doc=='xperiment') {
    mesr.dflt=cq(sig2,d1.c2,sigm,d2.c1,c1.c2,d1.p2,d2.p1,p1.p2,d2.scp1);
    mesr.plotdflt=mesr.ragdflt=cq(sig2,d1.c2,sigm,d2.c1);
    mesr.heatdflt=mesr.rocdflt=grep('scp',mesr.dflt,invert=T,value=T);
    mesr.order=mesr.dflt;
    n.mesr=length(mesr.dflt);
    col.mesr=c(colorRampPalette(RColorBrewer::brewer.pal(min(8,n.mesr-1),'Set1'))(n.mesr-1),
               'blue');
    ## use line widths, point cex to further discriminate measures
    ## sig2 is biggest. others gradually diminish. d2.scp1 is special - shouldn't be too small
    lwd.mesr=c(2,seq(1.5,0.75,len=n.mesr-2),1);
    cex.mesr=c(1,seq(0.9,0.5,len=n.mesr-2),0.75);
    lty.mesr=rep('solid',n.mesr);
    ## set names in all these lists
    ## CAUTION: have to use loop (not sapply) for scoping to work
    for (name in cq(col.mesr,lwd.mesr,cex.mesr,lty.mesr))
      assign(name,setNames(get(name),mesr.dflt));
  }
  ## at end, assign mesr parameters to global variables
  sapply(grep('mesr',ls(),value=T),function(what) assign(what,get(what),envir=.GlobalEnv));
  init.mesr<<-T;         # so dosmry will know init_mesr done
  invisible(T);
}

## initialize doc parameters
init_doc=function(
  subdoc=NULL,
  subdocx=
    if(doc=='xperiment') subdoc else if(is.null(subdoc)) NULL else match.arg(subdoc,cq(supp)),
  ## output directories. filename function ignores subdoc if NULL
  figdir=filename('figure',doc,subdocx,mdir), # directory for figures, eg, figure/repwr/m=1e4
  tbldir=filename('table',doc,subdocx,mdir),  # directory for tables, eg, table/repwr/m=1e4
  ## output modifiers
  outpfx=if(doc=='xperiment'|is.null(subdocx)) NULL else 'S',
  outsfx=if(doc=='xperiment'|is.null(subdocx)) NULL else letters, # used in figure and table blocks
  sectnum=if(doc=='xperiment'|is.null(subdocx)) NULL else T, # add section number to prefix eg, S1
  ## figures
  figpfx=outpfx,
  figsfx=outsfx,
  fignum=1,
  figblk=NULL,                  # index into figsfx if in figure block
  ## tables
  tblpfx=outpfx,
  tblsfx=outsfx,
  tblnum=1,
  tblblk=NULL,                  # index into tblsfx if in table block
  ## xtra figures - not included in document
  xfigpfx='X',
  xfigsfx=outsfx,
  xfignum=1,
  xfigblk=NULL,                 # index into xfigsfx if in figure block
  ## error cutoffs for plots
  fpr.cutoff=0.05,              # false positive rate cutoff for plots
  fnr.cutoff=0.20,              # false negative rate cutoff for plots
  ## clean, save
  save.out=T,
  save.fig=save.out,            # save figures (when called via dofig)
  save.tbl=save.out,            # save tables (when called via dotbl)
  save.txt.tbl=T,               # save txt tables. default T
  clean.out=F,
  clean.fig=clean.out,          # remove figdir
  clean.tbl=clean.out,          # remove tbldir
  ## plot control
  figscreen=if(doc=='readme') T else !save.fig,
                                 # plot figures on screen
  fignew=figscreen,              # plot each figure in new window
  figextra=T,                    # plot extra figures, too
  ## doc generation function
  docfun=get(paste(collapse='',c('doc_',doc,subdocx))),
  docsect=NULL,                  # all document sections. set by docfun
  end=NULL                       # placeholder for last parameter
  ) {
  subdoc=subdocx;                # to avoid confusion later
  ## assign parameters to global variables
  ## do it before calling any functions that rely on globals
  assign_global();
  ## init mesr parameters for doc
  init_docmesr();
  ## clean and create output directories
  outdir=c(figdir,tbldir);
  if (clean.fig) unlink(figdir,recursive=T);
  if (clean.tbl) unlink(tbldir,recursive=T);
  sapply(outdir,function(dir) dir.create(dir,recursive=TRUE,showWarnings=FALSE));
  invisible();
}
