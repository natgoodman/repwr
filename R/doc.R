#################################################################################
##
## Author:  Nat Goodman
## Created: 18-06-19
##          from repwr.R created 17-10-05
##          from swfdr.R created 17-10-02
##          from swfdr_base.R created 17-04-10
##          from fdr.R created 17-04-03
##          from pval_error.R created 16-11-01
## Restart: 18-02-15
##          from scope.R created 17-12-04
##
## Copyright (C) 2018 Nat Goodman.
## 
## Generate figures and tables for repwr.R blog post
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## --- Document Functions ---
## functions to make figures and tables for documents
## sect is which sections to run - for use during development
##   uses prefix matching and all matches run
## fignum is initial figure number
## fignew controls whether each figure drawn in new window
## figsave controls whether figures are saved
## docfun is document-specific function. default calculated from doc, eg, doc_repwr
dodoc=function(sect=NULL,fignum=1,fignew=F,figsave=T,docfun=get(paste(sep='_','doc',doc))) {
  docfun();
}
## --- Document Utilities ---
## utility functions to make figures and tables for documents
## run plot function, save if required, label result with function name
## CAUTION: ... interacts with partial argument matching to cause dofig args to be
##   matched by plot-function args. eg, 'doc' matched by 'd'. prepending with 'fig'
##   works only because no plot-function arg matches it
dofig=
  function(figfun,figname=NULL,figsect=parent(figsect,NULL),fignum=parent(fignum,1),
           fignew=parent(fignew,F),figsave=parent(figsave,save.fig),
           id=parent(id,NULL),
           ...) {
    figname=paste(collapse='_',c(figsect,figname));
    if (fignew) dev.new();
    
    dev=figfun(...,fignum=fignum);
    ## function may return multiple plots
    file=
      if (length(dev)==1) filename_fig(figname,fignum,id)
      else sapply(seq_along(dev), function(i) filename_fig(figname,fignum,id));
    for (i in seq_along(dev)) {
      if ((is.na(figsave)&!file.exists(file[i]))|(!is.na(figsave)&figsave))
        savePlot(file[i],device=dev[i]);
    }
    assign_parent(fignum,fignum+length(dev));
    setNames(dev,figname);
  }
## generate standard xdata for aggregated plots
xdata_doc=function(near=0,nx=2.5,n2.num=2,n1=seq(20,by=20,len=8)) {
  d2=round(seq(0,1,by=0.01),digits=5);
  do.call(rbind,lapply(n1,function(n1) {
    n2=seq(n1*nx,by=n1*nx,len=n2.num);
    xdata=expand.grid(n1=n1,n2=n2,d1=d,d2=d2);
    ## TODO: I don't thing the 1e-4 tolerance still needed
    subset(xdata,subset=(abs(d1-d2)<=(near+1e-4)));}))
}

## NOT YET PORTED
## run data-making function, save if required, store result in global workspace
dodata=function(what) {
  what=as.character(pryr::subs(what));
  f=get(what,envir=parent.frame(n=1),mode='function');
  data=f();
  what=sub('^do_','',what);
  if (save.data) save_data(data,what);
  ## fiddle with name to get the form we want
  assign(what,data,envir=.GlobalEnv);
  invisible(data);
}

