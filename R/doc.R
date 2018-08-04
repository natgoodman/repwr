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
## --- Generate Figures and Tables for Document ---
## sect is which sections to run - for use during development
##   uses prefix matching and all matches run
## fignum is initial figure number
## fignew controls whether each figure drawn in new window
## figsave controls whether figures are saved
##   if only one of fignew, figsave set, other is set to complement
## docfun is document-specific function. default calculated from doc, eg, doc_repwr
dodoc=
  function(sect=NULL,need.init=T,doc=parent(doc,'readme'),fignum=1,
           fignew=if (doc=='readme') T else F,
           docfun=get(paste(sep='_','doc',doc)),...) {
    if (need.init) {
      ## for sandbox runs, use doc-specific init
      if (doc=='xperiment') init_xperiment(doc=doc,...) else init(doc=doc,...);
    }
    docfun();
}
## --- Document Functions ---
## utility functions to make figures and tables for documents
## run plot function, save if required, label result with function name
## CAUTION: ... interacts with partial argument matching to cause dofig args to be
##   matched by plot-function args. eg, 'doc' matched by 'd'. prepending with 'fig'
##   works only because no plot-function arg matches it
dofig=
  function(figfun,figname=NULL,figsect=parent(figsect,NULL),fignum=parent(fignum,1),
           fignew=parent(fignew,F),id=parent(id,NULL),
           ...) {
    figname=paste(collapse='_',c(figsect,figname));
    if (fignew) dev.new();
    
    dev=figfun(...,fignum=fignum);
    ## function may return multiple plots
    file=
      if (length(dev)==1) filename_fig(figname,fignum,id)
      else sapply(seq_along(dev), function(i) filename_fig(figname,fignum,id));
    for (i in seq_along(dev)) {
      if ((is.na(save.fig)&!file.exists(file[i]))|(!is.na(save.fig)&save.fig))
        savePlot(file[i],device=dev[i]);
    }
    assign_parent(fignum,fignum+length(dev));
    setNames(dev,figname);
  }
## NOT YET PORTED
## run data-making function, save if required, store result in global workspace
## dodata=function(what) {
##   what=as.character(pryr::subs(what));
##   f=get(what,envir=parent.frame(n=1),mode='function');
##   data=f();
##   what=sub('^do_','',what);
##   if (save.data) save_data(data,what);
##   ## fiddle with name to get the form we want
##   assign(what,data,envir=.GlobalEnv);
##   invisible(data);
## }

