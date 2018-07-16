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
## General function to generate figures and tables for repwr.R documents
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## --- Document Functions ---
## general functions to make figures and tables for documents
## run plot function, save if required, label result with function name
## CAUTION: ... interacts with partial argument matching to cause dofig args to be
##   matched by plot-function args. eg, 'doc' matched by 'd'. prepending with 'fig'
##   works only because no plot-function arg matches it
dofig=
  function(figfun,figname=NULL,figsect=parent(figsect,NULL),fignum=parent(fignum,1),
           figdoc=parent(doc,'tech'),fignew=parent(fignew,F),figsave=parent(figsave,save.fig),
           id=parent(id,NULL),
           ...) {
    figname=paste(collapse='_',c(figsect,figname));
    if (fignew) dev.new();
    
    dev=figfun(...,fignum=fignum);
    ## function may return multiple plots
    file=
      if (length(dev)==1) filename_fig(figname,fignum,figdoc,id)
      else sapply(seq_along(dev), function(i) filename_fig(figname,fignum,figdoc,id));
    for (i in seq_along(dev)) {
      if ((is.na(figsave)&!file.exists(file[i]))|(!is.na(figsave)&figsave))
        savePlot(file[i],device=dev[i]);
    }
    assign_parent(fignum,fignum+length(dev));
    setNames(dev,figname);
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

