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
source('R/datman.R');
source('R/init.R');
source('R/plot.R');
source('R/select.R');
source('R/sim.R');
source('R/stats.R');
source('R/util.R');

## ---- run ----
## run the program
## parameters defined in init
run=function(...) {
  init(...);                     # process parameters & other initialization
  dopre();                       # precalculate or load global data
  dosim();                       # load saved simulations or do new ones
  dopost();                      # post-simulation pipeline
}
## --- Blog Functions ---
## make plots and tables for blog post.
## NOT YET PORTED
## run plot function, save if required, label result with function name
doplot=function(what,id=parent(id,NULL)) {
  what=as.character(pryr::subs(what));
  f=get(what,envir=parent.frame(n=1),mode='function');
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

