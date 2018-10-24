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
source('R/check.R');
source('R/datman.R');
source('R/doc.R');
source('R/doc_readme.R');
source('R/doc_repwr.R');
source('R/doc_resig.R');
source('R/doc_resigsupp.R');
source('R/doc_resig_fun.R');
source('R/doc_tech.R');
source('R/init.R');
source('R/plot.R');
source('R/select.R');
source('R/sim.R');
source('R/stats.R');
source('R/util.R');

## ---- run ----
## run the program
## parameters defined in init
run=function(need.init=T,...) {
  ## ## split ... args for init, dodata, dodoc. from stackoverflow.com/questions/4124900
  ## dots=list(...);
  ## init.args=dots[names(dots) %in% names(formals(init))];
  ## dodata.args=dots[names(dots) %in% names(formals(dodata))];
  ## dodoc.args=dots[names(dots) %in% c(names(formals(dodoc)),names(formals(init_doc)))];
  ## if (need.init) do.call('init',init.args);
  ## do.call('dodata',c(need.init=F,dodata.args)); # generate data - ie, run simulation
  ## do.call('dodoc',c(need.init=F,dodoc.args));   # generate figures, tables for doc
  if (need.init) wrap_fun(init);
  need.init=F;
  wrap_fun(dodata);                   # generate data - ie, run simulation
  wrap_fun(dodoc,init_doc); # generate figures, tables for doc
}
