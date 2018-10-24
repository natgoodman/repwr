#################################################################################
##
## Author:  Nat Goodman
## Created: 18-09-08
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
## Top level file for resig blog post.
## Specializes repwr.R for the very simple resig analysis, mainly by loading
##   sim_resig.R instead of sim.R
##
## This is a basic implementation meant to be run interactively.
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
source('R/repwr.R');
## overlay specialized files
source('R/sim_resig.R');

## ---- run ----
## run the program
## parameters defined in init
run=function(need.init=T,doc='resig',...) {
  if (need.init) wrap_fun(init);
  need.init=F;
  wrap_fun(dodata,...);               # generate data - ie, run simulation
  wrap_fun(dodoc,init_doc,...);       # generate figures, tables for doc
}
