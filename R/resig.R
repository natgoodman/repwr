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
source('R/check.R');
source('R/datman.R');
source('R/doc.R');
## source('R/doc_readme.R');
## source('R/doc_repwr.R');
source('R/doc_resig.R');
## source('R/doc_xperiment.R');
## source('R/doc_tech.R');
source('R/init.R');
source('R/plot.R');
source('R/select.R');
# source('R/sim.R');
source('R/sim_resig.R');
source('R/stats.R');
source('R/util.R');

## ---- run ----
## run the program
## parameters defined in init
run=function(doc='resig',...) {
  init(doc=doc,...);
  dodata(need.init=F,...);             # generate data - ie, run simulation
  dodoc(need.init=F,...);              # generate figures for doc
}
