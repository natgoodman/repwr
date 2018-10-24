#################################################################################
##
## Author:  Nat Goodman
## Created: 18-08-04
##          from repwr_uri_answer01.R created 18-08-02
##
## Copyright (C) 2018 Nat Goodman.
## 
## Experimental sandbox code related to replication power blog post
## Answers Uri's question in email Jul 28, 2018:
##    What's wrong with what I propose: having 80% power to reject d33%
##    if the truth is d=0?
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## source 'real code'. repwr.R get's 'em all
source('R/repwr.R');
## --- Eperimental sandbox functions ---
## run for this sandbox
run=function(need.init=T,...) {
  if (need.init) wrap_fun(init_xperiment);
  need.init=F;
  ## wrap_fun(dodata);                # use repwr data
  wrap_fun(dodoc,init_doc_xperiment,doc='xperiment'); # generate figures, tables for doc
}
## init for this sandbox. same parameters as rewpr
init_xperiment=function(doc='repwr',...) {
  ## call init with our arguments
  wrap_fun(init);
}
## init_doc for this sandbox
init_doc_xperiment=
  function(doc='xperiment',subdoc='uri_answer01',docfun=doc_xperiment,
           figdir=filename('figure',doc,subdoc,mdir),
           tbldir=filename('table',doc,subdoc,mdir),
           clean.out=T,figscreen=T,...) {
   doc<<-'xperiment';                    # bodily set doc to keep init_doc happy
   wrap_fun(init_doc);
 }

## make figures and tables for experimental sandbox code
## sect is which sections to run - for use during development
##   uses prefix matching and all matches run
doc_xperiment=function(sect=parent(sect,NULL)) {
  sect.all=c('20-50','20-500');
  if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all);
  mesr=cq(sig2,d2.scp1);
  sapply(sect,function(sect) {
##### 20-50: n1=20, n2=50
    if (sect=='20-50') {
      n1=20; n2=50; d1=d_sig(n1,0.05); d2=seq(0,1,0.1);
      d.33=uniroot(function(d0) p_d2t(n1,d=d1,d0,lower.tail=F)-0.33,interval=c(0,10))$root;
      d.scp=uniroot(function(d) p_d2t(n2,d,d0=d.33)-0.05,c(-10,10))$root;
      dofig(plotrate,'20-50_fpr',rate.rule='nonzro',,rate.type='error',mesr=mesr,
            d1=0,d2=d2,n1=n1,n2=n2,vline=c(d.scp,d.33),vhlty='dashed',legend='topright');
      dofig(plotrate,'20-50_fnr',rate.rule='nonzro',rate.type='error',mesr=mesr,
            d1=0.5,d2=d2,n1=n1,n2=n2,vline=c(d.scp,d.33),vhlty='dashed',legend='topright');
    }
##### 20-500: n1=20, n2=500
    if (sect=='20-500') {
      n1=20; n2=500; d1=d_sig(n1,0.05); d2=seq(0,1,0.1);
      d.33=uniroot(function(d0) p_d2t(n1,d=d1,d0,lower.tail=F)-0.33,interval=c(0,10))$root;
      d.scp=uniroot(function(d) p_d2t(n2,d,d0=d.33)-0.05,c(-10,10))$root;
      dofig(plotrate,'20-50_fpr',rate.rule='nonzro',,rate.type='error',mesr=mesr,
            d1=0,d2=d2,n1=n1,n2=n2,vline=c(d.scp,d.33),vhlty='dashed',legend='topright');
      dofig(plotrate,'20-50_fnr',rate.rule='nonzro',rate.type='error',mesr=mesr,
            d1=0.5,d2=d2,n1=n1,n2=n2,vline=c(d.scp,d.33),vhlty='dashed',legend='topright');
    }
    sect;
  })
}
