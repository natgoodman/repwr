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
run=function(save.fig=T,clean.fig=save.fig,...) {
  init_xperiment(save.fig=T,clean.fig=save.fig,...);
  ## dodata(need.init=F,...);              # use repwr data
  dodoc(need.init=F,docfun=doc_xperiment,...);  # generate figures for doc
}
## init for this sandbox
init_xperiment=init_uri_answer01=
  function(doc='xperiment',...) {
    subdoc='uri_answer01';
    init(doc='repwr',
         ## use repwr data and other params
         figdir=filename('figure','xperiment',subdoc,'m=1e4'),...);
  }
## make figures and tables for experimental sandbox code
## sect is which sections to run - for use during development
##   uses prefix matching and all matches run
doc_xperiment=doc_uri_answer01=
  function(sect=parent(sect,NULL),
           fignum=parent(fignum,1),fignew=parent(fignew,T),figsave=parent(figsave,F)) {
    if (missing(fignew)) fignew=!figsave else if (missing(figsave)) figsave=!fignew;
    sect.all=c('20-50','20-500');
    if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all);
    mesr=cq(sig2,d2.scp1);
##### 20-50: n1=20, n2=50
    if ((figsect='20-50') %in% sect) {
      n1=20; n2=50; d2=seq(0,1,0.1);
      d.33=uniroot(function(d0) p_d2t(n1,d=d1,d0,lower.tail=F)-0.33,interval=c(0,10))$root;
      d.scp=uniroot(function(d) p_d2t(n2,d,d0=d.33)-0.05,c(-10,10))$root;
      dofig(plotrate,'20-50_fpr',rate.rule='nonzro',,rate.type='error',mesr=mesr,
            d1=0,d2=d2,n1=n1,n2=n2,vline=c(d.scp,d.33),vhlty='dashed',legend='topright');
     dofig(plotrate,'20-50_fnr',rate.rule='nonzro',rate.type='error',mesr=mesr,
            d1=0.5,d2=d2,n1=n1,n2=n2,vline=c(d.scp,d.33),vhlty='dashed',legend='topright');
   }
##### 20-500: n1=20, n2=500
    if ((figsect='20-500') %in% sect) {
      n1=20; n2=500; d2=seq(0,1,0.1);
      d.33=uniroot(function(d0) p_d2t(n1,d=d1,d0,lower.tail=F)-0.33,interval=c(0,10))$root;
      d.scp=uniroot(function(d) p_d2t(n2,d,d0=d.33)-0.05,c(-10,10))$root;
      dofig(plotrate,'20-50_fpr',rate.rule='nonzro',,rate.type='error',mesr=mesr,
            d1=0,d2=d2,n1=n1,n2=n2,vline=c(d.scp,d.33),vhlty='dashed',legend='topright');
     dofig(plotrate,'20-50_fnr',rate.rule='nonzro',rate.type='error',mesr=mesr,
            d1=0.5,d2=d2,n1=n1,n2=n2,vline=c(d.scp,d.33),vhlty='dashed',legend='topright');
    }
    sect;
  }
