#################################################################################
##
## Author:  Nat Goodman
## Created: 18-08-04
##          from repwr_hack_proptrue.R created 18-07-25
##
## Copyright (C) 2018 Nat Goodman.
## 
## Experimental sandbox code related to replication power blog post
## Implements quick hack to determine whether adjusting prop.true
## has a big effect on repwr results. Short answer: no
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
  dodata(need.init=F,...); # generate data - ie, run simulation
  dodoc(need.init=F,docfun=doc_xperiment,...);  # generate figures for doc
}
## init for this sandbox
init_xperiment=init_hack_proptrue=
  function(doc='xperiment',clean=F,clean.fig=F,...) {
    subdoc='hack_proptrue';
    init(doc='xperiment',
         n=20*2^(0:4),d=c(0,0.2,0.5,0.8,1),m=1e3,mdir=paste_nv(m,m_pretty(m)),
         datadir=filename('data',docx,subdoc,mdir),figdir=filename('figure',docx,subdoc,mdir),
         clean=clean,clean.fig=clean.fig,...);
  }
## NG 18-07-21: quick hack to look at increasing number of false positives
##   for n=20,40 and d=0, runs dosim1 repeatedly to get 90% positives (default 900 out of 1000)
##   since there are 5 values of d, this sets FPR=0.9*m/5*m=0.18 vs. ~0.01 in standard sim
dosim=function() {
  ## do simulations! each value of n x d
  cases=expand.grid(n=n,d=d);
  apply(cases,1,function(case) {
    n=case['n']; d=case['d'];
    if (n<=40&d==0) {
      ## do 20x instances to get enough positives
      want.pos=round(0.9*m);
      sim=do.call(rbind,replicate(20,dosim1(n,d),simplify=F));
      sim.pos=subset(sim,subset=pval<=0.05);
      m.pos=nrow(sim.pos);
      if (m.pos>want.pos) sim.pos=sim.pos[sample.int(m.pos,want.pos),]
      else warning(paste(sep='','sim(',n,',',d,') only got ',m.pos,
                         ' positives. May not be enough'));
      sim.neg=subset(sim,subset=pval>0.05);
      sim.neg=sim.neg[sample.int(nrow(sim.neg),m-nrow(sim.pos)),];
      sim=rbind(sim.pos,sim.neg);
    } else sim=dosim1(n,d);             # do base simulation
    simr=dosimr(sim,n,d);               # add/substract columns to/from sim for results
  })
  invisible();
}
## make figures and tables for experimental sandbox code
## sect is which sections to run - for use during development
##   uses prefix matching and all matches run
doc_xperiment=doc_hack_proptrue=
  function(sect=parent(sect,NULL),
           fignum=parent(fignum,1),fignew=parent(fignew,T),figsave=parent(figsave,F)) {
    if (missing(fignew)) fignew=!figsave else if (missing(figsave)) figsave=!fignew;
    sect.all=cq(plotrate,heatrate,roc,multi_sig2,small_telescopes);
    if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all);
#################### nonzro
##### plotrate
    if ((figsect='plotrate') %in% sect) {
      doc='readme'; init(doc=doc); 
      dofig(plotrate,'nonzro_fpr',title.desc=doc,rate.rule='nonzro',mesr=mesr.heatdflt,
            d1=0,d2=seq(0.1,1,by=0.1),n1=20,n2=50,legend='topright')
      dofig(plotrate,'nonzro_fpr',title.desc=doc,rate.rule='nonzro',mesr=mesr.heatdflt,
            d1=0,d2=seq(0.1,1,by=0.1),n1=40,n2=100,legend='topright')
      doc='xperiment'; init(doc=doc); 
      dofig(plotrate,'nonzro_fpr',title.desc=doc,rate.rule='nonzro',mesr=mesr.heatdflt,
            d1=0,d2=seq(0.1,1,by=0.1),n1=20,n2=50,legend='topright')
      dofig(plotrate,'nonzro_fpr',title.desc=doc,rate.rule='nonzro',mesr=mesr.heatdflt,
            d1=0,d2=seq(0.1,1,by=0.1),n1=40,n2=100,legend='topright')
   }
##### heatrate
    if ((figsect='heatrate') %in% sect) {
      doc='readme'; init(doc=doc); 
      dofig(heatrate,'nonzro_fpr',title.desc=doc,rate.rule='nonzro',mesr=mesr.heatdflt,
            d1=0,d2=seq(0.1,1,by=0.1),n1=20,n2=50,smooth=F)
      dofig(heatrate,'nonzro_fpr',title.desc=doc,rate.rule='nonzro',mesr=mesr.heatdflt,
            d1=0,d2=seq(0.1,1,by=0.1),n1=40,n2=100,smooth=F)
      doc='xperiment'; init(doc=doc); 
      dofig(heatrate,'nonzro_fpr',title.desc=doc,rate.rule='nonzro',mesr=mesr.heatdflt,
            d1=0,d2=seq(0.1,1,by=0.1),n1=20,n2=50,smooth=F)
      dofig(heatrate,'nonzro_fpr',title.desc=doc,rate.rule='nonzro',mesr=mesr.heatdflt,
            d1=0,d2=seq(0.1,1,by=0.1),n1=40,n2=100,smooth=F)
    }
#####  plotroc
    if ((figsect='roc') %in% sect) {
      xdata=xdata_xperiment(near=0);
      doc='readme'; init(doc=doc); 
      dofig(plotroc,'exact',title.desc=paste('exact',doc),rate.rule='nonzro',
            xdata=xdata);
      doc='xperiment'; init(doc=doc); 
      dofig(plotroc,'exact',title.desc=paste('exact',doc),,rate.rule='nonzro',
            xdata=xdata);
      xdata=xdata_xperiment(near=1);
      doc='readme'; init(doc=doc); 
      dofig(plotroc,'inexact',title.desc=paste('inexact',doc),rate.rule='nonzro',
            xdata=xdata);
      doc='xperiment'; init(doc=doc); 
      dofig(plotroc,'inexact',title.desc=paste('inexact',doc),,rate.rule='nonzro',
            xdata=xdata);

    }
#####  plotrocm, plotragm
    if ((figsect='multi_sig2') %in% sect) {
      ## near exact. sig2
      near=round(c(0.01,0.05,0.1,0.2),digits=5);
      xdata=lapply(c(0,near,1),function(near) xdata=xdata_xperiment(near=near));
      names(xdata)=c('exact',paste(sep=' ','near',near),'inexact');
      doc='readme'; init(doc=doc); 
      dofig(plotrocm,'rocm',rate.rule='nonzro',xdata=xdata,title.desc=paste('near',doc),
            mesr='sig2');
      dofig(plotragm,'ragm',rate.rule='nonzro',xdata=xdata,title.desc=paste('near',doc),
             smooth='aspline',mesr='sig2');
      doc='xperiment'; init(doc=doc); 
      dofig(plotrocm,'rocm',rate.rule='nonzro',xdata=xdata,title.desc=paste('near',doc),
            mesr='sig2');
       dofig(plotragm,'ragm',rate.rule='nonzro',xdata=xdata,title.desc=paste('near',doc),
             smooth='aspline',mesr='sig2');
    }
##########
    if ((figsect='small_telescopes') %in% sect) {
      xdata=xdata_xperiment(near=0.1);
      doc='readme'; init(doc=doc); 
      dofig(plotroc,'roc',rate.rule='nonzro',xdata=xdata,
            title.desc=paste('near=0.1',doc),mesr=cq(sig2,d2.scp1));
      dofig(plotrag,'rag',rate.rule='nonzro',xdata=xdata,
            title.desc=paste('near=0.1',doc),smooth='aspline',mesr=cq(sig2,d2.scp1));
      doc='xperiment'; init(doc=doc); 
      dofig(plotroc,'roc',rate.rule='nonzro',xdata=xdata,
            title.desc=paste('near=0.1',doc),mesr=cq(sig2,d2.scp1));
      dofig(plotrag,'rag',rate.rule='nonzro',xdata=xdata,
            title.desc=paste('near=0.1',doc),smooth='aspline',mesr=cq(sig2,d2.scp1));
   }
    sect;
  }
## generate standard xdata for aggregated plots
xdata_xperiment=function(near=0,nx=2.5,n2.num=2,n1=c(20,40,60)) {
  d2=round(seq(0,1,by=0.01),digits=5);
  do.call(rbind,lapply(n1,function(n1) {
    n2=seq(n1*nx,by=n1*nx,len=n2.num);
    xdata=expand.grid(n1=n1,n2=n2,d1=d,d2=d2);
    ## TODO: I don't thing the 1e-4 tolerance still needed
    subset(xdata,subset=(abs(d1-d2)<=(near+1e-4)));}))
}

