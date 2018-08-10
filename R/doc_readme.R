#################################################################################
##
## Author:  Nat Goodman
## Created: 18-07-19
##          from doc_repwr.R created 18-06-20
##          from doc.R created 18-06-19
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
## Generate figures and tables for README
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## --- Generate Figures and Tables for repwr README ---
## make figures and tables for blog post
## sect is which sections to run - for use during development
##   uses prefix matching and all matches run
doc_readme=
  function(sect=parent(sect,NULL),
           fignum=parent(fignum,1),figscreen=parent(figscreen,T),fignew=parent(fignew,F)) {
    sect.all=cq(plotrate,heatrate,roc,rag,multi_sig2,small_telescopes);
    if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all);
#################### nonzro
##### plotrate
    if ((figsect='plotrate') %in% sect) {
      dofig(plotrate,'nonzro_fpr',rate.rule='nonzro',d=0,n1=20,n2=n,legend='topright');
      dofig(plotrate,'nonzro_fnr',rate.rule='nonzro',d1=0.5,d2=d,n1=20,n2=50,legend='topright');
      dofig(plotrate,'sameff_fpr',rate.rule='sameff',d1=0.5,d2=d[d!=0.5],n1=20,n2=50,
            legend='bottomright');
      dofig(plotrate,'sameff_fnr',rate.rule='sameff',d=d,n1=20,n2=50,legend='topright');
    }
##### heatrate
    if ((figsect='heatrate') %in% sect) {
      dofig(heatrate,'nonzro_fpr',rate.rule='nonzro',d=0,n1=20,n2=n);
      dofig(heatrate,'nonzro_fnr',rate.rule='nonzro',d1=0.5,d2=d,n1=20,n2=50);
      ## show several values of n2,d2 in one heatmap
      xdata=expand.grid(n1=20,n2=c(50,100,150),d1=0,d2=d);
      dofig(heatrate,'nonzro_fpr_multi',rate.rule='nonzro',xdata=xdata,vline=c(3.5,6.5,9.5,12.5),
            smooth=F);
      xdata=expand.grid(n1=20,n2=c(50,100,150),d1=0.5,d2=d);
      dofig(heatrate,'nonzro_fnr_multi',rate.rule='nonzro',xdata=xdata,vline=c(3.5,6.5,9.5,12.5),
            smooth=F);
    }
#####  plotroc
    if ((figsect='roc') %in% sect) {
      xdata=xdata_readme(near=0);
      dofig(plotroc,'exact',title.desc='exact replication',rate.rule='nonzro',xdata=xdata);
      xdata=xdata_readme(near=1);
      dofig(plotroc,'inexact',title.desc='inexact replication',rate.rule='nonzro',xdata=xdata);
    }
#####  plotrag
    if ((figsect='rag') %in% sect) {
      xdata=xdata_readme(near=0);
      dofig(plotrag,'exact',title.desc='exact replication',rate.rule='nonzro',xdata=xdata,
            smooth='aspline');
      xdata=xdata_readme(near=1);
      dofig(plotrag,'inexact',title.desc='inexact replication',rate.rule='nonzro',xdata=xdata,
            smooth='aspline');
    }
#####  plotrocm, plotragm
    if ((figsect='multi_sig2') %in% sect) {
      ## near exact. sig2
      near=round(c(0.01,0.05,0.1,0.2),digits=5);
      xdata=lapply(c(0,near,1),function(near) xdata=xdata_readme(near=near));
      names(xdata)=c('exact',paste(sep=' ','near',near),'inexact')
      dofig(plotrocm,'rocm',rate.rule='nonzro',xdata=xdata,title.desc='near exact',
            mesr='sig2');
      dofig(plotragm,'ragm',rate.rule='nonzro',xdata=xdata,title.desc='near exact',
            smooth='aspline',mesr='sig2');
    }
##########
    if ((figsect='small_telescopes') %in% sect) {
      xdata=xdata_readme(near=0.1);
      dofig(plotroc,'roc',rate.rule='nonzro',xdata=xdata,
            title.desc='near=0.1',mesr=cq(sig2,d2.scp1));
      dofig(plotrag,'rag',rate.rule='nonzro',xdata=xdata,
            title.desc='near=0.1',smooth='aspline',mesr=cq(sig2,d2.scp1));
    }
    sect;
  }
## generate standard xdata for aggregated plots
xdata_readme=function(near=0,nx=2.5,n2.num=2,n1=c(20,40,60)) {
  d2=round(seq(0,1,by=0.01),digits=5);
  do.call(rbind,lapply(n1,function(n1) {
    n2=seq(n1*nx,by=n1*nx,len=n2.num);
    xdata=expand.grid(n1=n1,n2=n2,d1=d,d2=d2);
    ## TODO: I don't thing the 1e-4 tolerance still needed
    subset(xdata,subset=(abs(d1-d2)<=(near+1e-4)));}))
}
