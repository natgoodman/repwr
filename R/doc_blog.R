#################################################################################
##
## Author:  Nat Goodman
## Created: 18-06-20
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
## Generate figures and tables for repwr.R blog post
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## --- Blog Post ---
## make figures and tables for blog post
## docs are blog, tech, readme
## sect is which sections to run - for use during development
##   uses prefix matching and all matches run
doc_blog=
  function(sect=NULL,fignonz=cq(nonzro,nonz1,nonz1or2,nonz1and2,nonz2),figtol=0,
           fignum=1,figdoc='blog',
           ## mesr=cq(sig2,sigm,d1.c2,d2.c1,c1.c2,d1.p2,d2.p1,p1.p2,d2.scp1),
           mesr=cq(sig2,d1.c2,sigm,d2.c1,c1.c2,d1.p2,d2.p1,p1.p2),
           fignew=F,figsave=T) {
    sect.all=cq(nonzro_exact_fpr,nonzro_exact_fnr,nonzro_exact_roc,
                nonzro_inexact_roc,nonzro_nearexact,sameff_nearexact,small_telescopes,
                last_placeholder);
    if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all);
    fignonz=match.arg(fignonz);
#################### nonzro
##### exact
   if ((figsect='nonzro_exact_fpr') %in% sect) {
      ## exact reps. false positives
      ## 2x2 panel for d=0
      dofig(plotrate,'n1=20_plot',rate.rule=fignonz,d=0,n1=20,n2=seq(50,by=50,len=10),
            mesr=mesr,legend='topright');
      dofig(plotrate,'n1=40_plot',rate.rule=fignonz,d=0,n1=40,n2=seq(100,by=100,len=10),
            mesr=mesr,legend='topright');
      dofig(heatrate,'n1=20_heat',rate.rule=fignonz,d=0,n1=20,n2=seq(50,by=50,len=10),
            smooth='none',mesr=mesr);
      dofig(heatrate,'n1=40_heat',rate.rule=fignonz,d=0,n1=40,n2=seq(100,by=100,len=10),
            smooth='none',mesr=mesr);
    }
    if ((figsect='nonzro_exact_fnr') %in% sect) {
      ## exact reps. false negatives. show several values of d in one heatmap
      xdata=subset(expand.grid(n1=20,n2=seq(50,by=100,len=5),d1=c(0.2,0.5,0.8),d2=c(0.2,0.5,0.8)),
                   subset=d1==d2);
      dofig(heatrate,'n1=20_heat',rate.rule=fignonz,xdata=xdata,smooth='none',mesr=mesr);
      ## draw vertical lines separating the 'panels'
      ## TODO: move this into heatrate
      abline(v=c(5.5,10.5));
      ## do it again for n=40
      xdata=subset(expand.grid(n1=40,n2=seq(100,by=200,len=5),d1=c(0.2,0.5,0.8),d2=c(0.2,0.5,0.8)),
                   subset=d1==d2);
      dofig(heatrate,'n1=40_heat',rate.rule=fignonz,xdata=xdata,smooth='none',mesr=mesr);
      abline(v=c(5.5,10.5));
     }
    if ((figsect='nonzro_exact_roc') %in% sect) {
      ## show same point with roc
      xdata=xdata_blog(near=0);
      dofig(plotroc,rate.rule=fignonz,xdata=xdata,mesr=mesr);
    }
##### inexact
    if ((figsect='nonzro_inexact_roc') %in% sect) {
      ## inexact reps. nonzro. key point: false positives are terrible
      xdata=xdata_blog(near=1);
      dofig(plotroc,rate.rule=fignonz,xdata=xdata,mesr=mesr);
   }
##### near exact
    if ((figsect='nonzro_nearexact') %in% sect) {
      ## near exact. sig1
      near=round(c(0.01,0.05,0.1,0.2),digits=5);
      xdata=lapply(c(0,near,1),function(near) xdata=xdata_blog(near=near));
      names(xdata)=c('exact',paste(sep=' ','near',near),'inexact')
      dofig(plotrocm,'rocm',rate.rule=fignonz,xdata=xdata,title.desc='near exact',mesr='sig2');
      dofig(plotragm,'ragm',rate.rule=fignonz,xdata=xdata,title.desc='near exact',x='n2',
            rate='fpr',mesr='sig2');         
    }
########## sameff
    if ((figsect='sameff_nearexact') %in% sect) {
      xdata=xdata_blog(near=1);
      dofig(plotroc,'unconstrained_delta=0.1',rate.rule='sameff',rate.tol=0.1,xdata=xdata,
            title.desc='d1, d2 unconstrained',mesr=mesr);
      dofig(plotroc,'unconstrained_delta=0.5',rate.rule='sameff',rate.tol=0.5,xdata=xdata,
            title.desc='d1, d2 unconstrained',mesr=mesr);

      xdata=subset(xdata,subset=d2<=d1);
      dofig(plotroc,'d2<=d1_delta=0.1',rate.rule='sameff',rate.tol=0.1,xdata=xdata,
            title.desc='d2 <= d1',mesr=mesr);
      dofig(plotroc,'d2<=d1_delta=0.5',rate.rule='sameff',rate.tol=0.5,xdata=xdata,
            title.desc='d2 <= d1',mesr=mesr);
    }
##########
    if ((figsect='small_telescopes') %in% sect) {
      xdata=xdata_blog(near=0.1);
      dofig(plotroc,'nonzro_roc',rate.rule='nonzro',xdata=xdata,
            title.desc='near=0.1',mesr=cq(sig2,d2.scp1));
      dofig(plotrag,'nonzro_rag',rate.rule='nonzro',xdata=xdata,
            title.desc='near=0.1',mesr=cq(sig2,d2.scp1),smooth='loess');
      xdata=xdata_blog(near=1);
      xdata=subset(xdata,subset=d2<=d1);
      dofig(plotroc,'sameff_roc',rate.rule='sameff',rate.tol=0.5,xdata=xdata,
            title.desc='d2 <= d1',mesr=cq(sig2,d2.scp1));
      dofig(plotrag,'sameff_rag',rate.rule='sameff',rate.tol=0.5,xdata=xdata,
            title.desc='d2 <= d1',mesr=cq(sig2,d2.scp1),smooth='loess');
    }
    sect;
  }
## generate standard xdata for aggregated plots
xdata_blog=function(near=0,nx=2.5,n2.num=2,n1=seq(20,by=20,len=8)) {
    d2=round(seq(0,1,by=0.01),digits=5);
    do.call(rbind,lapply(n1,function(n1) {
    n2=seq(n1*nx,by=n1*nx,len=n2.num);
    xdata=expand.grid(n1=n1,n2=n2,d1=d,d2=d2);
    ## TODO: I don't thing the 1e-4 tolerance still needed
    subset(xdata,subset=(abs(d1-d2)<=(near+1e-4)));}))
}
