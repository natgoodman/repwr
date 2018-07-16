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
## Generate figures and tables for repwr.R technical note
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## --- Technical Note ---
## make figures and tables for documents
## docs are blog, tech, readme
## sect is which sections to run - for use during development
##   uses prefix matching and all matches run
dotech=
  function(sect=NULL,fignum=1,doc='tech',mesr=cq(sig2,sigm,d1.c2,d2.c1,c1.c2,d1.p2,d2.p1,p1.p2),
           fignew=F) {
    sect.all=cq(exact_fpr,exact_fnr,exact_keypoint,
                inexact_nonzro_fpr,inexact_nonzro_fnr,inexact_nonzro_roc,
                inexact_sameff_fpr,inexact_sameff_fnr,inexact_sameff_roc,
                nearexact_nonzro_fpr,nearexact_nonzro_fnr,nearexact_nonzro_roc,
                last_placeholder);
    if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all);
 #################### exact
   if ((figsect='exact_fpr') %in% sect) {
      ## exact reps. false positives
      ## 2x2 panel for d=0
     dofig(plotrate,'n1=20_plot',rate.rule='nonzro',d=0,n1=20,n2=seq(50,by=50,len=10),
            mesr=mesr,legend='topright');
      dofig(heatrate,'n1=20_heat',rate.rule='nonzro',d=0,n1=20,n2=seq(50,by=50,len=10),
            mesr=mesr);
      dofig(plotrate,'n1=40_plot',rate.rule='nonzro',d=0,n1=40,n2=seq(100,by=100,len=10),
            mesr=mesr,legend='topright');
      dofig(heatrate,'n1=40_heat',rate.rule='nonzro',d=0,n1=40,n2=seq(100,by=100,len=10),
            mesr=mesr);
    }
    if ((figsect='exact_fnr') %in% sect) {
      ## exact reps. false negatives
      ## 2x2 panel for d=0.5.
      dofig(plotrate,'n1=20_plot',rate.rule='nonzro',d=0.5,n1=20,
            n2=seq(50,by=50,len=10),mesr=mesr,legend='right');
      dofig(heatrate,'n1=20_heat',rate.rule='nonzro',d=0.5,n1=20,
            n2=seq(50,by=50,len=10),mesr=mesr);
      dofig(plotrate,'n1=40_plot',rate.rule='nonzro',d=0.5,n1=40,
            n2=seq(100,by=100,len=10),mesr=mesr,legend='right');
      dofig(heatrate,'n1=40_heat',rate.rule='nonzro',d=0.5,n1=40,
            n2=seq(100,by=100,len=10),mesr=mesr);
    }
    if ((figsect='exact_key') %in% sect) {
      ## exact reps. key point: sig2 is all that works
      ## heat fpr & fnr side-by-side without smoothing
      dofig(heatrate,'d=05_heat',rate.rule='nonzro',d=0.5,n1=20,
            n2=seq(50,by=50,len=10),smooth='none',mesr=mesr);
      dofig(heatrate,'d=00_heat',rate.rule='nonzro',d=0,n1=20,
            n2=seq(50,by=50,len=10),smooth='none',mesr=mesr);
      dofig(heatrate,'d=08_heat',rate.rule='nonzro',d=0.8,n1=20,
            n2=seq(50,by=50,len=10),smooth='none',mesr=mesr);
    ## show same point with roc
    xdata=do.call(rbind,lapply(c(20,40,80,160),function(n1) {
      n2=seq(n1*2.5,by=n1*2.5,len=10);
      xdata=expand.grid(n1=n1,n2=n2,d1=d,d2=d);
      subset(xdata,subset=(d1==d2))}));
    dofig(plotroc,'roc',rate.rule='nonzro',xdata=xdata,x='n2',mesr=mesr);
    }
#################### inexact
########## nonzro
    if ((figsect='inexact_nonzro_fpr') %in% sect) {
      ## inexact reps. nonzro. false positives. key point: these are terrible
      dofig(plotrate,'n1=20_plot',rate.rule='nonzro',d1=0,d2=0.5,n1=20,n2=seq(50,by=50,len=10),
            mesr=mesr,legend='topright');
      dofig(plotrate,'n1=40_plot',rate.rule='nonzro',d1=0,d2=0.5,n1=40,n2=seq(100,by=100,len=10),
            mesr=mesr,legend='topright');
      dofig(plotrate,'n1=20_d2=d_plot',rate.rule='nonzro',d1=0,d2=d,n1=20,n2=50,
            mesr=mesr,legend='topright');
      dofig(plotrate,'n1=40_d2=d_plot',rate.rule='nonzro',d1=0,d2=d,n1=40,n2=100,
            mesr=mesr,legend='topright');
      }
    if ((figsect='inexact_nonzro_fnr') %in% sect) {
      ## false negatives. these are okay but doesn't matter since FPR so bad
      dofig(plotrate,'n1=20_plot',rate.rule='nonzro',d1=0.2,d2=0.5,n1=20,
            n2=seq(50,by=50,len=10),mesr=mesr,legend='topright');
      dofig(plotrate,'n1=40_plot',rate.rule='nonzro',d1=0.2,d2=0.5,n1=40,
            n2=seq(100,by=100,len=10),mesr=mesr,legend='topright');
      dofig(plotrate,'n1=20_d2=d_plot',rate.rule='nonzro',d1=0.2,d2=d,n1=20,
            n2=50,mesr=mesr,legend='topright');
      dofig(plotrate,'n1=40_d2=d_plot',rate.rule='nonzro',d1=0.2,d2=d,n1=40,
            n2=100,mesr=mesr,legend='topright');
      dofig(plotrate,'n1=20_d=d_plot',rate.rule='nonzro',
            d1=d[2:11],d2=d[2:11],n1=20,n2=100,mesr=mesr,legend='topright');
      dofig(plotrate,'n1=40_d=d_plot',rate.rule='nonzro',
            d1=d[2:11],d2=d[2:11],n1=40,n2=200,mesr=mesr,legend='topright');
    }
    if ((figsect='inexact_nonzro_roc') %in% sect) {
      ## look at it with ROC. fiddling with xdata doesn't change answer much
      xdata=do.call(rbind,lapply(c(20,40,80,160),function(n1) {
        n2=seq(n1*2.5,by=n1*2.5,len=10);
        xdata=expand.grid(n1=n1,n2=n2,d1=d,d2=d)}));
      dofig(plotroc,rate.rule='nonzro',xdata=xdata,x='n2',mesr=mesr)
      ## this one is better, but performance still unacceptable
      xdata=do.call(rbind,lapply(c(20,40,80,160),function(n1) {
        n2=seq(n1*2.5,by=n1*2.5,len=10);
        xdata=expand.grid(n1=n1,n2=n2,d1=d,d2=d);
        subset(xdata,subset=(d1<=d2))}));
      dofig(plotroc,rate.rule='nonzro',xdata=xdata,x='n2',mesr=mesr)
    }
#################### inexact
########## sameff
    if ((figsect='inexact_sameff_fpr') %in% sect) {
      ## false positives. use some subset of these. all bad
      dofig(plotrate,rate.rule='sameff',d1=0,d2=0.5,n1=20,n2=seq(50,by=50,len=10),
            mesr=mesr,legend='topright');
      dofig(plotrate,rate.rule='sameff',d1=0,d2=0.5,n1=40,n2=seq(100,by=100,len=10),
            mesr=mesr,legend='topright');
      dofig(plotrate,rate.rule='sameff',d1=0,d2=d[2:11],n1=20,n2=50,mesr=mesr,legend='topright');
      dofig(plotrate,rate.rule='sameff',d1=0,d2=d[2:11],n1=40,n2=100,mesr=mesr,legend='topright');

      dofig(plotrate,rate.rule='sameff',d1=d[1:10],d2=d[2:11],n1=20,n2=100,
            mesr=mesr,legend='topright');
      dofig(plotrate,rate.rule='sameff',d1=d[1:10],d2=d[2:11],n1=40,n2=100,
            mesr=mesr,legend='topright');
      dofig(plotrate,rate.rule='sameff',d1=d[1:9],d2=d[3:11],n1=40,n2=100,
            mesr=mesr,legend='topright');
      dofig(plotrate,rate.rule='sameff',d1=d[1:8],d2=d[4:11],n1=40,n2=100,
            mesr=mesr,legend='topright');
      ddd5=d[d!=0.5];
      dofig(plotrate,rate.rule='sameff',d1=ddd5,d2=0.5,n1=40,n2=100,
            mesr=mesr,legend='topright',xtitle='n');
      tol=.1;
      d2=0.5; d1=d[abs(d-d2)>tol];
      dofig(plotrate,rate.rule='nearff',rate.tol=tol,d1=d1,d2=d2,n1=40,n2=100,
            mesr=mesr,legend='topright',xtitle='n');
      d2=0.2; d1=d[abs(d-d2)>rate.tol];
      dofig(plotrate,rate.rule='nearff',rate.tol=tol,d1=d1,d2=d2,n1=40,n2=100,
            mesr=mesr,legend='topright',xtitle='n');
      }
    if ((figsect='inexact_sameff_fnr') %in% sect) {
      ## false negatives. okay but doesn't matter since FPR so bad
      dofig(plotrate,rate.rule='sameff',d1=d,d2=d,n1=20,n2=50,mesr=mesr,legend='topright');
      dofig(plotrate,rate.rule='sameff',d1=d,d2=d,n1=40,n2=100,mesr=mesr,legend='topright');
      tol=.1;
      d2=c(0.2,0.5,0.8);
      xdata=expand.grid(d1=d,d2=d2);
      xdata=subset(xdata,subset=abs(d1-d2)<=tol); xdata$n1=40; xdata$n2=100;
      dofig(plotrate,rate.rule='nearff',rate.tol=tol,xdata=xdata,mesr=mesr,legend='topright');
    }
    if ((figsect='inexact_sameff_roc') %in% sect) {
      ## look at it with ROC. fiddling with xdata doesn't change answer much
      xdata=do.call(rbind,lapply(c(20,40,80,160),function(n1) {n2=seq(n1*2.5,by=n1*2.5,len=10);
        xdata=expand.grid(n1=n1,n2=n2,d1=d,d2=d)}));
      dofig(plotroc,rate.rule='sameff',xdata=xdata,x='n2',mesr=mesr);
      xdata=do.call(rbind,lapply(c(20,40,80,160),function(n1) {n2=seq(n1*2.5,by=n1*2.5,len=10);
        xdata=expand.grid(n1=n1,n2=n2,d1=d,d2=d);
        subset(xdata,subset=(d1<=d2))}))
      dofig(plotroc,rate.rule='sameff',xdata=xdata,x='n2',mesr=mesr);
      ## using rate.tol helps a bit but still not acceptable
      xdata=do.call(rbind,lapply(c(20,40,80,160),function(n1) {n2=seq(n1*2.5,by=n1*2.5,len=10);
        xdata=expand.grid(n1=n1,n2=n2,d1=d,d2=d)}))
      dofig(plotroc,rate.rule='nearff',rate.tol=0.3,xdata=xdata,x='n2',mesr=mesr);
      dofig(plotroc,rate.rule='nearff',rate.tol=0.5,xdata=xdata,x='n2',mesr=mesr);
    }
#################### near exact
########## nonzro
    if ((figsect='nearexact_nonzro_fpr') %in% sect) {
      ## intuitively sig2 should still be good. actually, maybe not...
      ## some examples to get started
      dofig(plotrate,rate.rule='nonzro',d1=0,d2=0,n1=20,n2=seq(50,by=50,len=10),
            mesr=mesr,legend='topright');
      dofig(plotrate,rate.rule='nonzro',d1=0,d2=0,n1=20,n2=seq(50,by=50,len=10),
            mesr=mesr,legend='topright');
      dofig(plotrate,rate.rule='nonzro',d1=0,d2=0.1,n1=20,n2=seq(50,by=50,len=10),
            mesr=mesr,legend='topright');
      dofig(plotrate,rate.rule='nonzro',d1=0,d2=0.05,n1=20,n2=seq(50,by=50,len=10),
            mesr=mesr,legend='topright');
      dofig(plotrate,rate.rule='nonzro',d1=0,d2=0.1,n1=20,n2=seq(10,by=10,len=10),
            mesr=mesr,legend='topright');
      dofig(plotrate,rate.rule='nonzro',d1=0,d2=0.05,n1=20,n2=seq(10,by=10,len=10),
            mesr=mesr,legend='topright');
 
      dofig(plotrate,rate.rule='nonzro',d1=0,d2=seq(0,.1,by=.01),n1=20,n2=50,
            mesr=mesr,legend='topright');
      dofig(plotrate,rate.rule='nonzro',d1=0,d2=seq(0,.1,by=.01),n1=20,n2=70,
            mesr=mesr,legend='topright');
    }
    if ((figsect='nearexact_nonzro_roc') %in% sect) {
      ## this ROC shows sig2 almost works
      xdata=do.call(rbind,lapply(c(20,40,80,160),function(n1) {n2=seq(n1*2.5,by=n1*2.5,len=10);
        xdata=expand.grid(n1=n1,n2=n2,d1=d,d2=d);
        subset(xdata,subset=(abs(d1-d2)<=0.1))}));
      dofig(plotroc,rate.rule='nonzro',xdata=xdata,x='n2',mesr=mesr);
      ## this one shows sig2 works in some cases
      xdata=do.call(rbind,lapply(c(20,40,80,160),function(n1) {n2=seq(10,n1*5,len=10);
        xdata=expand.grid(n1=n1,n2=n2,d1=d,d2=d);
        subset(xdata,subset=(abs(d1-d2)<=0.1))}));
      dofig(plotroc,rate.rule='nonzro',xdata=xdata,x='n2',mesr=mesr);
      ## this one is a bit better
      xdata=do.call(rbind,lapply(c(20,40,80,160),function(n1) {n2=seq(n1,n1*2.5,len=10);
        xdata=expand.grid(n1=n1,n2=n2,d1=d,d2=d);
        subset(xdata,subset=(abs(d1-d2)<=0.1))}));
      dofig(plotroc,rate.rule='nonzro',xdata=xdata,x='n2',mesr=mesr);
      ## a bit more better
      xdata=do.call(rbind,lapply(c(80,160),function(n1) {n2=seq(n1,n1*3.5,len=10);
        xdata=expand.grid(n1=n1,n2=n2,d1=d,d2=d);
        subset(xdata,subset=(abs(d1-d2)<=0.1))}));
      dofig(plotroc,rate.rule='nonzro',xdata=xdata,x='n2',mesr=mesr);
      
    }
    sect;
  }
