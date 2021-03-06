################################################################################
##
## Author:  Nat Goodman
## Created: 18-10-18
##          from doc_resig.R created 18-09-05
##          by renaming doc_repwr.R created 18-06-20
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
## Generate figures and tables for resig.R blog post supplement
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
source('R/docfun_resigsupp.R');
## --- Generate Figures and Tables for resig Blog Post supplement ---
## sect is which sections to run - for use during development
##   uses prefix matching and runs all that match
doc_resigsupp=function(sect=parent(sect,NULL)) {
  sect.all=cq(exact_fpr,exact_fnr,inexact,nearexact);
  sect.desc.all=
    setNames(c('Exact replication','Exact replication',
               'Inexact replication',
               'Near exact replication'),sect.all);
  
  if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all);
  ## doc-wide variables
  d.nonzro=d[d!=0];
  n2=n[n>=50];                   # standard n2: exclude 20
  d2col=color_d();               # colors keyed by d
  n2col=color_n();               # colors keyed by n
  if (m==1e2) rate.empty<<-1;    # with m=1e2, many agg cases empty. set to 1 so don't bomb
  sapply(sect,function(sect) {
    if (!is.null(sectnum)) {
      ## compute section number. from stackoverflow.com/questions/5577727
      sectnum=which(sect==sect.all)[1];
      fignum<<-1; xfignum<<-1;
    }
    sect.desc=sect.desc.all[sect];
    extra=F;                         # start each section with extra off 
    figblk<<-NULL; xfigblk<<-NULL;   #   and blk off
    ## exact fpr
    if (sect=='exact_fpr') {
      ## these figures show the main point: FPR=sig.level/2 independent of n1, n2
      ## the factor of 2 is due to sdir
      xdata=xdata_exact(n1=n,n2=n2,d=0,by='n1');
      dofig(plotratm,'by_n1',xdata=xdata,x=cq(n2,d1,d2),col=n2col,smooth='spline',
            hline=c(fpr.cutoff/2,fpr.cutoff,fnr.cutoff),
            vhlty='dashed',vhlwd=c(1,0.5,0.5),vhcol=cq(red,black,black),plot.cutoff=F,
            title=title_resigsupp,title.rate='fpr',title.legend='n1',legend='topright');
      dofig(plotratm,'by_n1_nosdir',xdata=xdata,x=cq(n2,d1,d2),col=n2col,smooth='spline',
            hline=c(fpr.cutoff,fpr.cutoff,fnr.cutoff),posr.id='sig1_sig1',
            vhlty='dashed',vhlwd=c(1,0.5,0.5),vhcol=cq(red,black,black),plot.cutoff=F,
            title=title_resigsupp,title.rate='fpr',title.legend='n1',legend='topright');
      ## these tables further confirm the main point: FPR=sig.level/2 independent of n1, n2
      ## the factor of 2 is due to sdir
      drat.std=data_rate(xdata=expand.grid(n1=n,n2=n,d1=0,d2=0),mesr=cq(sig2))
      drat.nosdir=data_rate(xdata=expand.grid(n1=n,n2=n,d1=0,d2=0),mesr=cq(sig2),
                            posr.id='sig1_sig1');
      fpr.std=mean(drat.std$sig2);
      fpr.nosdir=mean(drat.nosdir$sig2);
      fpr=data.frame(fpr.std,fpr.nosdir);
      ## hmm... fpr.std=0.02583352. a bit high. what's going on?
      ## no correlation with n1, n2, mean(n1,n2)
      corfprstd=data.frame(cor.n1=with(drat.std,cor(n1,sig2)),
                           cor.n2=with(drat.std,cor(n2,sig2)),
                           cor.n1n2=cor(rowMeans(drat.std[,cq(n1,n2)]),drat.std$sig2));
      dotbl(fpr,corfprstd);
      
      ## remaining figures are extra
      extra=T;
      ## these two extra figures reprise Figure 1 of the blog post with n1=20 and n1=200
      xfigblk_start();
      dofig(plotrate,'n1=020',d=0,n1=20,n2=n2,smooth='spline',
            hline=c(fpr.cutoff/2,fpr.cutoff,fnr.cutoff),
            vhlty='dashed',vhlwd=c(1,0.5,0.5),vhcol=cq(red,black,black),plot.cutoff=F,
            title=title_resigsupp,title.rate='fpr',legend=NULL);
      dofig(plotrate,'fpr_n1=200',d=0,n1=200,n2=n2,smooth='spline',
            hline=c(fpr.cutoff/2,fpr.cutoff,fnr.cutoff),
            vhlty='dashed',vhlwd=c(1,0.5,0.5),vhcol=cq(red,black,black),plot.cutoff=F,
            title=title_resigsupp,title.rate='fpr',legend=NULL);
      ## these two extra figures repeat the analysis of the two 'real' figures with
      ## n1 and n2 depicted in opposite roles (n1 on the x-axis, n2 different lines)
      xfigblk_start();
      xdata=xdata_exact(n1=n,n2=n2,d=0,by='n2');
      dofig(plotratm,'by_n2',xdata=xdata,x=cq(n1,d1,d2),col=n2col,smooth='spline',
            hline=c(fpr.cutoff/2,fpr.cutoff,fnr.cutoff),
            vhlty='dashed',vhlwd=c(1,0.5,0.5),vhcol=cq(red,black,black),plot.cutoff=F,
            title=title_resigsupp,title.rate='fpr',title.legend='n2',legend='topright');
      dofig(plotratm,'by_n2_nosdir',xdata=xdata,x=cq(n1,d1,d2),col=n2col,smooth='spline',
            hline=c(fpr.cutoff,fpr.cutoff,fnr.cutoff),posr.id='sig1_sig1',
            vhlty='dashed',vhlwd=c(1,0.5,0.5),vhcol=cq(red,black,black),plot.cutoff=F,
            title=title_resigsupp,title.rate='fpr',title.legend='n2',legend='topright');
      ## boxplots also show no obvious correlation with n1, n2, or mean(n1,n2).
      xfigblk_start();
      boxlim=range(c(range(drat.std$sig2,na.rm=T),range(drat.nosdir$sig2,na.rm=T)));
      dofig(plotboxfpr_exact,'box',extra=T,drat=drat.std,posr.id='std',ylim=boxlim);
      dofig(plotboxfpr_exact,'box_nosdir',extra=T,drat=drat.nosdir,posr.id='sig1_sig1',
            ylim=boxlim);
    }
    ## exact fnr
    if (sect=='exact_fnr') {
      ## compute 1-power2 vs. n2 - specialized for plotfnr!
      power.n2=power_n2();
      ## construct xdata lists for n1=20 and n1=200
      xdata.020=lapply(d.nonzro,function(d) xdata=expand.grid(n1=20,n2=n2,d1=d,d2=d));
      names(xdata.020)=as.character(d.nonzro);
      xdata.200=lapply(d.nonzro,function(d) xdata=expand.grid(n1=200,n2=n2,d1=d,d2=d));
      names(xdata.200)=as.character(d.nonzro);
      ## plots using std posr - includes sdir
      figblk_start();
      dofig(plotfnr_exact,'n1=020',xdata=xdata.020,col=d2col,power.n2=power.n2)
      dofig(plotfnr_exact,'n1=200',xdata=xdata.200,col=d2col,power.n2=power.n2);
      ## plots using sig1_sig1 posr - omits sdir
      figblk_start();
      dofig(plotfnr_exact,'n1=020_nosdir',xdata=xdata.020,col=d2col,power.n2=power.n2,
            posr.id='sig1_sig1')
      dofig(plotfnr_exact,'n1=200_nosdir',xdata=xdata.200,col=d2col,power.n2=power.n2,
            posr.id='sig1_sig1');
      ## plot fnr vs n1, n2 aggregating over d
      figblk_start();
      xdata=xdata_exact(n1=n,n2=n2,d=d,by='n1');
      dofig(plotragm,'mean_d',xdata=xdata,x='n2',rate='fnr',col=n2col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate='fnr',title.desc='vs. n2,n1, agg d',
            title.legend='n1',legend='topright');
      dofig(plotfnr_exact_bias1,'mean_d_bias1',xdata=xdata,col=n2col);

      ## remaining figures are extra
      extra=T;
      ## These two extra figures demonstrate that sig1 bias is not due to sdir
      xfigblk_start();
      dofig(plotragm,'mean_d_nosdir',xdata=xdata,x='n2',rate='fnr',col=n2col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate='fnr',title.desc='vs. n2,n1, agg d',
            title.legend='n1',legend='topright',posr.id='sig1_sig1');
      dofig(plotfnr_exact_bias1,'mean_d_bias1_nosdir',xdata=xdata,col=n2col,posr.id='bias1');
      ## These two extra figures plot FNR vs. $1-power2$ across the entire dataset
      ## with same-direction switched on (first figure) and off (second figure). The
      ## dots are color-coded by power1, the power of the original study
      ## (smaller power has lighter dots). The points fall nicely along the
      ## diagonal in both cases, with less scatter without same-direction, as
      ## expected.
      xfigblk_start();
      xdata=expand.grid(n1=n,n2=n,d1=d.nonzro,d2=d);
      xdata=subset(xdata,subset=(d1==d2));
      dofig(plotfnrpwr_exact,'fnr_vs_power',xdata=xdata)
      dofig(plotfnrpwr_exact,'fnr_vs_power_nosdir',xdata=xdata,posr.id='sig1_sig1');
      ## The correlation of FNR and $1-power2$ across the entire dataset is
      ## 0.9997737 and 0.9999019 resp. with same-direction on and off.  When
      ## limiting to small values of $d2_{pop}$, $d2_{pop}\le0.2$, the
      ## correlations are 0.9989049 and 0.9995707.
      drat.std=dratfnrpwr_exact(xdata,posr.id='std');
      drat.nosdir=dratfnrpwr_exact(xdata,posr.id='sig1_sig1');
      corfnrpwr=do.call(rbind,lapply(d.nonzro,function(d) {
        drat.std=subset(drat.std,subset=(d2<=d));
        drat.nosdir=subset(drat.nosdir,subset=(d2<=d));
        cor.std=with(drat.std,cor(sig2,1-pwr2));
        cor.nosdir=with(drat.nosdir,cor(sig2,1-pwr2));
        data.frame(d2=d,cor.std,cor.nosdir);}))
      ## tables showing n2,d2 value that achieve 'cutoff' FNR
      fnr_n2byd2.020=cutoff_n2byd2(xdata.020);
      fnr_n2byd2.200=cutoff_n2byd2(xdata.200);
      fnr_n2byd2=merge(fnr_n2byd2.020,fnr_n2byd2.200,by=cq(cutoff,d2),suffixes=c('.020','.200'));
      fnr_n2byd2=fnr_n2byd2[with(fnr_n2byd2,order(cutoff,d2)),];
      ##
      fnr_d2byn2.020=cutoff_d2byn2(xdata.020);
      fnr_d2byn2.200=cutoff_d2byn2(xdata.200);
      fnr_d2byn2=merge(fnr_d2byn2.020,fnr_d2byn2.200,by=cq(cutoff,n2),suffixes=c('.020','.200'));
      fnr_d2byn2=fnr_d2byn2[with(fnr_d2byn2,order(cutoff,n2)),];
      dotbl(corfnrpwr,fnr_n2byd2,fnr_d2byn2);
      ## The next extra figure shows the values of $n2$ and $d2_{pop}$ that achieve FNR
      ## cutoffs of 0.05 and 0.20 for $n1=20$ and $n1=200$.
      xfigblk_end();                    # not xfigblk_start 'cuz there's only one of 'em
      dofig(plotfnrcutoff_exact,'fnr_cutoff',fnr_d2byn2=fnr_d2byn2);
    }
    ## inexact
    if (sect=='inexact') {
      ## fpr
      figblk_start();
      xdata.020=xdata_inexact(n1=20,d1=0);
      xdata.200=xdata_inexact(n1=200,d1=0);
      dofig(plotratm,'fpr_n1=020',xdata=xdata.020,x=cq(n1,n2,d1),col=d2col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate='fpr',title.legend='d2',legend='topright');
      dofig(plotratm,'fpr_n1=200',xdata=xdata.200,x=cq(n1,n2,d1),col=d2col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate='fpr',title.legend='d2',legend='topright');
      ## support statement: data are similar but not exactly the same.
      corfprdata=drat_cor(xdata.020,xdata.200);
      dotbl(corfprdata);
      ## extra plot shows similarity
      dofig(drat_plotsig2,'fpr_sigsig',xdata1=xdata.020,xdata2=xdata.200,extra=T);
      ## fnr
      xdata.020=xdata_inexact(n1=20,d1=0.5);
      xdata.200=xdata_inexact(n1=200,d1=0.5);
      dofig(plotratm,'fnr_n1=020',xdata=xdata.020,x=cq(n1,n2,d1),col=d2col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate='fnr',title.legend='d2',x.legend=8.45,y.legend=0.65);
      dofig(plotratm,'fnr_n1=200',xdata=xdata.020,x=cq(n1,n2,d1),col=d2col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate='fnr',title.legend='d2',x.legend=8.45,y.legend=0.65);
      ## support statement: data are similar but not exactly the same.
      corfnrdata=drat_cor(xdata.020,xdata.200);
      dotbl(corfnrdata);
      ## extra plot shows similarity
      dofig(drat_plotsig2,'fnr_sigsig',xdata1=xdata.020,xdata2=xdata.200,extra=T);
      ## fpr+fnr
      figblk_start();
      d2=c(0,0.1,0.2,0.5,1)
      xdata.fpr=xdata_inexact(n1=20,d1=0,d2=d2);
      xdata.fnr=xdata_inexact(n1=20,d1=0.5,d2=d2);
      xdata.020=xdata_rbind(xdata.fpr,xdata.fnr);
      dofig(plotragm,'fpr+fnr_n1=020',xdata=xdata.020,x=cq(n1,n2),col=d2col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate=cq(fpr,fnr),
            title.legend='d2',x.legend=9.1,y.legend=0.985,cex.legend=0.7);
      xdata.fpr=xdata_inexact(n1=200,d1=0,d2=d2);
      xdata.fnr=xdata_inexact(n1=200,d1=0.5,d2=d2);
      xdata.200=xdata_rbind(xdata.fpr,xdata.fnr);
      dofig(plotragm,'fpr+fnr_n1=200',xdata=xdata.200,x=cq(n1,n2),col=d2col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate=cq(fpr,fnr),
            title.legend='d2',x.legend=9.1,y.legend=0.985,cex.legend=0.7);
      ## extra plot shows similarity
      dofig(drat_plotsig2,'fpr+fnr_sigsig',xdata1=xdata.020,xdata2=xdata.200,extra=T);
      ## obvious rocm
      figblk_end();
      xdata=lapply(d,function(d2) xdata=expand.grid(n1=c(20,200),n2=n2,d1=d,d2=d2));
      names(xdata)=as.character(d);
      dofig(plotrocm,'rocm',xdata=xdata,x=cq(n1,n2),col=d2col,
            hline=c(fpr.cutoff,fnr.cutoff),vline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',
            vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,
            title.rate=NULL,title.desc='False negative vs. false positive rate',title.legend='d2');
    }
    ## nearexact
    if (sect=='nearexact') {
      ## sect.desc='Near exact replication';
      ## near=round(c(0,0.01,0.05,0.1,0.2),digits=5);
      near=seq(0.1,0.5,by=0.1);
      ## fpr
      figblk_start();
      xdata.020=xdata_near(n1=20,n2=n2,d1=0,near=near);
      xdata.200=xdata_near(n1=200,n2=n2,d1=0,near=near);
      dofig(plotragm,'fpr_n1=020',xdata=xdata.020,x=cq(n1,n2,d1),col=d2col,smooth='spline',
            rate='fpr',hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate='fpr',title.legend='near',legend='topright');
      dofig(plotragm,'fpr_n1=200',xdata=xdata.200,x=cq(n1,n2,d1),col=d2col,smooth='spline',
            rate='fpr',hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate='fpr',title.legend='near',legend='topright');
      ## support statement: data are similar but not exactly the same.
      corfprdata=drat_cor(xdata.020,xdata.200);
      dotbl(corfprdata);
      ## extra plot shows similarity
      dofig(drat_plotsig2,'fpr_sigsig',xdata1=xdata.020,xdata2=xdata.200,extra=T);
      ## fnr
      xdata.020=xdata_near(n1=20,n2=n2,d1=0.5,near=near);
      xdata.200=xdata_near(n1=200,n2=n2,d1=0.5,near=near);
      dofig(plotragm,'fnr_n1=020',xdata=xdata.020,x=cq(n1,n2,d1),col=d2col,smooth='spline',
            rate='fnr',hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate='fnr',title.legend='near',legend='topright');
      dofig(plotragm,'fnr_n1=200',xdata=xdata.200,x=cq(n1,n2,d1),col=d2col,smooth='spline',
            rate='fnr',hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate='fnr',title.legend='near',legend='topright');
       ## support statement: data are similar but not exactly the same.
      corfnrdata=drat_cor(xdata.020,xdata.200);
      dotbl(corfnrdata);
      ## extra plot shows similarity
      dofig(drat_plotsig2,'fnr_sigsig',xdata1=xdata.020,xdata2=xdata.200,extra=T);
      ## these extra figures plot fnr for single values of n1, n2, multiple d1
      extra=T; xfigblk_start();
      n1=20;
      sapply(c(100,200,300,400),function(n2x) {
        xdata=xdata_near(n1=n1,n2=n2x,d1=d.nonzro,near=near);
        figname=paste(sep='_','fnr',paste_nv(n1),paste_nv(n2,n2x));
        dofig(plotragm,figname,xdata=xdata,x=cq(n1,n2,d1),rate='fnr',col=d2col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate='fnr',title.legend='near',legend='topright');
      });
      extra=F;
      ## fpr+fnr
      figblk_start();
      xdata.fpr=xdata_near(n1=20,n2=n2,d1=0,near=near);
      xdata.fnr=xdata_near(n1=20,n2=n2,d1=0.5,near=near);
      xdata.020=xdata_rbind(xdata.fpr,xdata.fnr);
      dofig(plotragm,'fpr+fnr_n1=020',xdata=xdata.020,x=cq(n1,n2),col=d2col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate=cq(fpr,fnr),title.legend='near',legend='topright');
      xdata.fpr=xdata_near(n1=200,n2=n2,d1=0,near=near);
      xdata.fnr=xdata_near(n1=200,n2=n2,d1=0.5,near=near);
      xdata.200=xdata_rbind(xdata.fpr,xdata.fnr);
      dofig(plotragm,'fpr+fnr_n1=200',xdata=xdata.200,x=cq(n1,n2),col=d2col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate=cq(fpr,fnr),title.legend='near',legend='topright');
      ## extra plot shows similarity
      dofig(drat_plotsig2,'fpr+fnr_sigsig',xdata1=xdata.020,xdata2=xdata.200,extra=T);
      ## these extra figures do it systematically for all values of near, several values of n1, d1
      extra=T; xfigblk_start();
      near=d;
      n1=c(20,200);
      d1=c(0.2,0.5,0.8);
      ## CAUTION: must use loop, NOT sapply, for scoping of fignum to work!
      ## Hmmm... not sure it makes sense to do all these figures
      ## for (n1x in c(20,200)) {
      ##   for (d1x in c(0.2,0.5,0.8)) {
      sapply(c(20,200),function(n1x) {
        sapply(c(0.2,0.5,0.8),function(d1x) {
          xdata.fpr=xdata_near(n1=n1x,n2=n2,d1=0,near=near);
          xdata.fnr=xdata_near(n1=n1x,n2=n2,d1=d1x,near=near);
          xdata=xdata_rbind(xdata.fpr,xdata.fnr);
          figname=paste(sep='_','fpr+fnr',paste_nv(n1,n1x),paste_nv(d1,d1x));
          dofig(plotragm,figname,xdata=xdata,x=cq(n1,n2),col=d2col,smooth='spline',
                hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
                title=title_resigsupp,title.rate=cq(fpr,fnr),title.desc=paste_nv(d1,d1x),
                title.legend='near',legend='topright');
        })})
      extra=F;
      ## obvious rocm
      figblk_end();
      near=d;
      xdata=xdata_near(n1=c(20,200),n2=n2,d1=d,near=near);
      dofig(plotrocm,'rocm',xdata=xdata,x=cq(n1,n2),col=d2col,
            hline=c(fpr.cutoff,fnr.cutoff),vline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',
            vhlwd=0.5,plot.cutoff=F,
            title=title_resigsupp,title.rate=NULL,
            title.desc='False negative vs. false positive rate',title.legend='near');
    }
    sect;
  })
}
