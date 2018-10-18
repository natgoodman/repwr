#################################################################################
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
## --- Generate Figures and Tables for resig Blog Post supplement ___
## sect is which sections to run - for use during development
##   uses prefix matching and runs all that match
doc_resigsupp=function(sect=parent(sect,NULL)) {
  sect.all=cq(exact,inexact,nearexact);
  if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all);
  ## doc-wide variables
  d.nonzro=d[d!=0];
  n2=n[n>=50];                       # standard n2: exclude 20
  col=col_resig();                   # colors for plotratm, plotragm
  sapply(sect,function(sect) {
    ## compute section number. from stackoverflow.com/questions/5577727
    sectnum=if(sectnum) which(sect==sect.all)[1] else NULL;
    ## exact
    if (sect=='exact') {
      title.desc='Exact replication';
      dofig(plotrate,'fpr',d=0,n1=20,n2=n2,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fpr'),legend=NULL);
      xdata=lapply(d.nonzro,function(d) xdata=expand.grid(n1=20,n2=n2,d1=d,d2=d));
      names(xdata)=as.character(d.nonzro);
      dofig(plotratm,'fnr',xdata=xdata,x=cq(n1,n2),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fnr'),title.legend='d',legend='right');
      ## support statements: FPR=sig.level/2, FNR=1-power
      drat=data_rate(rate.rule='nonzro',xdata=expand.grid(n1=n,n2=n,d1=d,d2=d),mesr=cq(sig2));
      drat.exact=subset(drat,subset=d1==d2);
      drat.f=subset(drat.exact,subset=!true.dd);
      fpr_vs_siglevel=mean(drat.f$sig2-sig.level/2);
      drat.t=subset(drat.exact,subset=true.dd);
      drat.t$power=with(drat.t,power.t.test(n=n2,delta=d2)$power);
      fnr_vs_power=with(drat.t,cor(sig2,1-power));
      theory=t(setNames(c(fpr_vs_siglevel,fnr_vs_power),c('fpr-sig.level/2','cor(fnr,1-power)')));
      ## support final para: if d is small, n2 must be big...
      ## tables show n2,d2 value that achieve 'cutoff' FNR
      fnr_n2byd2=cutoff_n2byd2(xdata);
      fnr_d2byn2=cutoff_d2byn2(xdata);
      ## theoretical calculation of n2 value that achieve 'cutoff' FNR
      n2.20=power.t.test(delta=0.2,power=0.80)$n;
      n2.05=power.t.test(delta=0.2,power=0.95)$n;
      fnr_power=t(setNames(c(n2.20,n2.05),cq(n2.20,n2.05)));
      dotbl(theory,fnr_n2byd2,fnr_d2byn2,fnr_power);
    }
##### near exact
    if ((figsect='nearexact') %in% sect) {
      title.desc='Near exact replication';
      ## near=round(c(0,0.01,0.05,0.1,0.2),digits=5);
      near=seq(0,0.4,by=0.1);
      xdata=xdata_near(n1=20,n2=n2,d1=0,near=near);
      dofig(plotragm,'fpr',xdata=xdata,x=cq(n1,n2,d1),rate='fpr',col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fpr'),title.legend='near',legend='topright');
     
      xdata=xdata_near(n1=20,n2=n2,d1=0.5,near=near);
      dofig(plotragm,'fnr',xdata=xdata,x=cq(n1,n2,d1),rate='fnr',col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fnr'),title.legend='near',legend='topright');

      near=c(0.1,0.3);
      xdata.fpr=xdata_near(n1=20,n2=n2,d1=0,near=near);
      xdata.fnr=xdata_near(n1=20,n2=n2,d1=0.5,near=near);
      xdata=xdata_rbind(xdata.fpr,xdata.fnr);
      dofig(plotragm,'fpr+fnr',xdata=xdata,x=cq(n1,n2),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig(cq(fpr,fnr)),title.legend='near',legend='topright');
      ## support statements: near=0.1, n2=150 sweet spot with both error rates about 0.05
      ## For near=0.3 crossover point is n2=137 with error rates of about 0.15.
      n2crossover=crossover_n2bynear(xdata);
      ## error tables at end of Near Exact
      ## construct xdata for error tables
      near=c(0.1,0.3);
      xdata=xdata_near(n1=20,n2=n2,d1=d,near=near)
      err=nearexact_err(xdata);
      err_perl=err2perl(err);
      dotbl(n2crossover,err,err_perl);
      ## convert to almost the right thing in perl
      ## perl -n -e '@row=split; shift @row if /^\d+/; print "| ",join(" | ",@row)," |\n"' > foo
      ## copy-and-paste err_perl
      err<<-err;                        # for next section
    }
##### replication wise error rates
    if ((figsect='repwise') %in% sect) {
      title.desc='Replication-wise error';
      ## RW error tables
      prop.true=round(seq(0.05,0.95,by=0.05),digits=5);
      rwerr=nearexact_rwerr(err,prop.true);
      rwerr_perl=rwerr2perl(rwerr)
      dotbl(rwerr,rwerr_perl);
      ## convert to almost the right thing in perl
      ## perl -n -e '@row=split; shift @row if /^\d+/; print "| ",join(" | ",@row)," |\n"' > foo
      ## copy-and-paste reerr_perl
    }
##### supp_start - reset fignum to 1, set figpfx to 'S'. also tblnum, tblpfx
    if ((figsect='supp_start') %in% sect) {
      figpfx<<-tblpfx='S';
      fignum<<-1;
      tblnum<<-1;
    }
##### supp_exact
    if ((figsect='supp_exact') %in% sect) {
      figsect<<-sub('^supp_','',figsect);
      title.desc<<-'Exact replication';
      ## fpr
      dofig(plotrate,'fpr_n1=020',d=0,n1=20,n2=n2,smooth='spline',
            hline=c(fpr.cutoff/2,fpr.cutoff,fnr.cutoff),
            vhlty='dashed',vhlwd=c(1,0.5,0.5),vhcol=cq(red,black,black),plot.cutoff=F,
            title=title_resig('fpr'),legend=NULL);
      dofig(plotrate,'fpr_n1=200',d=0,n1=200,n2=n2,smooth='spline',
            hline=c(fpr.cutoff/2,fpr.cutoff,fnr.cutoff),
            vhlty='dashed',vhlwd=c(1,0.5,0.5),vhcol=cq(red,black,black),plot.cutoff=F,
            title=title_resig('fpr'),legend=NULL);
      ## confirm that FPR is sig.level/2 and that sdir is what causes FPR factor of 2
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
      ## boxplots also show no obvious correlation with n1, n2, or mean(n1,n2).
      ## TODO: decide whether to add as figures...
      ## boxplot(sig2~n1,data=drat.std,notch=F); abline(h=c(0.025,mean(drat.nosdir$sig2)),col='red',lty=cq(dashed,solid))
      ##  boxplot(sig2~n1,data=drat.nosdir,notch=F); abline(h=c(0.050,mean(drat.nosdir$sig2)),col='red',lty=cq(dashed,solid))
      ## fnr
      ## compute 1-power2 vs. n2 - specialized for plotfnr!
      power.n2=power_n2();
      ## construct xdata lists for n1=20 and n1=200
      xdata.020=lapply(d.nonzro,function(d) xdata=expand.grid(n1=20,n2=n2,d1=d,d2=d));
      names(xdata.020)=as.character(d.nonzro);
      xdata.200=lapply(d.nonzro,function(d) xdata=expand.grid(n1=200,n2=n2,d1=d,d2=d));
      names(xdata.200)=as.character(d.nonzro);
      ## plots using std posr - includes sdir
      dofig(plotfnr_exact,'fnr_n1=020',xdata=xdata.020,col=col)
      dofig(plotfnr_exact,'fnr_n1=200',xdata=xdata.200,col=col);
      ## plots using sig1_sig1 posr - omits sdir
      dofig(plotfnr_exact,'fnr_nosdir_n1=020',xdata=xdata.020,col=col,posr.id='sig1_sig1')
      dofig(plotfnr_exact,'fnr_nosdir_n1=200',xdata=xdata.200,col=col,posr.id='sig1_sig1');
      ## fnr vs 1-power
      xdata=expand.grid(n1=n,n2=n,d1=d.nonzro,d2=d);
      xdata=subset(xdata,subset=(d1==d2));
      dofig(plotfnrpwr_exact,'fnr_vs_power',xdata=xdata);
      dofig(plotfnrpwr_exact,'fnr_vs_power_nosdir',xdata=xdata,posr.id='sig1_sig1');
      ## correlation of fnr vs power2
      drat.std=dratfnrpwr_exact(xdata,posr.id='std');
      drat.nosdir=dratfnrpwr_exact(xdata,posr.id='sig1_sig1');
      corfnrpwr=do.call(rbind,lapply(d.nonzro,function(d) {
        drat.std=subset(drat.std,subset=(d2<=d));
        drat.nosdir=subset(drat.nosdir,subset=(d2<=d));
        cor.std=with(drat.std,cor(sig2,1-pwr2));
        cor.nosdir=with(drat.nosdir,cor(sig2,1-pwr2));
        data.frame(d2=d,cor.std,cor.nosdir);}))
      dotbl(corfnrpwr);
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
      dotbl(fnr_n2byd2,fnr_d2byn2);
      dofig(plotfnrcutoff_exact,'fnr_cutoff',fnr_d2byn2=fnr_d2byn2);
    }
##### supp_inexact
    if ((figsect='supp_inexact') %in% sect) {
      figsect=sub('^supp_','',figsect);
      title.desc='Inexact replication';
      BREAKPOINT()
      ## fpr
      xdata.020=xdata_inexact(n1=20,d1=0);
      xdata.200=xdata_inexact(n1=200,d1=0);
      dofig(plotratm,'fpr_n1=020',xdata=xdata.020,x=cq(n1,n2,d1),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fpr'),title.legend='d2',legend='topright');
      dofig(plotratm,'fpr_n1=200',xdata=xdata.200,x=cq(n1,n2,d1),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fpr'),title.legend='d2',legend='topright');
      ## support statement: data are similar but not exactly the same.
      corfprdata=drat_cor(xdata.020,xdata.200);
      dotbl(corfprdata);
      ## plot shows similarity. TODO: decide whether to add as figure
      ## drat_plotsig2(xdata.020,xdata.200)
      ## fnr
      xdata.020=xdata_inexact(n1=20,d1=0.5);
      xdata.200=xdata_inexact(n1=200,d1=0.5);
      dofig(plotratm,'fnr_n1=020',xdata=xdata.020,x=cq(n1,n2,d1),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fnr'),title.legend='d2',x.legend=8.45,y.legend=0.65);
      dofig(plotratm,'fnr_n1=200',xdata=xdata.020,x=cq(n1,n2,d1),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fnr'),title.legend='d2',x.legend=8.45,y.legend=0.65);
      ## support statement: data are similar but not exactly the same.
      corfnrdata=drat_cor(xdata.020,xdata.200);
      dotbl(corfnrdata);
      ## plot shows similarity. TODO: decide whether to add as figure
      ## drat_plotsig2(xdata.020,xdata.200)
      ## fpr+fnr
      d2=c(0,0.1,0.2,0.5,1)
      xdata.fpr=xdata_inexact(n1=20,d1=0,d2=d2);
      xdata.fnr=xdata_inexact(n1=20,d1=0.5,d2=d2);
      xdata.020=xdata_rbind(xdata.fpr,xdata.fnr);
      dofig(plotragm,'fpr+fnr_n1=020',xdata=xdata.020,x=cq(n1,n2),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig(cq(fpr,fnr)),
            title.legend='d2',x.legend=9.1,y.legend=0.985,cex.legend=0.7);
      xdata.fpr=xdata_inexact(n1=200,d1=0,d2=d2);
      xdata.fnr=xdata_inexact(n1=200,d1=0.5,d2=d2);
      xdata.200=xdata_rbind(xdata.fpr,xdata.fnr);
      dofig(plotragm,'fpr+fnr_n1=200',xdata=xdata.200,x=cq(n1,n2),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig(cq(fpr,fnr)),
            title.legend='d2',x.legend=9.1,y.legend=0.985,cex.legend=0.7);
      ## obvious rocm
      xdata=lapply(d,function(d2) xdata=expand.grid(n1=c(20,200),n2=n2,d1=d,d2=d2));
      names(xdata)=as.character(d);
      dofig(plotrocm,'rocm',xdata=xdata,x=cq(n1,n2),col=col,
            hline=c(fpr.cutoff,fnr.cutoff),vline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',
            vhlwd=0.5,plot.cutoff=F,
            title=title_resig(NULL,'False negative vs. false positive rate'),title.legend='d2');
    }
    ##### supp_nearexact
    if ((figsect='supp_nearexact') %in% sect) {
      title.desc='Near exact replication';
      BREAKPOINT()
      ## near=round(c(0,0.01,0.05,0.1,0.2),digits=5);
      near=seq(0.1,0.5,by=0.1);
      ## fpr
      xdata.020=xdata_near(n1=20,n2=n2,d1=0,near=near);
      xdata.200=xdata_near(n1=200,n2=n2,d1=0,near=near);
      dofig(plotragm,'fpr_n1=020',xdata=xdata.020,x=cq(n1,n2,d1),col=col,smooth='spline',
            rate='fpr',hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fpr'),title.legend='near',legend='topright');
      dofig(plotragm,'fpr_n1=200',xdata=xdata.200,x=cq(n1,n2,d1),col=col,smooth='spline',
            rate='fpr',hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fpr'),title.legend='near',legend='topright');
      ## support statement: data are similar but not exactly the same.
      corfprdata=drat_cor(xdata.020,xdata.200);
      dotbl(corfprdata);
      ## plot shows similarity. TODO: decide whether to add as figure
      ## drat_plotsig2(xdata.020,xdata.200)
      ## fnr
      xdata.020=xdata_near(n1=20,n2=n2,d1=0.5,near=near);
      xdata.200=xdata_near(n1=200,n2=n2,d1=0.5,near=near);
      dofig(plotragm,'fnr_n1=020',xdata=xdata.020,x=cq(n1,n2,d1),col=col,smooth='spline',
            rate='fnr',hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fnr'),title.legend='near',legend='topright');
      dofig(plotragm,'fnr_n1=200',xdata=xdata.200,x=cq(n1,n2,d1),col=col,smooth='spline',
            rate='fnr',hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fnr'),title.legend='near',legend='topright');
       ## support statement: data are similar but not exactly the same.
      corfnrdata=drat_cor(xdata.020,xdata.200);
      dotbl(corfnrdata);
      ## plot shows similarity. TODO: decide whether to add as figure
      ## drat_plotsig2(xdata.020,xdata.200)
      ## plot fnr for single values of n1, n2, multiple d1
      n1=20;
      ## CAUTION: must use loop, NOT sapply, for scoping of fignum to work!
      ## Hmmm... not sure it makes sense to do all these figures
      for (n2x in c(100,200,300,400)) {
        xdata=xdata_near(n1=n1,n2=n2x,d1=d.nonzro,near=near);
        figname=paste(sep='_','fnr',paste_nv(n1),paste_nv(n2,n2x));
        dofig(plotragm,figname,xdata=xdata,x=cq(n1,n2,d1),rate='fnr',col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fnr'),title.legend='near',legend='topright');
      }
      ## fpr+fnr
      xdata.fpr=xdata_near(n1=20,n2=n2,d1=0,near=near);
      xdata.fnr=xdata_near(n1=20,n2=n2,d1=0.5,near=near);
      xdata.020=xdata_rbind(xdata.fpr,xdata.fnr);
      dofig(plotragm,'fpr+fnr_n1=020',xdata=xdata.020,x=cq(n1,n2),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig(cq(fpr,fnr)),title.legend='near',legend='topright');
      xdata.fpr=xdata_near(n1=200,n2=n2,d1=0,near=near);
      xdata.fnr=xdata_near(n1=200,n2=n2,d1=0.5,near=near);
      xdata.200=xdata_rbind(xdata.fpr,xdata.fnr);
      dofig(plotragm,'fpr+fnr_n1=200',xdata=xdata.200,x=cq(n1,n2),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig(cq(fpr,fnr)),title.legend='near',legend='topright');
      ## do it systematically for all values of near, several values of n1, d1
      near=d;
      n1=c(20,200);
      d1=c(0.2,0.5,0.8);
      ## CAUTION: must use loop, NOT sapply, for scoping of fignum to work!
      ## Hmmm... not sure it makes sense to do all these figures
      for (n1x in c(20,200)) {
        for (d1x in c(0.2,0.5,0.8)) {
          xdata.fpr=xdata_near(n1=n1x,n2=n2,d1=0,near=near);
          xdata.fnr=xdata_near(n1=n1x,n2=n2,d1=d1x,near=near);
          xdata=xdata_rbind(xdata.fpr,xdata.fnr);
          figname=paste(sep='_','fpr+fnr',paste_nv(n1,n1x),paste_nv(d1,d1x));
          dofig(plotragm,figname,xdata=xdata,x=cq(n1,n2),col=col,smooth='spline',
                hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
                title=title_resig(cq(fpr,fnr),title.desc=paste_nv(d1,d1x)),
                title.legend='near',legend='topright');
          }}
      ## obvious rocm
      near=d;
      xdata=xdata_near(n1=c(20,200),n2=n2,d1=d,near=near);
      dofig(plotrocm,'rocm',xdata=xdata,x=cq(n1,n2),col=col,
            hline=c(fpr.cutoff,fnr.cutoff),vline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',
            vhlwd=0.5,plot.cutoff=F,
            title=title_resig(NULL,'False negative vs. false positive rate'),title.legend='near');
    }
    sect;
  })
}
