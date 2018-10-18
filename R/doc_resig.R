#################################################################################
##
## Author:  Nat Goodman
## Created: 18-09-05
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
## Generate figures and tables for resig.R blog post
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## --- Generate Figures and Tables for resig Blog Post ---
## sect is which sections to run - for use during development
##   uses prefix matching and runs all that match
doc_resig=function(sect=parent(sect,NULL)) {
  sect.all=cq(exact,nearexact,repwise);
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
    ## near exact
    if (sect=='nearexact') {
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
    ## replication wise error rates
    if (sect=='repwise') {
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
    sect;
  })
}
