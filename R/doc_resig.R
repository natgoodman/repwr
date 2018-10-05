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
## Generate figures and tables for repwr.R blog post
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## --- Generate Figures and Tables for resig Blog Post
## make figures and tables for blog post
## subdoc is blog, supplement, all. control which sections run
## sect is which sections to run - for use during development
##   uses prefix matching and all matches run
doc_resig=
  function(sect=parent(sect,NULL),
           subdoc='blog',subdocx=match.arg(subdoc,cq(blog,supplement,all)),
           figpfx=if((subdocx=='supplement')|(!is.null(sect)&all(grepl('^supp_',sect))))
                    'S' else NULL,
           tblpfx=if((subdocx=='supplement')|(!is.null(sect)&all(grepl('^supp_',sect))))
                    'S' else NULL,
           fignum=parent(fignum,1),tblnum=parent(tblnum,1),
           figscreen=parent(figscreen,T),fignew=parent(fignew,figscreen)) {
    subdoc=subdocx;                      # to avoid confusion later
    sect.blog=cq(exact,nearexact,repwise);
    sect.supp=cq(supp_start,supp_exact,supp_inexact,supp_nearexact);
    sect.all=c(sect.blog,sect.supp);
    if (is.null(sect)) {
      sect=if (subdoc=='blog') sect.blog else if (subdoc=='supplement') sect.supp else sect.all;
      }
    else sect=pmatch_choice(sect,sect.all);
    d.nonzro=d[d!=0];
    n2=n[n>=50];                       # standard n2: exclude 20
    col=col_resig();                   # colors for plotratm, plotragm
##### blog_start - reset fignum to 1, set figpfx to NULL
##### not necessary. included for stylistic consistency with 'supp'
    if ((figsect='blog_start') %in% sect) {
      figpfx=tblpfx=NULL;
      fignum=tblnum=1;
    }
##### exact
    if ((figsect='exact') %in% sect) {
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
      figpfx=tblpfx='S';
      fignum=tblnum=1;
    }
##### supp_exact
    if ((figsect='supp_exact') %in% sect) {
      figsect=sub('^supp_','',figsect);
      title.desc='Exact replication';
      ## fpr
      dofig(plotrate,'fpr_n1=020',d=0,n1=20,n2=n2,smooth='spline',
            hline=c(fpr.cutoff/2,fpr.cutoff,fnr.cutoff),
            vhlty='dashed',vhlwd=c(1,0.5,0.5),vhcol=cq(red,black,black),plot.cutoff=F,
            title=title_resig('fpr'),legend=NULL);
      dofig(plotrate,'fpr_n1=200',d=0,n1=200,n2=n2,smooth='spline',
            hline=c(fpr.cutoff/2,fpr.cutoff,fnr.cutoff),
            vhlty='dashed',vhlwd=c(1,0.5,0.5),vhcol=cq(red,black,black),plot.cutoff=F,
            title=title_resig('fpr'),legend=NULL);
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
      ## TODO: plot something from these tables
      dofig(plotfnrcutoff_exact,'fnr_cutoff',fnr_d2byn2=fnr_d2byn2);

    }
##### supp_inexact
    if ((figsect='supp_inexact') %in% sect) {
      figsect=sub('^supp_','',figsect);
      title.desc='Inexact replication';
      xdata=lapply(d,function(d)
        xdata=expand.grid(n1=20,n2=seq(50,by=50,len=10),d1=0,d2=d));
      names(xdata)=as.character(d);
      dofig(plotratm,'fpr',xdata=xdata,x=cq(n1,n2,d1),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fpr'),title.legend='d2',legend='topright');
      xdata=lapply(d,function(d)
        xdata=expand.grid(n1=20,n2=seq(50,by=50,len=10),d1=0.5,d2=d));
      names(xdata)=as.character(d);
      dofig(plotratm,'fnr',xdata=xdata,x=cq(n1,n2,d1),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fnr'),title.legend='d2',x.legend=8.45,y.legend=0.65);
   }
    if ((figsect='unused') %in% sect) {
      title.desc='Near exact replication';
      xdata=lapply(near,function(near) {
        do.call(rbind,lapply(d.nonzro,function(d1) {
          d2=round(seq(d1-near,d1+near,by=0.01),digits=5);
          ## trim d2 to [0,1]
          d2=d2[d2>=0&d2<=1];        
          expand.grid(n1=20,n2=150,d1=d1,d2=d2)}))});
      names(xdata)=as.character(near);
      dofig(plotragm,'fnr',xdata=xdata,x=cq(n1,n2,d1),rate='fnr',col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title.desc=title.desc,title.legend='near',legend='topright');
    }
    sect;
  }
## generate standard xdata for near-exact plots
xdata_near=function(n1,n2,d1,near,step=0.01) {
  d1=round(d1,digits=5);
  xdata=lapply(near,function(near) {
    do.call(rbind,lapply(d1,function(d1) {
      d2=round(seq(d1-near,d1+near,by=step),digits=5);
      d2=d2[d2>=0&d2<=1];             #trim d2 to [0,1]
      expand.grid(n1=20,n2=seq(50,by=50,len=10),d1=d1,d2=d2)}))});
  names(xdata)=as.character(near);
  xdata;
}
## merge compatible xdata lists for combined fpr,fnr near-exact plots
xdata_rbind=function(xdata1,xdata2) {
  name=unique(c(names(xdata1),names(xdata2)));
  xdata=lapply(name,function(name) rbind(xdata1[[name]],xdata2[[name]]));
  names(xdata)=name;
  xdata;
}
## palettes are RColorBrewer sequential palette names
col_resig=
  function(lo.brk=0.2,hi.brk=0.5,res=0.01,
           lo.palette='Reds',mid.palette='Blues',hi.palette='Greens') {
    lo.steps=(lo.brk/res)+1
    mid.steps=(hi.brk-lo.brk)/res;
    hi.steps=(1-hi.brk)/res;
    ## low end of RColorBrewer sequential palettes are too light
    lo.col=colorRampPalette(RColorBrewer::brewer.pal(9,lo.palette))(2*lo.steps)[-(1:lo.steps)];
    mid.col=colorRampPalette(RColorBrewer::brewer.pal(9,mid.palette))(2*mid.steps)[-(1:mid.steps)];
    hi.col=colorRampPalette(RColorBrewer::brewer.pal(9,hi.palette))(2*hi.steps)[-(1:hi.steps)];
    setNames(c(lo.col,mid.col,hi.col),seq(0,1,by=res));
  }
## generate title for doc_resig. simpler and shorter than general case
title_resig=
  function(rate.type=parent(rate.type,'error'),title.desc=parent(title.desc,NULL),
           fignum=parent(fignum,NULL),figpfx=parent(figpfx,NULL),
           posr.id=parent(posr.id,'std')) {
    if (is.null(rate.type)) rate.desc=NULL
    else {
      rate.desc=sapply(rate.type,function(rate.type) 
        switch(rate.type,
               pos='positive',neg='negative',error='error',correct='correct',
               fpr='false positive',fnr='false negative',
               tpr='true positive',tnr='true negative',
               roc='rate vs rate',rag='mean',ragm='mean'));
      rate.desc=paste(collapse=' and ',rate.desc);
      if (all(rate.type %notin% cq(roc,ragm))) rate.desc=paste(sep=' ',rate.desc,'rate');
    }
    posr.desc=if (posr.id=='std') NULL else 'w/o same-direction';
    fignum=paste(collapse='',c(figpfx,fignum));
    if (!is.null(fignum)) fignum=paste(sep='','Figure ',fignum,'.');
    paste(collapse=' ',c(fignum,title.desc,rate.desc,posr.desc));
  }
## these functions calculate crossover points
cutoff_n2byd2=function(xdata,cutoff=c(0.05,0.2)) {
  if (!is.data.frame(xdata)) xdata=do.call(rbind,xdata);
  drat=data_rate(xdata=xdata,x=cq(n2,d2),rate.type='error',mesr='sig2')
  byd2=split(drat,drat$d2);
  do.call(rbind,lapply(cutoff,function(cutoff) {
    do.call(rbind,lapply(byd2,function(drat) {
      n2sig2=function(n2) predict(smooth.spline(drat$n2,drat$sig2,spar=0.5),n2)$y;
      if (sign(n2sig2(50)-cutoff)==sign(n2sig2(500)-cutoff)) n2=NA
      else n2=uniroot(function(n2) n2sig2(n2)-cutoff,interval=c(50,500))$root;
      data.frame(cutoff=cutoff,d2=unique(drat$d2),n2=round(n2));
    }))}))}
cutoff_d2byn2=function(xdata,cutoff=c(0.05,0.2)) {
  if (!is.data.frame(xdata)) xdata=do.call(rbind,xdata);
  drat=data_rate(xdata=xdata,x=cq(n2,d2),rate.type='error',mesr='sig2')
  byn2=split(drat,drat$n2);
  do.call(rbind,lapply(cutoff,function(cutoff) {
    do.call(rbind,lapply(byn2,function(drat) {
      d2sig2=function(d2) predict(smooth.spline(drat$d2,drat$sig2,spar=0.5),d2)$y;
      if (sign(d2sig2(0)-cutoff)==sign(d2sig2(1)-cutoff)) d2=NA
      else d2=uniroot(function(d2) d2sig2(d2)-cutoff,interval=c(0,1))$root
      data.frame(cutoff=cutoff,n2=unique(drat$n2),d2=d2);
    }))}))}
crossover_n2bynear=function(xdata) {
  if (is.data.frame(xdata))
    stop('crossover_n2bynear needs a list of xdata data fromes, not a single data frame');
  near=names(xdata);
  do.call(rbind,lapply(near,function(near) {
    xdata=xdata[[near]];
    xdata.fpr=subset(xdata,subset=d1==0);
    drag.fpr=data_agg(xdata=xdata.fpr,x=cq(n2),rate='fpr',mesr='sig2');
    data.fpr=data.frame(n2=drag.fpr$byx,fpr=drag.fpr$fpr[,1]);
    xdata.fnr=subset(xdata,subset=d1!=0);
    d1=unique(xdata.fnr$d1);
    do.call(rbind,lapply(d1,function(d1) {
      xdata.fnr=xdata.fnr[xdata.fnr$d1==d1,,drop=F];
      drag.fnr=data_agg(xdata=xdata.fnr,x=cq(n2),rate='fnr',mesr='sig2');
      data.fnr=data.frame(n2=drag.fnr$byx,fnr=drag.fnr$fnr[,1]);
      n2fpr=function(n2,data=parent(data.fpr)) {
        predict(smooth.spline(data$n2,data$fpr,spar=0.5),n2)$y
      }
      n2fnr=function(n2,data=parent(data.fnr)) {
        predict(smooth.spline(data$n2,data$fnr,spar=0.5),n2)$y
      }
      if (sign(n2fpr(50)-n2fnr(50))==sign(n2fpr(500)-n2fnr(500))) 
        return(data.frame(near,d1,n2=NA,rate=NA))
      else
        n2=uniroot(function(n2) n2fpr(n2)-n2fnr(n2),interval=c(50,500))$root; rate=n2fpr(n2);
      data.frame(near,d1,n2=round(n2),rate)
    }));
  }))}
## error table for end of Near Exact
nearexact_err=function(xdata) {
  if (is.data.frame(xdata))
    stop('nearexact_err needs a list of xdata data fromes, not a single data frame');
  near=names(xdata);
  do.call(rbind,lapply(near,function(near) {
    xdata=xdata[[near]];
    xdata.fpr=subset(xdata,subset=d1==0);
    drag.fpr=data_agg(xdata=xdata.fpr,x=cq(n2),rate='fpr',mesr='sig2');
    data.fpr=data.frame(n2=drag.fpr$byx,fpr=drag.fpr$fpr[,1]);
    xdata.fnr=subset(xdata,subset=d1!=0);
    d1=unique(xdata.fnr$d1);
    do.call(rbind,lapply(d1,function(d1) {
      xdata.fnr=xdata.fnr[xdata.fnr$d1==d1,,drop=F];
      drag.fnr=data_agg(xdata=xdata.fnr,x=cq(n2),rate='fnr',mesr='sig2');
      data.fnr=data.frame(n2=drag.fnr$byx,fnr=drag.fnr$fnr[,1]);
      n2fpr=function(n2,data=parent(data.fpr)) {
        predict(smooth.spline(data$n2,data$fpr,spar=0.5),n2)$y
      }
      n2fnr=function(n2,data=parent(data.fnr)) {
        predict(smooth.spline(data$n2,data$fnr,spar=0.5),n2)$y
      }
      n2=sort(unique(xdata$n2));
      ## data.frame(near=as.numeric(near),d1,n2,fpr=n2fpr(n2),fnr=max(0,n2fnr(n2)));
      fpr=n2fpr(n2);
      fnr=sapply(n2fnr(n2), function(fnr) max(0,fnr));
      data.frame(near=as.numeric(near),d1,n2,fpr,fnr);
    }))}))}
## convert err table into format useful for perl
err2perl=function(tbl.err,n2=c(150,300,450),near=c(0.1,0.3),d1=c(0.2,0.5,0.8)) {
  tbl.err=tbl.err[tbl.err$n2%in%n2&tbl.err$near%in%near&tbl.err$d1%in%d1,];
  tbl.bycase=split(tbl.err,with(tbl.err,paste(n2,near)))
  fnr=do.call(rbind,lapply(tbl.bycase,function(tbl) t(tbl$fnr)));
  colnames(fnr)=cq(FNR_02,FNR_05,FNR_08)
  fpr=do.call(rbind,lapply(tbl.bycase,function(tbl) tbl[1,'fpr']))
  case=do.call(rbind,strsplit(names(tbl.bycase),' '))
  case=apply(case,1:2,as.numeric);
  colnames(case)=cq(n2,near)
  tbl=data.frame(case,FPR=fpr,fnr)
  rownames(tbl)=NULL
  tbl=round(tbl,digits=2)
  tbl=subset(tbl,subset=n2%in%c(150,300,450))
  tbl[with(tbl,order(n2,near)),];
}
## convert to almost the right thing in perl
## perl -n -e '@row=split; shift @row if /^\d+/; print "| ",join(" | ",@row)," |\n"' > foo
## copy-and-paste err_perl

## RW error table for Replication-wise error rates
nearexact_rwerr=function(err,prop.true) {
  do.call(rbind,lapply(prop.true,function(prop.true) 
    with(err,data.frame(n2,near,d1,prop.true,
                        rwfpr=rw_fpr(fpr,fnr,prop.true),rwfnr=rw_fnr(fpr,fnr,prop.true)))));
}
## convert rwerr table into format useful for perl
rwerr2perl=
  function(tbl.rwerr,n2=150,near=c(0.1,0.3),d1=c(0.2,0.5,0.8),prop.true=c(0.1,0.25,0.5,0.75,0.9)) {
    tbl.rwerr=tbl.rwerr[tbl.rwerr$n2%in%n2&tbl.rwerr$near%in%near&tbl.rwerr$d1%in%d1,];
    tbl.bycase=split(tbl.rwerr,with(tbl.rwerr,paste(n2,near,prop.true)));
    rwfpr=do.call(rbind,lapply(tbl.bycase,function(tbl) t(tbl$rwfpr)));
    rwfnr=do.call(rbind,lapply(tbl.bycase,function(tbl) t(tbl$rwfnr)));
    colnames(rwfpr)=cq(RWFPR_02,RWFPR_05,RWFPR_08);
    colnames(rwfnr)=cq(RWFNR_02,RWFNR_05,RWFNR_08);
    case=do.call(rbind,strsplit(names(tbl.bycase),' '));
    case=apply(case,1:2,as.numeric);
    colnames(case)=cq(n2,near,prop.true);
    tbl=data.frame(case,rwfpr,rwfnr);
    rownames(tbl)=NULL;
    tbl=round(tbl,digits=2);
    tbl=subset(tbl,subset=(n2==150&prop.true%in%c(0.1,0.25,0.5,0.75,0.9)));
    tbl[with(tbl,order(n2,near,prop.true)),];
  }
## convert to almost the right thing in perl
## perl -n -e '@row=split; shift @row if /^\d+/; print "| ",join(" | ",@row)," |\n"' > foo
## copy-and-paste rwerr_perl

## false discovery rate calculations. why does this give me SOOO much trouble??
rw_fpr=function(fpr=0.05,fnr=0.2,prop.true=.5) {
  prop.false=1-prop.true;
  prop.tp=(1-fnr)*prop.true;
  prop.fp=fpr*prop.false;
  prop.pos=prop.fp+prop.tp;
  rw.fpr=prop.fp/prop.pos;
  rw.fpr;
}
rw_fnr=function(fpr=0.05,fnr=0.2,prop.true=.5) {
  prop.false=1-prop.true;
  prop.tn=(1-fpr)*prop.false;
  prop.fn=fnr*prop.true;
  prop.neg=prop.fn+prop.tn;
  rw.fnr=prop.fn/prop.neg;
  rw.fnr;
}
####################
## functions for supplement
## plot fnr rates with 1-power overlaid
plotfnr_exact=
  function(xdata,posr.id='std',col=parent(col),power.n2=parent(power.n2),fignum=parent(fignum),
           cex.title=0.9) {
    ## plot simulated results
    plotratm(xdata=xdata,posr.id=posr.id,x=cq(n1,n2),col=col,smooth='spline',
             hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
             title=title_resig('fnr'),title.legend='d',x.legend=8.45,y.legend=0.65);
    ## plot theoretical values
    x.smooth=power.n2$n2.smooth/50;       # scale n2 to x-axis. CAUTION: not general
    y.smooth=power.n2$y.smooth;
    matlines(x.smooth,y.smooth,col=col[colnames(y.smooth)],lty='dashed')
  }
## compute 1-power2 vs. n2 - specialized for plotfnr!  
power_n2=function(n2=parent(n2),d.nonzro=parent(d.nonzro)) {
  y=do.call(cbind,lapply(d.nonzro,function(d2) 1-power.t.test(n=n2,d=d2)$power));
  n2.smooth=seq(min(n2),max(n2),len=100);
  y.smooth=splinem(n2,y,xout=n2.smooth);
  ## clamp y.smooth to [0,1]. interpolation can under- or over-shoot
  y.smooth=apply(y.smooth,1:2,function(y) if (!is.na(y)) max(min(y,1),0) else y);
  colnames(y.smooth)=d.nonzro;
  list(n2.smooth=n2.smooth,y.smooth=y.smooth);
}
## plot fnr vs 1-power2
plotfnrpwr_exact=function(xdata,posr.id='std',fignum=parent(fignum),cex.title=0.9) {
  drat=dratfnrpwr_exact(xdata,posr.id);
  ## n2col=setNames(colorRampPalette(cq(grey75,black))(11),n);
  col=colorRampPalette(cq(grey75,black))(101)[round(drat$pwr1*100)+1];
  ## plot(1-pwr2,drat$sig2,col=n2col[as.character(drat$n1)],pch=19,cex=0.75,
  plot(1-drat$pwr2,drat$sig2,col=col,pch=19,cex=0.75,
       xlab='1-power2 (theory)',ylab='false negative rate (simulated)',
       main=title_resig(NULL,'False negative rate vs. 1-power2'),cex.main=cex.title);
  abline(a=0,b=1,col='red');
  grid();
}
## construct drat for fnr vs 1-power analyses
dratfnrpwr_exact=function(xdata,posr.id='std') {
  drat=data_rate(rate.rule='nonzro',xdata=xdata,posr.id=posr.id,mesr=cq(sig2));
  drat$pwr1=with(drat,power.t.test(n=n1,delta=d1)$power);
  drat$pwr2=with(drat,power.t.test(n=n2,delta=d2)$power);
  drat;
}
## plot fnr cutoffs
plotfnrcutoff_exact=function(fnr_d2byn2,fignum=parent(fignum),cex.title=0.9) {
  data=merge(subset(fnr_d2byn2,subset=(cutoff==0.05)),subset(fnr_d2byn2,subset=(cutoff==0.20)),
             by='n2',suffixes=cq('.05','.20'));
  n2=fnr_d2byn2$n2;
  y=data[,grep('^d2',colnames(data),value=T)];
  n2.smooth=seq(min(n2),max(n2),len=100);
  y.smooth=splinem(n2,y,xout=n2.smooth);
  col=cq(red,red,blue,blue);
  lty=cq(solid,dashed,solid,dashed);
  matplot(n2.smooth,y.smooth,type='l',col=col,lty=lty,xlab='n2',ylab='d2',
          main=title_resig(NULL,'d2 vs. n2 for FNR cutoff=0.05, 0.20'),cex.main=cex.title);
  legend=sapply(strsplit(colnames(y.smooth),'\\.'),
                function(row) {
                  n2=as.numeric(row[2]);
                  cutoff=as.numeric(row[3])/100;
                  paste(sep=', ',n2,cutoff)})
  legend('topright',bty='n',col=col,lty=lty,cex=0.8,title='n1, cutoff',legend=legend);
  grid();
}
