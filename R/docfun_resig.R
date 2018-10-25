#################################################################################
##
## Author:  Nat Goodman
## Created: 18-10-25
##          by renaming doc_resig_fun.R created 18-10-18
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
## Functions used by doc_resig.R and doc_resigsupp.R
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## --- Generate and Manipulate xdata Data Frames ---
## generate standard xdata for exact plots
## CAUTION: sytlistically different from other xdata functions
##   I may change them all to this style if I like it
## generate standard xdata for near-exact plots
xdata_near=function(n1,n2=parent(n2),d1,near,step=0.01) {
  d1=round(d1,digits=5);
  xdata=lapply(near,function(near) {
    do.call(rbind,lapply(d1,function(d1) {
      d2=round(seq(d1-near,d1+near,by=step),digits=5);
      d2=d2[d2>=0&d2<=1];             #trim d2 to [0,1]
      expand.grid(n1=n1,n2=n2,d1=d1,d2=d2)}))});
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
## --- Plot-Related Utility Functions ---
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
           sect=parent(sect,NULL),sectnum=parent(sectnum,NULL),sect.desc=parent(sect.desc,NULL),
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
    fig=paste(sep='','Figure ',figname(name,sect,sectnum),'.');
    paste(collapse=' ',c(fig,sect.desc,rate.desc,title.desc,posr.desc));
  }
## --- Generate Analytic Tables ---
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
