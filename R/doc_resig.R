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
## sect is which sections to run - for use during development
##   uses prefix matching and all matches run
doc_resig=
  function(sect=parent(sect,NULL),
           fignum=parent(fignum,1),figscreen=parent(figscreen,T),fignew=parent(fignew,figscreen)) {
    sect.all=cq(exact,inexact,nearexact);
    if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all);
    d.nonzro=d[d!=0];
    col=col_resig();                    # colors for plotratm, plotragm
#################### nonzro
##### exact
    if ((figsect='exact') %in% sect) {
      title.desc='Exact replication';

      dofig(plotrate,'fpr',d=0,n1=20,n2=seq(50,by=50,len=10),smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fpr'),legend=NULL);
       xdata=lapply(d.nonzro,function(d)
        xdata=expand.grid(n1=20,n2=seq(50,by=50,len=10),d1=d,d2=d));
      names(xdata)=as.character(d.nonzro);
      dofig(plotratm,'fnr',xdata=xdata,x=cq(n1,n2),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fnr'),title.legend='d',legend='right');
      tbl.n2byd2=cutoff_n2byd2(xdata);
      tbl.d2byn2=cutoff_d2byn2(xdata);
      dotbl(fnr_n2byd2=tbl.n2byd2,fnr_d2byn2=tbl.d2byn2);
    }
##### inexact
    if ((figsect='inexact') %in% sect) {
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
##### near exact
    if ((figsect='nearexact') %in% sect) {
      title.desc='Near exact replication';
      ## near=round(c(0,0.01,0.05,0.1,0.2),digits=5);
      near=round(seq(0,0.4,by=0.1),digits=5);
      xdata=lapply(near,function(near) {
        d2=round(seq(0,near,by=0.01),digits=5);
        expand.grid(n1=20,n2=seq(50,by=50,len=10),d1=0,d2=d2)});
      names(xdata)=as.character(near);
      dofig(plotragm,'fpr',xdata=xdata,x=cq(n1,n2,d1),rate='fpr',col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fpr'),title.legend='near',legend='topright');
      ## these tables show n2, near values that achieve 'cutoff' FPR
      tbl.fpr=cutoff_n2bynear(xdata,rate='fpr');
      ## tbl.nearbyn2=cutoff_nearbyn2(xdata,rate='fpr');
      dotbl(fpr=tbl.fpr);
     
      xdata=lapply(near,function(near) {
        d1=0.5;
        d2=round(seq(d1-near,d1+near,by=0.01),digits=5);
        ## trim d2 to [0,1]
        d2=d2[d2>=0&d2<=1]; 
        xdata=expand.grid(n1=20,n2=seq(50,by=50,len=10),d1=d1,d2=d2)});
      names(xdata)=as.character(near);
      dofig(plotragm,'fnr',xdata=xdata,x=cq(n1,n2,d1),rate='fnr',col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig('fnr'),title.legend='near',legend='topright');

      near=round(c(0.1,0.3),digits=5);
      xdata=lapply(near,function(near) {
        d2=round(seq(0,near,by=0.01),digits=5);
        fpr=expand.grid(n1=20,n2=seq(50,by=50,len=10),d1=0,d2=d2);
        d1=0.5;
        d2=round(seq(d1-near,d1+near,by=0.01),digits=5);
        ## trim d2 to [0,1]
        d2=d2[d2>=0&d2<=1]; 
        fnr=expand.grid(n1=20,n2=seq(50,by=50,len=10),d1=d1,d2=d2);
        rbind(fpr,fnr);
      })
      names(xdata)=as.character(near);
      dofig(plotragm,'fpr+fnr',xdata=xdata,x=cq(n1,n2),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title=title_resig(cq(fpr,fnr)),title.legend='near',legend='topright');

      ## these tables show n2, near values that achieve 'cutoff' FNR
      tbl.fnr=cutoff_n2bynear(xdata,rate='fnr');
      ## tbl.nearbyn2=cutoff_nearbyn2(xdata,rate='fnr');
      dotbl(fpr=tbl.fpr);
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
           fignum=parent(fignum,NULL),posr.id=parent(posr.id,'std')) {
    rate.desc=sapply(rate.type,function(rate.type) 
      switch(rate.type,
             pos='positive',neg='negative',error='error',correct='correct',
             fpr='false positive',fnr='false negative',
             tpr='true positive',tnr='true negative',
             roc='rate vs rate',rag='mean',ragm='mean'));
    rate.desc=paste(collapse=' and ',rate.desc);
    if (all(rate.type %notin% cq(roc,ragm))) rate.desc=paste(sep=' ',rate.desc,'rate');
    posr.desc=if (posr.id=='std') NULL else paste_nv('posr',posr.id);
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
      d2=uniroot(function(d2) d2sig2(d2)-cutoff,interval=c(0,1))$root
      data.frame(cutoff=cutoff,n2=unique(drat$n2),d2=d2);
    }))}))}

cutoff_n2bynear=function(xdata,rate,cutoff=c(0.05,0.2)) {
  if (is.data.frame(xdata))
    stop('cutoff_n2bynear needs a list of xdata data fromes, not a single data frame');
  ## construct list of drags, one per xdata data frame
  drag=sapply(xdata,function(xdata) data_agg(xdata=xdata,x=cq(n2),rate=rate,mesr='sig2'),
              simplify=F);
  do.call(rbind,lapply(cutoff,function(cutoff) {
    do.call(rbind,lapply(seq_along(drag),function(i) {
    near=names(drag)[i]; drag=drag[[i]]; 
    drag=data.frame(n2=drag$byx,sig2=drag[[rate]]);
    n2sig2=function(n2) predict(smooth.spline(drag$n2,drag$sig2,spar=0.5),n2)$y;
    if (sign(n2sig2(50)-cutoff)==sign(n2sig2(500)-cutoff)) n2=NA
    else n2=uniroot(function(n2) n2sig2(n2)-cutoff,interval=c(50,500))$root;
    data.frame(cutoff=cutoff,near=near,n2=round(n2));
    }))}))}

cutoff_nearbyn2=function(xdata,rate,cutoff=c(0.05,0.2)) {
  if (is.data.frame(xdata))
    stop('cutoff_n2bynear needs a list of xdata data fromes, not a single data frame');
  ## construct list of drags, one per xdata data frame
  drag=sapply(xdata,function(xdata) data_agg(xdata=xdata,x=cq(n2),rate=rate,mesr='sig2'),
              simplify=F);
  drag=do.call(rbind,lapply(cutoff,function(cutoff) {
    do.call(rbind,lapply(seq_along(drag),function(i) {
      near=as.numeric(names(drag)[i]); drag=drag[[i]];
      data.frame(near=near,n2=drag$byx,sig2=drag[[rate]]);
    }))}));
  byn2=split(drag,drag$n2);
  do.call(rbind,lapply(cutoff,function(cutoff) {
    do.call(rbind,lapply(byn2,function(drag) {
      nearsig2=function(near) predict(smooth.spline(drag$near,drag$sig2,spar=0.5),near)$y;
      if (sign(nearsig2(min(drag$near))-cutoff)==sign(nearsig2(max(drag$near))-cutoff)) near=NA
      else near=uniroot(function(near) nearsig2(near)-cutoff,interval=range(drag$near))$root
      data.frame(cutoff=cutoff,n2=unique(drag$n2),near=near);
    }))}))}

n2fpr=function(n2,near=0.3) {
  data=data.fpr[[as.character(near)]];
  predict(smooth.spline(data$n2,data$fpr,spar=0.5),n2)$y
}
n2fnr=function(n2,d1=0.5,near=0.3) {
  data=data.fnr[[as.character(near)]];
  data=data[data$d1==d1,,drop=F];
  predict(smooth.spline(data$n2,data$fnr,spar=0.5),n2)$y
}

crossover=function(near=0.3,d1=0.5) {
  if (sign(n2fpr(50,near=near)-n2fnr(50,d1=d1,near=near))==
      sign(n2fpr(500,near=near)-n2fnr(500,d1=d1,near=near))) 
    return(data.frame(near,d1,n2=NA,rate=NA));
  n2=uniroot(function(n2) n2fpr(n2,near=near)-n2fnr(n2,d1=d1,near=near),interval=c(50,500))$root;
  rate=n2fpr(n2,near=near);
  data.frame(near,d1,n2=round(n2),rate)
}
