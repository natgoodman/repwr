#################################################################################
##
## Author:  Nat Goodman
## Created: 18-10-25
##          by copying doc_resig_fun.R created 18-10-18
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
xdata_exact=
  function(n1=parent(n1),n2=parent(n2),d1=parent(d1),d2=d1,d=d1,by=cq(d,n1,n2,d1,d2)) {
    d=round(d,digits=5);
    by=match.arg(by);
    if (by=='d') by='d1';
    xdata=expand.grid(n1=n1,n2=n2,d1=d);
    xdata$d2=xdata$d1;
    by=xdata[,by];
    split(xdata,by);
  }
## generate standard xdata for inexact plots
xdata_inexact=function(n1,n2=parent(n2),d1,d2=parent(d)) {
  d1=round(d1,digits=5);
  d2=round(d2,digits=5);
  xdata=lapply(d2,function(d2) xdata=expand.grid(n1=n1,n2=n2,d1=d1,d2=d2));
  names(xdata)=as.character(d2);
  xdata;
}
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
## colors keyed by d
color_d=
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
## colors keyed by n
color_n=
  function(lo.brk=100,hi.brk=350,
           lo.palette='Reds',mid.palette='Blues',hi.palette='Greens') {
    lo.steps=length(which(n<=lo.brk));
    mid.steps=length(which(n>lo.brk&n<=hi.brk));
    hi.steps=length(which(n>hi.brk));
    ## low end of RColorBrewer sequential palettes are too light
    lo.col=colorRampPalette(RColorBrewer::brewer.pal(9,lo.palette))(2*lo.steps)[-(1:lo.steps)];
    mid.col=colorRampPalette(RColorBrewer::brewer.pal(9,mid.palette))(2*mid.steps)[-(1:mid.steps)];
    hi.col=colorRampPalette(RColorBrewer::brewer.pal(9,hi.palette))(2*hi.steps)[-(1:hi.steps)];
    setNames(c(lo.col,mid.col,hi.col),n);
  }
## generate title for doc_resigsupp. simpler and shorter than general case
## like title_resig but puts Figure XXX on separate line
title_resigsupp=
  function(title.rate=parent(title.rate,'error'),title.desc=parent(title.desc,NULL),
           extra=parent(extra,F),
           sect=parent(sect,NULL),sectnum=parent(sectnum,NULL),sect.desc=parent(sect.desc,NULL),
           posr.id=parent(posr.id,'std')) {
    if (is.null(title.rate)) rate.desc=NULL
    else {
      rate.desc=sapply(title.rate,function(title.rate) 
        switch(title.rate,
               pos='positive',neg='negative',error='error',correct='correct',
               fpr='false positive',fnr='false negative',
               tpr='true positive',tnr='true negative',
               roc='rate vs rate',rag='mean',ragm='mean'));
      rate.desc=paste(collapse=' and ',rate.desc);
      if (all(title.rate %notin% cq(roc,ragm))) rate.desc=paste(sep=' ',rate.desc,'rate');
    }
    posr.desc=if (posr.id=='std') NULL else 'w/o same-direction';
    fig=paste(sep='','Figure ',figname(name,sect,sectnum),"\n");
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
####################
## --- Plot Functions for Supplement ---
## plot boxplots of fpr vs.n1
plotboxfpr_exact=
  function(drat,posr.id='std',cex.title=0.9,ylim=c(0,0.1),
           extra=parent(extra,F),
           sect=parent(sect,NULL),sectnum=parent(sectnum,NULL),sect.desc=parent(sect.desc,NULL)) {
    boxplot(sig2~n1,data=drat,notch=F,xlab='n1',ylab='false positive rate',ylim=ylim,
            main=title_resigsupp('fpr','vs. n1'),
            cex.main=cex.title);
    cutoff=if(posr.id=='std') 0.025 else 0.05;
    abline(h=c(cutoff,mean(drat$sig2)),col='red',lty=cq(dashed,solid));
    grid();
  }
    
## plot fnr rates with 1-power overlaid
plotfnr_exact=
  function(xdata,posr.id='std',d2col=parent(d2col),n2col=parent(n2col),
           power.n2=parent(power.n2),cex.title=0.9,
           extra=parent(extra,F),
           sect=parent(sect,NULL),sectnum=parent(sectnum,NULL),sect.desc=parent(sect.desc,NULL)) {
    ## plot simulated results
    plotratm(xdata=xdata,posr.id=posr.id,x=cq(n1,n2),col=col,smooth='spline',
             hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
             title=title_resigsupp('fnr'),title.legend='d',x.legend=8.45,y.legend=0.65);
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
plotfnrpwr_exact=
  function(xdata,posr.id='std',cex.title=0.9,
           sect=parent(sect,NULL),sectnum=parent(sectnum,NULL),sect.desc=parent(sect.desc,NULL)) {
    drat=dratfnrpwr_exact(xdata,posr.id);
    ## n2col=setNames(colorRampPalette(cq(grey75,black))(11),n);
    col=colorRampPalette(cq(grey75,black))(101)[round(drat$pwr1*100)+1];
    ## plot(1-pwr2,drat$sig2,col=n2col[as.character(drat$n1)],pch=19,cex=0.75,
    plot(1-drat$pwr2,drat$sig2,col=col,pch=19,cex=0.75,
         xlab='1-power2 (theory)',ylab='false negative rate (simulated)',
         main=title_resigsupp('fnr','vs. 1-power2'),cex.main=cex.title);
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
plotfnrcutoff_exact=
  function(fnr_d2byn2,cex.title=0.9,
           sect=parent(sect,NULL),sectnum=parent(sectnum,NULL),sect.desc=parent(sect.desc,NULL)) {
    data=merge(subset(fnr_d2byn2,subset=(cutoff==0.05)),subset(fnr_d2byn2,subset=(cutoff==0.20)),
               by='n2',suffixes=cq('.05','.20'));
    n2=fnr_d2byn2$n2;
    y=data[,grep('^d2',colnames(data),value=T)];
    n2.smooth=seq(min(n2),max(n2),len=100);
    y.smooth=splinem(n2,y,xout=n2.smooth);
    col=cq(red,red,blue,blue);
    lty=cq(solid,dashed,solid,dashed);
    matplot(n2.smooth,y.smooth,type='l',col=col,lty=lty,xlab='n2',ylab='d2',
            main=title_resigsupp(NULL,'d2 vs. n2 for FNR cutoff=0.05, 0.20'),cex.main=cex.title);
    legend=sapply(strsplit(colnames(y.smooth),'\\.'),
                  function(row) {
                    n2=as.numeric(row[2]);
                    cutoff=as.numeric(row[3])/100;
                    paste(sep=', ',n2,cutoff)})
    legend('topright',bty='n',col=col,lty=lty,cex=0.8,title='n1, cutoff',legend=legend);
    grid();
  }
## --- Miscellaneous Functions for Exploration ---
## compare drats, typically generated with different values of n1
## xdata may be data frames or lists of xdata data fromes
## x are 'same' columns between the xdatas
drat_cor=function(xdata1,xdata2,x=cq(n2,d1,d2)) {
  if (!is.data.frame(xdata1)) xdata1=do.call(rbind,xdata1);
  if (!is.data.frame(xdata2)) xdata2=do.call(rbind,xdata2);
  if (any(dim(xdata1)!=dim(xdata2)))
    stop(paste(sep='','incompatible xdata args: ',
               'nrow(xdata1)=',nrow(xdata1),'; ','nrow(xdata2)=',nrow(xdata2)));
  drat1=data_rate(xdata=xdata1,mesr='sig2');
  drat2=data_rate(xdata=xdata2,mesr='sig2');
  drat=merge(drat1,drat2,by=x,suffixes=c('1','2'))
  x1=drat$sig21;
  x2=drat$sig22;
  ## CAUTION: cor gives NA if either sig2 vector is constant. jiggle one element
  if (length(unique(x1))==1) x1[1]=x1[1]+.Machine$double.eps;
  if (length(unique(x2))==1) x2[1]=x2[1]+.Machine$double.eps;
  cor(x1,x2);
}
## plot sig2 vs sig2 for 2 drats
## TODO: this is very crude...
drat_plotsig2=function(xdata1,xdata2,x=cq(n2,d1,d2)) {
  if (!is.data.frame(xdata1)) xdata1=do.call(rbind,xdata1);
  if (!is.data.frame(xdata2)) xdata2=do.call(rbind,xdata2);
  if (any(dim(xdata1)!=dim(xdata2)))
    stop(paste(sep='','incompatible xdata args: ',
               'nrow(xdata1)=',nrow(xdata1),'; ','nrow(xdata2)=',nrow(xdata2)));
  drat1=data_rate(xdata=xdata1,mesr='sig2');
  drat2=data_rate(xdata=xdata2,mesr='sig2');
  drat=merge(drat1,drat2,by=x,suffixes=c('1','2'))
  x1=drat$sig21;
  x2=drat$sig22;
  plot(x1,x2,pch=16,cex=0.5);
  grid();
  abline(a=0,b=1,col='red');
}
