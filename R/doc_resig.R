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
## CAUTION: doesn't work to pass version via dodoc. knonw bug. see Outliner
doc_resig=
  function(sect=parent(sect,NULL),version=cq(original,shorter),
           fignum=parent(fignum,1),figscreen=parent(figscreen,T),fignew=parent(fignew,figscreen),
           tblnum=parent(tblnum,NULL)) {
    version=match.arg(version);
    sect.all=if(version=='original') cq(exact,inexact,nearexact,repwise)
             else cq(exact,nearexact,repwise);
    if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all);
    d.nonzro=d[d!=0];
    col=col_resig();                    # colors for plotratm, plotragm
#################### nonzro
##### exact
    if ((figsect='exact') %in% sect) {
      title.desc='Exact replication';
      dofig(plotrate,'fpr',d=0,n1=20,n2=seq(50,by=50,len=10),smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title.desc=title.desc,legend=NULL);
      xdata=lapply(d.nonzro,function(d)
        xdata=expand.grid(n1=20,n2=seq(50,by=50,len=10),d1=d,d2=d));
      names(xdata)=as.character(d.nonzro);
      dofig(plotratm,'fnr',xdata=xdata,x=cq(n1,n2),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title.desc=title.desc,title.legend='d',legend='right');
      ## tables shows n2,d2 value that achieve 'cutoff' FNR
      n2byd2=cutoff_n2byd2(xdata);
      d2byn2=cutoff_d2byn2(xdata);
      ## theoretical calculation of n2 value that achieve 'cutoff' FNR
      n2.20=power.t.test(delta=0.2,power=0.80)$n;
      n2.05=power.t.test(delta=0.2,power=0.95)$n;
      dotbl(fnr_n2byd2=n2byd2,fnr_d2byn2=d2byn2,fnr_power=c(n2.20,n2.05));
    }
##### inexact
    if ((figsect='inexact') %in% sect) {
      title.desc='Inexact replication';
      xdata=lapply(d,function(d)
        xdata=expand.grid(n1=20,n2=seq(50,by=50,len=10),d1=0,d2=d));
      names(xdata)=as.character(d);
      dofig(plotratm,'fpr',xdata=xdata,x=cq(n1,n2,d1),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title.desc=title.desc,title.legend='d2',legend='topright');
      xdata=lapply(d,function(d)
        xdata=expand.grid(n1=20,n2=seq(50,by=50,len=10),d1=0.5,d2=d));
      names(xdata)=as.character(d);
      dofig(plotratm,'fnr',xdata=xdata,x=cq(n1,n2,d1),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title.desc=title.desc,title.legend='d2',x.legend=8.45,y.legend=0.65);
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
            title.desc=title.desc,title.legend='near',legend='topright');
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
            title.desc=title.desc,title.legend='near',legend='topright');

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
            title.desc=title.desc,title.legend='near',x.legend=425/50,y.legend=1);

      ## construct xdata for cutoff and crossover tables
      near=d.nonzro;
      xdata=lapply(near,function(near) {
        do.call(rbind,lapply(d,function(d1) {
          d2=round(seq(d1-near,d1+near,by=0.01),digits=5);
          d2=d2[d2>=0&d2<=1];           # trim d2 to [0,1]
          expand.grid(n1=20,n2=seq(50,by=50,len=10),d1=d1,d2=d2);
        }))})
      names(xdata)=near;
      ## n2, near values that achieve 'cutoff' rate for xdata w/ range of d1
      n2cutoff=cutoff_n2bynear(xdata);
      n2crossover=crossover_n2bynear(xdata);
      dotbl(n2cutoff,n2crossover)
      ## error tables at end of Near Exact
      err=nearexact_err(xdata);
      err_perl=err2perl(err);
      dotbl(err,err_perl);
      ## convert to almost the right thing in perl
      ## perl -n -e '@row=split; shift @row if /^\d+/; print "| ",join(" | ",@row)," |\n"' > foo
      ## copy-and-paste err_perl
      err<<-err;                        # for next section
    }
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
## these functions calculate cutoff and crossover points
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

cutoff_n2bynear=function(xdata,rate=cq(fpr,fnr),cutoff=c(0.05,0.2)) {
  if (is.data.frame(xdata))
    stop('cutoff_n2bynear needs a list of xdata data fromes, not a single data frame');
  near=names(xdata);
  do.call(rbind,lapply(near,function(near) {
    xdata=xdata[[near]];
    d1=unique(xdata$d1);
    xdata.fpr=subset(xdata,subset=d1==0);
    do.call(rbind,lapply(d1,function(d1) {
      if ('fnr'%in%rate&d1==0) return();
      xdata.fnr=xdata[xdata$d1==d1,,drop=F];
      xdata=rbind(xdata.fpr,xdata.fnr);
      drag=data_agg(xdata=xdata,x=cq(n2),rate=rate,mesr='sig2');
      drag=data.frame(n2=drag$byx,drag[rate]);
      colnames(drag)=c('n2',rate);
      do.call(rbind,lapply(cutoff,function(cutoff) {
        n2=sapply(rate,function(rate) {
          n2rate=function(n2) predict(smooth.spline(drag$n2,drag[,rate],spar=0.5),n2)$y;
          if (sign(n2rate(50)-cutoff)==sign(n2rate(500)-cutoff)) n2=NA
          else n2=uniroot(function(n2) n2rate(n2)-cutoff,interval=c(50,500))$root;
        })
        ## data.frame(cutoff=cutoff,near=near,d1=d1,n2=t(round(n2)));
        data.frame(cutoff=cutoff,near=near,d1=d1,n2=t(round(n2)));
      }))}))}))}

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

## RE error table for Replication-wise error rates
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
