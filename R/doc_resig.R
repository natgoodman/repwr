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
            title.desc=title.desc,legend=NULL);
      xdata=lapply(d.nonzro,function(d)
        xdata=expand.grid(n1=20,n2=seq(50,by=50,len=10),d1=d,d2=d));
      names(xdata)=as.character(d.nonzro);
      dofig(plotratm,'fnr',xdata=xdata,x=cq(n1,n2),col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title.desc=title.desc,title.legend='d',legend='right');
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
      near=round(c(0,0.05,seq(0.1,0.5,by=0.1)),digits=5);
      xdata=lapply(near,function(near) {
        d2=round(seq(0,near,by=0.01),digits=5);
        expand.grid(n1=20,n2=seq(50,by=50,len=10),d1=0,d2=d2)});
      names(xdata)=as.character(near);
      dofig(plotragm,'fpr',xdata=xdata,x=cq(n1,n2,d1),rate='fpr',col=col,smooth='spline',
            hline=c(fpr.cutoff,fnr.cutoff),vhlty='dashed',vhlwd=0.5,plot.cutoff=F,
            title.desc=title.desc,title.legend='near',legend='topright');

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
  function(lo.brk=0.1,hi.brk=0.3,res=0.01,
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
