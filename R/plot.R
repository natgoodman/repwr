#################################################################################
##
## Author:  Nat Goodman
## Created: 18-05-03
##          from repwr.R restarted 18-02-15, created 17-10-05
##          from swfdr.R created 17-10-02
##          from swfdr_base.R created 17-04-10
##          from fdr.R created 17-04-03
##          from pval_error.R created 16-11-01
## Restart: 18-02-15
##          from scope.R created 17-12-04
##
## Copyright (C) 2018 Nat Goodman.
## 
## Plotting code for repwr.R
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################

## ---- Plot Functions ----

## plot rate or draw heat map vs. any of n1, n2, d1, d2
## specify drat by explicit parameter or get it from posr
## specify posr by explicit parameter, or id or from, relto smry types
## specify query (aka filter) by n1,n2,d1,d2 or xdata
##   n1,n2,d1,d2 - query is row-by-row combination
##   xdata - data.frame given desired combinations
##   nx is n multiplier: n2=nx*n1
##   d is synonym for d1, usually used when d1==d2
## rate.rule is keyword (eg, nonzro) or function that maps posr,rate.tol to logical vector
##   specify by keyword or function that maps posr,rate.tol to logical vector
## rate.type tells whether we want raw pos or neg rates, or error or correct rates
##   pos is no-op - uses posr as is; neg negates (really complements) everything
##   error negates  'true' cases to get false negative rate
##                  'false' cases already are false positive rate
##   correct negates 'false' cases to get true negative rate
##                   'true' cases already are true positive rate 
## truedd.multi tells what to do if rate implies mix of 'true' and 'false' cases
##   'error' means not allowed, 'asis' means ok but don't reorder,
##   'false.first','true.first' means move false (resp. true) cases first if x='auto'
## x tells which x variable drives the plot, or use x vars as given, or select automatically
## fpr.cutoff, fnr.cutoff tell where to draw dashed line for plot
##   or switch colors from 'lo' to 'hi' for heatmap
##   code sets cutoff to 1-cutoff for 'true' rates
## cutoff is default for pos, neg, and mixed plots
## plot.cutoff tells whether to plot cutoff at all
##   default, set in code, is T for 'error' or 'correct' rates, F for 'pos' or 'neg'
## xtitle tells whether to move single-valued vars from xdata and put in title
##   default 'no'. clearer in figures
## fignum is figure number. if not NULL "Figure fignum" prepended to title
## title.desc is additional text added to title
## xaxt (plot only) tells whether to let R compute x-axis and x-grid:
##   n means we do it; r,s,R are synonyms and mean R does it
##   CAUTION: auto works for reasonable params but not in general...
## smooth for plot tells whether to smooth data to make plot prettier
##   can be aspline, loess, none, T, F. default is aspline. T means aspline. F means none
## smooth for heatmap tells whether to reorder mesrs (y-axis) to make heatmap prettier
##   can be auto, none, 0, 1; 0 (resp. 1) mean sort toward 0 (resp. 1)
## plot.points tell whether to plot points on top of lines
## cex.points is cex for those extra points
## cex.single is cex for points drawn when there's only one x value
## vline,hline are vectors of x or y positions for extra vertical or horizontal lines
## vhlty is lty for these extra lines
##   quick hack to implement 'panels' in doc_repwr nonzro_exact_fnr section
plotrate=
  function(drat=NULL,posr=NULL,posr.id='std',
           rate.rule=cq(raw,nonzro,nonz1,nonz1or2,nonz1and2,nonz2,sameff,farzro,nearff,uni),
           rate.type=cq(error,pos,neg,correct),rate.tol=0,
           n1=n[1:5],n2=nx*n1,nx=2,d1=d,d2=d1,d=0.5,xdata=NULL,mesr=mesr.plotdflt,
           truedd.multi=c(cq(false.first,true.first,asis,error),FALSE),
           x=cq(auto,n1,n2,d1,d2,asis),xtitle=cq(none,auto,n,d,n1,n2,d1,d2),
           fpr.cutoff=parent(fpr.cutoff,.05),fnr.cutoff=parent(fnr.cutoff,0.20),cutoff=0.05,
           plot.cutoff=T,
           xaxt=cq(auto,n,r,s,R),xaxt.max=11,
           smooth=c(cq(aspline,loess,none),TRUE,FALSE),
           plot.points=F,cex.points=0.75,cex.single=1,
           vline=NULL,hline=NULL,vhlty='solid',
           title=NULL,fignum=NULL,title.desc=NULL,cex.title=0.9,ylab=NULL,xlab=NULL,
           legend.where='bottomright',x.legend=NULL,y.legend=NULL,cex.legend=0.8,
           xlim=NULL,ylim=c(0,1)) {
    ## init(must.exist=T);            # make sure environment initialized
    rate.rule=match.arg(rate.rule);
    rate.type=match.arg(rate.type);
    truedd.multi=as.character(truedd.multi);
    truedd.multi=match.arg(truedd.multi); if (truedd.multi=='FALSE') truedd.multi='error';
    x=match.arg(x);
    xtitle=match.arg(xtitle);
    xaxt=match.arg(xaxt);
    if (is.logical(smooth)) smooth=if (smooth) 'aspline' else 'none' else smooth=match.arg(smooth);
    check_mesr();
    ## if (missing(xlim)) xlim=range(x);
    if (is.null(drat)) drat=drat_order(data_rate(posr));
    xdata=drat$xdata; ydata=drat$ydata; true.dd=drat$true.dd; x=drat$x; rate.type=drat$rate.type;
    cutoff=rate_cutoff();               # set cutoff based on rate.type
    if (is.null(ylab)) ylab=ylab_rate();
    ## if (is.null(xlab)) xlab=x;
    if (x=='asis') {
      x=seq_len(nrow(xdata));
      xtitle='none';
      xaxt='n';
    } else x=xdata[,x];
    if (is.null(title)) {
      xdt=xdata_xtitle(xdata,xtitle);
      xdata=xdt$xdata; xtitle=xdt$xtitle;
      if (!is.null(fignum)) fignum=paste(sep=' ','Figure',fignum);
      title=paste(collapse="\n",c(fignum,title_rate(),title.desc,xtitle));
    }
    if (xaxt=='auto') if (nrow(xdata)<=xaxt.max) xaxt='n' else xaxt='s';
    if (xaxt=='n') {
      matplot(x,ydata,type='n',xlab=NA,ylab=ylab,main=title,cex.main=cex.title,
              xlim=xlim,ylim=ylim,xaxt='n');
      xaxis(at=x,labels=xdata,xlab=xlab);
      grid(nx=NA,ny=NULL);
      abline(v=x,lty='dotted',col='lightgray');
    } else {
      matplot(x,ydata,type='n',xlab=drat$x,ylab=ylab,main=title,cex.main=cex.title,
              xlim=xlim,ylim=ylim,xaxt='s');
      grid();
    }
    col=mesr2col(ydata);
    lwd=mesr2lwd(ydata);
    lty=mesr2lty(ydata);

    run.len=rle(true.dd)$lengths;
    run.end=cumsum(run.len);
    run.end=cumsum(run.len);
    run.start=c(1,head(run.end,-1)+1);
    run.intvl=cbind(run.start,run.end);
    apply(run.intvl,1,function(intvl) {
      i=intvl[1]:intvl[2];
      x=x[i]; ydata=ydata[i,,drop=FALSE];
      if (length(x)>1) {
        if (smooth=='none') {
          matlines(x,ydata,col=col,lty=lty,lwd=lwd);
        } else {
          ## smooth ydata so the plot will look nicer
          x.smooth=seq(min(x),max(x),len=100);
          if (smooth=='aspline') 
            y.smooth=asplinem(x,ydata,xout=x.smooth,method='improved')
          else y.smooth=loessm(x,ydata,xout=x.smooth);
          ## clamp y.smooth to [0,1]. interpolation can under- or over-shoot
          y.smooth=apply(y.smooth,c(1,2),function(y) if (!is.na(y)) max(min(y,1),0) else y);
          matlines(x.smooth,y.smooth,col=col,lty=lty,lwd=lwd);
        }
        if (plot.points) matpoints(x,ydata,col=col,pch=16,cex=cex.points);
      } else {
        if (length(x)==0) return();
        if (length(x)==1)
          ## only one x-point. draw points instead of lines. nothing to smooth obviously
          matpoints(x,ydata,col=col,pch=19,cex=cex.single);
      }});
    ## mesr_legend(ydata,where=legend.where,x=x.legend,y=y.legend);
    mesr_legend(mesr,where=legend.where,x=x.legend,y=y.legend);
    if (plot.cutoff) abline(h=cutoff,lty='dashed',lwd=0.5);
    ## plot extra lines if desired. nop if vline, hline NULL
    abline(v=vline,h=hline,lty=vhlty); 
    dev.cur();
  }
heatrate=
  function(drat=NULL,posr=NULL,posr.id='std',
           rate.rule=cq(raw,nonzro,nonz1,nonz1or2,nonz1and2,nonz2,sameff,farzro,nearff,uni),
           rate.type=cq(error,pos,neg,correct),rate.tol=0,
           n1=n[1:5],n2=nx*n1,nx=2,d1=d,d2=d1,d=0.5,xdata=NULL,mesr=mesr.heatdflt,
           ## truedd.multi=cq(false.first,true.first,asis,error),
           truedd.multi=c(cq(false.first,true.first,asis,error),FALSE),
           x=cq(auto,n1,n2,d1,d2,asis),xtitle=cq(none,auto,n,d,n1,n2,d1,d2),
           fpr.cutoff=parent(fpr.cutoff,.05),fnr.cutoff=parent(fnr.cutoff,0.20),cutoff=0.05,
           smooth=c(cq(auto,none),0,1,FALSE),
           title=NULL,fignum=NULL,title.desc=NULL,cex.title=0.9,ylab=NULL,xlab=NULL,
           vline=NULL,hline=NULL,vhlty='solid',
           legend.where='farright',x.legend=NULL,y.legend=NULL,cex.legend=0.75) {
    ## init(must.exist=T);            # make sure environment initialized
    rate.rule=match.arg(rate.rule);
    rate.type=match.arg(rate.type);
    truedd.multi=as.character(truedd.multi);
    truedd.multi=match.arg(truedd.multi); if (truedd.multi=='FALSE') truedd.multi='error';
    x=match.arg(x);
    xtitle=match.arg(xtitle);
    smooth=as.character(smooth);
    smooth=match.arg(smooth); if (smooth=='FALSE') smooth='none';
    check_mesr();
    if (is.null(drat)) drat=drat_order(data_rate(posr));
    xdata=drat$xdata; ydata=drat$ydata; true.dd=drat$true.dd; x=drat$x; rate.type=drat$rate.type;
    cutoff=rate_cutoff();               # set cutoff based on rate.type
    if (is.null(title)) {
      xdt=xdata_xtitle(xdata,xtitle);
      xdata=xdt$xdata; xtitle=xdt$xtitle;
      if (!is.null(fignum)) fignum=paste(sep=' ','Figure',fignum);
      title=paste(collapse="\n",c(fignum,title_rate(),title.desc,xtitle));
    }
    if (smooth!='none') {
      ## sort measures by distance from 0 or 1 to make nicer picture
      if (smooth=='auto') {
        if (rate.type %in% cq(neg,error,fpr,fnr)) yzero=rbind(zero=0,t(ydata))
        else yzero=rbind(zero=1,t(ydata));
      } else if (smooth=='0') yzero=rbind(zero=0,t(ydata))
      else yzero=rbind(zero=1,t(ydata))
      ydist=as.matrix(dist(yzero));
      ycol=rownames(yzero)[order(ydist['zero',])][2:ncol(ydist)];
      ydata=ydata[,ycol,drop=F];
    } else ycol=colnames(ydata);
    ## rownames(ydata)=x;
    ## log transform ydata so colors will work better
    ylog=ifelse(ydata<1e-2,2,-log10(ydata));
    x=seq_len(nrow(xdata));
    y=seq_len(ncol(ydata));
    ## setup colors and breakpoints for heatmap
    heat=heat_setup();
    brk=heat$brk; col=heat$col;
    ## initialize plot area then add heatmap. do it this way so we can add legend later
    legend.coord=heatlegend_coord(x,y);
    plot(x=c(min(x),legend.coord$x1),y=range(y),type='n',
         axes=F,main=title,cex.main=cex.title,xlab=NA,ylab=NA);
    ## bodily change max x in par('usr') to circumvent R's clever axis calculation
    par('usr'=c(par('usr')[1],legend.coord$x1,par('usr')[3:4]));
    image(x,y,ylog,add=T,axes=F,breaks=brk,col=col);

    xaxis(at=x,labels=xdata,xlab=xlab,lwd.axis=0);
    axis(2,at=y,labels=ycol,cex.axis=0.75,las=1,lwd=0,line=-0.5);
    box();
    abline(h=y+0.5,lty='dotted',col='lightgray',lwd=0.75);
    abline(v=x+0.5,lty='dotted',col='lightgray',lwd=0.75);
    heat_legend(legend.coord);
    ## plot extra lines if desired. nop if vline, hline NULL
    abline(v=vline,h=hline,lty=vhlty); 
    dev.cur();
  }
## 
## plot ROC-like graph
## plotroc plots rate vs rate for multiple measures across one data swath
## plotrocm plots rate vs rate for single measure across multiple data swaths
## specify drat by explicit parameter or get it from posr
## specify posr by explicit parameter, or id or from, relto smry types
## xrate is rate to plot on x-axis. default 'fpr'
## yrate is rate to plot on y-axis. default 'fnr'
## for plotroc, specify query (aka filter) by n1,n2,d1,d2 or xdata
##   n1,n2,d1,d2 - query is row-by-row combination
##   xdata - data.frame given desired combinations
##   nx is n multiplier: n2=nx*n1
##   d is synonym for d1, usually used when d1==d2
## for plotrocm, specify query (aka filter) by xdata
##   xdata - list of data.frames given desired combinations
## xrate is rate plotted on x-axis. default 'fpr'
## yrate is rate plotted on y-axis. default 'fnr'
## x tells which x variable to use for grouping when computing fpr,tpr points. default 'n2'
## fignum is figure number. if not NULL "Figure fignum" prepended to title
## title.desc is additional text added to title
## plot.points tells whether to plot points
## plot.lines tells whether to plot lines
## mesr is vector of mesrs for plotroc, single measure for plotrocm
plotroc=
  function(drat=NULL,posr=NULL,posr.id='std',
           rate.rule=cq(raw,nonzro,nonz1,nonz1or2,nonz1and2,nonz2,sameff,farzro,nearff,uni),
           rate.tol=0,
           n1=n[1:5],n2=nx*n1,nx=2,d1=d,d2=d1,d=0.5,mesr=mesr.rocdflt,
           xdata=expand.grid(n1=n1,n2=n2,d1=d1,d2=d2),
           ## x=cq(n2,n1,d1,d2),
           x=cq(n1,n2),
           xrate='fpr',yrate='fnr',
           fpr.cutoff=parent(fpr.cutoff,0.05),fnr.cutoff=parent(fnr.cutoff,0.20),
           tpr.cutoff=parent(tpr.cutoff,1-fnr.cutoff),tnr.cutoff=parent(tnr.cutoff,1-fpr.cutoff),
           plot.points=T,plot.lines=F,
           title=NULL,fignum=NULL,title.desc=NULL,cex.title=0.9,
           legend.where='topright',x.legend=NULL,y.legend=NULL,cex.legend=0.8,
           xlim=c(0,1),ylim=c(0,1)) {
    ## init(must.exist=T);            # make sure environment initialized
    rate.rule=match.arg(rate.rule);
    ## x=match.arg(x);
    ## TODO: change to x=pmatch_arg(x);
    check_mesr();
    ## if (is.null(drat)) drat=data_rate(posr,rate.type='pos',truedd.multi='asis');
    ## true.dd=drat$true.dd; drat=drat$drat;
    if (is.null(drat)) drat=data_rate(posr,rate.type='pos');
    ## drat.byx=split(drat,drat[,x]);
    drat.byx=split(drat,apply(drat[,x,drop=F],1,function(row) paste(collapse=' ',row)));
    if (is.null(title)) {
      if (!is.null(fignum)) fignum=paste(sep=' ','Figure',fignum);
      title=paste(collapse="\n",c(fignum,title_rate(rate.type='roc'),title.desc));
    }
    plot(x=NULL,y=NULL,type='n',xlab=rate2lab(xrate),ylab=rate2lab(yrate),
         main=title,cex.main=cex.title,xlim=xlim,ylim=ylim);
    col=col.mesr[mesr];
    cex=cex.mesr[mesr];
    x=rate2val(xrate); y=rate2val(yrate);
    if (plot.points) matpoints(x,y,col=col,cex=cex,pch=16);
    if (plot.lines) {
      lwd=lwd.mesr[mesr];
      lty=lty.mesr[mesr];
      matlines(x,y,col=col,lty=lty,lwd=lwd);
    }
    grid();
    abline(v=rate_cutoff(xrate),lty='dashed',lwd=0.5);
    abline(h=rate_cutoff(yrate),lty='dashed',lwd=0.5);
    ## abline(a=0,b=1,lty='dashed',lwd=0.5);
    mesr_legend(mesr,where=legend.where,x=x.legend,y=y.legend,plot.points=T);
    dev.cur();
  }
plotrocm=
  function(posr=NULL,posr.id='std',
           rate.rule=cq(raw,nonzro,nonz1,nonz1or2,nonz1and2,nonz2,sameff,farzro,nearff,uni),
           rate.tol=0,mesr='sig2',
           xdata,
           ## x=cq(n2,n1,d1,d2),
           x=cq(n1,n2),
           xrate='fpr',yrate='fnr',
           fpr.cutoff=parent(fpr.cutoff,0.05),fnr.cutoff=parent(fnr.cutoff,0.20),
           tpr.cutoff=parent(tpr.cutoff,1-fnr.cutoff),tnr.cutoff=parent(tnr.cutoff,1-fpr.cutoff),
           plot.points=T,plot.lines=F,
           title=NULL,fignum=NULL,title.desc=NULL,cex.title=0.9,
           legend.where='topright',x.legend=NULL,y.legend=NULL,cex.legend=0.8,
           xlim=c(0,1),ylim=c(0,1)) {
    ## init(must.exist=T);            # make sure environment initialized
    rate.rule=match.arg(rate.rule);
    ## x=match.arg(x);
    check_mesr();
    if (length(mesr)>1)
      stop(paste(sep=' ','plotrocm only plots a single measure, not',paste(collapse=', ',mesr)));
    if (is.data.frame(xdata))
      stop('plotrocm needs a list of xdata data fromes, not a single data frame');
    if (is.null(title)) {
      if (!is.null(fignum)) fignum=paste(sep=' ','Figure',fignum);
      title=paste(collapse="\n",
                  c(fignum,paste(sep=' ',title_rate(rate.type='roc'),'for',mesr),title.desc));
    }
    plot(x=NULL,y=NULL,type='n',xlab=rate2lab(xrate),ylab=rate2lab(yrate),
         main=title,cex.main=cex.title,xlim=xlim,ylim=ylim);
    n.xdata=length(xdata);
    col=colorRampPalette(c('red',RColorBrewer::brewer.pal(min(8,n.xdata-1),'Dark2')))(n.xdata);
    sapply(seq_len(n.xdata),function(i) {
      xdata=xdata[[i]];
      drat=data_rate(posr,rate.type='pos');
      ## drat.byx=split(drat,drat[,x]);
      drat.byx=split(drat,apply(drat[,x,drop=F],1,function(row) paste(collapse=' ',row)));
      x=rate2val(xrate); y=rate2val(yrate);
      if (plot.points) matpoints(x,y,col=col[i],pch=16);
      if (plot.lines) matlines(x,y,col=col[i]);
    });
    grid();
    abline(v=rate_cutoff(xrate),lty='dashed',lwd=0.5);
    abline(h=rate_cutoff(yrate),lty='dashed',lwd=0.5);
    ## abline(a=0,b=1,lty='dashed',lwd=0.5);
    rocm_legend(names(xdata),col=col,where=legend.where,x=x.legend,y=y.legend,plot.points=T);
    dev.cur();
  }
## 
## plot aggregated rates. aggregation is same as in plotroc
## plotrag plots aggregated rates for multiple measures across one data swath
## plotragm plots aggregated rates for single measure across multiple data swaths
## specify drat by explicit parameter or get it from posr
## specify posr by explicit parameter, or id or from, relto smry types
## x tells which x variable to use for grouping. default 'n2'
## rate is rates to plot on y-axis. default 'fpr','fnr'
## for plotrag, specify query (aka filter) by n1,n2,d1,d2 or xdata
##   n1,n2,d1,d2 - query is row-by-row combination
##   xdata - data.frame given desired combinations
##   nx is n multiplier: n2=nx*n1
##   d is synonym for d1, usually used when d1==d2
## for plotragm, specify query (aka filter) by xdata
##   xdata - list of data.frames giving desired combinations
## smooth tells whether to smooth data to make plot prettier
##   can be loess, aspline, none, T, F. default is loess. T means loess. F means none
##   CAUTION: loess makes prettier plots but suppresses 'waviness' caused by jumps
##            when n1 changes. to be safe, also try smooth='aspline' or plot.points
## fignum is figure number. if not NULL "Figure fignum" prepended to title
## title.desc is additional text added to title
## plot.points tells whether to plot points
## plot.lines tells whether to plot lines
## mesr is vector of mesrs for plotroc, single measure for plotrocm
plotrag=
  function(drat=NULL,posr=NULL,posr.id='std',
           rate.rule=cq(raw,nonzro,nonz1,nonz1or2,nonz1and2,nonz2,sameff,farzro,nearff,uni),
           rate.tol=0,
           n1=n[1:5],n2=nx*n1,nx=2,d1=d,d2=d1,d=0.5,mesr=mesr.ragdflt,
           xdata=expand.grid(n1=n1,n2=n2,d1=d1,d2=d2),
           ## x=cq(n2,n1,d1,d2),
           x=cq(n1,n2),
           rate=cq(fpr,fnr),
           fpr.cutoff=parent(fpr.cutoff,0.05),fnr.cutoff=parent(fnr.cutoff,0.20),
           tpr.cutoff=parent(tpr.cutoff,1-fnr.cutoff),tnr.cutoff=parent(tnr.cutoff,1-fpr.cutoff),
           smooth=c(cq(loess,aspline,none),TRUE,FALSE),
           plot.points=F,plot.lines=T,
           title=NULL,fignum=NULL,title.desc=NULL,cex.title=0.9,xlab=NULL,ylab='rate',
           legend.where='topright',x.legend=NULL,y.legend=NULL,cex.legend=0.8,
           xlim=NULL,ylim=c(0,1)) {
    ## init(must.exist=T);            # make sure environment initialized
    rate.rule=match.arg(rate.rule);
    ## x=match.arg(x);
    ## TODO: change to x=pmatch_arg(x);
    if (is.null(xlab)) xlab=NA;
    if (is.logical(smooth)) smooth=if (smooth) 'loess' else 'none' else smooth=match.arg(smooth);
    check_mesr();
    ## if (is.null(drat)) drat=data_rate(posr,rate.type='pos',truedd.multi='asis');
    ## true.dd=drat$true.dd; drat=drat$drat;
    if (is.null(drat)) drat=data_rate(posr,rate.type='pos');
    ## drat.byx=split(drat,drat[,x]);
    drat.byx=split(drat,apply(drat[,x,drop=F],1,function(row) paste(collapse=' ',row)));
    ## x=as.numeric(names(drat.byx));
    labels=as.data.frame(apply(do.call(rbind,strsplit(names(drat.byx),' ')),2,as.numeric));
    colnames(labels)=x;
    ## line below from StackExchange https://stackoverflow.com/questions/29482983/. Thx!!
    i=do.call(order,labels);
    drat.byx=drat.byx[i];
    labels=labels[i,];
    ## x=seq_along(names(drat.byx));
    x=seq_along(i);
    y=sapply(rate,function(rate) rate2val(rate),simplify=F);
    if (is.null(title)) {
      if (!is.null(fignum)) fignum=paste(sep=' ','Figure',fignum);
      title=paste(collapse="\n",c(fignum,title_rate(rate.type='rag'),title.desc));
    }
    if (is.null(xlim)) xlim=range(x);
    plot(x=NULL,y=NULL,type='n',xlab=xlab,ylab=ylab,main=title,cex.main=cex.title,
         xlim=xlim,ylim=ylim,xaxt='n');
    xaxis(at=x,labels=labels,xlab=xlab);
    col=col.mesr[mesr];
    cex=cex.mesr[mesr];
    lty=cq(solid,dashed,dotted,dotdash);
    sapply(seq_along(rate),function(i) {
      y=y[[i]];
      if (length(x)>1) {
        if (plot.lines) {
          lwd=lwd.mesr[mesr];
          lty=lty[i];
          if (smooth=='none') {
            matlines(x,y,col=col,lwd=lwd,lty=lty);
          } else {
          ## smooth ydata so the plot will look nicer
          x.smooth=seq(min(x),max(x),len=100);
            if (smooth=='aspline') 
              y.smooth=asplinem(x,y,xout=x.smooth,method='improved')
            else y.smooth=loessm(x,y,xout=x.smooth);
            ## clamp y.smooth to [0,1]. interpolation can under- or over-shoot
            y.smooth=apply(y.smooth,c(1,2),function(y) if (!is.na(y)) max(min(y,1),0) else y);
            matlines(x.smooth,y.smooth,col=col,lwd=lwd,lty=lty);
          }}}
      if (plot.points|length(x)==1) matpoints(x,y,col=col,cex=cex,pch=16);
    })
    grid();
    rate.cutoff=sapply(rate,rate_cutoff);
    abline(h=rate.cutoff,lty='dashed',lwd=0.5);
    rag_legend(mesr,rate,rate.lty=lty,where=legend.where,x=x.legend,y=y.legend);
    dev.cur();
  }
plotragm=
  function(posr=NULL,posr.id='std',
           rate.rule=cq(raw,nonzro,nonz1,nonz1or2,nonz1and2,nonz2,sameff,farzro,nearff,uni),
           rate.tol=0,mesr='sig2',
           xdata,
           ## x=cq(n2,n1,d1,d2),
           x=cq(n1,n2),
           rate=cq(fpr,fnr),
           fpr.cutoff=parent(fpr.cutoff,0.05),fnr.cutoff=parent(fnr.cutoff,0.20),
           tpr.cutoff=parent(tpr.cutoff,1-fnr.cutoff),tnr.cutoff=parent(tnr.cutoff,1-fpr.cutoff),
           smooth=c(cq(loess,aspline,none),TRUE,FALSE),
           plot.points=F,plot.lines=T,lty=cq(solid,dotted,dotdash,longdash),lwd=2,
           title=NULL,fignum=NULL,title.desc=NULL,cex.title=0.9,xlab=NULL,
           ylab=if(length(rate)==1) rate2lab(rate) else 'rate',
           legend.where='topright',x.legend=NULL,y.legend=NULL,cex.legend=0.8,
           xlim=NULL,ylim=c(0,1)) {
    ## init(must.exist=T);            # make sure environment initialized
    rate.rule=match.arg(rate.rule);
    ## x=match.arg(x);
    if (is.null(xlab)) xlab=NA;
    if (is.logical(smooth)) smooth=if (smooth) 'loess' else 'none' else smooth=match.arg(smooth);
    check_mesr();
    if (length(mesr)>1)
      stop(paste(sep=' ','plotragm only plots a single measure, not',paste(collapse=', ',mesr)));
    if (is.data.frame(xdata))
      stop('plotragm needs a list of xdata data fromes, not a single data frame');
    if (length(rate)==1) rate.desc=rate2lab(rate);
    if (is.null(title)) {
      if (!is.null(fignum)) fignum=paste(sep=' ','Figure',fignum);
      if (length(rate)==1) title.rate=paste(sep=' ',title_rate(rate.type='ragm'),ylab)
      else title.rate=paste(sep=' ',title_rate(rate.type='ragm'),'rate');
      title=paste(collapse="\n",c(fignum,paste(sep=' ',title.rate,'for',mesr),title.desc));
    }
    n.xdata=length(xdata);
    ## collect all labels and arrange along x-axis
    labels=unique(do.call(c,lapply(xdata,function(xdata) 
      unique(apply(xdata[,x,drop=F],1,function(row) paste(collapse=' ',row))))));
    labels=as.data.frame(apply(do.call(rbind,strsplit(labels,' ')),2,as.numeric));
    colnames(labels)=x;
    ## line below from StackExchange https://stackoverflow.com/questions/29482983/. Thx!!
    i=do.call(order,labels);
    labels=labels[i,,drop=F];
    if (is.null(xlim)) xlim=range(i);
    plot(x=NULL,y=NULL,type='n',xlab=xlab,ylab=ylab,main=title,cex.main=cex.title,
         xlim=xlim,ylim=ylim,xaxt='n');
    xaxis(at=1:nrow(labels),labels=labels,xlab=xlab);
    n.xdata=length(xdata);
    col=colorRampPalette(c('red',RColorBrewer::brewer.pal(min(8,n.xdata-1),'Dark2')))(n.xdata);
    sapply(seq_len(n.xdata),function(i) {
      xdata=xdata[[i]];
      col=col[i];
      drat=data_rate(posr,rate.type='pos');
      ## drat.byx=split(drat,drat[,x]);
      drat.byx=split(drat,apply(drat[,x,drop=F],1,function(row) paste(collapse=' ',row)));
      labels=as.data.frame(apply(do.call(rbind,strsplit(names(drat.byx),' ')),2,as.numeric));
      colnames(labels)=x;
      ## line below from StackExchange https://stackoverflow.com/questions/29482983/. Thx!!
      i=do.call(order,labels);
      drat.byx=drat.byx[i];
      x=1:length(i);
      y=sapply(rate,function(rate) rate2val(rate),simplify=F);
      sapply(seq_along(rate),function(i) {
        y=y[[i]];
        if (length(x)>1) {
          if (plot.lines) {
            lty=lty[i];
            if (smooth=='none') {
              matlines(x,y,col=col,lty=lty,lwd=lwd);
            } else {
              ## smooth ydata so the plot will look nicer
              x.smooth=seq(min(x),max(x),len=100);
              if (smooth=='aspline') 
                y.smooth=asplinem(x,y,xout=x.smooth,method='improved')
              else y.smooth=loessm(x,y,xout=x.smooth);
              ## clamp y.smooth to [0,1]. interpolation can under- or over-shoot
              y.smooth=apply(y.smooth,c(1,2),function(y) if (!is.na(y)) max(min(y,1),0) else y);
              matlines(x.smooth,y.smooth,col=col,lty=lty,lwd=lwd);
            }}}
      if (plot.points|length(x)==1) matpoints(x,y,col=col,cex=cex,pch=16);
      })})
    grid();
    rate.cutoff=sapply(rate,rate_cutoff);
    abline(h=rate.cutoff,lty='dashed',lwd=0.5);
    ragm_legend(names(xdata),rate,col=col,lty=lty,where=legend.where,x=x.legend,y=y.legend);
    dev.cur();
  }
   
## remove single-valued xvars from xdata and put in title
xdata_xtitle=function(xdata,xtitle) {
  if (xtitle=='none') xtitle=NULL
  else {
    ## find single valued xvars we want
    xvars=colnames(xdata);          # start with all of 'em
    if (xtitle=='n') {xvars=grep('^n',xvars,value=T)}
    else if (xtitle=='d') {xvars=grep('^d',xvars,value=T)}
    else if (xtitle!='auto') {xvars=xtitle};
    x.single=do.call(
      c,sapply(xvars,simplify=F,
               function(x) {xval=unique(xdata[[x]]); if (length(xval)==1) xval}));
    ## remove single-valued vars from xdata and put in title
    if (!is.null(x.single)) {
      ## remove vars from xdata without converting one column data frame into vector (sigh...)
      ## from https://stackoverflow.com/questions/4605206/drop-data-frame-columns-by-name
      xdata=within(xdata,rm(list=names(x.single)));
      ## coalesce pairs with same values
      x.n=x.single[grep('^n',names(x.single),value=T)];
      x.d=x.single[grep('^d',names(x.single),value=T)];
      if (length(x.n)==2 && length(unique(x.n))==1) x.n=c(n=unique(x.n));
      if (length(x.d)==2 && length(unique(x.d))==1) x.d=c(d=unique(x.d));
      ## format values as we want them
      x.n=sapply(x.n,n_pretty);
      x.d=sapply(x.d,d_pretty);
      x.single=c(x.n,x.d);
      ## finally turn it all into a string
      xtitle=paste(collapse=', ',paste(sep='=',names(x.single),x.single));
    } else {
      xtitle=NULL;
    }}
  list(xdata=xdata,xtitle=xtitle);
}


## determine cutoff for rate type
rate_cutoff=
  function(type=parent(rate.type),
           fpr.cutoff=parent(fpr.cutoff,.05),fnr.cutoff=parent(fnr.cutoff,0.20),
           tpr.cutoff=parent(tpr.cutoff,1-fnr.cutoff),tnr.cutoff=parent(tnr.cutoff,1-fpr.cutoff),
           cutoff=parent(cutoff,.05)) {
    switch(type,
           neg=cutoff,pos=1-cutoff,error=cutoff,correct=1-cutoff,
           fpr=fpr.cutoff,fnr=fnr.cutoff,tpr=1-fnr.cutoff,tnr=tnr.cutoff);
  }
## construct colors and breaks for heatmap
heat_setup=function(type=parent(rate.type),cutoff=parent(cutoff),steps=100) {
  hi.brk=seq(0,-log10(cutoff),length.out=(steps)+1);
  lo.brk=seq(-log10(cutoff),4,length.out=(steps)+1);
  blues=colorRampPalette(RColorBrewer::brewer.pal(5,'Blues'))(steps);
  reds=colorRampPalette(RColorBrewer::brewer.pal(4,'Reds'))(steps);
  brk=unique(c(hi.brk,lo.brk));
  col=c(rev(blues),reds);
  ## if (type %in% cq(neg,error,fpr,fnr)) col=c(blues,reds) else col=c(reds,blues); 
  list(brk=brk,col=col);
}
## generate 'rate' part of title for plots
title_rate=
  function(posr.id=parent(posr.id,'std'),rate.rule=parent(rate.rule,'raw'),
           rate.type=parent(rate.type,'error'),rate.tol=parent(rate.tol,0)) {
    if (is.function(rate.rule)) desc='user-defined rule'
    else {
      if (rate.rule=='raw') desc='raw'
      else {
        if (rate.rule=='nonzro') rate.rule='nonz1';
        desc=switch(rate.rule,
                    nonz1=if (rate.tol==0) 'study1 non-zero'
                           else paste(sep=' ','study1 >',rate.tol),
                    nonz1or2=if (rate.tol==0) 'study1 or study2 non-zero'
                           else paste(sep=' ','study1 or study2 >',rate.tol),
                    nonz1and2=if (rate.tol==0) 'study1 and study2 non-zero'
                           else paste(sep=' ','study1 and study2 >',rate.tol),
                    nonz2=if (rate.tol==0) 'study2 non-zero'
                           else paste(sep=' ','study2 >',rate.tol),
                    sameff=if (rate.tol==0) 'both studies have same effect'
                           else paste(sep=' ','effects differ by no more than',rate.tol),
                    ## for backwards compatibility
                    farzro=paste(sep=' ','study1 at least',rate.tol),
                    nearff=paste(sep=' ','effects differ by no more than',rate.tol),
                    ## not used
                    uni=paste(sep=' ','study1 at least',rate.tol,
                              'and effects differ by no more than',rate.tol));
        rule=switch(rate.rule,
                    nonz1='non-zero-1',
                    nonz1or2='non-zero-1or2',
                    nonz1and2='non-zero-1and2',
                    nonz2='non-zero-2',
                    sameff='same-effect',
                    ## for backwards compatibility
                    rate.rule);
        desc=paste(sep='',rule,' (',desc,')');
      }}
    rate.desc=switch(rate.type,
                     pos='positive',neg='negative',error='error',correct='correct',
                     fpr='false positive',fnr='false negative',
                     tpr='true positive',tnr='true negative',
                     roc='rate vs rate',rag='mean',ragm='mean');
    if (rate.type %notin% cq(roc,ragm)) rate.desc=paste(sep=' ',rate.desc,'rate');
    posr.desc=if (posr.id=='std') NULL else paste_nv('posr',posr.id);
    paste(collapse=' ',c(desc,rate.desc,posr.desc));
  }
## generate 'rate' part of ylab for plots
ylab_rate=function(rate.rule=parent(rate.rule),rate.type=parent(rate.type)) {
  if (is.function(rate.rule)) rate.rule='user-defined';
  if (rate.type %notin% cq(fpr,fnr,tpr,tnr)) rate.type=paste(sep=' ',rate.type,'rate');
  paste(sep=' ',rate.rule,rate.type);
}
## generate labels and x, y values (rates) for plotroc, plotrocm
rate2lab=function(rate)
  switch(rate,
         fpr='false positive rate',tpr='true positive rate',
         fnr='false negative rate',tnr='true negative rate',
         stop(paste(sep='','Unknown rate ',rate,'. Should be one of ',
                    paste(collapse=', ',cq(fpr,fnr,tpr,tnr)))));

rate2val=function(rate,drat.byx=parent(drat.byx),mesr=parent(mesr))
  switch(rate,
         fpr=do.call(rbind,
                     lapply(drat.byx,
                            function(drat) colMeans(drat[!drat$true.dd,mesr,drop=F],na.rm=T))),
         tnr=do.call(rbind,
                     lapply(drat.byx,
                            function(drat) colMeans(1-drat[!drat$true.dd,mesr,drop=F],na.rm=T))),
         tpr=do.call(rbind,
                     lapply(drat.byx,
                            function(drat) colMeans(drat[drat$true.dd,mesr,drop=F],na.rm=T))),
         fnr=do.call(rbind,
                     lapply(drat.byx,
                            function(drat) colMeans(1-drat[drat$true.dd,mesr,drop=F],na.rm=T))),
         stop(paste(sep='','Unknown rate ',rate,'. Should be one of ',
                    paste(collapse=', ',cq(fpr,fnr,tpr,tnr)))));
 
## get line properties for measures in matrix
mesr2col=function(matrix) col.mesr[colnames(matrix)];
mesr2lwd=function(matrix) lwd.mesr[colnames(matrix)];
mesr2lty=function(matrix) lty.mesr[colnames(matrix)];
mesr2cex=function(matrix) cex.mesr[colnames(matrix)];
## construct legend for measures in plotrate
mesr_legend=
  function(mesr,col=col.mesr[mesr],lwd=lwd.mesr[mesr],lty=lty.mesr[mesr],pt.cex=cex.mesr[mesr],
           where='bottomright',x=NULL,y=NULL,title='measure',
           bty='n',pch=16,cex=parent(cex.legend,0.8),title.col='black',
           plot.lines=T,plot.points=F) {
    if (missing(plot.lines)&!missing(plot.points)) {plot.lines=!plot.points}
    else {if (!missing(plot.lines)&missing(plot.points)) plot.points=!plot.lines;}
    if (!plot.lines) lty=NA;
    if (!plot.points) {pch=NA; pt.cex=NA}
    ## legend.col=col.mesr[legend.text];
    if (is.null(x)) x=where;
    seg.len=if (plot.lines) 6 else 2;
    legend(x,y,bty=bty,legend=mesr,cex=cex,
           title=title,title.col=title.col,
           pch=pch,pt.cex=pt.cex,col=col,lwd=lwd,lty=lty,seg.len=seg.len);
  }
## construct legend for rate colors in heatrate
heat_legend=
  function(coord,cutoff=parent(cutoff,0.05),heat=parent(heat),
           cex=parent(cex.legend,0.75),steps=100) {
    x0=coord$x0; y0=coord$y0; width=coord$width; height=coord$height;
    brk=heat$brk; col=heat$col;
    x=c(x0,x0+width);
    y=seq(y0,y0+height,length.out=2*steps+1)[1:(2*steps)];
    z=t(as.matrix(rev(head(brk,-1))));
    image(x,y,z,add=T,breaks=brk,col=col);
    ## add legend text
    text(x0+width/2,y0+height,"rates",pos=3,cex=cex);
    text(x0,y0,'0 ',adj=1,cex=cex);
    text(x0,y0+height/2,cutoff,adj=1,cex=cex);
    text(x0,y0+height,'1 ',adj=1,cex=cex);
  }
## compute coordinates for heat_legend
## return list of
##   x0,y0 - starting coords of legend image
##   width,height of legend image
##   x1 - end coord of main plat
heatlegend_coord=
  function(x,y,cutoff=parent(cutoff,0.05),width=0.1*max(x),
           yq=if(length(y)>1) quantile(y,prob=c(0.1,0.9)) else c(0.75,1.25),
           height=diff(yq),cex=parent(cex.legend,0.75)) {
    ## call plot to setup coordinate systme. needed by strwidth
    plot(x=c(min(x),1.2*max(x)),y=range(y),type='n',axes=F,xlab=NA,ylab=NA);
    ## add space to labels to get a smidge more room
    x0=max(x)+0.5+max(strwidth(c('0','1',cutoff),cex=cex))+strwidth(' ',cex=cex);
    y0=yq[1];
    x1=x0+width+strwidth(' ',cex=cex);
    list(x0=x0,y0=y0,width=width,height=height,x1=x1);
  }
## construct legend for plotrocm
rocm_legend=
  function(label,col,where='bottomright',x=NULL,y=NULL,title='replication type',
           bty='n',pch=16,cex=parent(cex.legend,0.8),title.col='black',
           lty='solid',lwd=1,pt.cex=1,
           plot.lines=T,plot.points=F) {
    if (missing(plot.lines)&!missing(plot.points)) {plot.lines=!plot.points}
    else {if (!missing(plot.lines)&missing(plot.points)) plot.points=!plot.lines;}
    if (is.null(x)) x=where;
    if (!plot.lines) lty=NA;
    seg.len=if (plot.lines) 4 else 2;
    if (!plot.points) pt.cex=NA;
    legend(x,y,bty=bty,legend=label,cex=cex,
           title=title,title.col=title.col,
           pch=pch,pt.cex=pt.cex,col=col,lwd=lwd,lty=lty,seg.len=seg.len);
  }
## construct legend for plotrag
rag_legend=
  function(mesr,rate,rate.lty=rate.lty,
           col=col.mesr[mesr],lwd=lwd.mesr[mesr],mesr.lty='solid',pt.cex=cex.mesr[mesr],
           where='bottomright',x=NULL,y=NULL,title='measure',
           bty='n',pch=16,cex=parent(cex.legend,0.8),title.col='black',
           plot.lines=T,plot.points=F,rate.legend=length(rate)>1) {
    if (missing(plot.lines)&!missing(plot.points)) {plot.lines=!plot.points}
    else {if (!missing(plot.lines)&missing(plot.points)) plot.points=!plot.lines;}
    if (!plot.lines) lty=NA;
    if (!plot.points) {pch=NA; pt.cex=NA}
    where.next=mesr_legend(
      mesr,col=col,lwd=lwd,lty=mesr.lty,pt.cex=pt.cex,where=where,x=x,y=y,title=title,
      bty=bty,pch=pch,cex=cex,title.col=title.col,plot.lines=plot.lines,plot.points=plot.points);
    if (!rate.legend) return(where.next);
    x=where.next$rect$left;
    y=where.next$rect$top-where.next$rect$h;
    seg.len=if (plot.lines) 6 else 2;
    legend(x,y,bty=bty,legend=rate,cex=cex,title='rate type',title.col=title.col,
           pch=pch,pt.cex=pt.cex,col='black',lwd=lwd,lty=rate.lty,seg.len=seg.len);
  }
## construct legend for plotragm
ragm_legend=
  function(label,rate,col,lty,where='bottomright',x=NULL,y=NULL,title='replication type',
           bty='n',pch=16,cex=parent(cex.legend,0.8),title.col='black',
           lwd=2,pt.cex=1,
           plot.lines=T,plot.points=F,rate.legend=length(rate)>1) {
    if (missing(plot.lines)&!missing(plot.points)) {plot.lines=!plot.points}
    else {if (!missing(plot.lines)&missing(plot.points)) plot.points=!plot.lines;}
    where.next=rocm_legend(
      label=label,col=col,where=where,x=x,y=y,title=title,
      bty=bty,pch=pch,cex=cex,title.col=title.col,lty='solid',lwd=lwd,pt.cex=pt.cex,
      plot.lines=plot.lines,plot.points=plot.points);
    if (!rate.legend) return(where.next);
    x=where.next$rect$left;
    y=where.next$rect$top-where.next$rect$h;
    if (!plot.lines) lty=NA;
    if (!plot.points) {pch=NA; pt.cex=NA}
    seg.len=if (plot.lines) 6 else 2;
    legend(x,y,bty=bty,legend=rate,cex=cex,title='rate type',,title.col=title.col,
           pch=pch,pt.cex=pt.cex,col='black',lwd=lwd,lty=lty,seg.len=seg.len);
  }
   
## annotate x axis
##  at is vector of x positions for ticks and labels. if NULL, no ticks or label
##  labels is vector or matrix. if matrix, columns are lines of annotion
##  xlab is overall label for axis
##  names is vector of names for annotation lines. if T, use colnames of labels
## 
xaxis=function(at,labels=NULL,xlab=NULL,tick=T,names=T,
           tcl=-0.4,
           line.labels=seq(.4,by=0.6,len=NCOL(labels)),
           line.xlab=min(max(line.labels)+1.2,par('mar')[1]-1.2),
           cex=1,col=par('fg'),lty='solid',lwd=1,cex.labels=0.75*cex,cex.xlab=cex,
           col.axis=col,col.labels=col,col.xlab=col,col.ticks=col,
           lwd.axis=lwd,lwd.ticks=lwd) {
    labels.axis=if (is.logical(labels)) labels else NA;
    axis(side=1,at=at,labels=labels.axis,col=col.axis,lty=lty,lwd=lwd.axis,
         tick=tick,tcl=tcl,col.ticks=col.ticks,lwd.ticks=lwd.ticks);
    if (is.na(labels.axis)) {
      ## emit each label line 
      if (is.vector(labels)) {
        mtext(labels,at=at,line=line.labels[1],cex=cex.labels,col=col.labels,side=1)
      } else {
        if (names&!is.null(colnames(labels))) {
          ## place annotation labels in left margin
          ## formula computes (min.plot-0.05*range of plot. empirically looks good
          at=c((1.05*par('usr')[1]-0.05*par('usr')[2]),at);
          ## at=c(par('usr')[1]-0.4,at); 
          labels=rbind(colnames(labels),labels);
        }
        for (j in seq_len(ncol(labels))) 
          mtext(labels[,j],at=at,line=line.labels[j],cex=cex.labels,col=col.labels,side=1);
      }
    }
    ## add overall label (hopefully) at bottom
    if (!is.null(xlab)) mtext(xlab,line=line.xlab,cex=cex.xlab,col=col.xlab,side=1);
  }

