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

## plot rate vs. any of n1, n2, d1, d2
## operates on full smry unless smryx set
##   default for n1 10, 20, 40, 80, 160
##   nx is n multiplier: n2=nx*n1
##   x tells which x variable drives the plot, or to use x vars as given
##     x='asis' also sets xtitle='none' and xaxt='n' so xaxis display will make sense
##   rate is the kind of rate we're plotting, eg, 'nonz1' or 'sameff'
##   error tells whether to plot error or correct rates
##   yesno.multiok tells whether okay to mix 'yes' and 'no' cases on same plot
##     eg, false positives and false negatives
##   fpr.cutoff, fnr.cutoff tell where to draw dashed line for 'good' rate, if plot.cutoff=T
##     code sets cutoff to 1-cutoff for 'correct' rates
##   cutoff is default for pos, neg, and mixed plots
##   plot.cutoff tells whether to plot cutoff at all
##     default, set in code, is T for 'error' or 'correct' rates, F for 'pos' or 'neg'
##   xtitle tells whether to move single-valued vars from xdata and put in title
##   xaxt tells whether to let R compute x-axis and x-grid:
##     n means we do it; r,s,R are synonyms and mean R does it
##     CAUTION: auto works for reasonable params but not in general...
##   smooth tells whether to smooth data to make plots prettier
##     can be T, F, aspline, loess, none. default (T) uses aspline
##   plot.points tell whether to plot points on top of lines
##   cex.points is cex for those extra points
##   cex.single is cex for points drawn when there's only one x value
##   d is synonym for d1, usually used when d1==d2
##   from.type, relto.type are smry types. single values or vectors keyed by mesr
##     from.type is smry source row
##     relto.type provides denominator in rate calculation
##     default is 'standard' interpretation
##   relto.multiok tells whether okay to plot mesrs with different relto.types
##     CAUTION: rates not comparable in this case!
plotrate_nndd=
  function(smryx,n1=n[1:5],n2=nx*n1,nx=2,d1=d,d2=d1,d=0.5,mesr=mesr.plotdflt,
           rate=cq(pos,neg,nonzro,nonz1,nonz1or2,nonz1and2,sameff,uni),
           error.rate=T,yesno.multiok=T,
           x=cq(auto,n1,n2,d1,d2,asis),xtitle=cq(auto,none,n,d,n1,n2,d1,d2),
           fpr.cutoff=parent(fpr.cutoff,.05),fnr.cutoff=parent(fnr.cutoff,0.20),cutoff=0.05,
           plot.cutoff=T,
           xaxt=cq(auto,n,r,s,R),xaxt.max=11,smooth=T,plot.points=F,cex.points=0.75,cex.single=1,
           from.type=parent(mesr.fromtype,'bsln'),relto.type=parent(mesr.reltotype,'sig1'),
           relto.multiok=F,
           title=NULL,cex.title=0.9,ylab=NULL,xlab=NULL,
           legend.where='bottomright',x.legend=NULL,y.legend=NULL,cex.legend=0.8,
           xlim=NULL,ylim=c(0,1)) {
    init(must.exist=T);            # make sure environment initialized
    rate=match.arg(rate);
    x=match.arg(x);
    xtitle=match.arg(xtitle);
    xaxt=match.arg(xaxt);
    if (is.logical(smooth)) smooth=if (smooth) 'aspline' else 'none'
    else smooth=match.arg(smooth,cq(aspline,loess,none));
    ## make sure measures & types legal and limit to type we need
    from.type=check_type(mesr,from.type,multiok=T);
    relto.type=check_type(mesr,relto.type,multiok=relto.multiok);
    check_mesr();
    ## get the smry subset for the parameters unless smryx passed as paramter
    if (missing(smryx)) {
      smry=get_data(smry);
      smryx=smry_select(smry,n1,n2,d1,d2,x=x);
    }
    smry=smryx$smry;
    ## if (missing(xlim)) xlim=range(x);
    data=data_rate(smry,mesr,rate,error.rate,from.type,relto.type);
    xdata=data$x; ydata=data$y; yes=data$yes;
    if (length(unique(yes))>1&&!yesno.multiok)
      stop("Cannot mix 'yes' and 'no' cases (eg, false positives and false negatives) on same plot unless yesno.multiok is TRUE");
    ## get rate type, eg, false positive or negative, and set cutoff
    rate.type=rate_type();
    cutoff=rate_cutoff();
    if (!is.null(yes)) {
      ## deal with 'yes' vs. 'no' cases. want 'no' cases first - these are false positives
      yes.run=rle(yes)$lengths;
      if (length(yes.run)==1|(length(yes.run)==2&yes[1]==F))
        ## in correct order already. use x ordering from smryx
        x=smryx$x
      else {
        ## sort 'no' cases first if allowed (ie, x='auto')
        if (x=='auto') {
          i=order(yes);
          xdata=xdata[i,]; ydata=ydata[i,]; yes=yes[i];
          x='asis'
        } else {
          stop('In line plots, all false positive cases have to be first unless x="auto"');
        }}
    } else {
      ## it's pos or neg - yes/no not relevant. use x ordering from smryx
      yes=rep(T,nrow(xdata));
      x=smryx$x;
    }
    if (is.null(ylab)) ylab=ylab.rate();
    ## if (is.null(xlab)) xlab=x;
    if (x=='asis') {
      x=seq_len(nrow(xdata));
      xtitle='none';
      xaxt='n';
    } else x=xdata[,x];
    if (is.null(title)) {
      xdt=xdata_xtitle(xdata,xtitle);
      xdata=xdt$xdata; xtitle=xdt$xtitle;
      title=paste(collapse="\n",c(title_rate(),xtitle));
    }
    if (xaxt=='auto') if (nrow(xdata)<=xaxt.max) xaxt='n' else xaxt='s';
    if (xaxt=='n') {
      matplot(x,ydata,type='n',xlab=NA,ylab=ylab,main=title,cex.main=cex.title,
              xlim=xlim,ylim=ylim,xaxt='n');
      xaxis(at=x,labels=xdata,xlab=xlab);
      grid(nx=NA,ny=NULL);
      abline(v=x,lty='dotted',col='lightgray');
    } else {
      matplot(x,ydata,type='n',xlab=xlab,ylab=ylab,main=title,cex.main=cex.title,
              xlim=xlim,ylim=ylim,xaxt='s');
      grid();
    }
    col=mesr2col(ydata);
    lwd=mesr2lwd(ydata);
    lty=mesr2lty(ydata);
    sapply(c(F,T),function(want) {
      ## if (want %notin% yes) return();
      x=x[yes==want]; ydata=ydata[yes==want,,drop=FALSE];
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
          y.smooth=apply(y.smooth,c(1,2),function(y) max(min(y,1),0));
          matlines(x.smooth,y.smooth,col=col,lty=lty,lwd=lwd);
        }
        if (plot.points) matpoints(x,ydata,col=col,pch=16,cex=cex.points);
      } else {
        if (length(x)==0) return();
        if (length(x)==1)
          ## only one x-point. draw points instead of lines. nothing to smooth obviously
          matpoints(x,ydata,col=col,pch=19,cex=cex.single);
      }});
    mesr_legend(ydata,where=legend.where,x=x.legend,y=y.legend);
    if (plot.cutoff) abline(h=cutoff,lty='dashed',lwd=0.5)
    dev.cur();
  }
## wrapper for plotrate_nndd which gets params from smry
plotrate_smry=function(smry,...) plotrate_nndd(n1=smry$n1,n2=smry$n2,d1=smry$d1,d2=smry$d2,...)

## draw heat map of rate vs. any of n1, n2, d1, d2
## operates on full smry unless smryx set
##   defaults set n1,n2,d1 to fixed values and d2 to d
##   nx is n multiplier: n2=nx*n1
##   rate is the kind of rate we're plotting, eg, 'nonz1' or 'sameff'
##   error tells whether to plot error or correct rates
##   yesno.multiok tells whether okay to mix 'yes' and 'no' cases on same plot
##     eg, false positives and false negatives
##   fpr.cutoff, fnr.cutoff tell where to switch colors from 'lo' to 'hi'
##     code sets cutoff to 1-cutoff for 'correct' rates
##   cutoff is default for pos, neg, and mixed plots
##   xtitle tells whether to move single-valued vars from xdata and put in title
##   d is synonym for d1, usually used when d1==d2
##   from.type, relto.type are smry types. single values or vectors keyed by mesr
##     from.type is smry source row
##     relto.type provides denominator in rate calculation
##     default is 'standard' interpretation
##   relto.multiok tells whether okay to plot mesrs with different relto.types
##     CAUTION: rates not comparable in this case!
heatrate_nndd=
  function(smryx,n1=n[1:5],n2=nx*n1,nx=2,d1=d,d2=d1,d=0.5,mesr=mesr.heatdflt,
           rate=cq(pos,neg,nonzro,nonz1,nonz1or2,nonz1and2,sameff,uni),
           error.rate=T,yesno.multiok=T,
           x=cq(asis,auto,n1,n2,d1,d2),xtitle=cq(auto,none,n,d,n1,n2,d1,d2),y=cq(auto,asis),
           fpr.cutoff=parent(fpr.cutoff,.05),fnr.cutoff=parent(fnr.cutoff,0.20),cutoff=0.05,
           ## xaxt=cq(n,auto,r,s,R),xaxt.max=11,
           from.type=parent(mesr.fromtype,'bsln'),relto.type=parent(mesr.reltotype,'sig1'),
           relto.multiok=F,
           title=NULL,cex.title=0.9,xlab=NULL,
           legend.where='farright',x.legend=NULL,y.legend=NULL,cex.legend=0.75) {
    init(must.exist=T);            # make sure environment initialized
    rate=match.arg(rate);
    x=match.arg(x);
    xtitle=match.arg(xtitle);
    y=match.arg(y);
    ## xaxt=match.arg(xaxt);
     ## make sure measures & types legal and limit to type we need
    from.type=check_type(mesr,from.type,multiok=T);
    relto.type=check_type(mesr,relto.type,multiok=relto.multiok);
    check_mesr();
    ## get the smry subset for the parameters unless smryx passed as paramter
    if (missing(smryx)) {
      smry=get_data(smry);
      smryx=smry_select(smry,n1,n2,d1,d2,x=x);
    }
    smry=smryx$smry; 
    ## if (missing(xlim)) xlim=range(x);
    data=data_rate(smry,mesr,rate,error.rate,from.type,relto.type);
    xdata=data$x; ydata=data$y; yes=data$yes;
    if (length(unique(yes))>1&&!yesno.multiok)
      stop("Cannot mix 'yes' and 'no' cases (eg, false positives and false negatives) on same heatmap unless yesno.multiok is TRUE");
    rate.type=rate_type();
    cutoff=rate_cutoff();
    if (!is.null(yes)) {
      ## deal with 'yes' vs. 'no' cases. want 'no' cases first - these are false positives
      yes.run=rle(yes)$lengths;
      if (length(yes.run)==1|(length(yes.run==2)&yes[1]==F))
        ## in correct order already. use x ordering from smryx
        x=smryx$x
      else {
        ## sort 'no' cases first if allowed (ie, x='auto')
        if (x=='auto') {
          i=order(yes);
          xdata=xdata[i,]; ydata=ydata[i,]; yes=yes[i];
          x='asis'
        }}
    } else {
      ## it's pos or neg - yes/no not relevant. use x ordering from smryx
      yes=rep(T,nrow(xdata));
      x=smryx$x;
    }
    ## if (is.null(xlab)) xlab=x;
    if (is.null(title)) {
      xdt=xdata_xtitle(xdata,xtitle);
      xdata=xdt$xdata; xtitle=xdt$xtitle;
      title=paste(collapse="\n",c(title_rate(),xtitle));
    }
    ## if (x=='asis') x=seq_len(nrow(data))
    ## else x=xdata[,x];
    ## if (xaxt=='auto') if (nrow(xdata)<=xaxt.max) xaxt='n' else xaxt='s';
    ## if (xaxt=='n') {
    if (y=='auto') {
      ## sort measures by distance from 0 or 1 to make nicer picture
      if (error.rate||rate=='neg') yzero=rbind(zero=0,t(ydata))
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
    dev.cur();
  }

## remove single-valued xvars from xdata and put in title
xdata_xtitle=function(xdata,xtitle) {
  if (xtitle=='none') xtitle=NULL
  else {
    ## find single valued xvars we want
    xvars=colnames(xdata);          # start with all of 'em
    if (xtitle=='n') {xvars=grep('^n',xvars)}
    else if (xtitle=='d') {xvars=grep('^d',xvars)}
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


## determine the type of rate we're plotting, eg, 'fpr' or 'fnr'
rate_type=function(rate=parent(rate),error.rate=parent(error.rate),yes=parent(yes)) {
  if (rate %notin% cq(pos,neg)) {
    if (error.rate&all(yes)) 'fnr'
    else if (error.rate&all(!yes)) 'fpr'
    else if (!error.rate&all(yes)) 'tpr'
    else if (!error.rate&all(!yes)) 'tnr'
    else if (error.rate) 'err'
    else 'corr';
  } else rate;
}
## determine cutoff for rate type
rate_cutoff=
  function(type=parent(rate.type),yes=parent(yes),
           fpr.cutoff=parent(fpr.cutoff,.05),fnr.cutoff=parent(fnr.cutoff,0.20),
           cutoff=parent(cutoff,.05)) {
    switch(type,
           neg=cutoff,pos=1-cutoff,err=cutoff,corr=1-cutoff,
           fpr=fpr.cutoff,fnr=fnr.cutoff,tpr=1-fnr.cutoff,tnr=1-fpr.cutoff);
  }
## construct colors and breaks for heatmap
heat_setup=function(cutoff=parent(cutoff),steps=100) {
  hi.brk=seq(0,-log10(cutoff),length.out=(steps)+1);
  lo.brk=seq(-log10(cutoff),4,length.out=(steps)+1);
  blues=colorRampPalette(RColorBrewer::brewer.pal(5,'Blues'))(steps);
  reds=colorRampPalette(RColorBrewer::brewer.pal(4,'Reds'))(steps);
  brk=unique(c(hi.brk,lo.brk));
  col=c(blues,reds);
  list(brk=brk,col=col);
  }
## generate 'rate' part of title for plots
title_rate=
  function(rate=parent(rate),type=parent(rate.type),relto.type=parent(relto.type)) {
    desc=switch(rate,
                pos='positive',
                neg='negative',
                nonzro='study1 non-zero',
                nonz1='study1 non-zero',
                nonz1or2='at least one study non-zero',
                nonz1and2='both studies non-zero',
                sameff='both studies have same effect',
                uni='both studies have same non-zero effect');
    desc=paste(sep='','(',desc,')');
    type=switch(type,
                err='error',
                corr='correct',
                fpr='false positive',
                fnr='false negative',
                tpr='true positive',
                tnr='true negative',
                none=NULL);
    if (length(unique(relto.type))==1&&relto.type!='sig1')
      relto=paste(sep=' ','-- relative to',unique(relto.type))
    else relto=NULL;
    paste(collapse=' ',c(rate,desc,type,'rate',relto));
}
## generate 'rate' part of ylab for plots
ylab.rate=function(rate=parent(rate),error.rate=parent(error.rate)) {
  err.or.cor=if (error.rate) 'error' else 'correct';
  paste(sep=' ',rate,err.or.cor,'rate');
}

## get line properties for measures in matrix
mesr2col=function(matrix) col.mesr[colnames(matrix)];
mesr2lwd=function(matrix) lwd.mesr[colnames(matrix)];
mesr2lty=function(matrix) lty.mesr[colnames(matrix)];
## construct legend for measures in plotrate
mesr_legend=
  function(ydata,where='bottomright',x=NULL,y=NULL,title='measure',
           bty='n',pch=19,cex=parent(cex.legend,0.8),title.col='black') {
  mesr=colnames(ydata);
  mesr=mesr[order(match(mesr,mesr.order))];
  ydata=ydata[,mesr,drop=FALSE];
  ## legend.col=col.mesr[legend.text];
  if (is.null(x)) x=where;
  col=mesr2col(ydata);
  lwd=mesr2lwd(ydata);
  lty=mesr2lty(ydata);
  legend(x,y,bty=bty,legend=mesr,cex=cex,
         title=title,title.col=title.col,
         col=col,lwd=lwd,lty=lty,seg.len=6);
  ## legend(x,y,bty=bty,legend=mesr,pch=pch,cex=cex,
  ##        title=title,title.col=title.col,
  ##        col=col,lwd=lwd,lty=lty,seg.len=8);
  ## col=col,lwd=lwd,lty=lty);
  ## lwd=lwd.mesr[legend.text])
}
## construct legend for rate colors in heatrate
heat_legend=
  function(coord,cutoff=parent(cutoff,0.05),cex=parent(cex.legend,0.75),steps=100) {
    x0=coord$x0; y0=coord$y0; width=coord$width; height=coord$height;
    heat=heat_setup();
    brk=heat$brk; col=heat$col;
 
    ## TODO: refactor color & breaks computation
    hi.brk=seq(0,-log10(cutoff),length.out=(steps)+1);
    lo.brk=seq(-log10(cutoff),4,length.out=(steps)+1);
    blues=colorRampPalette(RColorBrewer::brewer.pal(5,'Blues'))(steps);
    reds=colorRampPalette(RColorBrewer::brewer.pal(4,'Reds'))(steps);
    brk=unique(c(hi.brk,lo.brk));
    col=c(blues,reds);
    ## define legend grid
    ##   for some reason, doesn't work to let y range from y0 to y0+height
    ##   have to stop one unit smaller
    x=c(x0,x0+width);
    y=seq(y0,y0+height,length.out=2*steps+1)[1:(2*steps)];
    z=t(as.matrix(rev(c(head(hi.brk,-1),lo.brk[-1]))));
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
           yq=quantile(y,prob=c(0.1,0.9)),height=diff(yq),cex=parent(cex.legend,0.75)) {
    ## call plot to setup coordinate systme. needed by strwidth
    plot(x=c(min(x),1.2*max(x)),y=range(y),type='n',axes=F,xlab=NA,ylab=NA);
    ## add space to labels to get a smidge more room
    x0=max(x)+0.5+max(strwidth(c('0','1',cutoff),cex=cex))+strwidth(' ',cex=cex);
    y0=yq[1];
    x1=x0+width+strwidth(' ',cex=cex);
    list(x0=x0,y0=y0,width=width,height=height,x1=x1);
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

