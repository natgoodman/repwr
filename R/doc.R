#################################################################################
##
## Author:  Nat Goodman
## Created: 18-06-19
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
## --- Generate Figures and Tables for Document ---
## sect is which sections to run - for use during development
##   uses prefix matching and all matches run
## other params passed to init or init_doc. typically
##   subdoc - NULL or 'supp'. default: NULL
##   fignum, tblnum - initial figure, table number. default: 1
##   figpfx, tblpfx - prefix prepended to figure, table number, eg, 'S' for supp
##   figsfx, tblsfx - suffices appended to figure, table number in 'blocks' eg, a,b,c,...
##   save.out - save figures, tables to files. default: T
##   figscreen - plot figures on screen. default: T for readme, F otherwsie
##   fignew - plot each figure in new window. default figscreen
##   docfun document-specific function. default calculated from doc,subdoc eg, doc_repwr
dodoc=
  function(sect=NULL,need.init=T,doc=parent(doc,'readme'),...) {
    ## split ... args for init and init_doc. from stackoverflow.com/questions/4124900
    dots=list(...);
    iargs=c(list(doc=doc),dots[names(dots) %in% names(formals(init))]);
    dargs=dots[names(dots) %in% names(formals(init_doc))];
    if (need.init) {
      ## for sandbox runs, use doc-specific init
      if (is.na(pmatch(doc,'xperiment'))) do.call('init',iargs)
      else do.call('init_xperiment',iargs);
    }
    do.call('init_doc',dargs);
    docfun(sect=sect);
  }

## --- Document Functions ---
## utility functions to make figures and tables for documents
## run plot function, save if required, label result with function name
## CAUTION: ... interacts with partial argument matching to cause dofig args to be
##   matched by plot-function args. eg, 'doc' matched by 'd'. prepending with 'fig'
##   works only because no plot-function arg matches it
dofig=
  function(figfun,name=NULL,sect=parent(sect,NULL),sectnum=parent(sectnum,NULL),...) {
    fig=figname(name,sect,sectnum);
    file=filename_fig(name,sect,sectnum);
    BREAKPOINT();
    plot.to.file=((is.na(save.fig)&!file.exists(file))|(!is.na(save.fig)&save.fig));
    plot.to.screen=figscreen;           # for stylistic consistency
    ## NG 18-08-10: new scheme for plotting to file
    ##   plot to screen and file: dev.new here, dev.copy later
    ##   plot to screen only: dev.new here, no dev.copy later
    ##   plot to file only: png here, dev.off later
    ## plot.to.file only doesn't work if figfun returns multiple figures
    ##   trash multi-figure capability. we don't use it now
    if (!(plot.to.file||plot.to.screen)) {
      msg=paste(sep=' ',paste_nv(figscreen),'and',paste_nv(save.fig));
      if (is.na(save.fig)) msg=paste(sep=' ',msg,'and figure file',file,'exists');
      msg=paste(sep=' ',msg,'which means there is no where to plot the figure');
      stop(msg);
    }
    ##   bg='white' needed else image copied with transparent bg; renders as grey
    if (plot.to.screen) {
      dev=dev.new(bg='white');
      dev=dev.cur();
    }
    if (plot.to.file&!plot.to.screen) {
      ## png parameters found by trial and error. look reasonable
      ## TODO: learn the right way to do this!
      png(filename=file,height=8,width=8,units='in',res=200,pointsize=12);
      dev.png=dev.cur();
      }
    ## draw the figure!
    figfun(...,fignum=fig);
    if (plot.to.file&plot.to.screen) 
      ## png parameters found by trial and error. look reasonable
      ## TODO: learn the right way to do this!
      dev.png=dev.copy(png,filename=file,height=8,width=8,units='in',res=200,pointsize=12);
    ## always close plot.to.file device
    if (plot.to.file)  dev.off(dev.png);
    ## close plot.to.screen device unless user wants each figure in new window
    if (plot.to.screen&&!fignew) dev.off(dev);
    ## assign_parent(fignum,fignum+1);
    fignum<<-fignum+1;
    figname;
  }
## save one or more tables.
dotbl=
  function(...,sect=parent(sect,NULL),sectnum=parent(sectnum,NULL)) {
    tbl=list(...);                           # evaluates dots
    dots=match.call(expand.dots=FALSE)$...;  # doesn't evaluate dots
    tblname=sapply(seq_along(tbl),function(i) {
      name=names(tbl)[i];
      ## test for empty name. CAUTION: do it carefully lest R complains when name is empty list
      empty.name=if(length(name)==0) T else if(nchar(name)==0) T else F;
      if (empty.name) {
        name=as.character(dots[[i]]);
        data=get(name,envir=parent.frame(n=4)); # n=4 empirically determined
      } else data=tbl[[i]];
      file=filename_tbl(name,sect,sectnum);
      save_tbl(name,data,file=file);
      ## write.table(tbl[[name]],file=file,sep='\t',quote=F,row.names=F);
      tblnum<<-tblnum+1;
      tblname;})
    ## assign_parent(tblnum,tblnum+length(tbl));
    tblname;
 }
