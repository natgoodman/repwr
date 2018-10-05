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
## Data management -- files and cached objects -- for repwr.R: 
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################

## ---- Save and Load ----
## call with file or software attempts to construct file name
## CAUTION: when called in a fresh workspace, directory variables not yet defined
##### for data keyed by n, d
## save data in RData and optionally txt formats
save_nd=function(data,n,d,id=NULL,file=NULL,what,save,save.txt=F,keep) {
  if (is.null(file)) base=basename_nd(n,d,id,what)
  else base=desuffix(file);
  file=filename(base=base,suffix='RData');
  if ((is.na(save)&!file.exists(file))|(!is.na(save)&save)) {
    save(data,file=file);
    if (save.txt)
      write.table(data,file=filename(base=base,suffix='txt'),sep='\t',quote=F,row.names=F);
  }
  if (keep) keep_nd(data,n,d,id,what=what);
}
## load data from file
load_nd=function(file=NULL,n,d,id=NULL,what) {
  if (is.null(file)) file=filename_nd(n,d,id,what);
  what=load(file=file);               # what is name of saved data
  get(what);                          # return it
}
## get data already in memory (eg, in sim.list) or read from file
##   fail if data does not exist unless must.exist is FALSE
get_nd=function(n,d,id=NULL,what,load,keep,must.exist=T) {
  check_m();                            # are 'm' and 'datadir' consistent?
  case=casename_nd(n,d,id,short=T);
  if (is.na(load)|load) {
    what.list=get(paste(sep='.',what,'list'))
    data=what.list[[case]];
    if (is.null(data)) {
      file=filename_nd(n,d,id,what);
      if (file.exists(file)) data=load_nd(file=file)
      else {
        if (must.exist)
          stop(paste(sep=' ',what,'for case',case,
                     'not in memory and file',file,'does not exist'));
        data=NULL;
      }
    if (keep) keep_nd(data,n,d,id,what=what);
    }}
    else {
      if (must.exist)
        stop(paste(sep=' ','must.exist is TRUE but load is FALSE for',what,'case',cse));
      data=NULL;
    }
  invisible(data);
}
keep_nd=function(data,n,d,id=NULL,what) {
  case=casename_nd(n,d,id,short=T);
  what.list=get(paste(sep='.',what,'list'))
  what.list[[case]]=data;
  assign(paste(sep='.',what,'list'),what.list,envir=.GlobalEnv);
}
##### for data keyed by  n1, n2, d1, d2
save_nndd=function(data,n1,n2,d1,d2,id=NULL,file=NULL,what,save,save.txt=F,keep) {
  if (is.null(file)) base=basename_nndd(n1,n2,d1,d2,id,what)
  else base=desuffix(file);
  file=filename(base=base,suffix='RData');
  if ((is.na(save)&!file.exists(file))|(!is.na(save)&save)) {
    save(data,file=file);
    if (save.txt)
      write.table(data,file=filename(base=base,suffix='txt'),sep='\t',quote=F,row.names=F);
  }
  if (keep) keep_nndd(data,n1,n2,d1,d2,what=what);
}
## load data from file
load_nndd=function(file=NULL,n1,n2,d1,d2,id=NULL,what) {
  if (is.null(file)) file=filename_nndd(n1,n2,d1,d2,id,what);
  what=load(file=file);               # what is name of saved detl
  get(what);                          # return it
}
## get data already in memory (eg, in detl.list) or read from file
##   fail if data does not exist unless must.exist is FALSE
get_nndd=function(n1,n2,d1,d2,id=NULL,what,load,keep,must.exist=T) {
  check_m();                            # are 'm' and 'datadir' consistent?
  case=casename_nndd(n1,n2,d1,d2,id,short=T);
  if (is.na(load)|load) {
    what.list=get(paste(sep='.',what,'list'))
    data=what.list[[case]];
    if (is.null(data)) {
      file=filename_nndd(n1,n2,d1,d2,id,what);
      if (file.exists(file)) data=load_nndd(file=file)
      else {
        if (must.exist)
          stop(paste(sep=' ',what,'for case',case,
                     'not in memory and file',file,'does not exist'))
        else data=NULL;
      }
      if (keep) keep_nndd(data,n1,n2,d1,d2,what=what);
    }}
    else {
      if (must.exist)
        stop(paste(sep=' ','must.exist is TRUE but load is FALSE for',what,'case',case));
      data=NULL;
    }
  invisible(data);
}
keep_nndd=function(data,n1,n2,d1,d2,id=NULL,what) {
  case=casename_nndd(n1,n2,d1,d2,id,short=T);
  what.list=get(paste(sep='.',what,'list'))
  what.list[[case]]=data;
  assign(paste(sep='.',what,'list'),what.list,envir=.GlobalEnv);
}
##### sim
save_sim=function(sim,n,d,id=NULL,file=NULL,save=save.sim,save.txt=save.txt.sim,keep=keep.sim)
  save_nd(sim,n,d,id,file,what='sim',save,save.txt,keep);
load_sim=function(file=NULL,n,d,id=NULL) load_nd(file,n,d,id,what='sim');
get_sim=function(n,d,id=NULL,load=load.sim,keep=keep.sim,must.exist=T)
  get_nd(n,d,id,what='sim',load,keep,must.exist);
##### simr
save_simr=function(simr,n,d,id=NULL,file=NULL,save=save.simr,save.txt=save.txt.simr,keep=keep.simr)
  save_nd(simr,n,d,id,file,what='simr',save,save.txt,keep);
load_simr=function(file=NULL,n,d,id=NULL) load_nd(file,n,d,id,what='simr');
get_simr=function(n,d,id=NULL,load=load.simr,keep=keep.simr,must.exist=T)
  get_nd(n,d,id,what='simr',load,keep,must.exist);

##### si - s1 permutation indexes (mostly for testing)
save_si=function(si,n1,n2,d1,d2,id=NULL,file=NULL,save=save.si,save.txt=save.txt.si,keep=keep.si)
  save_nndd(si,n1,n2,d1,d2,id,file,what='si',save,save.txt,keep);
load_si=function(file=NULL,n1,n2,d1,d2,id=NULL)
  load_nndd(file,n1,n2,d1,d2,id,what='si');
get_si=function(n1,n2,d1,d2,id=NULL,load=load.si,keep=keep.si,must.exist=T)
  get_nndd(n1,n2,d1,d2,id,what='si',load,keep,must.exist);
##### detl - detailed measure results
save_detl=function(detl,n1,n2,d1,d2,id=NULL,file=NULL,
                   save=save.detl,save.txt=save.txt.detl,keep=keep.detl)
  save_nndd(detl,n1,n2,d1,d2,id,file,what='detl',save,save.txt,keep);
load_detl=function(file=NULL,n1,n2,d1,d2,id=NULL)
  load_nndd(file,n1,n2,d1,d2,id,what='detl');
get_detl=function(n1,n2,d1,d2,id=NULL,load=load.detl,keep=keep.detl,must.exist=T)
  get_nndd(n1,n2,d1,d2,id,what='detl',load,keep,must.exist);
##### smry - summary measure results
save_smry=function(smry,n1,n2,d1,d2,id=NULL,file=NULL,
                   save=save.smry,save.txt=save.txt.smry,keep=keep.smry)
  save_nndd(smry,n1,n2,d1,d2,id,file,what='smry',save,save.txt,keep);
load_smry=function(file=NULL,n1,n2,d1,d2,id=NULL)
  load_nndd(file,n1,n2,d1,d2,id,what='smry');
get_smry=function(n1,n2,d1,d2,id=NULL,load=load.smry,keep=keep.smry,must.exist=T)
  get_nndd(n1,n2,d1,d2,id,what='smry',load,keep,must.exist);
##### posr - positive rate results
save_posr=function(posr,id='std',file=NULL,
                   save=save.posr,save.txt=save.txt.posr,keep=keep.posr) {
  if (missing(file)) base=basename_posr(id)
  else base=desuffix(file);
  file=filename(base=base,suffix='RData');
  if ((is.na(save)&!file.exists(file))|(!is.na(save)&save)) {
    save(posr,file=file);
    if (save.txt)
      write.table(posr,file=filename(base=base,suffix='txt'),sep='\t',quote=F,row.names=F);
  }
  if (keep) keep_posr(posr,id);
  invisible(data);
} 
load_posr=
  function(file=NULL,id='std') {
    if (is.null(file)) file=filename_posr(id);
    what=load(file=file);               # what is name of saved data
    get(what);                          # return it
  }
get_posr=
  function(id='std',load=load.posr,keep=keep.posr,must.exist=T) {
    check_m();                            # check 'm' and 'datadir' consistent
    what='posr';
    case=casename_posr(id);
    if (is.na(load)|load) {
      data=posr.list[[case]];
      if (is.null(data)) {
        file=filename_posr(id);
        if (file.exists(file)) data=load_posr(file=file)
        else {
          if (must.exist)
            stop(paste(sep=' ',what,'for case',case,
                       'not in memory and file',file,'does not exist'));
          data=NULL;
        }
        if (keep) keep_posr(data,id);
      }}
    else {
      if (must.exist)
        stop(paste(sep=' ','must.exist is TRUE but load is FALSE for',what,'case',case));
      data=NULL;
    }
    invisible(data);
  }
keep_posr=function(data,id='std') {
  case=casename_posr(id);
  what=paste(sep='.','posr',case);
  assign(what,data,envir=.GlobalEnv); # assign globally
  posr.list[[case]]<<-data;             # assign to global list
}

##### data - top-level data saved in datadir
## save data in RData and txt formats
save_data=function(what,file=NULL,data=NULL,id=NULL,
                   save=save.data,save.txt=save.txt.data,keep=keep.data) {
  what=as.character(pryr::subs(what));
  if (missing(data) && exists(what,envir=parent.frame(n=1)))
    data=get(what,envir=parent.frame(n=1));
  if (is.null(data)) stop('Trying to save NULL object. Is "what" set correctly?');
  if (is.null(file)) base=basename_data(what,id)
  else base=desuffix(file);
  file=filename(base=base,suffix='RData');
  if ((is.na(save)&!file.exists(file))|(!is.na(save)&save)) {
    save(data,file=file);
    if (save.txt) {
      file=filename(base=base,suffix='txt');
      if (length(dim(data))==2) write.table(data,file=file,sep='\t',quote=F,row.names=F)
      else if (is.vector(data)) writeLines(data,file)
      else stop('Trying to save object with more than 2 dimensions as text. Is "what" set correctly?');
    }
  }
  if (keep) keep_data(name=what,data=data);
  invisible(data);
}
## load data from file
load_data=function(file=NULL,what=NULL,id=NULL) {
  if (is.null(file)&is.null(what))
    stop('Cannot load data unless file or what is set');
  if (is.null(file)) file=filename_data(what,id);
  what=load(file=file);               # what is name of saved data
  get(what);                          # return it
}
## get top-level data already in memory or read from file
##   fail if data does not exist unless must.exist is FALSE
get_data=function(what,id=NULL,load=load.data,keep=keep.data,must.exist=T,name=NULL) {
  check_m();                            # check 'm' and 'datadir' consistent
  if (is.na(load)|load) {
    if (!is.null(name)) what=name else what=as.character(pryr::subs(what));
    if (!is.null(id)) what=paste(sep='.',what,id);
    data=data.list[[what]];
    if (is.null(data)) {
      file=filename_data(what);
      if (file.exists(file)) {
        data=load_data(file=file);          # what is name of saved data
      } else {
        if (must.exist)
          stop(paste(sep=' ',what,'not in memory and file',file,'does not exist'));
        data=NULL;
      }
      if (keep&&!is.null(data)) keep_data(name=what,data=data);
    }}
    else {
      if (must.exist)
        stop(paste(sep=' ','must.exist is TRUE but load is FALSE for',what));
      data=NULL;
    }
  invisible(data);
}
## keep top-level data. assign globally and keep in global list
keep_data=function(what,id=NULL,data=NULL,name=NULL) {
 if (!is.null(name)) what=name else what=as.character(pryr::subs(what));
 if (missing(data) && exists(what,envir=parent.frame(n=1)))
   data=get(what,envir=parent.frame(n=1));
 if (is.null(data)) stop("Trying to keep NULL object. Is 'what' set correctly?");
 if (!is.null(id)) what=paste_id(what,id);
 assign(what,data,envir=.GlobalEnv); # assign globally
 data.list[[what]]<<-data;  # assign to global list
}

##### plot - save one or more plots - NOT USED - saving done in dofig
save_plot=function(dev,file=NULL,what=NULL,id=NULL,i=NULL) {
  if (is.null(file)&is.null(what)) stop('Cannot save plot unless file or what is set');
  if (is.null(file)) {
    file=
      if (length(dev)==1) filename_plot(what,id)
        else sapply(seq_along(dev), function(i) filename_plot(what,id=id,i=i));
  } else file=filename(desuffix(file),suffix='png');
  for (i in seq_along(dev)) savePlot(file[i],device=dev[i]);
}

##### table - saved in tbldir
save_tbl=function(what,file=NULL,data=NULL,
                  sect=parent(figsect,NULL),tblnum=parent(tblnum,NULL),id=parent(id,NULL),i=NULL,
                  save=save.tbl,save.txt=save.txt.tbl) {
  what=as.character(pryr::subs(what));
  if (missing(data) && exists(what,envir=parent.frame(n=1)))
    data=get(what,envir=parent.frame(n=1));
  if (is.null(data)) stop('Trying to save NULL object. Is "what" set correctly?');
  if (is.null(file)) base=basename_tbl(what,sect=sect,tblnum=tblnum,id=id,i=i)
  else base=desuffix(file);
  file=filename(base=base,suffix='RData');
  if ((is.na(save)&!file.exists(file))|(!is.na(save)&save)) {
    save(data,file=file);
    if (save.txt) {
      file=filename(base=base,suffix='txt');
      if (length(dim(data))==2) write.table(data,file=file,sep='\t',quote=F,row.names=F)
      else if (is.vector(data)) writeLines(as.character(data),file)
      else stop('Trying to save object with more than 2 dimensions as text. Is "what" set correctly?');
    }}
  invisible(data);
}

## ---- File Functions ----
##### names for data keyed by n, d
filename_nd=function(n,d,id=NULL,what,suffix='RData') 
  filename(basename_nd(n,d,id,what),suffix=suffix);
basename_nd=function(n,d,id=NULL,what)
  filename(get(paste(sep='',what,'dir')),base=what,tail=casename_nd(n,d,id))
casename_nd=function(n,d,id=NULL,short=F) {
  casename=
    if (short) paste(sep=',',n,d_pretty(d))
      else paste(sep=',',paste_nv(n),paste_nv(d,d_pretty(d)));
  paste_id(casename,id);
}
##### names for data keyed by n1, n2, d1, d2
filename_nndd=function(n1,n2,d1,d2,id=NULL,what,suffix='RData') 
  filename(basename_nndd(n1,n2,d1,d2,id,what),suffix=suffix);
basename_nndd=function(n1,n2,d1,d2,id=NULL,what)
  filename(get(paste(sep='',what,'dir')),base=what,tail=casename_nndd(n1,n2,d1,d2,id))
casename_nndd=function(n1,n2,d1,d2,id=NULL,short=F) {
  if (d1==d2) {
      casename=
        if (short) paste(sep=',',n1,n2,d_pretty(d1))
          else paste(sep=',',paste_nv(n1),paste_nv(n2),paste_nv('d1=d2',d_pretty(d1)));
    } else {
      casename=
        if (short) paste(sep=',',n1,n2,d_pretty(d1),d_pretty(d2))
        else paste(sep=',',paste_nv(n1),paste_nv(n2),
                   paste_nv(d1,d_pretty(d1)),paste_nv(d2,d_pretty(d2)));
    }
  paste_id(casename,id);
}
##### sim
filename_sim=function(n,d,id=NULL,suffix='RData') filename_nd(n,d,id,what='sim',suffix);
basename_sim=function(n,d,id=NULL) basename_nd(n,d,id,what='sim');
casename_sim=casename_nd;
##### simr
filename_simr=function(n,d,id=NULL,suffix='RData') filename_nd(n,d,id,what='simr',suffix);
basename_simr=function(n,d,id=NULL) basename_nd(n,d,id,what='simr');
casename_simr=casename_nd;
##### i1 (mostly for testing)
filename_i1=function(n1,n2,d1,d2,id=NULL,suffix='RData')
  filename_nndd(n1,n2,d1,d2,id,what='i1',suffix);
basename_i1=function(n1,n2,d1,d2,id=NULL) basename_nndd(n1,n2,d1,d2,id,what='i1');
casename_i1=casename_nndd;
##### detl
filename_detl=function(n1,n2,d1,d2,id=NULL,suffix='RData')
  filename_nndd(n1,n2,d1,d2,id,what='detl',suffix);
basename_detl=function(n1,n2,d1,d2,id=NULL) basename_nndd(n1,n2,d1,d2,id,what='detl');
casename_detl=casename_nndd;
##### smry
filename_smry=function(n1,n2,d1,d2,id=NULL,suffix='RData')
  filename(basename_smry(n1,n2,d1,d2,id),suffix=suffix);
basename_smry=function(n1,n2,d1,d2,id=NULL) basename_nndd(n1,n2,d1,d2,id,what='smry');
casename_smry=casename_nndd;
##### posr
filename_posr=function(id='std',suffix='RData')
  filename(base=basename_posr(id),suffix=suffix);
basename_posr=function(id='std')
  filename(posrdir,base='posr',tail=casename_posr(id));
casename_posr=function(id='std') id;
##### data - arbitrary objects saved in datadir
filename_data=function(what,id=NULL,suffix='RData')
  filename(basename_data(what,id),suffix=suffix);
basename_data=function(what,id=NULL) filename(datadir,base=paste_id(what,id));

##### figure - saved in figdir. may have numeric tail
filename_fig=function(figname,figpfx=NULL,fignum=NULL,id=NULL,i=NULL,suffix='png')
  filename(basename_fig(figname,figpfx,fignum,id,i),suffix=suffix);
basename_fig=function(figname,figpfx=NULL,fignum=NULL,id=NULL,i=NULL) {
  if (!is.null(i)) i=sprintf("%02i",i);
  if (!is.null(fignum)) fignum=c('figure',paste(collapse='',c(figpfx,sprintf("%03i",fignum))));
  base=paste(collapse='_',c(fignum,figname));
  basename=filename(figdir,base=base,tail=i);
  paste_id(basename,id);
}
##### table - saved in tbldir. may have numeric tail
filename_tbl=function(tblname,tblpfx=NULL,tblnum=NULL,id=NULL,i=NULL,suffix='png')
  filename(basename_tbl(tblname,tblpfx,tblnum,id,i),suffix=suffix);
basename_tbl=function(tblname,tblpfx=NULL,tblnum=NULL,id=NULL,i=NULL) {
  if (!is.null(i)) i=sprintf("%02i",i);
  if (!is.null(tblnum)) tblnum=c('table',paste(collapse='',c(tblpfx,sprintf("%03i",tblnum))));
  base=paste(collapse='_',c(tblnum,tblname));
  basename=filename(tbldir,base=base,tail=i);
  paste_id(basename,id);
}

## construct file or directory pathname from components
## wrapper for file.path with base, tail and suffix pasted on
##  base appended with '.'
##  tail components combined with '.'
##  suffix added unless already there
filename=function(...,base=NULL,tail=NULL,suffix=NULL) {
  if (!is.null(base)||!is.null(tail)) base=paste(collapse='.',c(base,tail));
  if (is.null(base)) file=file.path(...) else file=file.path(...,base);
  if (!is.null(suffix)) {
    ## remove leading '.' if present
    suffix=sub('^\\.','',suffix,perl=T);
    suffix.pattern=paste(collapse='|',paste(sep='','\\.',suffix,'$'));
    file=ifelse(grepl(suffix.pattern,file),file,paste(sep='.',file,suffix[1]));
  }
  file;
}
## remove suffix from filename
desuffix=function(file,suffix=c('RData','txt')) {
  if (!is.null(suffix)) {
    ## remove leading '.' if present
    suffix=sub('^\\.','',suffix,perl=T);
    suffix.pattern=paste(collapse='|',paste(sep='','\\.',suffix,'$'));
    file=sub(suffix.pattern,'',file);
  }
  file;
}
## filebasename same as filename but w/o suffix
filebasename=function(...) filename(...,suffix=NULL)
## construct directory pathname. synonym for filebasename
dirname=filebasename;

## clean specific data type. deletes directory, any top level files and in-memory list
cleanq=function(what,cleandir=T) {
  what=as.character(pryr::subs(what));
  ## delete top level files if exist, including ones with ids
  ## unlink(sapply(cq(RData,txt), function(suffix) filename_data(what,suffix=suffix)));
  unlink(filename(datadir,list.files(datadir,pattern=paste(sep='','^',what,'\\.'))));
  ## delete in-memory list
  whatlist=paste(sep='.',what,'list');
  if (exists(whatlist,envir=.GlobalEnv)) rm(list=whatlist,envir=.GlobalEnv);
  ## delete from top level data.list
  ## CAUTION: have to use loop (not sapply) for scoping to work
  if (exists('data.list',envir=.GlobalEnv)) {
    pat=paste(sep='','^',what,'(\\.|$)');
    for (name in grep(pat,names(data.list),value=T)) data.list[[name]]<<-NULL;
  }
  if (cleandir) {
    whatdir=paste(sep='',what,'dir');
    ## delete directory if exists
    if (exists(whatdir,envir=.GlobalEnv)) unlink(get(whatdir,envir=.GlobalEnv),recursive=T);
  }
}
