##juli_2_developed by Hyun-Tae Shin 
##version_20170518
##update: 20160907 inluding counter fusion in SGI output
##update: 20160929 changing known annotation for SGI output, 20170518 debugged line685,686 
##Command: Rscript script input_bam output_path sample_name thread
##Requirement: samtools, index file of bam, gap_hg19.txt, refGene.txt, CosmicFusionExport_V76.tsv, Pfam-A.full.human,HGNC_GeneName_UniProtID_160524.txt

#################
######path#######
#################
samtools.path='/usr/bin/samtools'
# samtools.path='/data/ngs_tools/samtools-0.1.19/samtools' # TODO
# gapdata='/data/ngs_ref/cs_v2.2/gap_hg19.txt' # TODO
# refgene='/data/ngs_ref/cs_v2.2/refGene.txt' # TODO
# cosmic='/data/ngs_ref/cs_v2.2/CosmicFusionExport_V76.tsv' # TODO
# pfampath='/data/ngs_ref/cs_v2.2/Pfam-A.full.human' # TODO
# gpmappath='/data/ngs_ref/cs_v2.2/HGNC_GeneName_UniProtID_160524.txt' # TODO

#################
#####args########
#################

args <- commandArgs(TRUE)
input.bam=args[1]
output=args[2]
Sample_ID=args[3]
thread=args[4]

gapdata=args[5]
refgene=args[6]
cosmic=args[7]
pfampath=args[8]
gpmappath=args[9]

panelgene=args[10]

if(!is.na(thread)){thread=as.numeric(thread)}else{thread=1} #default 1

output.bam=paste0(output,"/",Sample_ID)

#################
#####parameter###
#################

repeat.ran=3.62 ##mean repeat no from 10000 random sequence read
MR=0.8 ##missmatch ratio
t1=c('ALK','BCL2','BCR','BRAF','EGFR','FGFR1','FGFR2','FGFR3','NTRK1','PDGFRA','RARA','RET','ROS1')

##########################
#####library & option#####
##########################

library(data.table)
library(gtools)
library(foreach)
library(doMC)
registerDoMC(thread)
library(compiler)
enableJIT(3)
options(scipen = 999)
library(ggplot2)

####################
####data loading####
####################

gap=read.table('/mnt/c/Users/dlgkr/op/juli_ref/gap_hg19.txt', sep='\t',stringsAsFactors = F)
ref=read.table('/mnt/c/Users/dlgkr/op/juli_ref/refGene.txt',sep='\t',stringsAsFactors = F)
cos=read.table('/mnt/c/Users/dlgkr/op/juli_ref/CosmicFusionExport_V76.tsv',sep='\t',stringsAsFactors = F,header=T)
pfam=read.table('/mnt/c/Users/dlgkr/op/juli_ref/Pfam-A.full.human',sep='\t',stringsAsFactors = F,fill=T, comment.char = "%")
gpmap=read.table('/mnt/c/Users/dlgkr/op/juli_ref/HGNC_GeneName_UniProtID_160524.txt',sep='\t',stringsAsFactors = F,fill=T,header=T,quote = "")
input.bam = "/mnt/c/Users/dlgkr/op/git_upload/cancerscan/sample/CD_22_13185_BD_D_SCN_1.recal.sorted.bam"
Sample_ID = 'CD_22_13185_BD_D_SCN_1'
output.bam=paste0("/mnt/c/Users/dlgkr/op/git_upload/cancerscan/sample","/",Sample_ID)

#################
###function######
#################

leading = function(num, full){ ## function for png output numbering like 0001, 0002, ---
    string = unlist(strsplit(as.character(num), ''))
    if (length(string)<full){
        count = full-length(string)
        return(paste0(c(rep('0', count), string), collapse=''))
    }
    
}

endfun=function(cigar){
  out=sum(as.numeric(unlist(strsplit(gsub('[[:digit:]]+[SHIN]','',cigar),'[MD]'))))
  return(out)
}

gapfun=function(chr,pos){
  out=sum(gap$V1==chr & gap$V2 <= pos & pos<=gap$V3)
  return(out) 
}

ranfun=function(chr,ou){
  tab.chr=table(chr)
  bed.chr=unique(chr)
  out=c()
  for(i in 1:length(bed.chr)){
    no=ceiling(sum(chr==bed.chr[i])/ou)
    for(j in 1:no){
      if(j!=no){ran=paste(bed.chr[i],((j-1)*ou+1),(j*ou),sep='/')}
      if(j==no){ran=paste(bed.chr[i],((j-1)*ou+1),tab.chr[names(tab.chr)==bed.chr[i]],sep='/')}
      out=c(out,ran)
    }
  }
  return(out)
}

pbedfun=function(bed.ran){
  pbed.ran=unlist(strsplit(bed.ran,'/')) 
  pbed=bed[bed$V1==pbed.ran[1],]
  pbed=pbed[order(pbed$V2),]
  pbed=pbed[c(as.numeric(pbed.ran[2]):as.numeric(pbed.ran[3])),]
  s=e=pbed$V2
  s[pbed$V3==0]=s[pbed$V3==0]-readlen
  e[pbed$V3==0]=e[pbed$V3==0]+2*medis
  s[pbed$V3==1]=s[pbed$V3==1]-2*medis
  e[pbed$V3==1]=e[pbed$V3==1]+readlen
  s[s<0]=e[e<0]=1
  out=data.table(pbed.ran[1],s,e,pbed[,V1:=NULL])
  colnames(out)=paste0('V',seq(1:ncol(out)))
  return(out)
}

tagfun=function(satag){
  a=unlist(strsplit(satag,'\\t'))
  out=gsub('SA:Z:','',a[grepl('SA:Z',a)])
  return(out)
}

scbasefun.sm=function(cigar,base){
  if(grepl('S',cigar)){
    slen=as.numeric(gsub('S[[:alnum:]]+$','',cigar))
    pre <- paste0('[A-Z]{',nchar(base)-slen,'}$')
    a=gsub(pre,'',base)
    sbase=paste(rev(substring(a,1:nchar(a),1:nchar(a))),collapse="")
  }
  if(grepl('H',cigar)){sbase=''}
  return(sbase)
}

scbasefun.ms=function(cigar,base){
  if(grepl('S',cigar)){
    slen=as.numeric(gsub('S$','',gsub('^[[:alnum:]]+[INMD]','',cigar)))
    pre <- paste0('^[A-Z]{',nchar(base)-slen,'}')
    sbase=gsub(pre,'',base)
  }
  if(grepl('H',cigar)){sbase=''}
  return(sbase)
}

mbasefun.sm=function(cigar,base){
  if(!grepl('H',cigar)){
    if( grepl('[S].+M.+[S]',cigar)){
      ss=as.numeric(gsub('S[[:alnum:]]+$','',cigar))
      se=as.numeric(gsub('S$','',gsub('^[[:alnum:]]+[INMD]','',cigar)))
    }else{
      ss=as.numeric(gsub('S[[:alnum:]]+$','',cigar))
      se=0
    }
    pre <- paste0('^[A-Z]{',ss,'}')
    post <- paste0('[A-Z]{',se,'}$')
    sbase=gsub(post,'',gsub(pre,'',base))
  }
  if(grepl('H',cigar)){sbase=base}
  return(sbase)
}

mbasefun.ms=function(cigar,base){
  if(!grepl('H',cigar)){
    if( grepl('[S].+M.+[S]',cigar)){
      ss=as.numeric(gsub('S[[:alnum:]]+$','',cigar))
      se=as.numeric(gsub('S$','',gsub('^[[:alnum:]]+[INMD]','',cigar)))
    }else{
      ss=0
      se=as.numeric(gsub('S$','',gsub('^[[:alnum:]]+[INMD]','',cigar)))
    }
    pre <- paste0('^[A-Z]{',ss,'}')
    post <- paste0('[A-Z]{',se,'}$')
    a=gsub(post,'',gsub(pre,'',base))
  }  
  if(grepl('H',cigar)){
    a=base
  }
  sbase=paste(rev(substring(a,1:nchar(a),1:nchar(a))),collapse="")
  return(sbase)
}

repbasefun=function(scbase){
  no=seq(1:max(nchar(scbase)))            
  tmpfun=function(no){
    a=substr(scbase,no,no)
    a=a[a!=""]
    b=names(table(a))[table(a)==max(table(a))][1]
    return(b)
  }
  c=paste(sapply(no,tmpfun),collapse='')
  return(c)
}

slidefun=function(s1,s2){
  mat='n'
  if(!grepl(2,s1)){ # for hard clip
    len=nchar(s1)+nchar(s2)
    ss2=c(rep(NA,nchar(s1)),unlist(strsplit(s2,'')),rep(NA,nchar(s1)))
    com=rep(NA,len)
    for( n in round(nchar(s1)/2):round((nchar(s1))+(nchar(s2)/2))){ 
      ss1=c(rep(NA,n-1),unlist(strsplit(s1,'')),rep(NA,len+1-n))
      com[n]=sum(ss1==ss2,na.rm=T)
    }
    if( max(com,na.rm=T) > min(nchar(s1),nchar(s2))*0.8 ){mat='y'} ## 80% of short base agreement
  }else{mat='h'}
  return(mat)
}

repeat.fun=function(br.str){  ## 6 base permutation
  read=br.str
  max.nor=c()
  for( j in 1:6 ){
    per=permutations(4,j,v=c('A','T','G','C'),repeats.allowed=TRUE)
    nor.ev=c()
    for( k in 1:nrow(per) ){
      max.rep=max(nchar(unlist(strsplit(gsub('[ATGC]+',' ',gsub(paste(per[k,],collapse=''),'M',read)),' '))))
      nor.ev=c(nor.ev,(max.rep-1)*j)
    }
    max.nor=c(max.nor,max(nor.ev))
  }
  return(max(max.nor)/repeat.ran)
}

reffun=function(gn,pos){
  if( !grepl('Flanking',gn) & !grepl('Compl',gn) ){
    nm=as.numeric(gsub('^.+_','',ref$V2[ref$V13==gn & ref$V5 <= pos & pos <= ref$V6]))
    gnnm=ref$V2[ref$V13==gn & ref$V5 <= pos & pos <= ref$V6][nm==min(nm)]    
    gnref=ref[ref$V2==gnnm & ref$V5 <= pos & pos <= ref$V6, ][1,]    
    es=as.numeric(unlist(strsplit(gnref$V10,',')))
    ee=as.numeric(unlist(strsplit(gnref$V11,',')))
    is=ee[-length(ee)]-1
    ie=es[-1]+1
    cod=as.numeric(unlist(strsplit(gnref$V16,',')))
    if( sum(es <= pos & pos <= ee)==1 ){
      if( gnref$V4 == '+' ){
        gnfun=paste0('_Exon(',which(es <= pos & pos <= ee),'/',gnref$V9,')')
        cod=paste0('_Frame(',cod[which(es <= pos & pos <= ee)],',',
                   cod[which(es <= pos & pos <= ee)+1],')')
      }
      if( gnref$V4 == '-' ){
        gnfun=paste0('_Exon(',gnref$V9-which(es <= pos & pos <= ee)+1,'/',gnref$V9,')')
        cod=paste0('_Frame(',rev(cod)[gnref$V9-which(es <= pos & pos <= ee)+1],',',
                   rev(cod)[gnref$V9-which(es <= pos & pos <= ee)],')')
      }
    }
    if( sum(is <= pos & pos <= ie)==1 ){
      if( gnref$V4 == '+' ){
        gnfun=paste0('_Intron(',which(is <= pos & pos <= ie),'/',gnref$V9-1,')')
        cod=paste0('_Frame(',cod[which(is <= pos & pos <= ie)+1],')')
      }
      if( gnref$V4 == '-' ){
        gnfun=paste0('_Intron(',gnref$V9-which(is <= pos & pos <= ie),'/',gnref$V9-1,')')
        cod=paste0('_Frame(',rev(cod)[gnref$V9-which(is <= pos & pos <= ie)+1],')')
      }
    }
    if( gnref$V5 <= pos & pos < gnref$V7 & grepl('Exon',gnfun)){
      if( gnref$V4 == '+' ){gnfun=gsub('_Exon','_5pUTR',gnfun)}
      if( gnref$V4 == '-' ){gnfun=gsub('_Exon','_3pUTR',gnfun)}
      cod=''
    }
    if( gnref$V8 < pos & pos <= gnref$V6 & grepl('Exon',gnfun)){
      if( gnref$V4 == '+' ){gnfun=gsub('_Exon','_3pUTR',gnfun)}
      if( gnref$V4 == '-' ){gnfun=gsub('_Exon','_5pUTR',gnfun)}
      cod=''
    }
    info=paste0(gnnm,gnfun,cod)
  }else{
    info=NA
  }
  return(info)
}

mapfun=function(gn,info,bp){
  pref=ref[ref$V2==gsub('_[[:alnum:]]+[(p].+$','',info) & ref$V13==gn,][1,]
  oes=as.numeric(unlist(strsplit(pref$V10,",")))
  oee=as.numeric(unlist(strsplit(pref$V11,",")))
  st=en=Domain=rep(NA,length(oes))
  if(pref$V4=='+'){
    exon=oee-oes
    nexon=exon/sum(exon)*0.2
    intron=oes[-1]-oee[-length(oee)]
    nintron=intron/sum(intron)*0.2
    for(n in 1:length(exon)){
      if(n==1){
        st[n]=0.05
        en[n]=st[n]+nexon[n]
      }else{
        st[n]=en[n-1]+nintron[n-1]
        en[n]=st[n]+nexon[n]
      }
    }
    if(grepl('Intron',info)){
      no=as.numeric(gsub('^.+[(]','',gsub('/.+$','',info)))
      nbp=en[no]+nintron[no]*((bp-(oee)[no])/intron[no])
    }
    if(!grepl('Intron',info)){
      no=as.numeric(gsub('^.+[(]','',gsub('/.+$','',info)))
      nbp=st[no]+nexon[no]*((bp-(oes)[no])/exon[no])
    }
  }
  if(pref$V4=='-'){
    exon=rev(oee-oes)
    nexon=exon/sum(exon)*0.2
    intron=rev(oes[-1]-oee[-length(oee)])
    nintron=intron/sum(intron)*0.2
    for(n in 1:length(exon)){
      if(n==1){
        st[n]=0.05
        en[n]=st[n]+nexon[n]
      }else{
        st[n]=en[n-1]+nintron[n-1]
        en[n]=st[n]+nexon[n]
      }
    }
    if(grepl('Intron',info)){
      no=as.numeric(gsub('^.+[(]','',gsub('/.+$','',info)))
      nbp=en[no]+nintron[no]*((rev(oee)[no]-bp)/intron[no])
    }
    if(!grepl('Intron',info)){
      no=as.numeric(gsub('^.+[(]','',gsub('/.+$','',info)))
      nbp=st[no]+nexon[no]*((rev(oes)[no]-bp)/exon[no])
    }
  }
  if(sum(gpmap[,2]==gn)!=0){
    ppfam=pfam[grepl(gpmap[gpmap[,2]==gn,4][1],pfam$V1),]
    if( nrow(ppfam)!=0){
      index=((oes <= pref$V7 & pref$V7 < oee) | (pref$V7 < oes & oee < pref$V8) | (oes < pref$V8 & pref$V8 <=oee))
      es=oes[index]
      es[1]=pref$V7
      ee=oee[index]
      ee[sum(index)]=pref$V8
      Domain=rep(NA,length(es))
      if(pref$V4=='+'){exon=ee-es}else{exon=rev(ee-es)}
      cumexon=cumsum(exon)
      for(i in 1:nrow(ppfam)){
        dstr=unlist(strsplit(gsub('[[:space:]]+.+$','',gsub('^.+/','',ppfam$V1[i])),'-',))
        dst=as.numeric(dstr[1])*3-2
        den=as.numeric(dstr[2])*3
        Domain[(dst <= c(1,cumexon[-length(cumexon)]+1) & cumexon <= den )|(c(1,cumexon[-length(cumexon)]+1)<=dst & dst<=cumexon)|(c(1,cumexon[-length(cumexon)]+1)<=den & den<=cumexon)]=gsub('^.+[[:space:]]+','',ppfam$V2[i])
      }
      utr=which(index==FALSE)
      utr5=utr[utr<(which(index==TRUE)[1])]
      utr3=utr[which(index==TRUE)[sum(index)]<utr]
      Domain=c(rep(NA,length(utr5)),Domain,rep(NA,length(utr3)))
    }
  }
  out=data.frame(st,en,Domain,nbp,stringsAsFactors=F)
  return(out)
}

#######################################################
### Measuring of read length and median insert size ###
#######################################################

print(paste("Measuring of read length and median insert size",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

system(paste(samtools.path,"view -f 2","-@",thread, input.bam,"|","cut -f10","|","head",">",paste0(output.bam,".tmp0")))

readlen=median(as.numeric(sapply(fread(paste0(output.bam,".tmp0"),header=F),nchar)))

system(paste(samtools.path,"view -f 2","-@",thread,input.bam,"|","cut -f3,9",">",paste0(output.bam,".tmp1")))

dat=fread(paste0(output.bam,".tmp1"))
is=abs(as.numeric(dat$V2))
medis=median(is)
chr=unique(dat$V1)
chr=chr[chr!='chrM'] ##w/o chrM

pdf(paste0(output.bam,".pdf"))
hist(is[is<1500],breaks=1500,main =paste0("Distribution of insert size(","read length:",readlen,")"),xlab ="Insert size")
abline(v=medis,col="red")
dev.off()

dat=data.table(medis,readlen)
write.table(dat,paste0(output.bam,".info"),sep="\t",quote=F,row.names=F)

system(paste('rm -rf',paste0(output.bam,'*tmp*')))

###################################################################
### Fusion detection using information of discordant pair reads ###
###################################################################

print(paste("Fusion detection using information of discordant pair reads",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

######################
### BAM separation ###
######################

print(paste("BAM separation",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

pb <- txtProgressBar(min = 0, max = length(chr), style = 3)
out=foreach(c=1:length(chr)) %dopar% { 
  setTxtProgressBar(pb, c)
  system( paste(samtools.path,"view -bh -F 12",input.bam,chr[c],">",paste0(output.bam,'.',chr[c])) )
}
close(pb)

##########################################
### Identification of candidate breaks ###
##########################################

print(paste("Identification of candidate breaks",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

pb <- txtProgressBar(min = 0, max = length(chr), style = 3)
out=foreach(c=1:length(chr)) %dopar% { 
  setTxtProgressBar(pb, c)
  
  system(paste(samtools.path,"view -F 2",paste0(output.bam,'.',chr[c]),"|","awk '$6~/[S]/'","|","cut -f4,6",">",paste0(output.bam,'.',c,".tmp0")))
  system(paste(samtools.path,"view -f 256",paste0(output.bam,'.',chr[c]),"|","cut -f4,6",">>",paste0(output.bam,'.',c,".tmp0")))
  system(paste(samtools.path,"view",paste0(output.bam,'.',chr[c]),"|","awk '$6~/[SH]/'","|","grep 'SA:Z'","|","cut -f4,6",">>",paste0(output.bam,'.',c,".tmp0")))
  
  dat=unique(fread(paste0(output.bam,'.',c,".tmp0")))
  V3=(as.numeric(dat$V1)+as.numeric(sapply(dat$V2,endfun))-1)
  dat=data.table(dat,V3) 
  
  sta.br=data.table(table(dat$V1[grepl('[SH][[:digit:]]+M',dat$V2)]),0)
  end.br=data.table(table(dat$V3[grepl('M[[:digit:]]+[SH]',dat$V2)]),1)
  sh.br=rbind(sta.br,end.br)
  sh.br=sh.br[order(as.numeric(sh.br$V1)),]
  sh.br=data.table(chr[c],sh.br[,V1],sh.br[,V2],sh.br[,N])
  
  write.table(sh.br,paste0(output.bam,'.',c,".bed1.tmp"),sep="\t",quote=F,row.names=F,col.names=F)
}
close(pb)

for(c in 1:length(chr)){
  if(c==1){system(paste('cat',paste0(output.bam,'.',c,".bed1.tmp"),">",paste0(output.bam,".bed1")))}
  if(c!=1){system(paste('cat',paste0(output.bam,'.',c,".bed1.tmp"),">>",paste0(output.bam,".bed1")))}
}

system(paste('rm -rf',paste0(output.bam,'*tmp*')))

################################################
### Filtering breaks without a counter break ###
################################################

print(paste("Filtering breaks without a counter break",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

bed=fread(paste0(output.bam,".bed1"))
out=mapply(gapfun,bed$V1,bed$V2)
bed=bed[out==0,]
bed.ran=ranfun(bed$V1,10000)

pb <- txtProgressBar(min = 0, max = length(bed.ran), style = 3)
out=foreach(b=1:length(bed.ran)) %dopar% { 
  setTxtProgressBar(pb, b)
  
  bchr=gsub('/.+$','',bed.ran[b])
  pbed=pbedfun(bed.ran[b])
  write.table(pbed,paste0(output.bam,'.',b,".bed1.tmp"),sep="\t",quote=F,row.names=F,col.names=F)
  
  input.bed=paste0(output.bam,'.',b,".bed1.tmp")
  system(paste(samtools.path,"view -bh -F 2 -L",input.bed,paste0(output.bam,'.',bchr),">",paste0(output.bam,'.',b,".tmp0")))
  system(paste(samtools.path,"view -F 16",paste0(output.bam,'.',b,".tmp0"),"|","cut -f4,6,7,8",">",paste0(output.bam,'.',b,".tmp0f")))
  system(paste(samtools.path,"view -f 16",paste0(output.bam,'.',b,".tmp0"),"|","cut -f4,6,7,8",">",paste0(output.bam,'.',b,".tmp0r")))
  system(paste(samtools.path,"view -f 256 -L",input.bed,paste0(output.bam,'.',bchr),"|","cut -f4,6,7,8",">>",paste0(output.bam,'.',b,".tmp0f"))) 
  
  if(file.info(paste0(output.bam,'.',b,".tmp0f"))$size != 0){
    fdat=data.table(fread(paste0(output.bam,'.',b,'.tmp0f')),1)
    colnames(fdat)=paste0('V',c(1:ncol(fdat)))
  }else{fdat=c()}
  if(file.info(paste0(output.bam,'.',b,".tmp0r"))$size != 0){
    rdat=data.table(fread(paste0(output.bam,'.',b,'.tmp0r')),0)
    colnames(rdat)=paste0('V',c(1:ncol(rdat)))
  }else{rdat=c()}
  dat=rbind(fdat,rdat)
  
  system(paste(samtools.path,"view -L",input.bed,paste0(output.bam,'.',bchr),"|","awk '$6~/[SH]/'","|","grep 'SA:Z'",">",paste0(output.bam,'.',b,'.tmp1')))
  if(file.info(paste0(output.bam,'.',b,'.tmp1'))$size != 0){
    system(paste('cat',paste0(output.bam,'.',b,'.tmp1'),"|","cut -f4,6,7,8",">",paste0(output.bam,'.',b,'.tmp10')))
    system(paste('cat',paste0(output.bam,'.',b,'.tmp1'),"|","cut -f12-",">",paste0(output.bam,'.',b,'.tmp11')))
    
    pdat=fread(paste0(output.bam,'.',b,'.tmp10'))
    tag=fread(paste0(output.bam,'.',b,'.tmp11'),sep='$',header=F)
    satag=as.vector(sapply(tag,tagfun))
    V3=gsub(',[[:print:]]+$','',satag)
    V4=as.numeric(gsub(',[[:print:]]+$','',gsub('^[[:alnum:]]+,','',satag)))
    tag.dat=data.table(pdat[grepl('SA:Z',tag),c('V3','V4'):=NULL],V3,V4)
    V5=rep(0,nrow(tag.dat))
    tag.dat=data.table(tag.dat,V5)
  }else{tag.dat=c()}
  tdat=unique(rbind(dat,tag.dat))
  tdat$V5[grepl('M[[:digit:]]+[SH]',tdat$V2)]=1
  tdat$V5[grepl('[SH][[:digit:]]+M',tdat$V2)]=0
  tdat$V5[grepl('[SH].+M.+[SH]',tdat$V2)]=2
  tdat$V3[tdat$V3=='=']=bchr
  V6=as.numeric(tdat$V1)+as.numeric(sapply(tdat$V2,endfun))-1
  tdat=data.table(tdat,V6)
  
  pairfun=function(chr,pos,ori){
    if(ori==0){index= tdat$V5!=1 & ((!grepl('[SH]',tdat$V2) & (pos-2) < tdat$V1 & tdat$V1 < (pos+medis*2)) | (grepl('[SH]',tdat$V2) & (pos==tdat$V1))) }
    if(ori==1){index= tdat$V5!=0 & ((!grepl('[SH]',tdat$V2) & (pos-medis*2) < tdat$V6 & tdat$V6 < (pos+2)) | (grepl('[SH]',tdat$V2) & (pos==tdat$V6))) }
    pdat=tdat[index,]
    pdat=pdat[ !(chr==pdat$V3 & (pos-medis*2) < pdat$V4 & pdat$V4 < (pos+medis*2)),]
    pdat=pdat[!duplicated(paste(pdat$V1,pdat$V6)),]
    tc=table(pdat$V3)
    tc=names(tc)[tc>=2]
    out=0
    if(length(tc)!=0){
      pdat=pdat[order(pdat$V3,pdat$V4),]
      overfun=function(tc){
        cpdat=pdat[pdat$V3==tc,]
        ind.no=which(abs(c(cpdat$V4[-1],-(medis*2))-cpdat$V4) < medis*2)
        cpdat=cpdat[unique(c(ind.no,ind.no+1)),]
        out=0
        if(sum(grepl('[SH]',cpdat$V2))!=0){out=nrow(cpdat)}
        return(out)
      }
      out=sum(sapply(tc,overfun))
    }
    return(out)
  }
  
  out=mapply(pairfun,bchr,pbed$V4,pbed$V5)
  pbed=data.table(pbed[,c('V2','V3'):=NULL],out)
  write.table(pbed[out>=2,],paste0(output.bam,'.',b,".bed2.tmp"),sep="\t",quote=F,row.names=F,col.names=F)
}  
close(pb) 

for( b in 1:length(bed.ran)){
  if(b==1){system(paste('cat',paste0(output.bam,'.',b,".bed2.tmp"),">",paste0(output.bam,".bed2")))}
  if(b!=1){system(paste('cat',paste0(output.bam,'.',b,".bed2.tmp"),">>",paste0(output.bam,".bed2")))}         
}

system(paste('rm -rf',paste0(output.bam,'*tmp*'))) 

#############################################
### Generation of representative sequence ###
#############################################

print(paste("Generation of representative sequence",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

if(file.info(paste0(output.bam,".bed2"))$size != 0){
  bed=fread(paste0(output.bam,".bed2"))
  bed.ran=ranfun(bed$V1,10000)
  
  pb <- txtProgressBar(min = 0, max = length(bed.ran), style = 3)
  out=foreach(b=1:length(bed.ran)) %dopar% { 
    setTxtProgressBar(pb, b)
    
    bchr=gsub('/.+$','',bed.ran[b])
    pbed=pbedfun(bed.ran[b])
    write.table(pbed,paste0(output.bam,'.',b,".bed2.tmp"),sep="\t",quote=F,row.names=F,col.names=F)
    
    input.bed=paste0(output.bam,'.',b,".bed2.tmp")
    system(paste(samtools.path,"view -L",input.bed,paste0(output.bam,'.',bchr),"|","awk '$6~/[SH]/'","|","cut -f4,6,10",">",paste0(output.bam,'.',b,'.tmp')))
    
    dat=unique(fread(paste0(output.bam,'.',b,'.tmp')))
    V4=as.numeric(dat$V1)+as.numeric(sapply(dat$V2,endfun))-1
    dat=data.table(dat,V4)
    
    mmfun=function(pos,ori){
      out=NA
      if(ori==0){
        index=(pos==dat$V1)&grepl('[SH][[:digit:]]+M',dat$V2)
        scbase=mapply(scbasefun.sm,dat$V2[index],dat$V3[index])
      }
      if(ori==1){
        index=(pos==dat$V4)&grepl('M[[:digit:]]+[SH]',dat$V2)
        scbase=mapply(scbasefun.ms,dat$V2[index],dat$V3[index])
      }
      if( max(nchar(scbase)) >= 10 ){
        sccom=seq(1:max(nchar(scbase)))
        funseqcom=function(no){
          a=substr(scbase,no,no)
          a=a[a!=""]
          out=sum(a==names(table(a))[table(a)==max(table(a))][1])
          return(out)
        }
        if( sum(sapply(sccom,funseqcom))/sum(nchar(scbase)) > MR){ ## filtering MR
          if(ori==0){mbase=mapply(mbasefun.sm,dat$V2[index],dat$V3[index])}
          if(ori==1){mbase=mapply(mbasefun.ms,dat$V2[index],dat$V3[index])}
          sc.rep=repbasefun(scbase)
          m.rep=repbasefun(mbase)
          out=paste(gsub('[GC]',0,gsub('[AT]',1,sc.rep)),gsub('[GC]',0,gsub('[AT]',1,m.rep)),sep='/')
        }
      }
      if( max(nchar(scbase)) < 10 & sum(scbase=='')!=0){ ## for hard clipping
        if(ori==0){mbase=mapply(mbasefun.sm,dat$V2[index],dat$V3[index])}
        if(ori==1){mbase=mapply(mbasefun.ms,dat$V2[index],dat$V3[index])}
        m.rep=repbasefun(mbase)
        out=paste('2',gsub('[GC]',0,gsub('[AT]',1,m.rep)),sep='/')
      }
      return(out)
    }
    
    out=mapply(mmfun,pbed$V4,pbed$V5)
    pbed=data.table(pbed[,c('V2','V3'):=NULL],out)
    write.table(pbed[!is.na(out),],paste0(output.bam,'.',b,".bed3.tmp"),sep="\t",quote=F,row.names=F,col.names=F)
  }
  close(pb) 
  
  for( b in 1:length(bed.ran)){
    if(b==1){system(paste('cat',paste0(output.bam,'.',b,".bed3.tmp"),">",paste0(output.bam,".bed3")))}
    if(b!=1){system(paste('cat',paste0(output.bam,'.',b,".bed3.tmp"),">>",paste0(output.bam,".bed3")))}         
  }
  
  system(paste('rm -rf',paste0(output.bam,'*tmp*')))
}

###########################
### Gathering read name ###
###########################

print(paste("Gathering read name",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

if(file.info(paste0(output.bam,".bed3"))$size != 0){
  bed=fread(paste0(output.bam,".bed3"))
  bed.ran=ranfun(bed$V1,10000)
  
  pb <- txtProgressBar(min = 0, max = length(bed.ran), style = 3)
  out=foreach(b=1:length(bed.ran)) %dopar% { 
    setTxtProgressBar(pb, b)
    
    bchr=gsub('/.+$','',bed.ran[b])
    pbed=pbedfun(bed.ran[b])
    write.table(pbed,paste0(output.bam,'.',b,".bed3.tmp"),sep="\t",quote=F,row.names=F,col.names=F)
    
    input.bed=paste0(output.bam,'.',b,".bed3.tmp")
    system(paste(samtools.path,"view -bh -F 2 -L",input.bed,paste0(output.bam,'.',bchr),">",paste0(output.bam,'.',b,".tmp0")))
    system(paste(samtools.path,"view -F 16",paste0(output.bam,'.',b,".tmp0"),"|","cut -f1,4,6,7,8",">",paste0(output.bam,'.',b,".tmp0f")))
    system(paste(samtools.path,"view -f 16",paste0(output.bam,'.',b,".tmp0"),"|","cut -f1,4,6,7,8",">",paste0(output.bam,'.',b,".tmp0r")))
    system(paste(samtools.path,"view -f 256 -L",input.bed,paste0(output.bam,'.',bchr),"|","cut -f1,4,6,7,8",">>",paste0(output.bam,'.',b,".tmp0f")))     

    if(file.info(paste0(output.bam,'.',b,".tmp0f"))$size != 0){
      fdat=data.table(fread(paste0(output.bam,'.',b,'.tmp0f')),1)
      colnames(fdat)=paste0('V',c(1:ncol(fdat)))
    }else{fdat=c()}
    if(file.info(paste0(output.bam,'.',b,".tmp0r"))$size != 0){
      rdat=data.table(fread(paste0(output.bam,'.',b,'.tmp0r')),0)
      colnames(rdat)=paste0('V',c(1:ncol(rdat)))
    }else{rdat=c()}
    dat=rbind(fdat,rdat)
  
    system(paste(samtools.path,"view -L",input.bed,paste0(output.bam,'.',bchr),"|","awk '$6~/[SH]/'","|","grep 'SA:Z'",">",paste0(output.bam,'.',b,'.tmp1')))
    if(file.info(paste0(output.bam,'.',b,'.tmp1'))$size != 0){
      system(paste('cat',paste0(output.bam,'.',b,'.tmp1'),"|","cut -f1,4,6,7,8",">",paste0(output.bam,'.',b,'.tmp10')))
      system(paste('cat',paste0(output.bam,'.',b,'.tmp1'),"|","cut -f12-",">",paste0(output.bam,'.',b,'.tmp11')))
      
      pdat=fread(paste0(output.bam,'.',b,'.tmp10'))
      tag=fread(paste0(output.bam,'.',b,'.tmp11'),sep='$',header=F)
      satag=as.vector(sapply(tag,tagfun))
      V4=gsub(',[[:print:]]+$','',satag)
      V5=as.numeric(gsub(',[[:print:]]+$','',gsub('^[[:alnum:]]+,','',satag)))
      tag.dat=data.table(pdat[grepl('SA:Z',tag),c('V4','V5'):=NULL],V4,V5)
      V6=rep(0,nrow(tag.dat))
      tag.dat=data.table(tag.dat,V6)
    }else{tag.dat=c()}
    tdat=unique(rbind(dat,tag.dat)) 
    tdat$V6[grepl('M[[:digit:]]+[SH]',tdat$V3)]=1
    tdat$V6[grepl('[SH][[:digit:]]+M',tdat$V3)]=0 
    tdat$V6[grepl('[SH].+M.+[SH]',tdat$V3)]=2
    tdat$V4[tdat$V4=='=']=bchr
    V7=as.numeric(tdat$V2)+as.numeric(sapply(tdat$V3,endfun))-1
    tdat=data.table(tdat,V7)
  
    namefun=function(chr,pos,ori){
      if(ori==0){index= tdat$V6!=1 & ((!grepl('[SH]',tdat$V3) & (pos-2) < tdat$V2 & tdat$V2 < (pos+medis*2)) | (grepl('[SH]',tdat$V3) & (pos==tdat$V2))) }
      if(ori==1){index= tdat$V6!=0 & ((!grepl('[SH]',tdat$V3) & (pos-medis*2) < tdat$V7 & tdat$V7 < (pos+2)) | (grepl('[SH]',tdat$V3) & (pos==tdat$V7))) }
      pdat=tdat[index,]
      pdat=pdat[ !(chr==pdat$V4 & (pos-medis*2) < pdat$V5 & pdat$V5 < (pos+medis*2)),]
      pdat=pdat[order(match(pdat$V1,names(sort(table(pdat$V1),decreasing=T)))),]
      pdat=pdat[!duplicated(paste(pdat$V2,pdat$V7)),]
      tc=table(pdat$V4)
      tc=names(tc)[tc>=2]
      out=''
      if(length(tc)!=0){
        pdat=pdat[order(pdat$V4,pdat$V5),]
        overfun=function(tc){
          cpdat=pdat[pdat$V4==tc,]
          ind.no=which(abs(c(cpdat$V5[-1],-(medis*2))-cpdat$V5) < medis*2)
          name=unique(cpdat$V1[unique(c(ind.no,ind.no+1))])
          out=''
          if(length(name)>=2 & sum(grepl('[HS]',cpdat$V3[unique(c(ind.no,ind.no+1))]))!=0){
            out=paste(name,collapse='/')
          }
          return(out)
        }
        out=sapply(tc,overfun)
        out=paste(out[out!=''],collapse='/')
      }
      return(out)
    }
    
    out=mapply(namefun,bchr,pbed$V4,pbed$V5)
    bbed=data.table(pbed,out)
    bbed=bbed[out!='',]
    write.table(bbed[,c('V2','V3'):=NULL],paste0(output.bam,'.',b,".bed4.tmp1"),sep="\t",quote=F,row.names=F,col.names=F)
    
    if(sum(out!='')!=0){
      pbed=pbed[,c('V2','V3'):=NULL]
      pbed=pbed[as.logical(out!=''),] #degug 20170518
      out=out[as.logical(out!='')] #degug 20170518
      for(i in 1:nrow(pbed)){
        ppbed=data.table(pbed[i,],unlist(strsplit(out[i],'/')))
        if(i==1){nbed=ppbed}else{nbed=rbind(nbed,ppbed)}
      }
      V3=paste(nbed$V1,nbed$V4,nbed$V5,nbed$V8,nbed$V6+nbed$V7)
      nbed=data.table(nbed$V2,V3)
    }else{nbed=NULL}
    write.table(nbed,paste0(output.bam,'.',b,".bed4.tmp2"),sep="\t",quote=F,row.names=F,col.names=F)
  }
  close(pb) 
  
  for( b in 1:length(bed.ran)){
    if(b==1){
      system(paste('cat',paste0(output.bam,'.',b,".bed4.tmp1"),">",paste0(output.bam,".bed41")))
      system(paste('cat',paste0(output.bam,'.',b,".bed4.tmp2"),">",paste0(output.bam,".bed42")))
    }
    if(b!=1){
      system(paste('cat',paste0(output.bam,'.',b,".bed4.tmp1"),">>",paste0(output.bam,".bed41")))
      system(paste('cat',paste0(output.bam,'.',b,".bed4.tmp2"),">>",paste0(output.bam,".bed42")))
    }         
  }

  system(paste('rm -rf',paste0(output.bam,'*tmp*'))) 
}

##########################################
### Sequence comparison between breaks ###
##########################################

print(paste("Sequence comparison between breaks",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

if(file.info(paste0(output.bam,".bed41"))$size != 0){
  bed=fread(paste0(output.bam,".bed41"))
  
  pb <- txtProgressBar(min = 0, max = nrow(bed), style = 3)
  count.br=foreach(b=1:nrow(bed)) %dopar% {
    setTxtProgressBar(pb, b)
    count.br=NA
    
    name=unique(unlist(strsplit(bed$V7[bed$V1==bed$V1[b] & bed$V2==bed$V2[b]  & bed$V3==bed$V3[b]],'/')))
        
    bedfun=function(name){
      out=list(system(paste('cat',paste0(output.bam,".bed42"),"|",paste0("awk '$1==\"",name,"\"'"),"|","cut -f2"),intern = TRUE))
      return(out)
    }

    bed.can=unlist(sapply(name,bedfun))
    chr.str=gsub('[[:space:]].+','',bed.can)
    pos.str=as.numeric(gsub('^.+[[:space:]]','',gsub('[[:space:]][01][[:space:]].+$','',bed.can)))
    bed.can=bed.can[ !(chr.str==bed$V1[b] & (bed$V2[b]-medis*2)< pos.str & pos.str < (bed$V2[b]+medis*2)) ]
    
    if(sum(!is.na(bed.can))!=0 ){    
      if( max(table(bed.can)) >=2){
        basecomfun=function(bed.str){
          out=NA
          bed.split=unlist(strsplit(bed.str,' '))        
          rbase=unlist(strsplit(bed$V6[b],'/'))
          pbase=unlist(strsplit(bed.split[4],'/'))
          if(slidefun(rbase[1],pbase[2])!='n' & slidefun(pbase[1],rbase[2])!='n'){out=paste(bed.split[1],bed.split[2],bed.split[3],sep='/')}
          return(out)
        } 
        can=sapply(unique(bed.can[!is.na(bed.can)]),basecomfun)
        if( sum(!is.na(can))!=0 ){
          can=can[!is.na(can)]
          can=data.table(names(can),can)
          tab=data.table(table(bed.can))
          colnames(can)=colnames(tab)=c('V1','V2')
          can.tab=merge(can,tab,by='V1')
          colnames(can.tab)=c('V1','V2','V3')
          sd=as.numeric(gsub('^.+[[:space:]]','',can.tab$V1))
          count.br=can.tab$V2[ can.tab$V3==max(can.tab$V3) & sd==max(sd)][1]
        }
      }
    }
    count.br
  }
  close(pb)
  
  bed=data.table(bed[,c('V6','V7'):=NULL],unlist(count.br))
  bed=bed[!is.na(unlist(count.br)),]
  
  if(nrow(bed)!=0){write.table(bed,paste0(output.bam,".bed5"),sep="\t",quote=F,row.names=F,col.names=F)}
  
  system(paste('rm -rf',paste0(output.bam,'*tmp*'))) 
}

system(paste('rm -rf',paste0(output.bam,'*chr*')))

###############################################################
### Fusion detection using information of proper pair reads ###
###############################################################

print(paste("Fusion detection using information of proper pair reads",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

######################
### BAM separation ###
######################

print(paste("BAM separation",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

pb <- txtProgressBar(min = 0, max = length(chr), style = 3)
out=foreach(c=1:length(chr)) %dopar% { 
  setTxtProgressBar(pb, c)
  system( paste(samtools.path,"view -H",input.bam,">",paste0(output.bam,'.',chr[c],'.tmp')) )
  system( paste(samtools.path,"view -f 2",input.bam,chr[c],"|","awk '$6~/[S]/'",">>",paste0(output.bam,'.',chr[c],'.tmp')) )
  system( paste(samtools.path,"view -Sbh",paste0(output.bam,'.',chr[c],'.tmp'),">",paste0(output.bam,'.',chr[c])) )
}
close(pb)

system(paste('rm -rf',paste0(output.bam,'*tmp*'))) 

##########################################
### Identification of candidate breaks ###
##########################################

print(paste("Identification of candidate breaks",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

pb <- txtProgressBar(min = 0, max = length(chr), style = 3)
out=foreach(c=1:length(chr)) %dopar% { 
  setTxtProgressBar(pb, c)
  
  system(paste(samtools.path,"view",paste0(output.bam,'.',chr[c]),"|","cut -f4,6",">",paste0(output.bam,'.',c,".tmp0")))

  dat=unique(fread(paste0(output.bam,'.',c,".tmp0")))
  V3=(as.numeric(dat$V1)+as.numeric(sapply(dat$V2,endfun))-1)
  dat=data.table(dat,V3) 
  
  sta.br=data.table(table(dat$V1[grepl('[S][[:digit:]]+M',dat$V2)]),0)
  end.br=data.table(table(dat$V3[grepl('M[[:digit:]]+[S]',dat$V2)]),1)
  sh.br=rbind(sta.br,end.br)
  sh.br=sh.br[order(as.numeric(sh.br$V1)),]
  sh.br=data.table(chr[c],sh.br[,V1],sh.br[,V2],sh.br[,N])
  sh.br=sh.br[sh.br$V4>=3,]

  write.table(sh.br,paste0(output.bam,'.',c,".bed6.tmp"),sep="\t",quote=F,row.names=F,col.names=F)
}
close(pb)

for(c in 1:length(chr)){
  if(c==1){system(paste('cat',paste0(output.bam,'.',c,".bed6.tmp"),">",paste0(output.bam,".bed6")))}
  if(c!=1){system(paste('cat',paste0(output.bam,'.',c,".bed6.tmp"),">>",paste0(output.bam,".bed6")))}
}

system(paste('rm -rf',paste0(output.bam,'*tmp*')))

##############################
### Filtering duplications ###
##############################

print(paste("Filtering duplications",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

if(file.info(paste0(output.bam,".bed6"))$size != 0){
  bed=fread(paste0(output.bam,".bed6"))
  out=mapply(gapfun,bed$V1,bed$V2)
  bed=bed[out==0,]
  bed.ran=ranfun(bed$V1,10000)
  
  pb <- txtProgressBar(min = 0, max = length(bed.ran), style = 3)
  out=foreach(b=1:length(bed.ran)) %dopar% { 
    setTxtProgressBar(pb, b)
    
    bchr=gsub('/.+$','',bed.ran[b])
    pbed=pbedfun(bed.ran[b])
    write.table(pbed,paste0(output.bam,'.',b,".bed6.tmp"),sep="\t",quote=F,row.names=F,col.names=F)
    
    input.bed=paste0(output.bam,'.',b,".bed6.tmp")
    system(paste(samtools.path,"view -L",input.bed,paste0(output.bam,'.',bchr),"|","cut -f1,4,6,9",">",paste0(output.bam,'.',b,".tmp")))
 
    dat=unique(fread(paste0(output.bam,'.',b,".tmp")))
    V5=rep(0,nrow(dat))
    dat=data.table(dat,V5)
    dat$V5[grepl('M[[:digit:]]+[S]',dat$V3)]=1
    dat$V5[grepl('[S][[:digit:]]+M',dat$V3)]=0 
    dat$V5[grepl('[S].+M.+[S]',dat$V3)]=2
    V6=(as.numeric(dat$V2)+as.numeric(sapply(dat$V3,endfun))-1)
    tdat=data.table(dat,V6) 
    
    namefun=function(pos,ori){
      if(ori==0){index= tdat$V5!=1 & pos==tdat$V2 }
      if(ori==1){index= tdat$V5!=0 & pos==tdat$V6 }
      pdat=tdat[index,]
      pdat=pdat[order(match(pdat$V1,names(sort(table(pdat$V1),decreasing=T)))),]
      pdat=pdat[!duplicated(paste(pdat$V2,pdat$V6)),]
      out=length(unique(pdat$V1))
      return(out)
    }

    out=mapply(namefun,pbed$V4,pbed$V5)
    bbed=data.table(pbed,out)
    bbed=bbed[out>=3,]
    write.table(bbed[,c('V2','V3'):=NULL],paste0(output.bam,'.',b,".bed7.tmp"),sep="\t",quote=F,row.names=F,col.names=F)
  }
  close(pb) 
  
  for( b in 1:length(bed.ran)){
    if(b==1){system(paste('cat',paste0(output.bam,'.',b,".bed7.tmp"),">",paste0(output.bam,".bed7")))}
    if(b!=1){system(paste('cat',paste0(output.bam,'.',b,".bed7.tmp"),">>",paste0(output.bam,".bed7")))}         
  }
  
  system(paste('rm -rf',paste0(output.bam,'*tmp*'))) 
}

#############################################
### Generation of representative sequence ###
#############################################

print(paste("Generation of representative sequence",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

if(file.info(paste0(output.bam,".bed7"))$size != 0){
  bed=fread(paste0(output.bam,".bed7"))
  bed.ran=ranfun(bed$V1,10000)
  
  pb <- txtProgressBar(min = 0, max = length(bed.ran), style = 3)
  out=foreach(b=1:length(bed.ran)) %dopar% { 
    setTxtProgressBar(pb, b)
    
    bchr=gsub('/.+$','',bed.ran[b])
    pbed=pbedfun(bed.ran[b])
    write.table(pbed,paste0(output.bam,'.',b,".bed7.tmp"),sep="\t",quote=F,row.names=F,col.names=F)
    
    input.bed=paste0(output.bam,'.',b,".bed7.tmp")
    system(paste(samtools.path,"view -L",input.bed,paste0(output.bam,'.',bchr),"|","cut -f4,6,10",">",paste0(output.bam,'.',b,'.tmp')))
    
    dat=unique(fread(paste0(output.bam,'.',b,'.tmp')))
    V4=as.numeric(dat$V1)+as.numeric(sapply(dat$V2,endfun))-1
    dat=data.table(dat,V4)
    
    mmfun=function(pos,ori){
      out=NA
      if(ori==0){
        index=(pos==dat$V1)&grepl('[S][[:digit:]]+M',dat$V2)
        scbase=mapply(scbasefun.sm,dat$V2[index],dat$V3[index])
      }
      if(ori==1){
        index=(pos==dat$V4)&grepl('M[[:digit:]]+[S]',dat$V2)
        scbase=mapply(scbasefun.ms,dat$V2[index],dat$V3[index])
      }
      if( sum(nchar(scbase) >= 10)>=3 ){ # 3>= sc with >=10 base
        sccom=seq(1:max(nchar(scbase)))
        funseqcom=function(no){
          a=substr(scbase,no,no)
          a=a[a!=""]
          out=sum(a==names(table(a))[table(a)==max(table(a))][1])
          return(out)
        }
        if( sum(sapply(sccom,funseqcom))/sum(nchar(scbase)) > MR){ ## filtering MR
          if(ori==0){mbase=mapply(mbasefun.sm,dat$V2[index],dat$V3[index])}
          if(ori==1){mbase=mapply(mbasefun.ms,dat$V2[index],dat$V3[index])}
          sc.rep=repbasefun(scbase)
          m.rep=repbasefun(mbase)
          out=paste(gsub('[GC]',0,gsub('[AT]',1,sc.rep)),gsub('[GC]',0,gsub('[AT]',1,m.rep)),sep='/')
        }
      }
      return(out)
    }
    
    out=mapply(mmfun,pbed$V4,pbed$V5)
    pbed=data.table(pbed[,c('V2','V3'):=NULL],out)
    write.table(pbed[!is.na(out),],paste0(output.bam,'.',b,".bed8.tmp"),sep="\t",quote=F,row.names=F,col.names=F)
  }
  close(pb) 
  
  for( b in 1:length(bed.ran)){
    if(b==1){system(paste('cat',paste0(output.bam,'.',b,".bed8.tmp"),">",paste0(output.bam,".bed8")))}
    if(b!=1){system(paste('cat',paste0(output.bam,'.',b,".bed8.tmp"),">>",paste0(output.bam,".bed8")))}         
  }
  
  system(paste('rm -rf',paste0(output.bam,'*tmp*')))
}

##########################################
### Sequence comparison between breaks ###
##########################################

print(paste("sequence comparison between breaks",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

if(file.info(paste0(output.bam,".bed8"))$size != 0){
  bed=fread(paste0(output.bam,".bed8"))
  bed$V6=paste(bed$V1,bed$V2,bed$V3,bed$V6,bed$V5)
  
  pb <- txtProgressBar(min = 0, max = nrow(bed), style = 3)
  count.br=foreach(b=1:nrow(bed)) %dopar% {
    setTxtProgressBar(pb, b)
    count.br=NA

    bed.can=bed$V6[ bed$V1==bed$V1[b] & bed$V2[b]-medis*3 <= bed$V2 & bed$V2 <= bed$V2[b]+medis*3]
    bed.can=bed.can[abs(as.numeric(gsub('[[:space:]].+$','',gsub('^[[:alnum:]]+[[:space:]]','',bed.can)))-bed$V2[b]) >10] #length of sc >=10
   
    if(length(bed.can)!=0 ){    
        basecomfun=function(bed.str){
          out=NA
          rbed.split=unlist(strsplit(bed$V6[b],' ')) 
          pbed.split=unlist(strsplit(bed.str,' '))        
          rbase=unlist(strsplit(rbed.split[4],'/'))
          pbase=unlist(strsplit(pbed.split[4],'/'))
          if(slidefun(rbase[1],pbase[2])!='n' & slidefun(pbase[1],rbase[2])!='n'){out=paste(pbed.split[1],pbed.split[2],pbed.split[3],sep='/')}
          return(out)
        } 
        can=sapply(bed.can,basecomfun)
        if( sum(!is.na(can))!=0 ){
          can=can[!is.na(can)]
          no.can=as.numeric(gsub('^.+[[:space:]]','',names(can)))
          count.br=can[ no.can==max(no.can)][1]
        }
    }
    count.br
  }
  close(pb)
  
  bed=data.table(bed[,c('V6'):=NULL],unlist(count.br))
  bed=bed[!is.na(unlist(count.br)),]
  
  if(nrow(bed)!=0){write.table(bed,paste0(output.bam,".bed9"),sep="\t",quote=F,row.names=F,col.names=F)}
  
  system(paste('rm -rf',paste0(output.bam,'*tmp*'))) 
}
system(paste('rm -rf',paste0(output.bam,'*chr*')))

############################################
### Measuring uniqueness of sequence ######
############################################

print(paste(" Measuring uniqueness of sequence ",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

if(file.info(paste0(output.bam,".bed5"))$size != 0 |file.info(paste0(output.bam,".bed9"))$size != 0 ){
  if(file.info(paste0(output.bam,".bed5"))$size != 0){bed1=fread(paste0(output.bam,".bed5"))}else{bed1=c()}
  if(file.info(paste0(output.bam,".bed9"))$size != 0){bed2=fread(paste0(output.bam,".bed9"))}else{bed2=c()}
  bed=rbind(bed1,bed2)
  bed=bed[!duplicated(paste(bed$V1,bed$V2,bed$V3,bed$V6)),]
  bed.ran=ranfun(bed$V1,10000)
  
  pb <- txtProgressBar(min = 0, max = length(bed.ran), style = 3)
  out=foreach(b=1:length(bed.ran)) %dopar% { 
    setTxtProgressBar(pb, b)
    
    bchr=gsub('/.+$','',bed.ran[b])
    pbed=pbedfun(bed.ran[b])
    write.table(pbed,paste0(output.bam,'.',b,".bed9.tmp"),sep="\t",quote=F,row.names=F,col.names=F)
    
    input.bed=paste0(output.bam,'.',b,".bed9.tmp")
    system(paste(samtools.path,"view -L",input.bed,input.bam,"|","awk '$6~/[S]/'","|","cut -f4,6,10",">",paste0(output.bam,'.',b,'.tmp')))
    
    if(file.info(paste0(output.bam,'.',b,'.tmp'))$size != 0){
      dat=unique(fread(paste0(output.bam,'.',b,'.tmp')))
      V4=as.numeric(dat$V1)+as.numeric(sapply(dat$V2,endfun))-1
      dat=data.table(dat,V4)
      
      repfun=function(pos,ori){
        if(ori==0){
          index=(pos==dat$V1)&grepl('[S][[:digit:]]+M',dat$V2)
          scbase=mapply(scbasefun.sm,dat$V2[index],dat$V3[index])
        }
        if(ori==1){
          index=(pos==dat$V4)&grepl('M[[:digit:]]+[S]',dat$V2)
          scbase=mapply(scbasefun.ms,dat$V2[index],dat$V3[index])
        }
        out=NA
        if(sum(nchar(scbase))!=0){out=round(repeat.fun(repbasefun(scbase)),2)}
        return(out)
      }   
      out=mapply(repfun,pbed$V4,pbed$V5)
    }else{out=rep(NA,nrow(pbed))}
    pbed=data.table(pbed[,c('V2','V3'):=NULL],out)
    write.table(pbed,paste0(output.bam,'.',b,".bed10.tmp"),sep="\t",quote=F,row.names=F,col.names=F)
  }
  close(pb)
  
  for( b in 1:length(bed.ran)){
    if(b==1){system(paste('cat',paste0(output.bam,'.',b,".bed10.tmp"),">",paste0(output.bam,".bed10")))}
    if(b!=1){system(paste('cat',paste0(output.bam,'.',b,".bed10.tmp"),">>",paste0(output.bam,".bed10")))}         
  }
  
  system(paste('rm -rf',paste0(output.bam,'*tmp*')))
}

##################################
### annotation & visualization ###
##################################

print(paste("Annotation & visualization",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))

if(file.info(paste0(output.bam,".bed10"))$size != 0 ){
  
  #### vcf generation ####
  
  vcf=read.table(paste0(output.bam,".bed10"),sep='\t',stringsAsFactors = F)
  colnames(vcf)=c('chr','pos','ori','spl','dis','counter','rep')  ## column name assignment
  
  event=c()
  skip=c()
  
  for( i in 1:nrow(vcf)){
    if( sum(skip==i)==0 ){
      a=data.frame(vcf$chr[i],vcf$pos[i],vcf$ori[i],vcf$dis[i],vcf$spl[i],vcf$rep[i],stringsAsFactors = F)
      b=unlist(strsplit(vcf$counter[i],'/'))
      b.no=which(vcf$chr==b[1]&vcf$pos==b[2]&vcf$ori==b[3]&vcf$counter==paste(vcf$chr[i],vcf$pos[i],vcf$ori[i],sep='/'))
      if( length(b.no)!=0 & !(a[1]==b[1] & a[2]==b[2])){
        b=vcf[b.no,]
        skip=c(skip,b.no) ##skip a same combination 
        b=data.frame(b$chr,b$pos,b$ori,b$dis,b$spl,b$rep,stringsAsFactors = F)
        c=data.frame(a,b,stringsAsFactors = F)
        colnames(c)=paste0('V',c(1:ncol(c)))
        
        ##identification of fusion event
        
        if(c$V1!=c$V7){
          Event='Interchromosomal_translocation'
        }
        if(c$V1==c$V7){
          if(c$V3==c$V9){
            Event='Inversion'
          }
          if(c$V3!=c$V9){
            if( c(c$V3,c$V9)[c(c$V2,c$V8)==min(c(c$V2,c$V8))]==1 ){
              Event='Deletion'
            }else{
              Event='Tandem'
            }
          }
        }
        c=data.frame(c,Event,stringsAsFactors = F)
        
        ##identification of gene names
        
        if( length(unique(ref$V4[ref$V3==c$V1 & (ref$V5 <= c$V2 & c$V2 <= ref$V6)]))== 2 ){
          pref=ref[ref$V3==c$V1 & (ref$V5 <= c$V2 & c$V2 <= ref$V6),]
          prefp=names(sort(table(pref$V13[pref$V4=='+']),decreasing=T))[1]
          prefn=names(sort(table(pref$V13[pref$V4=='-']),decreasing=T))[1]
          GeneA=c(prefp,prefn)
          StrandA=c('+','-')
        }
        if( length(unique(ref$V4[ref$V3==c$V1 & (ref$V5 <= c$V2 & c$V2 <= ref$V6)]))== 1 ){
          GeneA=names(sort(table(ref$V13[ref$V3==c$V1 & (ref$V5 <= c$V2 & c$V2 <= ref$V6)]),decreasing=T))[1]
          GeneA=c(GeneA,paste0('Compl_',GeneA))
          str=unique(ref$V4[GeneA[1]==ref$V13])
          StrandA=c(str,c('-','+')[str!=c('-','+')]) 
        }
        if( sum(ref$V3==c$V1 & (ref$V5 <= c$V2 & c$V2 <= ref$V6) )== 0 ){
          pref=ref[ref$V3==c$V1,]
          mn=min(abs(pref$V5 - c$V2),abs(c$V2 - pref$V6))
          bp1=pref[(abs(pref$V5 - c$V2)==mn | abs(c$V2 - pref$V6)==mn),]
          GeneA=c(paste('Flanking',paste(unique(bp1$V13),collapse='&'),sep='_'),paste('Flanking',paste(unique(bp1$V13),collapse='&'),sep='_'))
          StrandA=c('+','-')
        }
        if( length(unique(ref$V4[ref$V3==c$V7 & (ref$V5 <= c$V8 & c$V8 <= ref$V6)]))==2 ){
          pref=ref[ref$V3==c$V7 & (ref$V5 <= c$V8 & c$V8 <= ref$V6),]
          prefp=names(sort(table(pref$V13[pref$V4=='+']),decreasing=T))[1]
          prefn=names(sort(table(pref$V13[pref$V4=='-']),decreasing=T))[1]
          GeneB=c(prefp,prefn)
          StrandB=c('+','-')
        }
        if( length(unique(ref$V4[ref$V3==c$V7 & (ref$V5 <= c$V8 & c$V8 <= ref$V6)]))==1  ){
          GeneB=names(sort(table(ref$V13[ref$V3==c$V7 & (ref$V5 <= c$V8 & c$V8 <= ref$V6)]),decreasing=T))[1]
          GeneB=c(GeneB,paste0('Compl_',GeneB))
          str=unique(ref$V4[GeneB[1]==ref$V13])
          StrandB=c(str,c('-','+')[str!=c('-','+')])      
        }
        if( sum(ref$V3==c$V7 & (ref$V5 <= c$V8 & c$V8 <= ref$V6) )==0  ){
          pref=ref[ref$V3==c$V7,]
          mn=min(abs(pref$V5 - c$V8 ),abs(c$V8  - pref$V6))
          bp1=pref[(abs(pref$V5 - c$V8)==mn | abs(c$V8  - pref$V6)==mn),]
          GeneB=c(paste('Flanking',paste(unique(bp1$V13),collapse='&'),sep='_'),paste('Flanking',paste(unique(bp1$V13),collapse='&'),sep='_'))
          StrandB=c('+','-')
        }
        
        cg1=data.frame(c,GeneA[1],StrandA[1],GeneB[1],StrandB[1],stringsAsFactors=F)
        cg2=data.frame(c,GeneA[2],StrandA[2],GeneB[1],StrandB[1],stringsAsFactors=F)
        cg3=data.frame(c,GeneA[1],StrandA[1],GeneB[2],StrandB[2],stringsAsFactors=F)
        cg4=data.frame(c,GeneA[2],StrandA[2],GeneB[2],StrandB[2],stringsAsFactors=F)
        colnames(cg1)=colnames(cg2)=colnames(cg3)=colnames(cg4)
        cg=rbind(cg1,cg2,cg3,cg4)
        
        colnames(cg)=paste0('V',c(1:17))
        event=rbind(event,cg)
      }
    }
  }
  
  dat=event
  filter=rep('x',nrow(event))
  
  for( i in 1:nrow(dat)){
    if( dat$V3[i]==dat$V9[i] & dat$V15[i]!=dat$V17[i] ){
      filter[i]="o"
    }
    if( dat$V3[i]!=dat$V9[i] & dat$V15[i]==dat$V17[i] ){
      filter[i]="o"
    }
  }
  
  dat=event[filter=='o',]
  write.table(dat,paste0(output,'/',paste0(Sample_ID,'.vcf')),sep='\t',col.names=T,row.names=F,quote=F)
  
  #### annotation ####
  
  for(i in 1:nrow(dat)){
    a=mapply(reffun,dat$V14[i],dat$V2[i])[1]
    b=mapply(reffun,dat$V16[i],dat$V8[i])[1]
    if( dat$V3[i] != dat$V9[i] ){
      if( dat$V15[i] =='+' ){
        gene5p=c(dat$V14[i],dat$V16[i])[c(dat$V3[i],dat$V9[i])==1]
        gene3p=c(dat$V14[i],dat$V16[i])[c(dat$V3[i],dat$V9[i])==0]
        info5p=c(a,b)[c(dat$V3[i],dat$V9[i])==1]
        info3p=c(a,b)[c(dat$V3[i],dat$V9[i])==0]
      }
      if( dat$V15[i] =='-' ){
        gene5p=c(dat$V14[i],dat$V16[i])[c(dat$V3[i],dat$V9[i])==0]
        gene3p=c(dat$V14[i],dat$V16[i])[c(dat$V3[i],dat$V9[i])==1]
        info5p=c(a,b)[c(dat$V3[i],dat$V9[i])==0]
        info3p=c(a,b)[c(dat$V3[i],dat$V9[i])==1]
      }
    }
    if( dat$V3[i] == dat$V9[i] ){
      if( dat$V3[i] ==1 ){
        gene5p=c(dat$V14[i],dat$V16[i])[c(dat$V15[i],dat$V17[i])=='+']
        gene3p=c(dat$V14[i],dat$V16[i])[c(dat$V15[i],dat$V17[i])=='-']
        info5p=c(a,b)[c(dat$V15[i],dat$V17[i])=='+']
        info3p=c(a,b)[c(dat$V15[i],dat$V17[i])=='-']
      }
      if( dat$V3[i] ==0 ){
        gene5p=c(dat$V14[i],dat$V16[i])[c(dat$V15[i],dat$V17[i])=='-']
        gene3p=c(dat$V14[i],dat$V16[i])[c(dat$V15[i],dat$V17[i])=='+']
        info5p=c(a,b)[c(dat$V15[i],dat$V17[i])=='-']
        info3p=c(a,b)[c(dat$V15[i],dat$V17[i])=='+']
      }
    }
    c=paste0(gene5p,'->',gene3p)
    d=NA
    if(grepl('Intron',a) & grepl('Intron',b)){
      if( gsub('^.+Frame','',a)==gsub('^.+Frame','',b) ){
        d='Inframe'
      }else{
        d='Outframe'
      }
    }
    if(sum(grepl('Exon',c(a,b)))==1){
      if( grepl('Intron',info5p) ){
        if( gsub('^.+Frame','',info5p)==gsub('[[:digit:]]+,','',gsub('^.+Frame','',info3p)) ){
          d='Possible_Inframe'
        }else{
          d='Possible_Outframe'
        }
      }
      if( grepl('Exon',info5p) ){
        d='Possible_Outframe'
      }
    }
    e=NA
    if( dat$V14[i]!=dat$V16[i] ){
      if( sum(grepl(dat$V14[i],cos$Translocation.Name) & grepl(dat$V16[i],cos$Translocation.Name))!=0 ){
        g=cos$Translocation.Name[grepl(dat$V14[i],cos$Translocation.Name) & grepl(dat$V16[i],cos$Translocation.Name)]
        g1=names(sort(table(gsub('[{].+$','',g)),decreasing=T))[1]
        if( gsub('->.+$','',c)==g1 ){
          e1=sort(table(cos$Primary.site[grepl(dat$V14[i],cos$Translocation.Name) & grepl(dat$V16[i],cos$Translocation.Name)]),decreasing=T)
          e=paste(paste0(names(e1),'(',e1,')'),collapse ='_')
        }
      }
    }
    
    t='T3'
    if( sum(dat$V14[i]==t1 | dat$V16[i]==t1)!=0 ){
      t='T2'
      if(!(is.na(e))){
        t='T1'
      }
    }
    
    f=data.frame(t,a,b,c,d,e,stringsAsFactors = F)
    
    if(i==1){
      ff=f
    }else{
      ff=rbind(ff,f)
    }
  }
    
  t.dat=data.frame(Sample_ID,ff[,1],dat,ff[,-1],stringsAsFactors = F)
  colnames(t.dat)=c('Sample_ID','Tier','ChrA','BreakA','OriA','DisA','SplitA','RepeatA','ChrB','BreakB','OriB','DisB','SplitB','RepeatB','Event','GeneA','StrGeneA','GeneB','StrGeneB','InfoA','InfoB','Direction','Frame','Cosmic')
  
  tt.dat=t.dat[order(t.dat$DisA+t.dat$DisB,decreasing=T),]
  ft.dat=tt.dat[ !grepl('Compl',tt.dat$GeneA) & !grepl('Compl',tt.dat$GeneB) & !grepl('Flanking',tt.dat$GeneA) & !grepl('Flanking',tt.dat$GeneB) & grepl('NM_',tt.dat$InfoA) & grepl('NM_',tt.dat$InfoB),] #gene-gene fusion
  
  write.table(tt.dat,paste0(output,'/',paste0(Sample_ID,'.annotated.vcf')),sep='\t',col.names=T,row.names=F,quote=F)
  write.table(ft.dat,paste0(output,'/',paste0(Sample_ID,'.annotated.vcf.filtered')),sep='\t',col.names=T,row.names=F,quote=F)
  
  ##Column name change to V1 for SGI
  
  p.dat=ft.dat  
  p.dat$Event=paste(p.dat$Event,p.dat$Frame)
  p.dat=p.dat[,-c(8,14,23)]
  Total_reads=p.dat$DisA+p.dat$DisB
  Totalreads=p.dat$SplitA+p.dat$SplitB
  p.dat=data.frame(p.dat,Totalreads,Total_reads,stringsAsFactors=F)

  gene_data <- read.table(panelgene, sep="\t", header=T)
  v2gene <- as.character(gene_data$Entrez_accession)
  
  csfun=function(GeneA,GeneB){if(sum(GeneA==v2gene)!=0 | sum(GeneB==v2gene)!=0){out='o'}else{out='x'}}
  pp.dat=p.dat=p.dat[mapply(csfun,p.dat$GeneA,p.dat$GeneB)=='o',] ##filtering non target gene
  
  intra=quantile(p.dat$Total_reads[p.dat$GeneA==p.dat$GeneB])[4]+1.5*IQR(p.dat$Total_reads[p.dat$GeneA==p.dat$GeneB])
  inter=quantile(p.dat$Total_reads[p.dat$GeneA!=p.dat$GeneB])[4]+1.5*IQR(p.dat$Total_reads[p.dat$GeneA!=p.dat$GeneB])  
  p.dat=p.dat[p.dat$Tier=='T1' | !is.na(p.dat$Cosmic) | (p.dat$GeneA==p.dat$GeneB & p.dat$Total_reads >= intra) | (p.dat$GeneA!=p.dat$GeneB & p.dat$Total_reads >= inter),] ##filtering under Q3+1.5IQG and Tier2/3 and Novel
  
  p.dat=p.dat[p.dat$GeneA!=p.dat$GeneB,] ##filtering intrageneic events

  ##for couner fusion
  pp.dat=pp.dat[pp.dat$GeneA!=pp.dat$GeneB,]
  rc=paste0(gsub('^.+->','',p.dat$Direction),'->',gsub('->.+$','',p.dat$Direction))
  rc.dat=pp.dat[pp.dat$Direction[pp.dat$GeneA!=pp.dat$GeneB] %in% rc & pp.dat$Total_reads > quantile(pp.dat$Total_reads)[3],] ##more than Q2
  p.dat=rbind(p.dat,rc.dat)
  
  colnames(p.dat)=c("ID","Tier","ChrA","BreakpointA","ReadA","Cover_ReadA","cnt_ReadA","ChrB","BreakpointB",
                     "ReadB","Cover_ReadB","cnt_ReadB","isrealfusion","GeneA","OriA","GeneB","OriB","Read.posA","Read.posB",
                     "Direction","Known","Total_reads","Totalreads")
  p.dat=data.frame(p.dat$ID,p.dat$Tier,p.dat$GeneA,p.dat$ReadA,p.dat$cnt_ReadA,p.dat$GeneB,p.dat$ReadB,
                    p.dat$cnt_ReadB,p.dat$Total_reads,p.dat$OriA,p.dat$OriB,p.dat$Cover_ReadA,p.dat$ChrA,
                    p.dat$BreakpointA,p.dat$Read.posA,p.dat$Cover_ReadB,p.dat$ChrB,p.dat$BreakpointB,
                    p.dat$Read.posB,p.dat$Direction,p.dat$Known,p.dat$isrealfusion,p.dat$Totalreads,
                    stringsAsFactors=F)
  colnames(p.dat)=gsub('p.dat.','',colnames(p.dat))

  p.dat$Known[!is.na(p.dat$Known)]='Cosmic'
  p.dat$Known[is.na(p.dat$Known)]='Novel'

  write.table(p.dat,paste0(output,'/',paste0(Sample_ID,'.SV.final.txt')),sep='\t',col.names=T,row.names=F,quote=F)
  
  #### visualization #### 
  
  if(nrow(ft.dat)!=0 ){
    plots <- list()
    for(p in 1:nrow(ft.dat)){
      g5=unlist(strsplit(ft.dat$Direction[p],'->'))[1]
      g3=unlist(strsplit(ft.dat$Direction[p],'->'))[2]
      if(g5!=g3){
        info5=c(ft.dat$InfoA[p],ft.dat$InfoB[p])[c(ft.dat$GeneA[p],ft.dat$GeneB[p])==g5][1]
        info3=c(ft.dat$InfoA[p],ft.dat$InfoB[p])[c(ft.dat$GeneA[p],ft.dat$GeneB[p])==g3][1]
        bp5=c(ft.dat$BreakA[p],ft.dat$BreakB[p])[c(ft.dat$GeneA[p],ft.dat$GeneB[p])==g5][1]
        bp3=c(ft.dat$BreakA[p],ft.dat$BreakB[p])[c(ft.dat$GeneA[p],ft.dat$GeneB[p])==g3][1]
      }
      if(g5==g3){
        if(ft.dat$StrGeneA[p]=='+'){
          info5=c(ft.dat$InfoA[p],ft.dat$InfoB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==min(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]
          info3=c(ft.dat$InfoA[p],ft.dat$InfoB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==max(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]
          bp5=c(ft.dat$BreakA[p],ft.dat$BreakB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==min(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]
          bp3=c(ft.dat$BreakA[p],ft.dat$BreakB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==max(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]
        }
        if(ft.dat$StrGeneA[p]=='-'){
          info5=c(ft.dat$InfoA[p],ft.dat$InfoB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==max(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]
          info3=c(ft.dat$InfoA[p],ft.dat$InfoB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==min(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]
          bp5=c(ft.dat$BreakA[p],ft.dat$BreakB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==max(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]
          bp3=c(ft.dat$BreakA[p],ft.dat$BreakB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==min(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]        
        }
      }
      
      gd5=data.frame(0.55,mapfun(g5,info5,bp5),stringsAsFactors=F)
      gd3=data.frame(0.55,mapfun(g3,info3,bp3),stringsAsFactors=F)
      gd3$st=gd3$st+0.5
      gd3$en=gd3$en+0.5
      gd3$nbp=gd3$nbp+0.5
      colnames(gd5)[1]=colnames(gd3)[1]='V1'
      
      fgd5=gd5[c(1:as.numeric(gsub('^.+[(]','',gsub('/.+$','',info5)))),]
      if(!grepl('Intron',info5)){
        fgd5$en[as.numeric(gsub('^.+[(]','',gsub('/.+$','',info5)))]=fgd5$nbp[1]
        no5=0.5-fgd5$en[nrow(fgd5)]
        fgd5$st=fgd5$st+no5
        fgd5$en=fgd5$en+no5
      }
      if(grepl('Intron',info5)){
        no5=(0.5-(fgd5$nbp[1]-fgd5$en[as.numeric(gsub('^.+[(]','',gsub('/.+$','',info5)))]))-fgd5$en[nrow(fgd5)]
        fgd5$st=fgd5$st+no5
        fgd5$en=fgd5$en+no5
      }
      fgd3=gd3[c(as.numeric(gsub('^.+[(]','',gsub('/.+$','',info3))):nrow(gd3)),]
      if(!grepl('Intron',info3)){
        fgd3$st[1]=fgd3$nbp[1]
        no3=fgd3$st[1]-0.5
        fgd3$st=fgd3$st-no3
        fgd3$en=fgd3$en-no3
      }
      if(grepl('Intron',info3)){
        fgd3=fgd3[-1,]
        no3=fgd3$nbp[1]-0.5
        fgd3$st=fgd3$st-no3
        fgd3$en=fgd3$en-no3
      }      
      fgd=rbind(fgd5,fgd3)
      fgd$V1=0.35
      
      gd=rbind(gd5,gd3,fgd)  
      if(sum(!is.na(gd$Domain))==0){gd$Domain='NA'}
      gd$Domain=factor(gd$Domain)
      arr1=c(fgd5$nbp[1],0.498,0.52,0.39)
      arr2=c(fgd3$nbp[1],0.502,0.52,0.39)
      arr.dat=data.frame(rbind(arr1,arr2),stringsAsFactors = F)

      text0=c(0.05,0.98,paste('ID:',Sample_ID))     
      text1=c(0.05,0.94,paste('Fusion:',ft.dat$Direction[p]))
      text2=c(0.05,0.90,paste('Frame:',ft.dat$Frame[p]))
      text3=c(0.05,0.86,paste('Event type:',ft.dat$Event[p]))
      text4=c(0.05,0.82,paste('Supporting reads:',ft.dat$DisA[p]+ft.dat$DisB[p]))
      #text5=c(0.05,0.75,paste('Cosmic:',ft.dat$Cosmic[p]))
      text6=c(0.05,0.78,'GeneName/Breakpoint/Strand/FusionInfo:')
      InA=paste0(' ',ft.dat$GeneA[p],'/',ft.dat$ChrA[p],':',ft.dat$BreakA[p],'/',paste0('(',ft.dat$StrGeneA[p],')'),'/',ft.dat$InfoA[p])
      InB=paste0(' ',ft.dat$GeneB[p],'/',ft.dat$ChrB[p],':',ft.dat$BreakB[p],'/',paste0('(',ft.dat$StrGeneB[p],')'),'/',ft.dat$InfoB[p])
      if(g5!=g3){
        text7=c(0.05,0.74,c(InA,InB)[c(ft.dat$InfoA[p],ft.dat$InfoB[p])==info5][1])
        text8=c(0.05,0.7,c(InA,InB)[c(ft.dat$InfoA[p],ft.dat$InfoB[p])==info3][1])
      }
      if(g5==g3){
        if(ft.dat$StrGeneA[p]=='+'){
          text7=c(0.05,0.74,c(InA,InB)[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==min(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))])
          text8=c(0.05,0.7,c(InA,InB)[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==max(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))])
        }
        if(ft.dat$StrGeneA[p]=='-'){
          text7=c(0.05,0.74,c(InA,InB)[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==max(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))])
          text8=c(0.05,0.7,c(InA,InB)[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==min(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))])
        }
      }
      text9=c(0.05,0.1,paste0(p,'/',nrow(ft.dat)))
      text10=c(0.05,0.6,paste("5'",g5))
      text11=c(0.55,0.6,paste("5'",g3))
      text12=c(0.45,0.31,ft.dat$Direction[p])
      
      if (ft.dat[p, 'Tier'] == 'T1') {
        text.dat2=data.frame(rbind(text10, text11, text12), stringsAsFactors = F)
        text.dat2$X1=as.numeric(text.dat2$X1)
        text.dat2$X2=as.numeric(text.dat2$X2)
        plots[[p]]=ggplot() +
          geom_segment(data=gd,aes(x=st, xend=en, y=V1, yend=V1,colour=Domain), size=5) +
          geom_segment(data=arr.dat,aes(x=X1, y=X3, xend=X2, yend=X4),
                       arrow = arrow(length = unit(0.02, "npc"))) +
          geom_text(data=text.dat2,aes(x=X1, y=X2, label=X3),hjust=0, size=2)+
          ylim(c(0,1)) + 
          guides(col = guide_legend(ncol = 6, byrow = FALSE), 
                 shape = guide_legend(override.aes = list(size = 0.5)),
                 color = guide_legend(override.aes = list(size = 0.5))) +
          theme(axis.ticks=element_blank(),axis.title=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
                panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
                legend.title = element_text(size = 6), legend.text = element_text(size = 6),legend.position=c(0.5, 0.25),
                legend.key = element_blank(),legend.key.size = unit(0.02, "in"), legend.key.width = unit(0.01, "in"),
                legend.background =element_blank())
          
        num = leading(p, 4)
        png(paste0(output,'/',paste0(Sample_ID,'_output.T1.',num,'.png')), 
            width = 4, height=5, unit = "in", res=400)
        invisible(lapply(plots, print))
        dev.off()
      }

      # text.dat=data.frame(rbind(text0,text1,text2,text3,text4,text6,text7,text8,text9,text10,text11,text12),stringsAsFactors = F) #text5 delete
      # text.dat$X1=as.numeric(text.dat$X1)
      # text.dat$X2=as.numeric(text.dat$X2)
        
      # plots[[p]]=ggplot() +
      #   geom_segment(data=gd,aes(x=st, xend=en, y=V1, yend=V1,colour=Domain), size=5) +
      #   geom_segment(data=arr.dat,aes(x = X1, y = X3, xend = X2, yend = X4),
      #                arrow = arrow(length = unit(0.02, "npc"))) +
      #   geom_text(data=text.dat,aes(x=X1, y=X2, label=X3), hjust=0, size=2) +
      #   ylim(c(0,1)) +
      #   guides(col = guide_legend(ncol = 6, byrow = FALSE), 
      #          shape = guide_legend(override.aes = list(size = 0.5)),
      #          color = guide_legend(override.aes = list(size = 0.5))) +
      #   theme(axis.ticks=element_blank(),axis.title=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
      #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
      #         legend.title = element_text(size = 6), legend.text = element_text(size = 6),legend.position=c(0.5, 0.25),
      #         legend.key.size = unit(0.02, "in"), legend.key.width = unit(0.01, "in"),
      #         legend.key = element_blank(),legend.background =element_blank())

      # num = leading(p, 4)
      # png(paste0(output,'/',paste0(Sample_ID,'_output.',num,'.png')), 
      #     width = 4, height=5, unit = "in", res=400)
      # invisible(lapply(plots, print))
      # dev.off()      
    }
    
    # pdf(paste0(output,'/',paste0(Sample_ID,'_output.pdf')),paper="a4")
    # invisible(lapply(plots, print))
    # dev.off()
  }  
}else{
  print(paste("No result",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))
}

system(paste('rm -rf',paste0(output.bam,'*bed*')))
print(paste("Finish",format(Sys.time(),"%Y-%b-%d %H:%M:%S")))
