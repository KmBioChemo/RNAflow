suppressWarnings(suppressMessages({ library(pheatmap); library(matrixStats); library(RColorBrewer); library(grid) }))
OUT <- "paper/figure3_rebuild/panels"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
clean_term <- function(x){ x<-sub("^HALLMARK_","",x); gsub("_"," ",x) }
D <- readRDS("paper/.panel_cache.rds")
render_hm <- function(n_top,w,h,tag){
  m<-as.matrix(D$gv); v<-matrixStats::rowVars(m,na.rm=TRUE)
  keep<-order(v,decreasing=TRUE)[seq_len(min(n_top,nrow(m)))]; m<-m[keep,,drop=FALSE]; rownames(m)<-clean_term(rownames(m))
  md<-D$tm[match(colnames(m),D$tm[[colnames(D$tm)[1]]]),,drop=FALSE]
  ann<-data.frame(row.names=colnames(m), cancer_type=as.character(md[["cancer_type"]]), check.names=FALSE)
  png(file.path(OUT,paste0(tag,".png")),width=w,height=h,units="in",res=300,bg="white")
  pheatmap(m,scale="row",annotation_col=ann,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=FALSE,
    color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),border_color=NA,fontsize=9,main="")
  dev.off(); cat("wrote",tag,w,"x",h,"\n") }
render_hm(25,15.0,4.2,"a_gsva_b25")
