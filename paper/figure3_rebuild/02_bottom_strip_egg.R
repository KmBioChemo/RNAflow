suppressWarnings(suppressMessages({ library(ggplot2); library(egg); library(grid) }))
OUT <- "paper/figure3_rebuild/panels"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FONT <- "Liberation Sans"
`%||%` <- function(a,b) if(is.null(a)) b else a
clean_term <- function(x){ x<-sub("^HALLMARK_","",x); gsub("_"," ",x) }
wrap_label <- function(x, width=42) vapply(x, function(s) paste(strwrap(s,width=width),collapse="\n"), character(1), USE.NAMES=FALSE)
OI <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#7F7F7F")
source("R/fig_theme.R"); source("R/fig_palettes.R"); source("R/fig_wgcna.R"); source("R/module_enrichment.R")
theme_bare <- function(b=12) theme_publication() + theme(
  text=element_text(family=FONT), plot.title=element_blank(),
  axis.title=element_text(size=b), axis.text=element_text(size=b-2.5),
  legend.title=element_text(size=b-2,face="bold"), legend.text=element_text(size=b-3),
  legend.position="bottom", legend.key.size=unit(10,"pt"),
  legend.margin=margin(0,0,0,0), legend.box.spacing=unit(2,"pt"), plot.margin=margin(3,6,3,6))
D <- readRDS("paper/.panel_cache.rds")
sc <- merge(D$pca$scores, D$tm, by="sample")
b <- ggplot(sc, aes(PC1,PC2,colour=cancer_type)) + geom_point(size=2.2,alpha=.92,stroke=0) +
  scale_colour_manual(values=OI,name="Cancer type") + guides(colour=guide_legend(nrow=2,override.aes=list(size=2.4))) +
  labs(x=sprintf("PC1 (%.1f%%)",D$pca$pct[1]),y=sprintf("PC2 (%.1f%%)",D$pca$pct[2])) + theme_bare()
c <- fig_soft_threshold(D$sft, mode="publication")
c$layers <- Filter(function(L) !inherits(L$geom,"GeomText"), c$layers)
c$facet$params$labeller <- as_labeller(c("Scale-free fit (signed R2)"="Scale-free fit R²","Mean connectivity"="Connectivity"))
c <- c + theme_bare() + theme(strip.text=element_text(size=8))
d <- fig_module_trait(D$mt, mode="publication", text_size=2.6) +
  scale_fill_gradient2(low="#2980B9",mid="white",high="#C0392B",midpoint=0,limits=c(-1,1),
    breaks=c(-1,-0.5,0,0.5,1), name="Correlation") + theme_bare() +
  guides(fill=guide_colourbar(barwidth=unit(2.6,"cm"),barheight=unit(0.3,"cm"),
    title.position="top",title.hjust=0.5,ticks.colour="grey40")) +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8),axis.text.y=element_text(size=8),
    legend.text=element_text(size=8),legend.title=element_text(size=9,face="bold"),
    plot.caption=element_text(size=8.5,colour="grey25",hjust=0.5,margin=margin(t=3)))
e <- fig_module_enrichment(D$mod_enrich, max_terms=12, mode="publication") + theme_bare(11) +
  theme(legend.box="vertical", legend.box.just="left", legend.spacing.y=unit(1,"pt")) +
  guides(size=guide_legend(title="Gene count",order=1), colour=guide_colourbar(title="-log10 FDR",order=2,barwidth=unit(2.4,"cm"),barheight=unit(0.28,"cm")))
png(file.path(OUT,"fig3_strip.png"), width=15, height=4.2, units="in", res=300, bg="white")
egg::ggarrange(b,c,d,e, nrow=1, widths=c(1,1,1.15,1), labels=c("b","c","d","e"),
  label.args=list(gp=grid::gpar(fontface="bold", fontsize=16), hjust=-0.2, vjust=1.2))
invisible(dev.off()); cat("wrote fig3_strip (egg, letters only)\n")
