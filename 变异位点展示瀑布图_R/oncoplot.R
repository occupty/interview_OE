library("ComplexHeatmap")
setwd("E:/interview/scripts/瀑布图/")
snp <- read.table("SNP_genelist.txt", sep = "\t", header = T, row.names = 1)
indel <- read.table("Indel_genelist.txt", sep = "\t", header = T, row.names = 1)
data <- matrix(rep("0",44*19), nrow = 44)
dt <- data.frame(data)
for (i in seq(1,44)){
    for (j in seq(1,19)){
        dt[i,j] = paste(snp[i,j],indel[i,j],sep = ";")
    }
}
for (i in seq(1,44)){
     for (j in seq(1,19)){
          if (dt[i,j] == ";"){dt[i,j] = ""}
  }
}
for (i in seq(1,44)){
    for (j in seq(1,19)){
        if (dt[i,j] == ";indel"){dt[i,j] = "indel;"}
    }
}
rownames(dt) = rownames(snp)
colnames(dt) = colnames(snp)
col = c("missense_SNV" = "#836fff", "UTR5" = "#ff0000", "synonymous_SNV" = "#228b22", "UTR3" = "#00bfff", "splicing" = "#ffa500", "indel" = "#000000")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "white", col = NA))
  },
  missense_SNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["missense_SNV"], col = NA))
  },
  UTR5 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["UTR5"], col = NA))
  },
  synonymous_SNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["synonymous_SNV"], col = NA))
  },
  UTR3 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["UTR3"], col = NA))
  },
  splicing = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["splicing"], col = NA))
  },
  indel = function(x, y, w, h) {
    x1 = x - w/2 + unit(0.25, "mm")
    x2 = x + w/2 - unit(0.25, "mm")
    x3 = x + w/2 - unit(0.25, "mm")
    y1 = y + h/2 - unit(0.25, "mm")
    y2 = y + h/2 - unit(0.25, "mm")
    y3 = y - h/2 + unit(0.25, "mm")
    grid.polygon(unit.c(x1, x2, x3), unit.c(y1, y2, y3), id = NULL,
                 gp = gpar(fill = col["indel"], col = NA))
  }
)
heatmap_legend_param = list(title = "test_indel", at = c("missense_SNV", "UTR5", "synonymous_SNV", "UTR3", "splicing", "indel"), 
                            labels = c("missense SNV", "UTR5", "synonymous SNV", "UTR3", "splicing", "Indel"))
pdf("test.pdf", height = 7, width = 10)
oncoPrint(dt, 
          alter_fun = alter_fun, col = col, 
          column_title = "TITLE", heatmap_legend_param = heatmap_legend_param, show_column_names = TRUE,
          left_annotation =  rowAnnotation(
            rbar = anno_oncoprint_barplot(
              axis_param = list(direction = "reverse")
            )),
          right_annotation = NULL)
dev.off()
