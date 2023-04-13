library(circlize)
#library(grid)
#(graphics)
library(ComplexHeatmap)
dir <- "E:/interview/scripts/GOcircos/"
rawdata <- "MvsF_ALL.GO_enrichment_result.xls"
setwd(dir)
raw.dt <- read.table(rawdata, head = T, sep = "\t")
plot.dt <- head(raw.dt[c("GO_accession", "Term_type", "Corrected_pValue", "DEG_item", "Bg_item")], 20)
plot.dt$Corrected_pValue <- -log10(plot.dt$Corrected_pValue)
# plot.dt$Corrected_pValue <- runif(20, min = 0, max = 25)
num.convert250 <- function(n){
  if (n <= 2){
    return(n)
  } else {
    return (2 + num.convert250((n - 2)/5))
  }
}

batch_pro250 <- function(n){
  for (i in 1:length(n)){
    if (n[i] > 2){
      n[i] = (n[i] - 2)/4
      n[i] = 2 + num.convert250(n[i])
    }
  }
  return (n)
}

# 输入一个数字
# 如果n > 10
#   tmp = (n - 10)/9
#   使用tmp进行回归运算 
#     如果 tmp <= 10 无需回归运算
#       输出 （tmp/10）
#     如果 tmp > 10 需要多层回归
#       输出 （10 + 回归((tmp -10)/10)）
# out = 10 + 回归输出

num.convert10000 <- function(n){
  if (n <= 10){
    return (n)
  } else {
    return(10 + num.convert10000(((n - 10)/10)))
  }
}
batch_pro10000 <- function(n){
  for (i in 1:length(n)){
    if (n[i] > 10){
      n[i] = (n[i] - 10)/9
      n[i] = 10 + num.convert10000(n[i])
    }
  }
  return (n)
}

axis_conf <- function(n, ...){
  if (max(n, ...) <= 10){
    sct_len <<- 10
    sct_mj <<- c(2, 4, 6, 8)
    sct_lab <<- c(2, 4, 6, 8)
    return (n)
  } else if (max(n, ...) <= 1250){
    sct_len <<- 10
    sct_mj <<- c(0, 2, 4, 6, 8, 10)
    sct_lab <<- c(0, 2, 10, 50, 250, NA)
    return (batch_pro250(n))
  } else if (max(n, ...) > 1250){
    sct_len <<- 45
    sct_mj <<- c(0, 10, 20, 30, 40, 45)
    sct_lab <<- c(0, 10, 100, 1000, 10000, NA)
    return (batch_pro10000(n))
  }
}

col.generator <- function(a, b){
  col = c(rep(b[1], length(a)))
  # col = factor(rep(b[1], length(a)), levels = b)
  for (i in 1:length(a)){
    col[i] = b[a[i]]
  }
  return (col)
}
qval_col <- c("#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#99000d")
names(qval_col) <- c("(0,1.3]", "(1.3,2]","(2,5]","(5,10]","(10,15]","(15,20]",">=20")
col.qValue <- function(a){
  col = c(rep(NA, length(a)))
  for (i in 1:length(a)){
    if (a[i] <= 1.3){
      col[i] = "#fee5d9"
    } else if (a[i] < 2){
      col[i] = "#fcbba1"
    } else if (a[i] < 5){
      col[i] = "#fc9272"
    } else if (a[i] < 10){
      col[i] = "#fb6a4a"
    } else if (a[i] < 15){
      col[i] = "#ef3b2c"
    } else if (a[i] < 20){
      col[i] = "#cb181d"
    } else if (a[i] >= 20){
      col[i] = "#99000d"
    }
  }
  return (col)
}
plot.dt$gene_num <- axis_conf(plot.dt$DEG_item, plot.dt$Bg_item)
plot.dt$all_num <- axis_conf(plot.dt$Bg_item)
plot.dt$rich_factor <- plot.dt$DEG_item/plot.dt$Bg_item
# plot.dt$rich_factor <- runif(10, min = 0, max = 1)
plot.dt$Bg_item <- as.character(plot.dt$Bg_item)
plot.dt$DEG_item <- as.character(plot.dt$DEG_item)
plot.dt <- plot.dt[order(plot.dt$Term_type, -plot.dt$rich_factor),]
plot.dt$start <- 0
plot.dt$end <- sct_len
plot_data <- plot.dt[c("GO_accession", "start", "end")]
type_col <- c(molecular_function = "#F7CB16", biological_process = "#F93887", cellular_component = "#954276")
qvalue_col <- col.qValue(plot.dt$Corrected_pValue)

color1 <- col.generator(plot.dt$Term_type, type_col)
pdf("GOcircos.pdf", width = 8, height = 8)
circos.par(gap.degree = 0.5, start.degree = 90)
circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1) 
circos.track(
  ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = color1,  
  panel.fun = function(x, y) {
    ylim = get.cell.meta.data('ycenter')  
    xlim = get.cell.meta.data('xcenter')
    sector.name = get.cell.meta.data('sector.index')  
    circos.axis(h = 'top', 
                major.at = sct_mj, labels = sct_lab, minor.ticks = 0, 
                labels.col = "black", labels.facing = "clockwise", labels.niceFacing = TRUE, 
                col = "grey", labels.cex = 0.6) 
    circos.text(xlim, ylim, sector.name, cex = 0.6, col = "black", facing = "bending", font = 2, niceFacing = TRUE)  
  })

plot_data <- plot.dt[c("GO_accession", "start", "all_num", "Bg_item")]

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08,  bg.border = NA, stack = TRUE,  
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, col = qvalue_col[get.cell.meta.data('sector.numeric.index')], border = NA, ...)  
    ylim = get.cell.meta.data('ycenter')  
    xlim = (region[[1]] + region[[2]])/2
    circos.text(xlim, ylim, value[[1]], cex = 0.4, facing = "bending", niceFacing = TRUE)  
  } )


plot_data <- plot.dt[c("GO_accession", "start", "gene_num", "DEG_item")]

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08,  bg.border = NA, stack = TRUE,  
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, col = "#69115d", border = NA, ...)  
    ylim = get.cell.meta.data('ycenter')  
    xlim = (region[[1]] + region[[2]])/2
    circos.text(xlim, ylim, value[[1]], cex = 0.4, col = "white", facing = "bending", niceFacing = TRUE)  
  } )

plot_data <- plot.dt[c("GO_accession", "start", "end", "rich_factor")]

circos.genomicTrack(
  plot_data, ylim = c(0, 1), track.height = 0.3, bg.col = 'gray95', bg.border = NA,  
  panel.fun = function(region, value, ...) {
    sector.name = get.cell.meta.data('sector.index')
    for (i in (1:9)*0.1){
      circos.lines(c(0, max(region)), c(i, i), col = 'gray', lwd = 0.3)
    }
    circos.genomicRect(region, value, col = color1[get.cell.meta.data('sector.numeric.index')], border = NA, ytop.column = 1, ybottom = 0, ...)
  } )

type_lgd <- Legend(title = "Term Type", title_gp = gpar(fontsize = 8, fontface = "bold"), 
                   at = names(type_col), labels_gp = gpar(fontsize = 6),
                   legend_gp = gpar(fill = type_col), 
             grid_height = unit(3, 'mm'), grid_width = unit(3, 'mm'), 
             row_gap = unit(0.8, "mm"))

gene_lgd <- Legend(title = "Gene Number",  title_gp = gpar(fontsize = 8, fontface = "bold"), 
                   at = c("Background Number", "Input Number"), labels_gp = gpar(fontsize = 6),
                   legend_gp = gpar(fill = c("#fb6a4a", "#69115d")), 
                   grid_height = unit(3, 'mm'), grid_width = unit(3, 'mm'), 
                   row_gap = unit(0.8, "mm"))

qval_lgd <- Legend(
  labels = names(qval_col), labels_gp = gpar(fontsize = 6),
  title = "-log10(Qvalue)", title_gp = gpar(fontsize = 8, fontface = "bold"),
  graphics = list(
    function(x, y, w, h) grid.points(x, y, gp = gpar(col = qval_col[1]), pch = 16),
    function(x, y, w, h) grid.points(x, y, gp = gpar(col = qval_col[2]), pch = 16),
    function(x, y, w, h) grid.points(x, y, gp = gpar(col = qval_col[3]), pch = 16),
    function(x, y, w, h) grid.points(x, y, gp = gpar(col = qval_col[4]), pch = 16),
    function(x, y, w, h) grid.points(x, y, gp = gpar(col = qval_col[5]), pch = 16),
    function(x, y, w, h) grid.points(x, y, gp = gpar(col = qval_col[6]), pch = 16),
    function(x, y, w, h) grid.points(x, y, gp = gpar(col = qval_col[7]), pch = 16)), 
  row_gap = unit(1, "mm"))

lgd_pkg <- packLegend(type_lgd, gene_lgd, qval_lgd, max_height = unit(50, "mm"))
draw(lgd_pkg)

dev.off()
circos.clear()

