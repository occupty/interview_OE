setwd("E:/interview/scripts/GO气泡组合图/")
go_raw <- read.table("SUvsCK_ALL.GO_enrichment_result.xls", head = T, sep = "\t")
diff_raw <- read.table("SUvsCK.Differential_analysis_results.xls", head = T, sep = "\t")
go_dt <- head(go_raw[c(1,2,3,5,10)],30)
diff_dt <- diff_raw[c(1,4)]
colnames(go_dt)<-c( 'ID', 'term','category','adj_pval','genes')
colnames(diff_dt)<-c("ID","logFC")
# pdf("GOBubble.pdf", height = 10, width = 15)
# GOBubble(circ, table.legend = T, table.col = T)
# dev.off()
# pdf("GOBubble_multiple.pdf", height = 10, width = 40)
# GOBubble(circ, display = "multiple")
# dev.off()
library(ggplot2)
library(ggrepel)
go_dt$count <- head(go_raw$DEG_item, 30)
bub_col <- c("#9cbf01", "#74abf9", "#efb134")
names(bub_col) <- c("cellular_component", "biological_process", "molecular_function")
go_updown <- read.table("SUvsCK.GO_enrichment_result_up_down.xls", head = T, sep = "\t")
go_updown$zscore <- (go_updown$Up - go_updown$Down)/sqrt(go_updown$DEG_item)
go_dt$zscore <- head(go_updown$zscore, 30)
go_dt$adj_pval <- -log10(go_dt$adj_pval)

multi_bubble <- ggplot(go_dt) + geom_point(aes(zscore, adj_pval, colour = category, size = count), alpha = 0.5) +
  scale_size_continuous(range = c(5,30)) + 
  geom_text_repel(aes(zscore, adj_pval), label = go_dt$ID) + 
  theme(legend.title = element_text(face = "bold"),
        legend.key=element_rect(fill= NA ,color= NA ),
        legend.direction="horizontal",
        legend.position="bottom") + guides(size=FALSE) + 
  theme(panel.background=element_rect(fill="white",color="grey50"),
        panel.grid.major=element_line(size=0.1,linetype =1,color="grey70"),
        panel.grid.minor=element_line(size=0.1,linetype =1,color="grey70")) + 
  theme(axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold")) + 
  xlab(label = "z-score") + ylab(label = "-log(adj P-value)") +
  scale_color_manual(values = bub_col)

library(ggpubr)
go_tab <- go_dt[c(1,2)]
color_c <- c(rep(NA,30))
for (i in 1:30){
  color_c[i] = bub_col[go_dt[i,3]]
}
table_out <- ggtexttable(go_tab, rows = NULL,
            theme = ttheme(
              colnames.style = colnames_style(color = "black", fill = "grey70", face = "bold"),
              tbody.style = tbody_style(color = "black", fill = color_c)
            )
)

p1 <- ggarrange(multi_bubble, table_out, nrow = 1)

facet_bubble <- ggplot(go_dt) + geom_point(aes(zscore, adj_pval, colour = category, size = count), alpha = 0.5) +
  scale_size_continuous(range = c(5,30)) + 
  geom_text_repel(aes(zscore, adj_pval), label = go_dt$ID) + 
  theme(legend.title = element_text(face = "bold"),
        legend.key=element_rect(fill= NA ,color= NA ),
        legend.direction="horizontal",
        legend.position="bottom") + guides(size=FALSE) + 
  theme(panel.background=element_rect(fill="white",color="grey50"),
        panel.grid.major=element_line(size=0.1,linetype =1,color="grey70"),
        panel.grid.minor=element_line(size=0.1,linetype =1,color="grey70")) + 
  theme(axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold")) + 
  xlab(label = "z-score") + ylab(label = "-log(adj P-value)") +
  scale_color_manual(values = bub_col) + 
  facet_wrap(~category, scales = "free")

p2 <- ggarrange(p1, facet_bubble, ncol = 1, labels = c("a", "b"))

pdf("GO_Bubble.pdf", width = 20, height = 20)
p2
dev.off()
# write.table(circ, "GO_Bubble.xls", row.names = F, sep = "\t")
