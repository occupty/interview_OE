setwd("E:/interview/scripts/多组通路矩阵点图/KEGG/")
library(ggplot2)
for (i in c("all", "up", "down")){
  term <- read.table(paste("reconsitution_", i, ".xls", sep = ""), head = T, sep = "\t")
  term <- head(term, 120)
  term$KEGG_Pathway <- factor(term$KEGG_Pathway, levels = rev(unique(term$KEGG_Pathway)))
  term$GeneNumber = as.numeric(ifelse(term$GeneNumber == "-", NA, term$GeneNumbe))
  term$Pvalue = as.numeric(ifelse(term$Pvalue == "-", max(term$Pvalue), term$Pvalue))
  ggplot(term, aes(x = Compare, y = KEGG_Pathway)) + 
    geom_point(aes(color = Pvalue, size = GeneNumber)) + 
    scale_color_gradient2(midpoint = mean(c(min(term$Pvalue), max(term$Pvalue))), 
                          high = "#4172ad", mid = "#f4eeba", low = "#d63228") +
    xlab("Group") + ylab("KEGG Pathway") + 
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank(),
          legend.key = element_rect(fill = NA),
          axis.line.x = element_line(colour = "black"),
          axis.ticks.x = element_line(colour = "black"),
          axis.line.y  = element_line(colour = "black"),
          axis.line.y.right = element_line(colour = "black"),
          axis.ticks.y = element_line(colour = "black"),
          axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(colour = "black"),
          ) 
  ggsave(paste("KEGG_", i, ".pdf", sep = ""), 
         height = ifelse(i == "all", 6, 8), width = ifelse(i == "all", 6, 8))
}
