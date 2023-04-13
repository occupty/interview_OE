setwd("E:/interview/scripts/多组通路矩阵点图/GO/")
type <- c("BP", "CC", "MF")
names(type) <- c("biological_process", "cellular_component", "molecular_function")
for (i in c("all", "up", "down")){
  term <- read.table(paste("reconsitution_", i, ".xls", sep = ""), head = T, sep = "\t")
  term <- head(term, 120)
  term$GO_Term <- factor(term$GO_Term, levels = rev(unique(term$GO_Term)))
  term$Term_Type <- type[term$Term_Type]
  library(ggplot2)
  ggplot(term, aes(x = Compare, y = GO_Term)) + 
    geom_point(aes(color = Pvalue, size = GeneNumber)) + 
    facet_grid(Term_Type~., scales = "free", space = "free") + 
    scale_color_gradient2(midpoint = mean(c(min(term$Pvalue), max(term$Pvalue))), 
                          high = "#4172ad", mid = "#f4eeba", low = "#d63228") +
    # scale_size_continuous(breaks = c(100, 1000, 4000, 8000)) +
    xlab("Group") + ylab("GO  Term") + 
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          strip.background = element_rect(fill = NA, color = "black"),
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
  ggsave(paste("GO_", i, ".pdf", sep = ""), height = 8, width = ifelse(i == "down", 8.5, 10))
}
