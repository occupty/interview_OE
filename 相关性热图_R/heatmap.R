# 三个参数：fpkm数据，代谢物表达量数据，输出文件名
args <- commandArgs(TRUE)
library("psych")
library("pheatmap")
library("tidyverse")
fpkm <- read.table(args[1], header = T, row.names = 1, sep = "\t") # 导入fpkm数据，数据格式：gene_id\tsample1\tsample2\t ....
compd <- read.table(args[2], header = T, row.names = 1, sep = "\t") # 导入代谢物数据，数据格式：compound\tsample1\tsample2\t ....
fpkm.t <- t(fpkm) # 转置用于corr.test函数，corr.test()识别每列为一条信息
compd.t <- t(compd) # 同上
cortest_psy <- corr.test(compd.t, fpkm.t, method = "spearman") # corr.test()相关性分析：method 可以选择“pearson”, “kendall”, “spearman”三种
pheatmap(cortest_psy$r,number_color = "black", # 在每个cell中填写内容的颜色设置
	color = colorRampPalette(c("#4a3896","white","#d61f26"))(256), # 使用colorRamPalette()创建热图所需的颜色参数
	treeheight_col = 50, treeheight_row = 50, # 热图中聚类树高度设置，默认50
	cellwidth = 18, cellheight = 12, # 热图每个cell长宽设置
	display_numbers = matrix(ifelse(cortest_psy$p<0.010,"**",ifelse(cortest_psy$p>0.010 & cortest_psy$p<0.050,"*","")), nrow(cortest_psy$p)), # 在每个cell中填写内容的设置，0.01<p<0.05为"*", p<0.01为"**"
	scale = "row", # 按行进行标准化，可选参数'none', 'row' or 'column'
	border_color = NA, # 去除cell分割线
	filename = paste0(args[3], ".pdf")) # 设置文件名
psy_r <- data.frame(Compound = rownames(cortest_psy$r), cortest_psy$r)
psy_p <- data.frame(Compound = rownames(cortest_psy$p), cortest_psy$p)
psy_r <- gather(psy_r, key = "Gene", value = "R", -Compound);
psy_p <- gather(psy_p, key = "Gene", value = "Pvalue", -Compound);
psy_out <- merge(psy_r, psy_p, by = c("Compound", "Gene"))
write.table(psy_out,paste0(args[3],"_R_P.xls"), sep="\t", quote = F, row.names = F)
