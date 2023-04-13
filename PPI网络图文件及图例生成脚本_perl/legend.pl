use Cwd;
$colordef = shift;
open DEF, "<$colordef";
$colordef =~ s/_legendColor.txt//;
while (<DEF>){
	chomp;
	@line = split /\t/, $_;
	$set{$line[0]} = $line[1];
	$setstr .= "'".$line[0]."'='".$line[1]."', ";
}
$setstr =~ s/, $//;
$pwd = getcwd();
#print $pwd;
$Rfile =<<END;
setwd("$pwd");
library(ComplexHeatmap)
nodeColor <- c($setstr);
if (length(nodeColor) > 20){ncol = 2}else{ncol = 1}
lgd_node <- Legend(title = "Gene", ncol = 2, labels = names(nodeColor), legend_gp = gpar(fill = nodeColor))
#lgd_edge <- Legend(labels = c('positive', 'negative'), title = "Correlation", type = "lines", legend_gp = gpar(col = c('red', 'blue')), grid_width = unit(1, "cm"), grid_height = unit(2, "mm"))
#lgd_pkg <- packLegend(lgd_node, lgd_edge, max_height=unit(50, "mm"))
pdf("${colordef}_Legend.pdf")
draw(lgd_node)
dev.off()
END

open R, "|/allwegene8/Individual_software/R/R-4.0.3/install/bin/R --vanilla --slave";
print R $Rfile;
