$genedes = "/allwegene3/database/ref/Ensembl_animals_release100/mus_musculus/ensemble_bioMart/gene_description.txt";
$cmp = shift;
$ppi = $cmp.".ppi";
$gene2ptn = $cmp."_unique.xls";
$colordef = "colordef";
open PI, "<$ppi";
<PI>;
while (<PI>){
    chomp;
    @line = split /\t/, $_; 
    ${$edge{$line[0]}}{$line[1]} = $line[2];
}
foreach $center (keys %edge){
	$length{$center} = keys %{$edge{$center}};
}
open GP, "<$gene2ptn";
while (<GP>){
    chomp;
    @line = split /\t/, $_; 
	push @{$ptn{$line[1]}}, $line[0];
}

open GD, "<$genedes";
<GD>;
while (<GD>){
    chomp;
    @line = split /\t/, $_; 
    $des{$line[0]} = $line[1];
}

foreach $protein (keys %ptn){
	@{$ptn2name{$protein}} = map {$des{$_} ? $des{$_} : $_} @{$ptn{$protein}};
}

open CL, "<$colordef";
open CLOUT, ">${cmp}_legendColor.txt";

foreach $center (sort {$length{$b} <=> $length{$a}} keys %length) {
	#unless (defined @{$ptn{$center}}){print $center." not defined."}
	#print join("|", @{$ptn{$center}})."\t$length{$center}\n";
	$colorline = <CL>;
	chomp $colorline;
	@{$top10ct{$center}} = split /\t/, $colorline;
	@{${$top10ct{$center}}[1]} = split /,/, ${$top10ct{$center}}[1];
	print CLOUT join("|", @{$ptn2name{$center}})."\t${$top10ct{$center}}[0]\n";
	last if (++$n == 10);
}
close CLOUT;

system("perl ./legend.pl ${cmp}_legendColor.txt");

sub color_merge{
    $old_color = $_[0];
    $new_color = $_[1];
    @out_color = map {int(($$old_color[$_] + $$new_color[$_])/2)} (0,1,2);
    return @out_color;
}

$edgedef = "edgedef>node1 VARCHAR,node2 VARCHAR,weight DOUBLE\n";
$nodedef = "nodedef>name VARCHAR,label VARCHAR,color VARCHAR,visible BOOLEAN,labelvisible BOOLEAN\n";
foreach $ct (keys %top10ct){
    $nodedef .= $ct.",'".join("|", @{$ptn2name{$ct}})."','".join(",", @{${$top10ct{$ct}}[1]})."',true,true\n";
    foreach $nd (keys %{$edge{$ct}}){
        $edgedef .= "$ct,$nd,${$edge{$ct}}{$nd}\n";
		#$nodedef .= $nd.",,true,false\n" unless (exists $edge{$nd})
		next if (exists $top10ct{$nd});
        @{$color{$nd}} = ($color{$nd}) ? &color_merge($color{$nd}, ${$top10ct{$ct}}[1]) : @{${$top10ct{$ct}}[1]};
    }   
}
foreach $nd (keys %color){
	$nodedef .= $nd.",,'".join(",", @{$color{$nd}})."',true,false\n";
}

open EDGE, ">${cmp}.gdf";
print EDGE $nodedef.$edgedef;
