#!/usr/bin/perl -w
# Copyright (C) 2016  Kai Tang <tangkai.ustc@gmail.com>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


use strict;
use File::Spec;
use Carp;
use Text::NSP::Measures::2D::Fisher2::twotailed;

my $debug = 0;
my $delete_tmp = 1;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $bin_size     = 200;
my $sliding_size = 50;
my $allowed_gap  = 100;

my $covered_num = int ($bin_size / $sliding_size);

my $DMC_cutoff  = 10;
my $depth_cutoff = 4;# 4;
my $p_value_cutoff = 0.05;
my %chr_len = ("chr1"=>30427671, "chr2"=>19698289, "chr3"=>23459830, "chr4"=>18585056, "chr5"=>26975502);

my $R_file_name = "p_adjust_BH_v0.3.R";
my ($volume,$dir_base,$perl_file) =  File::Spec->splitpath( $0 );
my $R_script = File::Spec-> catfile($dir_base, $R_file_name);

croak($R_file_name . " should be in the same directory with step2_DMR.pl") unless (-e $R_script);


my %bin_last_index;
# bin start from 0
for my $chr (sort keys %chr_len){
	my $len = $chr_len{$chr};
	my $last = int( ($len - 1) / $sliding_size);
	$bin_last_index{$chr} = $last;
}
if($debug){
	foreach my $chr (sort keys %bin_last_index){
		print STDERR join("\t", ($chr, $bin_last_index{$chr})), "\n";
	}
	print STDERR "\n\n";
	#exit;
}

my $usage = "$0 \n <WT_isMeth> <mut_isMeth> <outdir> <outpre>\n\n";
croak $usage unless (@ARGV == 4);

my $wt_in  = shift or croak;
my $mut_in = shift or croak;
my $outdir = shift or croak;
my $outpre = shift or croak;

croak unless (-d $outdir);
croak $wt_in  unless (-e $wt_in);
croak $mut_in unless (-e $mut_in);

my $output_hyper = File::Spec->catfile($outdir, $outpre."_hyper_list.txt");
my $output_hypo  = File::Spec->catfile($outdir, $outpre."_hypo_list.txt");

croak "\n output_hyper exits \n\n" if (-e $output_hyper);
croak "\n output_hypo exits \n\n" if (-e $output_hypo);


my $temp_p_file =File::Spec->catfile($outdir, $outpre."_BinSize". $bin_size . "_sliding" . $sliding_size . "_BH_adj_p_value". $p_value_cutoff . "_depth". $depth_cutoff. "_p.txt"); 

my $temp_p_adj_file =File::Spec->catfile($outdir, $outpre."_BinSize". $bin_size . "_sliding" . $sliding_size . "_BH_adj_p_value". $p_value_cutoff . "_depth". $depth_cutoff. "_p_adj.txt"); 

croak "\n\n$temp_p_file exists!!\n\n\n" if(-e $temp_p_file);
croak "\n\n$temp_p_adj_file exists!!\n\n\n" if(-e $temp_p_adj_file);

if ($debug){
	
	print STDERR "output:\n";
	print STDERR join("\n", ($output_hyper, $output_hypo)), "\n\n";
	
	print STDERR "\n\nOK\n\n";
	exit;
}

my (%CG_num, %CHG_num, %CHH_num); # number of CXX that meet the requirement that both depth >= $depth_cutoff
my (%CG_depth, %CHG_depth, %CHH_depth);#depth sum
my (%mCG, %mCHG, %mCHH);               #mCXX sum


open (WT,  $wt_in )  or croak "cannot open  $wt_in: $!";
open (MUT, $mut_in ) or croak "cannot open $mut_in: $!";


croak if(-e $temp_p_file);
open (P, ">$temp_p_file ") or croak "cannot open $temp_p_file: $!";
croak if(-e $temp_p_adj_file);
print P join("\t", ("chr", "start", "end", "p_CG","CG_hyper",  "p_CHG", "CHG_hyper", "p_CHH", "CHH_hyper")), "\n";


my (%dep_num_wt, %dep_num_mut);
my (%mC_num_wt,  %mC_num_mut);

my (%hyper_DMCs, %hypo_DMCs);

my @raw_list_hyper;
my @raw_list_hypo;


my ($l_wt, $l_mut);

print STDERR "\n\nstart read isMeth files...\t";

my $h_wt  = <WT>;
my $h_mut = <MUT>;

#my $line = 0;

while($l_wt = <WT>, $l_mut = <MUT>){
	chomp $l_wt;
	chomp $l_mut;
	my @a_wt  = split "\t", $l_wt;
	my @a_mut = split "\t", $l_mut;
	
	croak unless ($a_wt[1] == $a_mut[1]);
	my $chr = $a_wt[0];
	my $pos = $a_wt[1];
	my $type = $a_wt[3];

	
	my $dep_wt = $a_wt[5];
	my $dep_mut = $a_mut[5];
	
	my $mC_wt = $a_wt[4];
	my $mC_mut = $a_mut[4];
	
	my $per_wt  = $a_wt[6];
	my $per_mut = $a_mut[6];
	

	if($per_wt == 0){
		if($mC_mut >= 2){
			$hyper_DMCs{$chr}->[$pos] = $type;
		}
	}elsif($per_mut / $per_wt >= 2){
		$hyper_DMCs{$chr}->[$pos] = $type;
	}
	
	if($per_mut == 0){
		if($mC_wt >= 2){
			$hypo_DMCs{$chr}->[$pos] = $type;
		}
	}elsif($per_wt / $per_mut >= 2){
		$hypo_DMCs{$chr}->[$pos] = $type;
	}

		
	my $base_index = int( ($pos - 1) / $sliding_size);
		
	for my $i (($base_index - $covered_num)..($base_index + $covered_num)){
		if($i >= 0){
			my $interval_start = $i * $sliding_size + 1;
			my $interval_end   = $interval_start + $bin_size - 1;
			if($pos >= $interval_start and $pos <= $interval_end){
				#$ref_den->{$chr}->[$i]++;
				$dep_num_wt{$type}->{$chr}->[$i]  += $dep_wt ;
				$dep_num_mut{$type}->{$chr}->[$i] += $dep_mut;
				
				$mC_num_wt {$type}->{$chr}->[$i] += $mC_wt;
				$mC_num_mut{$type}->{$chr}->[$i] += $mC_mut;
			#	$ref_C_num {$type}->{$chr}->[$i] ++;
			}
		}
	}
}

close(WT);
close(MUT);

print STDERR "Done\n";

my @types = ("CG", "CHG", "CHH");

print STDERR "start cal p-value...\t";


foreach my $chr(sort keys %chr_len){
	my $last = $bin_last_index{$chr};
	
	for (my $i = 0; $i <= $bin_last_index{$chr}; $i++){
		
		my $start = $i * $sliding_size + 1;;
		my $end   = $start + $bin_size - 1;
		
		my ($total_C_num, $total_mC_wt, $total_dep_wt) = (0) x 3;
		my ( $total_mC_mut, $total_dep_mut ) = (0) x 2;
 		my @p_array = ();
 		
		foreach my $type (@types){
		#	my $ref_num = 0; 
			my ($mC_wt , $dep_wt ) = (0) x 2;
			my ($mC_mut, $dep_mut) = (0) x 2;
			if(defined $dep_num_wt{$type}->{$chr}->[$i]) { $dep_wt = $dep_num_wt{$type}->{$chr}->[$i] }
			if(defined $dep_num_mut{$type}->{$chr}->[$i]){ $dep_mut = $dep_num_mut{$type}->{$chr}->[$i] }
			if(defined $mC_num_wt{$type}->{$chr}->[$i])  { $mC_wt = $mC_num_wt{$type}->{$chr}->[$i] }
			if(defined $mC_num_mut{$type}->{$chr}->[$i]) { $mC_mut = $mC_num_mut{$type}->{$chr}->[$i] }
			
			my $n11 = $mC_wt;
			my $n1p = $mC_wt + $mC_mut;
			my $np1 = $dep_wt;
			my $npp = $dep_wt + $dep_mut;
			my $p_val = 1;
			$p_val = cal_Fisher_exact_test($n11, $n1p, $np1, $npp);
			if($p_val eq "error"){
				$p_val = 1;
			}
			
			my $hyper_label = "n";
			if($p_val != 1){
				my $wt_per  = $mC_wt / $dep_wt;
				my $mut_per = $mC_mut / $dep_mut;
				if($mut_per > $wt_per) {$hyper_label = "h"}
				if($mut_per < $wt_per) {$hyper_label = "l"}
			}
			
			push @p_array, $p_val, $hyper_label;
			
		}
		print P join ("\t", ($chr, $start, $end, @p_array)), "\n";
	}	
}
close(P);

%dep_num_wt  = ();
%dep_num_mut = ();
%mC_num_wt = ();
%mC_num_mut = ();

print STDERR "Done\n";

#####################
#	R cmd
#####################
#	    1   2        3        4     5    6
my $R_cmd = "R --slave --vanilla --args $temp_p_file $temp_p_adj_file < $R_script ";

print STDERR $R_cmd, "\n";
`$R_cmd`;

croak unless (-e $temp_p_adj_file);

#print STDERR "Done\n\n";


#my ($last_chr, $last_start, $last_end ) = ("chr0", 0, 0);

print STDERR "read adj p-value...\t";


open(ADJ, $temp_p_adj_file) or croak "cannot open $temp_p_adj_file";
my $h = <ADJ>;
while(<ADJ>){
	chomp;
	my @a = split "\t";
	my ($chr, $start, $end) = @a[0..2];
	my ( $p_CG, $label_CG,  $p_CHG, $label_CHG,  $p_CHH , $label_CHH ) = @a[3..8];
	if(  $p_CG <= $p_value_cutoff or  $p_CHG <= $p_value_cutoff or $p_CHH <= $p_value_cutoff ){
		
		my $flag_hyper = 0;
		
		if($p_CG <= $p_value_cutoff and $label_CG eq "h"){
			$flag_hyper = 1;
		}
		if($p_CHG <= $p_value_cutoff and $label_CHG eq "h"){
			$flag_hyper = 1;
		}
		if($p_CHH <= $p_value_cutoff and $label_CHH eq "h"){
			$flag_hyper = 1;
		}
		if($flag_hyper){
			push @raw_list_hyper, join("_", ( $chr, $start, $end )) ;
		}
		my $flag_hypo = 0;
		
		if($p_CG <= $p_value_cutoff and $label_CG eq "l"){
			$flag_hypo = 1;
		}
		if($p_CHG <= $p_value_cutoff and $label_CHG eq "l"){
			$flag_hypo = 1;
		}
		if($p_CHH <= $p_value_cutoff and $label_CHH eq "l"){
			$flag_hypo = 1;
		}
		if($flag_hypo){
			push @raw_list_hypo, join("_", ( $chr, $start, $end )) ;
		}
		
	} 
}
close(ADJ);

print STDERR "Done\n\n";




##################
# hyper
##################
my @list_hyper_merged;
my @hyper_list_border;
my %region_in_hyper;
my (%mC_list_hyper, %dep_list_hyper, %percentage_list_hyper );
my @final_list_hyper;


my %C_number_hyper;

print STDERR "handle hyper list...\n";

#print STDERR "merge_raw_list...\n";
merge_list(\@raw_list_hyper, \@list_hyper_merged, $allowed_gap);

print STDERR "adjust border and filter DMC...\n";
adjust_list_border(\@list_hyper_merged, \@hyper_list_border, \%hyper_DMCs, $DMC_cutoff);

%hyper_DMCs = ();

record_region_pos(\@hyper_list_border, \%region_in_hyper);


print STDERR "read_isMeth_get_mC_dep\n";

#hyper_list
#read_isMeth_get_mC_dep($mut_in, $wt_in, \%region_in_hyper, \%mC_list_hyper, \%dep_list_hyper, \%percentage_list_hyper);
read_isMeth_get_mC_dep($mut_in, $wt_in, \%region_in_hyper, \%mC_list_hyper, \%dep_list_hyper, \%percentage_list_hyper, \%C_number_hyper);
print STDERR "Done\n";

print STDERR "begin_last_filter\n";

#last_filter(\@hyper_list_border, \%mC_list_hyper, \%dep_list_hyper, \%percentage_list_hyper, \@final_list_hyper);
last_filter(\@hyper_list_border, \%mC_list_hyper, \%dep_list_hyper, \%percentage_list_hyper, \@final_list_hyper, \%C_number_hyper);

print STDERR "output...\n";
#$output_hyper
output_hyper_list($output_hyper, \@final_list_hyper);

print STDERR "Hyper Done\n\n\n";


@list_hyper_merged = ();
@hyper_list_border = ();
%region_in_hyper   = ();
%mC_list_hyper = ();
%dep_list_hyper  = ();
%percentage_list_hyper = ();
@final_list_hyper = ();
%C_number_hyper = ();


##################
# hypo
##################

my @list_hypo_merged;	
my @hypo_list_border;
my %region_in_hypo;
my (%mC_list_hypo, %dep_list_hypo, %percentage_list_hypo);
my @final_list_hypo;

my %C_number_hypo ;

print STDERR "handle hypo list...\n";
merge_list(\@raw_list_hypo,  \@list_hypo_merged, $allowed_gap);

print STDERR "adjust border and filter DMC...\n";
adjust_list_border(\@list_hypo_merged, \@hypo_list_border, \%hypo_DMCs, $DMC_cutoff);

%hypo_DMCs = ();

record_region_pos(\@hypo_list_border, \%region_in_hypo);


#hypo_list
print STDERR "read_isMeth_get_mC_dep\n";
#read_isMeth_get_mC_dep($wt_in, $mut_in, \%region_in_hypo, \%mC_list_hypo, \%dep_list_hypo, \%percentage_list_hypo);
read_isMeth_get_mC_dep($wt_in, $mut_in, \%region_in_hypo, \%mC_list_hypo, \%dep_list_hypo, \%percentage_list_hypo, \%C_number_hypo);

print STDERR "Done\n";

print STDERR "begin_last_filter\n";
#last_filter(\@hypo_list_border, \%mC_list_hypo, \%dep_list_hypo, \%percentage_list_hypo, \@final_list_hypo);
last_filter(\@hypo_list_border, \%mC_list_hypo, \%dep_list_hypo, \%percentage_list_hypo, \@final_list_hypo , \%C_number_hypo);
print STDERR "output...\n";
output_hypo_list ($output_hypo, \@final_list_hypo);

if($delete_tmp){
	my $cmd_delete = "rm -f $temp_p_file $temp_p_adj_file";
}

print STDERR "Done\n\n";

exit;


sub output_hyper_list{
	my ($file, $ref) = @_;
	croak if(-e $file);
	
	open(OUT, ">>$file") or croak;
	
				
	print OUT join("\t", ( "chr", "start", "end")), "\n";
				
	my $last_index = scalar(@{$ref}) - 1;
	
	for my $i(0..$last_index){
		my @tmp_out =  @{$ref->[$i]};

		print OUT join("\t", @tmp_out[0..2]), "\n";
	}
	
	close(OUT);	
}


sub output_hypo_list{
	my ($file, $ref) = @_;
	croak if(-e $file);
	
	open(OUT, ">>$file") or croak;
	print OUT join("\t", ( "chr", "start", "end", "DMC_num", "DMC_type", "CG_num", "CHG_num", "CHH_num", 
							"wt_CG", "wt_CG_per",	"wt_CHG", "wt_CHG_per","wt_CHH", "wt_CHH_per",
						   "mut_CG", "mut_CG_per","mut_CHG", "mut_CHG_per","mut_CHH", "mut_CHH_per",
			  			   "avge_meth_level_wt"  , "avge_meth_level_mut")), "\n";
				
	my $last_index = scalar(@{$ref}) - 1;
	
	for my $i(0..$last_index){
		
		print OUT join("\t",  @{$ref->[$i]}), "\n";
	}
	
	close(OUT);	
	
}

sub last_filter{
	my ($list_ref, $mC_list_ref, $dep_list_ref, $percentage_list_ref, $final_list_ref, $C_num_ref) = @_;
	
	my @types = ("CG", "CHG", "CHH");
	
	
	
	my $last_index = scalar(@{$list_ref}) - 1;
	
	foreach my $i (0..$last_index){
	#	my ($chr, $start, $end, $num) = split "_", $list_ref->[$i];
		my @a = split "_", $list_ref->[$i];
		
		my ($chr, $start, $end, $num) = @a[0..3];
		
		my (%divide_h, %divide_l); # record like "$mC_mut/$dep_mut="
		my (%per_h, %per_l);
		
		foreach my $type(@types){
			my ($mC_h , $dep_h ) = (0) x 2;
			my ($mC_l, $dep_l) = (0) x 2;
			
			if(defined $mC_list_ref->{"hyper"}->{$type}->[$i]) { $mC_h =  $mC_list_ref->{"hyper"}->{$type}->[$i]}
			if(defined $mC_list_ref->{"hypo"} ->{$type}->[$i]) { $mC_l =  $mC_list_ref->{"hypo"}->{$type}->[$i]}
		
			if(defined $dep_list_ref->{"hyper"}->{$type}->[$i]) { $dep_h =  $dep_list_ref->{"hyper"}->{$type}->[$i]}
			if(defined $dep_list_ref->{"hypo"}->{$type}->[$i]) {  $dep_l =  $dep_list_ref->{"hypo"}->{$type}->[$i]}
			
			$divide_h{$type} = "$mC_h/$dep_h=";
			$divide_l{$type} = "$mC_l/$dep_l=";
			
			if($dep_h != 0) { $per_h{$type} = sprintf("%.4f",  100 * $mC_h/$dep_h);		}
			else            { $per_h{$type} = "NA";  }	

			if($dep_l != 0)  { $per_l{$type}  = sprintf("%.4f",  100 * $mC_l/$dep_l); }
			else			 { $per_l{$type}  = "NA";}	
		}
		
		
		my $avge_h  =  cal_mean(\@{$percentage_list_ref->{"hyper"}->[$i]});
		my $avge_l  =  cal_mean(\@{$percentage_list_ref->{"hypo"}->[$i]});
		
		my ($cg_num, $chg_num, $chh_num) = (0) x 3;
		if(defined $C_num_ref->{$i}->{"CG"})  {$cg_num   = $C_num_ref->{$i}->{"CG"} }
		if(defined $C_num_ref->{$i}->{"CHG"}) {$chg_num  = $C_num_ref->{$i}->{"CHG"} }
		if(defined $C_num_ref->{$i}->{"CHH"}) {$chh_num  = $C_num_ref->{$i}->{"CHH"} }
		
		
		if($avge_l eq "NONE" or $avge_l == 0){
			if ($avge_h > 0){
		#		push @{$final_list_ref}, [$chr, $start, $end, $num, $cg_num, $chg_num, $chh_num,
				push @{$final_list_ref}, [ @a, $cg_num, $chg_num, $chh_num,
							 $divide_h{"CG"}, $per_h{"CG"},  $divide_h{"CHG"}, $per_h{"CHG"}, $divide_h{"CHH"}, $per_h{"CHH"},
					 		 $divide_l{"CG"}, $per_l{"CG"},    $divide_l{"CHG"}, $per_l{"CHG"},   $divide_l{"CHH"}, $per_l{"CHH"} , $avge_h ,$avge_l];
			}
		}elsif($avge_h / $avge_l >=2 ){
	#		push @{$final_list_ref}, [$chr, $start, $end, $num, $cg_num, $chg_num, $chh_num,
			push @{$final_list_ref}, [@a, $cg_num, $chg_num, $chh_num,
							 $divide_h{"CG"}, $per_h{"CG"},  $divide_h{"CHG"}, $per_h{"CHG"}, $divide_h{"CHH"}, $per_h{"CHH"},
					 		 $divide_l{"CG"}, $per_l{"CG"},    $divide_l{"CHG"}, $per_l{"CHG"},   $divide_l{"CHH"}, $per_l{"CHH"} , $avge_h ,$avge_l];
		}
			
	}
	print STDERR "last_filter\t" , scalar(@{$final_list_ref}), "\n";
}


sub read_isMeth_get_mC_dep{
	my ($hyper_file, $hypo_file, $region_in_list_ref, 
	$mC_list_ref, $dep_list_ref, $percentage_list_ref, $C_num_ref) = @_;
	
	croak unless (-e $hyper_file);
	croak unless (-e $hypo_file);
	
	open(HYPER, $hyper_file) or croak;
	open(HYPO,  $hypo_file)  or croak;
	
	my $h_hyper = <HYPER>;
	my $h_hypo  = <HYPO> ;
	
	my ($l_h, $l_l) ;# hyper_high_h ; hypo_low_l;
	
	while ($l_h =  <HYPER>, $l_l = <HYPO>){
		chomp $l_l ;
		chomp $l_h ;
		
		my @a_h = split "\t", $l_h ;
		my @a_l = split "\t", $l_l ;

		my ($chr, $pos) = @a_h[0..1];
		
		next unless (defined $region_in_list_ref->{$chr}->[$pos] );
		
		my $type = $a_h[3];
		
		my $index = $region_in_list_ref->{$chr}->[$pos];
		
		######
		#
		############
		$C_num_ref->{$index}->{$type} ++;
		#######################
		
		my ($mC_h , $dep_h, $per_h ) = @a_h[4..6]; 
		my ($mC_l , $dep_l, $per_l ) = @a_l[4..6]; 
		
		$mC_list_ref->{"hyper"}->{$type}->[$index] += $mC_h;
		$mC_list_ref->{"hypo"}->{$type}->[$index]  += $mC_l;
		
		$dep_list_ref->{"hyper"}->{$type}->[$index] += $dep_h;
		$dep_list_ref->{"hypo"}->{$type}->[$index]  += $dep_l;
		
		push @{$percentage_list_ref->{"hyper"}->[$index]} , $per_h ;
		push @{$percentage_list_ref->{"hypo"}->[$index]}  , $per_l ;
	}
	
	
	close(HYPER);
	close(HYPO);
}


# adjust_list_border(\@list_hyper_merged, \@hyper_list_border, \%DMC, DMC_cutoff);
sub adjust_list_border{
	my ($ref_merged_list, $ref_border_list, $ref_DMC, $DMC_cutoff_sub) = @_;
	
	my $last_index = scalar(@{$ref_merged_list}) - 1;
	
	foreach my $i (0..$last_index){
	#	my ($chr, $start, $end) = split "_", $ref_merged_list->[$i];
		my @pts = split "_", $ref_merged_list->[$i];
		my ($chr, $start, $end) = @pts[0..2];
		my $DMC_num = 0;
		
		my %type_nums = ("CG" => 0, "CHG" => 0, "CHH" => 0 );
		
		for my $j ($start..$end){
			if(defined $ref_DMC->{$chr}->[$j]){
				$DMC_num++;
				my $type_temp = $ref_DMC->{$chr}->[$j];
				$type_nums{$type_temp}++;
			}
		}
	
		
		if($DMC_num >= $DMC_cutoff_sub ){
			my ($real_start, $real_end) = ($start, $end) ;
			for my $j($start..$end){
				if(defined $ref_DMC->{$chr}->[$j]){
					$real_start = $j;
					last;
				}
			}
			for (my $j = $end; $j >= $real_start; $j--){
				if(defined $ref_DMC->{$chr}->[$j]){
					$real_end = $j;
					last;
				}
			}
	#		push @{$ref_border_list}, join("_" , ($chr, $real_start, $real_end, $DMC_num));
	
			my 	$type_label = join("+", ( $type_nums{CG}, $type_nums{CHG}, $type_nums{CHH} ));
			push @{$ref_border_list}, join("_" , ($chr, $real_start, $real_end, $DMC_num, $type_label));
		}
	}

	print STDERR "after_adjusted\t", scalar( @{$ref_border_list} ), "\n";
}

# adjust_list_border(\@list_hyper_merged, \@hyper_list_border, \%DMC, DMC_cutoff);
# record_region_pos(\@list, \%hash)
sub record_region_pos{
	my ($ref_list, $ref_h) = @_;
	my $last_index = scalar(@{$ref_list}) - 1;
	
	foreach my $i (0..$last_index){
	#	my ($chr, $start, $end, $num) = split "_", $ref_list->[$i];
		my @pts = split "_", $ref_list->[$i];
		my ($chr, $start, $end) = @pts[0..2];
		
		for my $j($start..$end){
			$ref_h->{$chr}->[$j] = $i;
		}
	}
}


# merge_list(\@raw_list_hyper, \@list_hyper_merged, $allowed_gap);
sub merge_list{
	my ($ref_raw, $ref_new, $gap) = @_;
	my ($last_chr, $last_start, $last_end ) = ("chr0", 0, 0);
	my $last_index = scalar( @{$ref_raw} ) - 1;
	
	print STDERR "before_merge\t", $last_index + 1, "\n";
	#print STDERR "last_index = $last_index\n\n";
	
	for my $i(0..$last_index){
#		my ($chr, $start, $end) = split "_", $ref_raw->[$i];
		my @pts = split "_", $ref_raw->[$i];
		my ($chr, $start, $end) = @pts[0..2];
		if($chr ne $last_chr || $start > $last_end + $gap + 1 ){
			if($last_chr ne "chr0") {
				push @{$ref_new}, join("_", ( $last_chr, $last_start, $last_end)) ;
			}
			($last_chr, $last_start, $last_end ) =  ($chr, $start, $end);
		}
		else{
			$last_end = $end;
		}
	
	}
	if($last_chr ne "chr0") {
		push @{$ref_new}, join("_", ( $last_chr, $last_start, $last_end)) ;
	}
	
	print STDERR "after_merge\t", scalar(@{$ref_new}), "\n";

}


sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too
}

#		           WT	    mut
#		           word2   ~word2
#mC	     word1    n11      n12 | n1p
#No.T	~word1    n21      n22 | n2p
#      	          --------------
#		          np1      np2   npp

#$n11 = $wt_mC;
#$n1p = $wt_mC + $mut_mC;
#$np1 = $wt_depth;
#$npp = $wt_depth + $mut_depth;

sub cal_Fisher_exact_test{
	my($n11, $n1p ,$np1, $npp) = @_;
	my  $p_value = calculateStatistic( n11=>$n11,
					   n1p=>$n1p,
					   np1=>$np1,
					   npp=>$npp);
	my $errorCode;
	if( ($errorCode = getErrorCode())){
#		print STDERR $errorCode." - ".getErrorMessage(),"\n\n";
  		return "error";
	}
  	else{
		return $p_value;
	}
}

sub cal_mean{
	my ($ref) = @_;
	my $sum = 0;
	my $num = 0;
	
	foreach my $item (@{$ref}){
		$sum+=$item;
		$num++;
	}
	if($num == 0){
		return "NONE";
	}else{
		return ( sprintf ( "%.4f", 100 * $sum / $num));
	}
}
