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

my $debug = 0;

my $depth_cutoff = 4;

my $usage = "$0 \n <isMeht_file> <sample_label> <bed_like_file> <output> \n\n";
die $usage unless (@ARGV == 4);

my $isMeth_file 	= shift or die;
my $sample_label 	= shift or die;
my $bed_file		= shift or die;



my $output = shift or die;

die unless (-e $bed_file);
die unless (-e $isMeth_file);


die if (-e $output);




my @bed_list;
my %pos;


my $last_one_index = read_bed_file($bed_file, \@bed_list);
# head is index 0 and "\n" is chompped.

print STDERR "last: line_num: ", $last_one_index, "\n\n";

record_region_pos(\@bed_list, \%pos);


my ( %C_type_nums,  %called_mC_nums, %meth_level_sum,  %seq_mC_sum,  %seq_depth_sum ); #Weighted methylation level

# %{wt/mut}->{1/2}->{CG/CHG/CHH/total}

Initialization_type_num(\%C_type_nums, $last_one_index);

Initialization(\%called_mC_nums, $last_one_index);
Initialization(\%meth_level_sum, $last_one_index);
Initialization(\%seq_mC_sum, $last_one_index);
Initialization(\%seq_depth_sum, $last_one_index);

#read_isMeth_get_info( $isMeth_file_WT,  $isMeth_file_mut,  \%pos,
#			\%C_type_nums,  \%called_mC_nums, \%meth_level_sum,  \%seq_mC_sum,  \%seq_depth_sum,
#			$depth_cutoff		
#);

read_isMeth_get_info_0_1( 
			$isMeth_file,
			\%pos,
			\%C_type_nums,  
		#	\%called_mC_nums, 
			\%meth_level_sum,  
			\%seq_mC_sum,  
			\%seq_depth_sum,
			$depth_cutoff		
);


check_type_num(\%C_type_nums, $last_one_index );


check(\%called_mC_nums, $last_one_index );
check(\%meth_level_sum, $last_one_index );
check(\%seq_mC_sum, $last_one_index );
check(\%seq_depth_sum, $last_one_index );

output_list ($output, \@bed_list, 
			  \%C_type_nums,  \%called_mC_nums, \%meth_level_sum,  \%seq_mC_sum,  \%seq_depth_sum,
			  $sample_label
			  );


exit;
####################
# add Nov16, 2013
###################
sub simple_chr{
	my ($chr) = @_;
	if( $chr =~ /chr/i){
		$chr =~  s/chr//i;
	}
	if($chr eq "M" ){
		$chr = "Mt";
	}elsif( $chr eq "C"){
		$chr = "Pt";
	}
	return $chr;
}
#####################

sub cal_meth{
	my ($first_ref, $second_ref, $ref_a) = @_;
	
	my @types_sub = ("CG", "CHG", "CHH", "total");
	my $i = -1;
	foreach my $type_sub (@types_sub){
		$i++;
		my ($first, $second) = ($first_ref->{$type_sub}  , $second_ref->{$type_sub} ); 
		if( $second != 0){
			$ref_a->[$i] = eval sprintf( "%.3f", 100 * $first / $second );
		}
	}
	
}


sub assign{
	my ($first_ref, $ref_a) = @_;
	
	my @types_sub = ("CG", "CHG", "CHH", "total");
	my $i = -1;
	foreach my $type_sub (@types_sub){
		$i++;
		$ref_a->[$i] = $first_ref->{$type_sub};
	}
}



#done
sub check_type_num{
	my ($ref, $last_index_sub) = @_;
	#my @labs_sub = ("wt", "mut");
	my @types_sub = ("CG", "CHG", "CHH");
	
		for my $i (1..$last_index_sub){
			
			my $sum = 0;
			foreach my $type_sub (@types_sub){
				$sum  +=  $ref->{$i}->{$type_sub} ;
			}
			
			my $t = $ref->{$i}->{"total"};
			if ( abs( $sum - $t ) > 0.01 ){
				print STDERR "maybe wrong  \n", $i, "\n\n";
				print STDERR "sum = $sum != $t \n\n";
			}
			
			
		}
	
}


#modified
sub check{
	my ($ref, $last_index_sub) = @_;
	#my @labs_sub = ("wt", "mut");
	my @labs_sub = ( "mut");

	my @types_sub = ("CG", "CHG", "CHH");
	
	foreach my $lab_sub (@labs_sub){
		for my $i (1..$last_index_sub){
			
			my $sum = 0;
			foreach my $type_sub (@types_sub){
				$sum  +=  $ref->{$lab_sub}->{$i}->{$type_sub} ;
			}
			
			my $t = $ref->{$lab_sub}->{$i}->{"total"};
			if ( abs( $sum - $t ) > 0.01 ){
				print STDERR "maybe wrong  \n", join("\t", ( $lab_sub, $i )), "\n\n";
				print STDERR "sum = $sum != $t \n\n";
			}
			
			
		}
	}
}


#modified
sub Initialization_type_num{
	my ($ref, $last_index_sub) = @_;
#	my @labs_sub = ("wt", "mut");
	my @types_sub = ("CG", "CHG", "CHH", "total");
	for my $i (1..$last_index_sub){
		foreach my $type_sub (@types_sub){
			$ref->{$i}->{$type_sub} = 0;
			
		}
	}
}


# Initialization(\%seq_depth_sum, $last_one_index);

# modified
sub Initialization{
	my ($ref, $last_index_sub) = @_;
	#my @labs_sub = ("wt", "mut");
	my @labs_sub = ( "mut");
	
	my @types_sub = ("CG", "CHG", "CHH", "total");
	foreach my $i (1..$last_index_sub){
		foreach my $type_sub (@types_sub){
			foreach my $lab_sub (@labs_sub){
				$ref->{$lab_sub}->{$i}->{$type_sub} = 0;
			}
		}
	}
}

sub read_bed_file{
	my ($file, $ref) = @_;
	die unless (-e $file);
	
	open(IN, $file) or die;
	my $i = -1;
	while(<IN>){
		$i++;
		chomp;
		$ref->[$i] = $_;
	}
	close(IN);
	
	return $i;
}


sub record_region_pos{
	my ($ref_list, $ref_h) = @_;
	my $last_index = scalar(@{$ref_list}) - 1;
	
	foreach my $i (1..$last_index){
		my @a = split "\t", $ref_list->[$i];
		my ($chr, $start, $end ) = @a[0..2];
		#######
		$chr = simple_chr($chr);
		#########
		for my $j($start..$end){
			$ref_h->{$chr}->[$j] = $i; #record index
		}
	}
}



sub output_list{
	my ($file, $ref_list, 	
	    $C_type_nums_ref,  $called_mC_nums_ref, $meth_level_sum_ref,  $seq_mC_sum_ref,  $seq_depth_sum_ref,
	    $sample_label_sub
	) = @_;
	
	my @types = ("CG", "CHG", "CHH", "total");

	die if ( -e $file);
	open(OUT, ">>$file") or die;
	my $last_index = scalar(@{$ref_list}) - 1;
	my $head = $ref_list->[0];
	
	
	print OUT join("\t", ($head, 
			      "mCG_" . $sample_label_sub, "mCHG_" . $sample_label_sub, "mCHH_" . $sample_label_sub, "mC_" . $sample_label_sub,
	
	)), "\n";
	
	foreach my $i (1..$last_index){

		
		my @weighted_level_mut = ("NA") x 4;
		

		my @C_type_nums_a    = (0) x 4;
		
		assign(\%{$C_type_nums_ref -> {$i}}, \@C_type_nums_a);
		

		cal_meth( \% {$seq_mC_sum_ref->{"mut"}->{$i} }, \% {$seq_depth_sum_ref->{"mut"}->{$i} }, \@weighted_level_mut );
		
		
		print OUT join("\t", ( $ref_list->[$i],
				      @weighted_level_mut
				)), "\n";		
	}	
	close(OUT);	
}



sub read_isMeth_get_info_0_1{
	my ($isMeth_file_mut_sub,
		$pos_ref,
		$C_type_nums_ref,  
	#	$called_mC_nums_ref, 
		$meth_level_sum_ref,  
		$seq_mC_sum_ref, 
		 $seq_depth_sum_ref,
		$depth_cutoff_sub) = @_;
	
	die unless (-e $isMeth_file_mut_sub) ;
	
	open(MUT, $isMeth_file_mut_sub) or die;
	
	my $h_mut  = <MUT>;
	
	my  $l_mut;
	
	while( $l_mut = <MUT>){
		chomp $l_mut;
		my @a_mut = split "\t", $l_mut;
		my ($chr, $pos) = @a_mut[0..1];
		##############
		$chr = simple_chr($chr);
		 
		###############
		
		next unless (defined $pos_ref->{$chr}->[$pos] );
		
		my ($mC_mut, $dep_mut, $per_mut ) = @a_mut[4..6]; 
		
		if( $dep_mut >= $depth_cutoff_sub ){
			
			my $index = $pos_ref->{$chr}->[$pos];
			my $type = $a_mut[3];

			$C_type_nums_ref-> {$index}->{$type}++;
			$C_type_nums_ref-> {$index}->{"total"}++;
				
			$meth_level_sum_ref->{"mut"}->{$index}->{$type} += $per_mut;
			$meth_level_sum_ref->{"mut"}->{$index}->{ "total" } += $per_mut;
			$seq_mC_sum_ref->{"mut"}->{$index}->{$type} += $mC_mut;
			$seq_mC_sum_ref->{"mut"}->{$index}->{"total" } += $mC_mut;
			$seq_depth_sum_ref->{"mut"}->{$index}->{$type} += $dep_mut;
			$seq_depth_sum_ref->{"mut"}->{$index}->{"total" } += $dep_mut;
		}
	}
	
	close(MUT);
}
