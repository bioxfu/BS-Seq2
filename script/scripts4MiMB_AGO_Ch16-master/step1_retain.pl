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

my $debug = 0;

my $usage = "$0 \n<outdir> <cutoff> <number_of_files> <isMeth> <label> [<isMeth> <label> ...] \n\n";
croak $usage unless(@ARGV >= 7);

my $outdir = shift or croak $usage;
my $dep_cutoff = shift or croak $usage;
my $num_of_files = shift or croak $usage;
my $last_index = $num_of_files - 1;

croak($outdir . "does not exist") unless (-d $outdir);

croak ("input file number and labels is NOT equal to " . $num_of_files) unless (@ARGV == 2 * $num_of_files);

my @labels;
my @files;
my @outputs;

for my $i(0..$last_index){
	$files[$i]  = shift or croak;
	$labels[$i] = shift or croak;
	$outputs[$i] = File::Spec->catfile($outdir,  $labels[$i] . "_dep" . $dep_cutoff .  "_Meth.txt" );
	croak ("file " . $files[$i] . "does not exist" )unless (-e $files[$i]);
	croak("output " . $outputs[$i] . " already exists") if (-e $outputs[$i]);
}

print STDERR "depth cutoff is $dep_cutoff \n\n";
if($debug){
	print STDERR join("\n", @labels),  "\n\n";
	print STDERR join("\n", @files),   "\n\n";
	print STDERR join("\n", @outputs), "\n\n";
	exit;
}

my @fhr;
my @fhw;
my $head = join("\t", ("chr", "pos", "strand", "type", "num_C", "depth", "meth_level"));
for my $i(0..$last_index){
	open($fhr[$i], "<" , $files[$i]) or croak "$files[$i]";
	open($fhw[$i], ">" , $outputs[$i]) or croak $outputs[$i], ": $!";
	
	print {$fhw[$i]} $head, "\n";
}

my $l;
my $total_num = 0;

while($l = readline($fhr[0]) ){
	my $flag_print = 0;
	my @out = ();
	my $i = 0;

	@{$out[$i]} = handle_line($l);
	my ($chr, $pos, $strand, $type, $num_C, $depth, $meth_level) = @{$out[$i]};

	if($depth < $dep_cutoff){
		$flag_print = 1;
	}
	
	for $i(1..$last_index){
		$l = readline($fhr[$i]);
		@{$out[$i]} = handle_line($l);
		my ($tmp_chr, $tmp_pos, $tmp_strand, $tmp_type, $tmp_num_C, $tmp_depth, $tmp_meth_level) = @{$out[$i]};
		croak("lines are NOT consistent") unless ($chr eq $tmp_chr and $tmp_pos eq $pos);
		if($tmp_depth < $dep_cutoff){
			$flag_print = 1;
		}
	}
	
	if($flag_print == 0){
		next if ($chr =~ /chloroplast|mitochondria|chrC|chrM|Pt|Mt/i); # skip chloroplast and mitochondria;
		$total_num++;
		for my $j(0..$last_index){
			print {$fhw[$j]} join("\t", @{$out[$j]}) , "\n";
		}
	}
}

for my $i(0..$last_index){
	close($fhr[$i]);
	close($fhw[$i]);
}

print STDERR "\nnumber of cytosines meet the requirement all depth >= $dep_cutoff:\n";
print STDERR $total_num, "\n\n";

exit;

sub round {
    my ($number) = shift;
    return int($number + 0.5 * ($number <=> 0)); # take care of negative numbers too
}

sub handle_line{
	my ($line) = shift;
	chomp $line;
	my ($chr, $pos, $end, $total, $meth_level, $strand) = split "\t", $line; 
	$chr = lc($chr);
	$pos++;#in acgt-count output,base count in a reference starts with 0 
	my ($type, $depth) = split ":", $total;
	my $num_C = round($depth * $meth_level);
	if($type eq "CpG"){
		$type = "CG";
	}
	$meth_level = eval sprintf("%.4f", $meth_level);
	return ($chr, $pos, $strand, $type, $num_C, $depth, $meth_level);
}
