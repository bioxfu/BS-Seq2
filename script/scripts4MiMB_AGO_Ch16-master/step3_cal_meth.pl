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

my $debug = 0 ;

my $script_base = "cal_meth_level_single.pl";
my ($volume,$dir_base,$perl_file) =  File::Spec->splitpath( $0 );
my $script = File::Spec-> catfile($dir_base, $script_base);
croak($script . " should be in the same directory with step3_cal_meth.pl") unless (-e $script);


my $usage = "$0 \n <isMeth_dir> <input_bed> <outdir> <output_name>\n\n";
croak $usage unless(@ARGV == 4);

my $dir_isMeth = shift or croak;
my $orignal_file = shift or croak;
my $outdir = shift or croak;
my $output_name = shift or croak;


croak unless (-e $orignal_file);
croak unless (-d $dir_isMeth);
croak unless (-d $outdir);

my $output = File::Spec->catfile($outdir, $output_name);
croak if(-e $output);

opendir(DIR, $dir_isMeth);
my @files_raw = grep /\_Meth.txt$/, sort readdir DIR ;
closedir DIR;

my @files;
my @labels;

my $k = -1;
foreach my $file (@files_raw){
	my $curr_label = get_label($file);
	$k++;
	$files[$k] = $file;
	$labels[$k] = $curr_label;
}



if($debug ){
	print STDERR join("\n", @files), "\n\n\n";
}

#my $i = -1;

my $f0 = $orignal_file . ".temp0";

my $cp_cmd = "cp $orignal_file $f0";
print STDERR $cp_cmd, "\n\n";
if(!$debug){
	`$cp_cmd`;
}

#foreach my $file (@files){
for my $i (0..$#files){
	my $file = $files[$i];
	my $label = $labels[$i];
	
	my $bench_file =  $orignal_file . ".temp" . $i ;
	if(!$debug){
		croak unless (-e $bench_file);
	}
	my $isMeth_file  =  File::Spec->catfile($dir_isMeth, $file );
	croak unless (-e $isMeth_file);
	my $out_file =  $orignal_file . ".temp" . ($i+1) ;
	
	croak if(-e $out_file);
	#<isMeht_file> <sample_label> <bed_like_file> <output> 	
	my $cmd = "$script $isMeth_file $label $bench_file  $out_file";
	
	print STDERR $cmd, "\n\n";
	if(!$debug){
		`$cmd`;
	}
}

my $final = $orignal_file . ".temp" . ($#files+1);

my $cp_out_cmd = "cp $final $output";
print STDERR $cp_out_cmd, "\n\n";
if(!$debug){
	`$cp_out_cmd`;
}

#my $rm_cmd = "rm -f $bed_dir/*txt.temp* ";
my $rm_cmd = "rm -f " . $orignal_file . ".temp* ";
print STDERR $rm_cmd, "\n\n";
if(!$debug){
	`$rm_cmd`;
}


exit;

#get_label($file);
sub get_label{
	my ($f) = @_;
	if ($f =~ /(\S+)_Meth/){
		return $1;
	}else{
		croak $f, "file names should look like *_Meth.txt \n\n";
	}
}
