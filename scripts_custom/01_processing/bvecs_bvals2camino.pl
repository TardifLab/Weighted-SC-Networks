#! /usr/bin/env perl

###################
### This function will write diffusion directions from bvecs bvals to a .scheme file 
###  needed for camino (and AMICO)
###################

require 5.001;
use Getopt::Tabular;
use MNI::Startup qw(nocputimes);
use MNI::Spawn;
use MNI::FileUtilities qw(check_output_dirs);
use File::Basename;

if($0 =~ /[\/A-Za-z_0-9]+\/([A-Za-z0-9_]+)\.pl/) {$program = $1;}	#program name without path
$Usage = <<USAGE;

Usage: $program.pl -vec <bvecs from dcm2nii> -val <bvals from dcm2nii> -o <.scheme>
-help for options

### This function will write diffusion directions from bvecs bvals to a .scheme file 
###  needed for camino (and AMICO)


USAGE

@args_table = (["-vec","string",1,\$bvec,"b vector file from dcm2nii"],
	       ["-val","string",1,\$bval,"b values file from dcm2nii"],
	       ["-o","string",1,\$scheme,"name of scheme file for camino"]
	       );

Getopt::Tabular::SetHelp ($Usage, '');

GetOptions(\@args_table, \@ARGV,\@args)|| exit 1;

open(SCHEME,"> $scheme") or die "can't create $scheme: $!";

## open vector and bvalue files
open(BVEC,"< $bvec") or die "can't open $bvec: $!";
@bvecs = <BVEC>; #slurp all files, each line at an index
open(BVAL,"< $bval") or die "can't open $bval: $!";
@bval = <BVAL>; #slurp all files, each line at an index
#print "n_dir:$n_dir\n\n";

#print "$bval[0]\n $bval[1]\n $bvecs[0]\n\n $bvecs[1]\n\n $bvecs[2]\n\n";

#### Some monkeying around to get the right format for the list

print SCHEME "VERSION: BVECTOR\n";

@b_values = split(" ", $bval[0]);
@x_direction = split(" ", $bvecs[0]);
@y_direction = split(" ", $bvecs[1]);
@z_direction = split(" ", $bvecs[2]);
for ($i=0; $i<scalar(@b_values); $i++){
	print SCHEME "@x_direction[$i] @y_direction[$i] @z_direction[$i] @b_values[$i]\n"; 
}

