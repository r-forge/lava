#!/usr/bin/perl

$dir = $ARGV[0];
$cmd = "grep Version ".$dir."/DESCRIPTION";
$val = `$cmd`;
@comp = split(/ /,$val);
$version = $comp[1];
chop($version);

@filename = split(/\//,$dir);
$n = @filename;
$file = $filename[$n-1];

if ($#ARGV>0) {
  print($version . "\n");
} else {
  print($file . "_" . $version . ".tar.gz\n")
}



