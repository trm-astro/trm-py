#!/usr/bin/perl -w

use strict;

use PGPLOT;

(@ARGV == 1) or die "usage: device\n";

my $device = shift;

my (@xt, @yt, $n);

# speed of light

my $C = 299792.6;

my $x1 = 845;
my $x2 = 870;

pgopen($device);
pgslw(2);
pgscf(2);
pgsch(1.5);
pgscr(3,0,0.5,0);
pgenv($x1, $x2, -100, 800, 0, 0);
pglab('Wavelength (nm)', "Flux density", " ");

pgsci(3);
#pgsls(4);

@xt = ();
@yt = ();
$n = 0;

open(SPEC, "8498.fit");
while(<SPEC>){
    ($xt[$n],$yt[$n++]) = split(' ');
}
close(SPEC);

foreach (@xt) {
    $_ /= $C;
    $_ *= 849.802;
    $_ += 849.802;
}

pgmove($x1, 0);
pgdraw($xt[0], $yt[0]);

pgline(scalar(@xt), \@xt, \@yt);

my $xsave = pop @xt;
my $ysave = pop @yt;

@xt = ();
@yt = ();
$n = 0;


open(SPEC, "8542.fit");
while(<SPEC>){
    ($xt[$n],$yt[$n++]) = split(' ');
}
close(SPEC);

foreach (@xt) {
    $_ /= $C;
    $_ *= 854.209;
    $_ += 854.209;
}


pgmove($xsave, $ysave);
pgdraw($xt[0], $yt[0]);

pgline(scalar(@xt), \@xt, \@yt);

$xsave = pop @xt;
$ysave = pop @yt;

@xt = ();
@yt = ();
$n = 0;

open(SPEC, "8662.fit");
while(<SPEC>){
    ($xt[$n],$yt[$n++]) = split(' ');
}
close(SPEC);

foreach (@xt) {
    $_ /= $C;
    $_ *= 866.214;
    $_ += 866.214;
}

pgmove($xsave, $ysave);
pgdraw($xt[0], $yt[0]);

pgline(scalar(@xt), \@xt, \@yt);

$xsave = pop @xt;
$ysave = pop @yt;

pgmove($xsave, $ysave);
pgdraw($x2, 0);

# finally plot the data

@xt = ();
@yt = ();
my @et = ();
$n = 0;

open(SPEC, "triplet.dat");
$n = 0;
while(<SPEC>){
    ($xt[$n],$yt[$n],$et[$n++]) = split(' ');
}
close(SPEC);

foreach (@et) {
#    $_ *= 3;
}
foreach (@xt) {
    $_ /= 10.;
}

pgsci(2);
pgerrb(2, scalar(@xt), \@xt, \@yt, \@et, 0);
pgerrb(4, scalar(@xt), \@xt, \@yt, \@et, 0);
pgsci(1);
pgsch(0.8);
pgpt(scalar(@xt), \@xt, \@yt, 17);

pgclos();

exit;
