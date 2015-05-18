#!/usr/bin/perl

@dirs = qw (modules src);

foreach $dir (@dirs) {

    chdir $dir;

    opendir DIR, ".";
    @files = grep /\.F90$/, readdir DIR;
    closedir DIR;

    foreach $file (@files) {
	$new = $file;
	$new =~ s/\.F90/\.f90/;
	rename $file, $new;
    }

    change_makefile ();

    chdir "..";
}

change_configure ();

exit;

sub change_configure {

    rename "configure.bak", "configure";
}

sub change_makefile {

    my $makefile = "Makefile";

    rename "Makefile.bak", "Makefile";
}
