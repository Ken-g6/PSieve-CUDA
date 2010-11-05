#!/usr/bin/perl
# Convert OpenCL code to a C-style string.
# I, the author, Ken Brazier, place this Perl code in the Public Domain.
use strict;

# Load the code in C++ only (to avoid app.c) and start the string!
print <<EOF;
#ifdef __cplusplus
#ifndef _APPCL_CL
#define _APPCL_CL
const char* appcl= \\
EOF

my $in_ifdef = 0;
my $in_blockquote = 0;
while(<>) {
	# Chomp.
	s/[\r\n]*$//;
	# Remove C++ comments - can be dangerous but makes smaller code!
	# The unless makes it safer, but possibly less effective.
	s/\/\/.*$// unless /".*\/\/.*"/;

	# Remove a few C comments as well.
	s/\/\*.*\*\///g;

	# Remove C block comments.  Dangerous, but perhaps not as dangerous as not removing them!
	if($in_blockquote) {
		if(s/^.*\*\///) { $in_blockquote = 0; }
		else { next; }
	}
	if(s/\/\*.*$//) {
		++$in_blockquote;
	}

	# Fix \'s not at the ends of lines.
	s/\\(.)/\\\\$1/g;

	# Remove (some) whitespace.
	s/^[ 	]*//;
	s/[ 	]*$//;
	next if(/^$/);

	# Remove ifdefs and everything in them.
	# Dangerous if you #define anything!
	# Also Dangerous if you use #else.
	#if(/^#if/) {
	#$in_ifdef++;
	#next;
	#}
	#if($in_ifdef > 0) {
	#--$in_ifdef if(/^#endif$/);
	#next;
	#}

	# Escape quotes - vital!
	s/"/\\"/g;

	# Print the line in a string.  Adjacent strings are concatenated.
	if(s/\\$//) {
		# If string ended in an escaped newline, just print the string without a newline.
		print '"', $_, "\" \\\n";
	} else {
		print '"', $_, "\\n\" \\\n";
	}
}
print ";\n#endif\n#endif\n";
