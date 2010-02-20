Proth Prime Search
Searches for Proth primes, for the Proth Prime Search Extended sieve.

If you are running Linux, just run `./ppsieve.sh'.  You can probably ignore the
rest of this readme.  Minimum requirements for binaries are kernel 2.6.15,
GlibC 2.3.4.  If you don't meet these, compile from source (see the end of this
file.)

An example command line to run the Proth Prime Search Extended sieve:
ppsieve-######### -p####G -P####G -ippse_37TE1.txt -ffpps2_####G-####G.txt -t4 -q 

A similar command is found in ppsieve.bat for Windows users.

===============================================================

ppsieve-######### 	The executable.  Pick the right one for your OS and CPU.

-p####G  sieve start value (#### 000 000 000)

-P####G  sieve end value (#### 000 000 000)

-i       input sieve filename

-f       output factorfile (this one you shall upload when range is completed)

-t4      Only multicore systems: `-t --threads N'
         will spread the work of the sieve over N cores. Although it is a little
         LESS productive, it's easier to manage. (-t4 on a quad, -t2 on a duo,
	 not used on single core - ok to delete)
         
-q       quiet mode

-a	 Alternate algorithm, faster on AMD processors.  `-a yes' turns it on;
	 `-a no' forces it off.  Default is an educated guess at what's best.

Default is not to change process priority.
-z ##   sets priority (idle, low, normal, or give a number from 0 to 19)

This is far from finished. Run `./ppsieve -h' for a list of options.

Options can be given default values in a configuration file and overridden
on the command line. E.g. to start 2 threads add this line to ppconfig.txt:

  threads=2


Integer arguments to config file or command line options can be given using
this shorthand: K=10^3, M=10^6, G=10^9, T=10^12, P=10^15, or k=2^10, m=2^20,
g=2^30, t=2^40, p=2^50. So for example the integer 1000000 can be given as
any of: 1M, 1000K, 1000000, 1e6, 10e5, etc.


For best results, adjust the time spend sieving candidates p by specifying
the prime sieve depth q (e.g. sieve p for factors up to q=10e6):

  ./ppsieve -i FILE -p 7700e12 -P 7701e12 -Q 10e6

"This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version."

Source code should be found at http://home.att.net/~k.brazier/programs/
