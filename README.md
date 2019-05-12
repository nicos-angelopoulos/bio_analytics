Computational biology data analytics.

Collects a number of biological data analytics tasks.
This library provides tools for the bio_db served data,
empowering downstream analyses of user's experimental data.

Installation and loading: 
==
?- pack_install(bio_analytics).  
      % also installs pack(lib). other dependencies are at load time via lib
?- use_module(library(lib)).
?- lib(bio_analytics).
==

The library comes with one experimental dataset: data/silac/bt.csv which
is used in the example files in directory examples/ .

@author nicos angelopoulos
@version  0.1 2019/4/22
@version  0.2 2019/5/08
@version  0.3 2019/5/11
@licence  MIT

See 
[bio_analytics](http://stoics.org.uk/~nicos/sware/bio_analytics)
[bio_db](http://stoics.org.uk/~nicos/sware/bio_db)
