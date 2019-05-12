
:- use_module(library(lib)).

:- set_prolog_flag(lib_suggests_warns,debug).
% use the above to see all suggested code dependencies. This also installs R dependencies.
% the one below installs all suggested code
% :- initialization(set_prolog_flag(lib_suggests_warns,install), now).

:- lib(bio_analytics).
:- debug(bt_go_string).
:- debug(exp_go_over_string_graphs).
:- debug(exp_gene_family_string_graph).

/** bt_go_string.

Run GO term ever-representation (exp_go_over/3) 
analysis on the data/silac/bt.csv data and then 
plot graphs for 

At load time the code sets flag that will 
auto-load R libraries that are suggested by bio_analytics.

The following will create a bunch of svg images in relative 
directory bt_go_graphs/ . 

==
?- [pack('bio_analytics/examples/bt_over_go_string')].
?- bt_over_go_string.
==
@author nicos angelopoulos
@version  0.1 2019/5/5
@see exp_go_over/3

*/
bt_over_go_string :-
    absolute_file_name( pack('bio_analytics/data/silac/bt.csv'), ExpF ),
    Mess1 = 'Starting gene ontology over-representation analysis and string graph plotting for bt.csv',
    debug( bt_go_string, Mess1, [] ),
    Dir = bt_over_go_graphs,
    exp_go_over_string_graphs( ExpF, _Gov, Dir, [go_over_pv_cut(0.01)] ),
    Mess2 = 'String graphs for gene ontology over-represented terms in dir: ~p',
    debug( bt_go_string, Mess2, [Dir] ).
