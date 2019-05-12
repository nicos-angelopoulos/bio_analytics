
:- use_module(library(lib)).

% :- set_prolog_flag(lib_suggests_warns,debug).
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
plot graphs for the bt.csv experiment. 
This version uses the stored over representation csv in 
data/silac/bt_go_vers.csv.

The following will create a bunch of svg images in relative 
directory bt_go_graphs/ . 

==
?- [pack('bio_analytics/examples/bt_go_string')].
?- bt_go_string.
==
@author nicos angelopoulos
@version  0.1 2019/5/5
@see exp_go_over/3

*/
bt_go_string :-
    absolute_file_name( pack('bio_analytics/data/silac/bt.csv'), ExpF ),
    Mess1 = 'Starting gene ontology over-representation analysis and string graph plotting for bt.csv',
    debug( bt_go_string, Mess1, [] ),
    Dir = bt_go_graphs,
    absolute_file_name( pack('bio_analytics/data/silac/bt_go_over.csv'), Gov ),
    % exp_go_over_string_graphs( ExpF, Gov, Dir, [go_over_pv_cut(0.01)] ),
    exp_go_over_string_graphs( ExpF, Gov, Dir, [] ),
    Mess2 = 'String graphs for gene ontology over-represented terms are in dir: ~p',
    debug( bt_go_string, Mess2, [Dir] ).
