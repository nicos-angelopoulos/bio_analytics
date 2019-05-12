
:- use_module(library(lib)).

% :- set_prolog_flag(lib_suggests_warns,debug).
% use the above to see all suggested code dependencies.
% the one below installs all suggested code
:- initialization(set_prolog_flag(lib_suggests_warns,install), now).

:- lib(bio_analytics).
:- debug(bt_go).

/** bt_go.

Run GO term ever-representation (exp_go_over/3) 
analysis on the data/silac/bt.csv data.

At load time the code sets flag that will 
auto-load R libraries that are suggested by bio_analytics.

@author nicos angelopoulos
@version  0.1 2019/5/5
@see exp_go_over/3

*/

bt_go :-
    absolute_file_name( pack('bio_analytics/data/silac/bt.csv'), DatF ),
    tmp_file( go_over, TmpF ),
    Mess1 = 'Starting gene ontology over-representation analysis for bt.csv',
    debug( bt_go, Mess1, [] ),
    exp_go_over( DatF, TmpF, [to_file(true)] ),
    Mess2 = 'Gene ontology over-represented terms in file: ~p',
    debug( bt_go, Mess2, [TmpF] ).
