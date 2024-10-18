
:- lib(bio_gid_ints/2).

bio_gid_std_defaults([]).

/** bio_gid_std( +Id, +GidS, -StdS ).

Map GidS of type ID to start term represention. 

Typically this is to map GidS (one Id or a list of ids) that might be atomic from 
an operation such as getting them from R, to integers.

Currently only ncbi ids are converted to integers. Everything else is passed through.

Examples
==
?- bio_gid_std( ncbi, '373904', Std ).

?- bio_gid_std( ncbi, ['373904'], Std ).
==

@author nicos angelopoulos
@version  0.1 2024/10/18
@tbd other Id types

*/

bio_gid_std( ncbi, GidS, StdS ) :-
     !,
     ( is_list(GidS) ->
          bio_gid_ints( GidS, StdS )
          ;
          bio_gid_ints( [GidS], [StdS] )
     ).

bio_gid_std( _, GidS, StdS ) :-
     GidS = StdS.
