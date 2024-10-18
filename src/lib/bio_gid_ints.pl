
:- use_module(library(lib)).

bio_gid_ints_defaults([]).

/** bio_gid_ints( AtomsOr, Ints ).

Map a list of atoms or integers to their integer values.

Opts
  * debug(Dbg=false)
    informational, progress messages

Examples
==
?- bio_gid_ints(['373904'], Ints).
==

@author nicos angelopoulos
@version  0.1 2024/10/18

*/

bio_gid_ints( [], [] ).
bio_gid_ints( [E|Es], [I|Is] ) :-
     ( integer(E) -> I = E
                   ; atom_number(E,I),
                     integer(I)
     ),
     bio_gid_ints( Es, Is ).
