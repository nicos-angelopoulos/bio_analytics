/** go_org_symbol( +Org,  +Go, -Symb ).

Get symbols of Go id for Organism.

Go can be either a GO: atom or an integer.

==
[debug]  ?- go_org_symbol( mouse, 2, Symbs ), write(Symbs), nl, fail.
Akt3
Mef2a
Mgme1
Mpv17
Mrpl15
Mrpl17
Mrpl39
Msto1
Opa1
Slc25a33
Slc25a36
Tymp
false.

[debug]  ?- go_org_symbol( hs, 'GO:0000002', Symbs ), write(Symbs), nl, fail.
AKT3
LONP1
MEF2A
MGME1
MPV17
MSTO1
OPA1
PIF1
SLC25A33
SLC25A36
SLC25A4
TYMP
false.
==

@author nicos angelopoulos
@version  0:1 2019/04/07
*/
go_org_symbol( Alias, Go, Symb ) :-
    bio_db_organism( Alias, Org ),
    go_id( Go, _, Gi ),
    go_org_symbol_1( Org, Gi, Symb ).
    
go_org_symbol_1( gallus, Gi, Symb ) :-
    map_gont_gallus_gont_symb( Gi, _Rel, _Evid, Symb ).
go_org_symbol_1( hs, Gi, Symb ) :-
    map_gont_gont_symb( Gi, _Evid, Symb ).
go_org_symbol_1( mouse, Gi, Symb ) :-
    map_gont_mouse_gont_symb( Gi, _Evid, Symb ).
