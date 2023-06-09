
/** bio_org.

Bio analytics supports a number of organism specific internal predicates.

This is a doc predicate listing them:
  * bio_conductor_annot_dbi_org/3
  * bio_symbols_map/4
  * go_org_symbol/3

Examples
==
?- bio_org.
==

These predicates are collated in a single source file src/bio_org.pl so that 
in future it is easier to add new organisms and more extensive support for 
existing organisms.

@author nicos angelopoulos
@version  0.1 2023/06/07

*/

bio_org.


/** bio_symbols_map(+Org, +FromID, +Values, -Symbols).

Convert Org-anism specific values of FromID ids to Symbols.

==
?- bio_analytics:bio_symbols_map(human, ensg, [114783], Hymbs).


?- bio_symbols_map(pig, ensg, [], Pymbs).

?- bio_symbols_map(mouse, symb, [], Mymbs).
==

@author nicos angelopoulos
@version  0:2 2023/06/09, added pig, moved to this file and added doc.
@see bio_symbols/3

*/
bio_symbols_map( _Org, symb, Values, Symbs ) :-
     !,
     Values = Symbs.
bio_symbols_map( gallus, FromID, Values, Symbs ) :-
     ( Values = [_-_|_] -> 
          maplist( bio_symbol_paired_chicken(FromID), Values, Symbs )
          ;
          maplist( bio_symbol_map_chicken(FromID), Values, Symbs )
     ).
bio_symbols_map( human, FromID, Values, Symbs ) :-
     ( Values = [_-_|_] -> 
          maplist( bio_symbol_paired_human(FromID), Values, Symbs )
          ;
          maplist( bio_symbol_map_human(FromID), Values, Symbs )
     ).
bio_symbols_map( mouse, FromID, Values, Symbs ) :-
     ( Values = [_-_|_] -> 
          maplist( bio_symbol_paired_mouse(FromID), Values, Symbs )
          ;
          maplist( bio_symbol_map_mouse(FromID), Values, Symbs )
     ).
bio_symbols_map( pig, FromID, Values, Symbs ) :-
     ( Values = [_-_|_] -> 
          maplist( bio_symbol_paired_pig(FromID), Values, Symbs )
          ;
          maplist( bio_symbol_map_pig(FromID), Values, Symbs )
     ).

bio_symbol_paired_chicken( FromID, Id-Val, Symb-Val ) :-
     bio_symbol_map_chicken( FromID, Id, Symb ).
bio_symbol_paired_human( FromID, Id-Val, Symb-Val ) :-
     bio_symbol_map_human( FromID, Id, Symb ).
bio_symbol_paired_mouse( FromID, Id-Val, Symb-Val ) :-
     bio_symbol_map_mouse( FromID, Id, Symb ).
bio_symbol_paired_pig( FromID, Id-Val, Symb-Val ) :-
     bio_symbol_map_pig( FromID, Id, Symb ).

bio_symbol_map_chicken( ncbi, Ncbi, Symb ) :-
     cgnc_galg_cgnc_ncbi( Cgnc, Ncbi ),
     cgnc_galg_cgnc_symb( Cgnc, Symb ).

/*
23.06.09

?- findall( EnsG, (hgnc_homs_hgnc_ensg(Hgnc,EnsG), hgnc_homs_hgnc_symb(Hgnc,Symb),\+ ense_homs_ensg_symb(EnsG,Symb)), EnsGs ).
EnsGs = ['ENSG00000277656', 'ENSG00000227442', 'ENSG00000196101', 'ENSG00000227357', 'ENSG00000227099', 'ENSG00000233697', 'ENSG00000276387', 'ENSG00000275434', 'ENSG00000272546'|...].

?- findall( EnsG, (hgnc_homs_hgnc_ensg(Hgnc,EnsG), hgnc_homs_hgnc_symb(Hgnc,Symb),\+ ense_homs_ensg_symb(EnsG,Symb)), EnsGs ), length(EnsGs,LenInHGNCOnly).
EnsGs = ['ENSG00000277656', 'ENSG00000227442', 'ENSG00000196101', 'ENSG00000227357', 'ENSG00000227099', 'ENSG00000233697', 'ENSG00000276387', 'ENSG00000275434', 'ENSG00000272546'|...],
LenInHGNCOnly = 74.

?- findall( EnsG, (ense_homs_ensg_symb(EnsG,Symb), \+ hgnc_homs_hgnc_ensg(Hgnc,EnsG)), EnsGs ), length(EnsGs,InEnseOnly).
EnsGs = ['ENSG00000232995', 'ENSG00000252969', 'ENSG00000252448', 'ENSG00000201944', 'ENSG00000281394', 'ENSG00000281859', 'ENSG00000199934', 'ENSG00000201898', 'ENSG00000280498'|...],
InEnseOnly = 360.

There is a single mismatch, prefer hgnc
?- hgnc_homs_hgnc_ensg(Hgnc,EnsG), hgnc_homs_hgnc_symb(Hgnc,Symb),\+ ense_homs_ensg_symb(EnsG,Symb).
Hgnc = 4641,
EnsG = 'ENSG00000277656',
Symb = 'GSTT1' .
*/
bio_symbol_map_human( ensg, EnsG, Symb ) :-
     % prefer Hgnc, if it exists
     ( hgnc_homs_ensg_hgnc(EnsG,Hgnc) ->
          hgnc_homs_hgnc_symb( Hgnc, Symb )
          ;
          ( ense_homs_ensg_symb(EnsG,Symb) ->
               true
               ;
               Symb = EnsG
          )
     ).

/* 23.06.09: Ncbi IDs....

?- findall( Ncbi, (hgnc_homs_ncbi_symb(Ncbi,Symb), ncbi_homs_ensg_ncbi(EnsG,Ncbi), \+ bio_analytics:bio_symbol_map_human(ensg, EnsG, Symb2)), Ncbis), length(Ncbis,OnlyInHgnc).
Ncbis = [1590, 2630, 2753, 2910, 3137, 4500, 5369, 5703, 6952|...],
OnlyInHgnc = 428.

?- findall( Ncbi, (hgnc_homs_ncbi_symb(Ncbi,Symb), ncbi_homs_ensg_ncbi(EnsG,Ncbi), \+ bio_analytics:bio_symbol_map_human(ensg, EnsG, Symb)), Ncbis), length(Ncbis,OnlyInHgnc).
Ncbis = [1590, 2630, 2753, 2910, 3137, 4207, 4500, 5369, 5703|...],
OnlyInHgnc = 489.

?- findall( Ncbi, (ncbi_homs_ensg_ncbi(EnsG,Ncbi), bio_analytics:bio_symbol_map_human(ensg, EnsG, Symb), \+ hgnc_homs_ncbi_symb(Ncbi,Symb2)), Ncbis ), length(Ncbis,OnlyinNcbiOrDiff).
Ncbis = [105373378, 101180901, 101927446, 401312, 105369632, 100128966, 101927057, 105369209, 100130283|...],
OnlyinNcbiOrDiff = 174.

?- findall( Ncbi, (ncbi_homs_ensg_ncbi(EnsG,Ncbi), bio_analytics:bio_symbol_map_human(ensg, EnsG, Symb), \+ hgnc_homs_ncbi_symb(Ncbi,Symb)), Ncbis ), length(Ncbis,OnlyinNcbi).
Ncbis = [93655, 105373378, 128854680, 101180901, 104797536, 101927446, 401312, 105369632, 63914|...],
OnlyinNcbi = 235.

?- findall( Symb1-Symb2, (ncbi_homs_ensg_ncbi(EnsG,Ncbi), bio_analytics:bio_symbol_map_human(ensg, EnsG, Symb1),hgnc_homs_ncbi_symb(Ncbi,Symb2),Symb1 \== Symb2), Pairs ), length(Pairs,Diffs).
Pairs = ['ST7'-'ST7-OT3', 'DUSP13B'-'DUSP13A', 'RHOXF1'-'LINC01402', 'SMIM8'-'LINC01590', 'PDE10A'-'LINC00473', 'USP9Y'-'TTTY15', 'ASB3'-'GPR75-ASB3', 'FAM98A'-'RASGRP3-AS1', ... - ...|...],
Diffs = 61.
*/

bio_symbol_map_human( ncbi, Ncbi, Symb ) :-
     % prefer HGNC, there are some clashes...
     ( hgnc_homs_ncbi_symb(Ncbi,Symb) ->
          true
          ;
          ncbi_homs_ensg_ncbi( EnsG, Ncbi ), 
          bio_symbol_map_human( ensg, EnsG, Symb )
     ).
bio_symbol_map_human( hgnc, Hgnc, Symb ) :-
     hgnc_homs_hgnc_ncbi( Hgnc, Symb ).

bio_symbols_map_mouse( ncbi, Ncbi, Symb ) :-
     % fixme: additionals ?
     mgim_musm_mgim_ncbi( Mgim, Ncbi ),
     mgim_musm_mgim_symb( Mgim, Symb ).

bio_symbols_map_pig( ensg, EnsG, Symb ) :-
     % fixme: check if there are alternatives+additionals ?
     ense_suss_ensg_symb( EnsG, Symb ).
 

/** bio_conductor_annot_dbi_org(+Org, -DbiToken, -DbiOrg).

The Annotation DBI token and organism strings for bio_db Org-anism.

==
?- bio_conductor_annot_dbi_org(human, Tkn, Dorg).
Tkn = "Hs",
Dorg = "Homo sapiens".

?- bio_db_organism(gallus,Ggallus), bio_conductor_annot_dbi_org(Ggallus, Tkn, Dorg).
Ggallus = chicken,
Tkn = "Gg",
Dorg = "Gallus gallus".
==

@author nicos angelopoulos
@version  0:1 2023/06/07

*/
bio_conductor_annot_dbi_org(chicken, "Gg", "Gallus gallus").
bio_conductor_annot_dbi_org(  human, "Hs", "Homo sapiens").
bio_conductor_annot_dbi_org(  mouse, "Mm", "Mus musculus").
bio_conductor_annot_dbi_org(    pig, "Ss", "Sus scrofa").

% maybe we should move the pred below to a bio_conductor_annot_dbi.pl source file ?

/** bio_conductor_annot_dbi_org_lib(+DbiToken, +EnsLoaded, -Lib).

From an annotation dbi organism token, construct Dbi library string and possibly load it.

When EnsLoaded is grounded to =|true\= the library is loaded. 

==
?- bio_conductor_annot_dbi_org(pig, PigTkn, PigOrg), 
   bio_conductor_annot_dbi_org_lib(PigTkn, false, PigDbiLib ).
      
PigTkn = "Ss",
PigOrg = "Sus scrofa",
PigDbiLib = "org.Ss.eg.db".
==

@author nicos angelopoulos
@version  0:1 2023/06/07
@see bio_conductor_annot_dbi_org/3
@tbd ensure loaded rather than load everytime (most likely the message says loading only)

*/
bio_conductor_annot_dbi_org_lib( DbiTkn, Load, DbiLib ) :-
     string_concat( "org.", DbiTkn, DbiLibPfx ),
     string_concat( DbiLibPfx, ".eg.db", DbiLib ),
     ( Load == true -> 
               lib( bioc(DbiLib) )      % fixme: ensure loaded, instead
               ;
               true
     ).

/** go_org_symbol(+Org,  +Go, -Symb).

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
@version  0:2 2023/06/09,  added pig, moved to bio_org.pl
*/
go_org_symbol( Alias, Go, Symb ) :-
    bio_db_organism( Alias, Org ),
    go_id( Go, _, Gi ),
    go_org_symbol_1( Org, Gi, Symb ).
    
go_org_symbol_1( gallus, Gi, Symb ) :-
    gont_galg_gont_symb( Gi, _Rel, _Evid, Symb ).
go_org_symbol_1(  human, Gi, Symb ) :-
    gont_homs_gont_symb( Gi, _Evid, Symb ).
go_org_symbol_1(  mouse, Gi, Symb ) :-
    gont_musm_gont_symb( Gi, _Evid, Symb ).
go_org_symbol_1(    pig, Gi, Symb ) :-
    gont_suss_symb_gont( Gi, _Rel, _Evid, Symb ).
