
/** bio_org.

Bio analytics supports a number of organism specific internal predicates.

This is a doc predicate listing them:
  * bio_conductor_annot_dbi_org/3
  * bio_symbols_map/4

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


/** bio_symbols_map( +Org, +IdTkn, +Values, -Symbols ).

Convert Org-anism specific values of IdTkn ids to Symbols.

==
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
bio_symbols_map( gallus, ExpId, Values, Symbs ) :-
     bio_symbols_map_gallus( ExpId, Values, Symbs ).
bio_symbols_map( hs, ExpId, Values, Symbs ) :-
     bio_symbols_map_hs( ExpId, Values, Symbs ).
bio_symbols_map( mouse, ExpId, Values, Symbs ) :-
     bio_symbols_map_mouse( ExpId, Values, Symbs ).
bio_symbols_map( pig, ExpId, Values, Symbs ) :-
     bio_symbols_map_pig( ExpId, Values, Symbs ).

bio_symbols_map_gallus( ncbi, Values, Symbs ) :-
     % fixme: additionals ?
     findall( SymbTerm, ( 
                      member(Valu,Values),
                      ( Valu = Ncbi-V -> SymbTerm=Symb-V; Valu = Ncbi, SymbTerm=Symb ),
                      cgnc_galg_cgnc_ncbi(Cgnc,Ncbi),
                      cgnc_galg_cgnc_symb(Cgnc,Symb)
                    ),
                         Symbs ).

bio_symbols_map_hs( ncbi, Values, Symbs ) :-
     % fixme: additionals ?
     findall( SymbTerm, ( member(Valu,Values),
                      ( Valu = Ncbi-V -> SymbTerm=Symb-V; Valu = Ncbi, SymbTerm=Symb ),
                      hgnc_homs_hgnc_ncbi(Cgnc,Ncbi),
                      hgnc_homs_hgnc_symb(Cgnc,Symb)
                    ),
                         Symbs ).
bio_symbols_map_mouse( ncbi, Values, Symbs ) :-
     % fixme: additionals ?
     findall( SymbTerm, ( member(Valu,Values),
                      ( Valu = Ncbi-V -> SymbTerm=Symb-V; Valu = Ncbi, SymbTerm=Symb ),
                      mgim_musm_mgim_ncbi(Cgnc,Ncbi),
                      mgim_musm_mgim_symb(Cgnc,Symb)
                    ),
                         Symbs ).
bio_symbols_map_pig( ensg, EnsGs, Symbs ) :-
     % fixme: check if there are alternatives+additionals ?
     % fixme: can we remove the -KV malarky from all 4 organisms- see upstream ?
     findall( SymbTerm, ( member(Valu,EnsGs),
                          (Valu = EnsG-V -> SymbTerm=Symb-V; Valu=EnsG, SymbTerm=Symb ),
                          ense_suss_ensg_symb(EnsG,Symb)
                        ), Symbs ).
 

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


