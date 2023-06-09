
/** go_org_symbols( +GoT, +Org, -Symbols ).

Gets the symbols belonging to a GO term.

==
?- go_org_symbols( 'GO:0000375', human, Symbs ), length( Symbs, Len ).
Symbs = ['SCAF11', 'SLU7', 'SRSF10', 'WDR83', 'DBR1', 'MPHOSPH10', 'PRPF4', 'PRPF6', 'SF3B1'|...],
Len = 25.

?- go_org_symbols( 'GO:0000375', mouse, Symbs ), length( Symbs, Len ).
Symbs = ['Rnu4atac', 'Rnu6atac', 'Srrm1', 'Sf3a2', 'Srsf10', 'Dbr1', 'Scaf11', 'Slu7', 'Srsf10'|...],
Len = 11.

?- go_org_symbols( 375, human, Symbs ), length( Symbs, Len ).
Symbs = ['SCAF11', 'SLU7', 'SRSF10', 'WDR83', 'DBR1', 'MPHOSPH10', 'PRPF4', 'PRPF6', 'SF3B1'|...],
Len = 25.

?- go_org_symbols( 375, chicken, Symbs ), length( Symbs, Len ).
Symbs = ['DBR1'],
Len = 1.

?- go_org_symbols( 375, pig, Symbs ), length( Symbs, Len ).
==

@author nicos angelopoulos
@version  0.1 2019/4/7
@see go_org_symb/3, just a findall of
@see go_symbols_reach/3 for a version that travells the ontology

*/
go_org_symbols( Go, Org, Symbs ) :-
	findall( GoSymb, go_org_symbol(Org,Go,GoSymb), Symbs ).
