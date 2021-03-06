
/** go_org_symbols( +GoT, +Org, -Symbols ).

Gets the symbols belonging to a GO term.

==
?- go_org_symbols( 'GO:0000375', human, Symbs ), length( Symbs, Len ).
Symbs = ['BCAS2', 'DBR1', 'DDX23', 'GEMIN2', 'KHSRP'|...],
Len = 26.

?- go_org_symbols( 'GO:0000375', mouse, Symbs ), length( Symbs, Len ).
Symbs = ['Dbr1', 'Rbm17', 'Rnu4atac', 'Rnu6atac', 'Scaf11', 'Sf3a2', 'Slu7', 'Srrm1', 'Srsf10'|...],
Len = 10.

?- go_org_symbols( 375, human, Symbs ), length( Symbs, Len ).
Symbs = ['BCAS2', 'DBR1', 'DDX23', 'GEMIN2', 'KHSRP', 'LSM1', 'MPHOSPH10', 'PRPF3', 'PRPF4'|...],
Len = 26.
==

@author nicos angelopoulos
@version  0.1 2019/4/7
@see go_org_symb/3, just a findall of
@see go_symbols_reach/3 for a version that travells the ontology

*/
go_org_symbols( Go, Org, Symbs ) :-
	findall( GoSymb, go_org_symbol(Org,Go,GoSymb), Symbs ).
