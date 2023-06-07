
go_symbols_reach_defaults( Defs ) :-
	Defs = [ descent(true), 
             org(hs),
	         as_child_includes(true),
		     as_child_consists_of(true),
		     as_child_regulates(false),
		     as_child_positively_regulates(false),
		     as_child_negatively_regulates(false)
		  ].

/** go_symbols_reach( +GoT, -Symbols, +Opts ).

Gets the symbols belonging to a GO term. Descents to GO child relations, 
which by default are includes (reverse of is_a) and consists_of (reverse of part_of)
to pick up Symbols recursively.

Opts 
  * org(Org=hs)
    should be recognised by 1st arg of bio_db_organism/2.
  * go(GO=lib_org)
    which implementation of GO to follow (v0.3 this changes default behaviour). 
    Some of the options below are not functional with =|GO=lib_org|=.
    Alterantive: =\bio_db\= (which used to be the default, before the option was introduced).

  * descent(Desc=true)
    whether to collect symbols from descendant GO terms
  * as_child_includes(Inc=true)
    collect from edge_gont_include/2
  * as_child_consists_of(Cns=true)
    collect from edge_gont_consists_of/2
  * as_child_regulates(Reg=false)
    collect from edge_gont_regulates/2
  * as_child_negatively_regulates(Reg=false)
    collect from edge_gont_negatively_regulates/2
  * as_child_positively_regulates(Reg=false)
    collect from edge_gont_positively_regulates/2
  * debug(Dbg=false)
    see options_append/3

Listens to debug(go_symbols_reach).

==
?- go_symbols_reach( 'GO:0000375', Symbs, [] ), length( Symbs, Len ).
Symbs = ['AAR2', 'ALYREF', 'AQR', 'ARC', 'BCAS2', 'BUD13', 'BUD31', 'CACTIN', 'CASC3'|...],
Len = 293.

?- go_symbols_reach( 'GO:0000375', Symbs, org(mouse) ), length( Symbs, Len ).
Symbs = ['4930595M18Rik', 'Aar2', 'Aqr', 'Bud13', 'Bud31', 'Casc3', 'Cdc40', 'Cdc5l', 'Cdk13'|...].
Len = 190.

?- go_symbols_reach( 375, Symbs, true ), length( Symbs, Len ).
Symbs = ['AAR2', 'ALYREF', 'AQR', 'ARC', 'BCAS2', 'BUD13', 'BUD31', 'CACTIN', 'CASC3'|...],
Len = 293.
==

@author nicos angelopoulos
@version  0.1 2015/7/26
@version  0.2 2019/4/7,          added org, moved to new pack
@version  0.3 2023/6/7,          added option go(GoImpl), which changes default behaviour in comparison to past

*/
go_symbols_reach( GO, Symbs, Args ) :-
	options_append( go_symbols_reach, Args, Opts ),
	options( descent(Desc), Opts ),
	options( org(Org), Opts ),
	go_symbs_descent_term( Desc, GoT, Child, Term, Opts ),
	go_symbs( [GO], Org, GoT, Child, Term, [], [], Symbs ).

go_symbs( [], _Org, _GoT, _Ch, _Term, _, Symbs, Symbs ).
go_symbs( [GoIn|GOs], Org, GoT, ChGo, FTerm, GoSeen, Seen, Symbs ) :-
    go_id( GoIn, _GoAtm, GoId ),
	ord_add_element( GoSeen, GoId, NxGoSeen ),
	% findall( GoSymb, gont_homs_gont_symb(GO,GoSymb), GoSymbs ),
	findall( GoSymb, go_org_symbol(Org,GoId,GoSymb), GoSymbs ),
	sort( GoSymbs, OSymbs ),
	debug( go_symbols_reach, '~w, Symbols: ~w', [GoIn,OSymbs] ),
	ord_union( OSymbs, Seen, NxSeen ),
    Dbg = debug( go_symbols_reach, 'calling, ~w', [FTerm] ),
	findall( ChGo, (GoT=GoId,Dbg,FTerm), ChGos ),
	sort( ChGos, ChGosOrd ),
	debug( go_symbols_reach, '~w, Children: -~w', [GoIn,ChGosOrd] ),
	ord_subtract( ChGosOrd, GoSeen, ChGosAdd ),
	ord_union( GOs, ChGosAdd, NxGos ),
	go_symbs( NxGos, Org, GoT, ChGo, FTerm, NxGoSeen, NxSeen, Symbs ).

go_symbs_descent_term( false, _, false, false, _ ).
go_symbs_descent_term( true, GoT, Child, Term, Opts ) :-
	Rships = [ includes, consists_of, 
                regulates, positively_regulates, negatively_regulates
			 ],
	findall( Pname, (
				        member(Rship,Rships),
				        atom_concat( as_child_, Rship, Oname ),
				        Opt =.. [Oname,true],
				        options(Opt,Opts),
				        atom_concat(edge_gont_,Rship,Pname)
	                ),
				      Pnames ),
	go_symbs_descent_disjunction( Pnames, GoT, Child, Term ).

go_symbs_descent_disjunction( [], _, false, false ).
go_symbs_descent_disjunction( [H|T], Term, Child, Goal ):-
	G =.. [H,Term,Child],
	go_symbs_descent_disjunction_1( T, Term, Child, G, Goal ).


go_symbs_descent_disjunction_1( [], _Term, _Child, Goal,  Goal ).
go_symbs_descent_disjunction_1( [H|T], Term, Child, Left, Goal ) :-
	G =.. [H,Term,Child],
	go_symbs_descent_disjunction_1( T, Term, Child, (Left;G), Goal ).

% bio_db script no longer produce this: 22.12.21
edge_gont_includes( A, B ) :-
     bio_db:edge_gont_is_a( B, A ).

edge_gont_consists_of( A, B ) :-
     bio_db:edge_gont_part_of( B, A ).
