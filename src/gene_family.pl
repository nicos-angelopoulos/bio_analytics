
%% gene_family( +Family, +Org, -Symbols ).
%
% If Family is a known family alias for organism Org, it is expanded to a list of its constituent gene symbols.
% 
% Currently 2 organisms are supported: human and mouse.<br>
% Family can be a gene ontology term (atom of the form GO:XXXXXX).<br>
% If family is a list of numbers or atoms that map to numbers, then they are taken to be Entrez ids which
% are converted to gene symbols in Symbols.
%
% If family is a list of symbols, it is passed on to Symbols.
% 
% Listens to debug(gene_family).
%
% Known families:
% * autophagy/human
%  from: http://autophagy.lu/clustering/index.html
% 
%==
% Located family as the bio_analytics gene collection: autophagy
% ?- gene_family( autophagy, human, Auto ), length( Auto, Len ).
% Auto = ['AMBRA1', 'APOL1', 'ARNT', 'ARSA', 'ARSB', 'ATF4', 'ATF6', 'ATG10', 'ATG12'|...],
% Len = 232.
%
% ?- debug( gene_family ).
% ?- gene_family( 375, hs, Symbs ), length( Symbs, Len ).
% % Located GO term as the family identifier.
% Symbs = ['BCAS2', 'DBR1', 'DDX23', 'GEMIN2', 'KHSRP', 'LSM1', 'MPHOSPH10', 'PRPF3', 'PRPF4'|...],
% Len = 26.
% 
% ?- gene_family( 'GO:0000375', hs, Symbs ), length( Symbs, Len ).
% % Located GO term as the family identifier.
% Symbs = ['BCAS2', 'DBR1', 'DDX23', 'GEMIN2', 'KHSRP', 'LSM1', 'MPHOSPH10', 'PRPF3', 'PRPF4'|...],
% Len = 26.
%
% ?- gene_family( [55626, 8542, 405], hs, Auto ).
% Converted input from Entrezes to Symbols.
% Auto = ['AMBRA1', 'APOL1', 'ARNT'].
%
% ?- gene_family( unknown, hs, Auto ).
% ERROR: Unhandled exception: gene_family(cannot_find_input_family_in_the_known_ones(unknown,[autophagy]))
%==
% Family datasets are in pack(bio_analytics/data/families).
%
% @author nicos angelopoulos
% @version  0:1 2019/3/5,  from old code
% @tbd error reporting via print_message/2
%
gene_family( Family, OrgPrv, Symbols ) :-
    ( bio_db_organism(OrgPrv,Org) -> true; Org = OrgPrv ),
	gene_family_known( Family, Org, Symbols ),
	!,
	debug( gene_family, 'Located family as the bio_analytics gene collection: ~w', [Family] ).
gene_family( GoTerm, Org, Symbols ) :-
	atomic( GoTerm ),
	( integer(GoTerm) -> GoTerm = UpGoTerm; upcase_atom(GoTerm,UpGoTerm) ),
    go_id( UpGoTerm, _GoAtom, GoInt ),
	% atom_concat( 'GO:', GoAtmInt, UpGOTerm ),
    % atom_number( GoAtmInt, GoInt ),
	!,
    % gene_family_org_go_term( Org, GoInt, Symbols ),
    go_symbols_reach( GoInt, Symbols, [organism(Org)] ),
	debug( gene_family, 'Located GO term as the family identifier.', [] ).
gene_family( Family, Org, Symbols ) :-
	is_list( Family ),
    maplist( is_symbol(Org), Family ),
	% maplist( hgnc_homs_hgnc_symb, _, Family ),
	!,
	Symbols = Family,
	debug( gene_family, 'Assertained input as list of symbols- passing it through.', [] ).
gene_family( Family, Org, Symbols ) :-
	is_list( Family ),
    maplist( entz_symb(Org), Family, Symbols ),
	!,
	debug( gene_family, 'Converted input from Entrezes to Symbols.', [] ).
gene_family( Family, _Org, _Symbols ) :-
	is_list( Family ),
	!,
	length( Family, Len ),
	throw( gene_family(family_inpput_not_symbols_or_entrez_ids(length(Len))) ).
gene_family( Family, Org, _Symbols ) :-
	findall( Fam, gene_family_known(Fam,Org,_), Families ),
	throw( gene_family(cannot_find_input_family_in_the_known_ones(Family,Org,Families)) ).

/*
gene_family_org_go_term( hs, GoInt, Symbols ) :-
    gene_family_org_go_term( human, GoInt, Symbols ).
gene_family_org_go_term( human, GoInt, Symbols ) :-
	findall( Symb, gont_homs_gont_symb(GoInt,_Evid,Symb), Symbols ).
gene_family_org_go_term( mouse, GoInt, Symbols ) :-
	findall( Symb, gont_musm_gont_symb(GoInt,_Evid,Symb), Symbols ).
    */

gene_family_known( autophagy, hs, Auto ) :-
	mtx_column( pack('bio_analytics/data/families/autophagy'), 'Symbol', Auto ).
