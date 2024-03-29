gene_family_defaults([org(hs)]).

%% gene_family( +Family, -Symbols, +Opts ).
%
% If Family is a known family alias it is expanded to a list of its constituent gene symbols.
% 
% Family can be a gene ontology term (atom of the form GO:XXXXXX).
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
% Opts (v0.2)
%  * org(Org=hs)
%    organism for the symbols
%
% Opts are also passed to: go_symbols_reach/3.
%
%==
% Located family as the bio_analytics gene collection: autophagy
% ?- gene_family(autophagy, Auto, org(human)), length(Auto, Len).
% Auto = ['AMBRA1', 'APOL1', 'ARNT', 'ARSA', 'ARSB', 'ATF4', 'ATF6', 'ATG10', 'ATG12'|...],
% Len = 232.
%
% ?- debug(gene_family).
% ?- gene_family(375, Symbs, []), length(Symbs, Len).
% % Located GO term as the family identifier.
% Symbs = ['BCAS2', 'DBR1', 'DDX23', 'GEMIN2', 'KHSRP', 'LSM1', 'MPHOSPH10', 'PRPF3', 'PRPF4'|...],
% Len = 26.
% 
% ?- gene_family('GO:0000375', Symbs, []), length(Symbs, Len).
% % Located GO term as the family identifier.
% Symbs = ['BCAS2', 'DBR1', 'DDX23', 'GEMIN2', 'KHSRP', 'LSM1', 'MPHOSPH10', 'PRPF3', 'PRPF4'|...],
% Len = 26.
%
% ?- gene_family([55626, 8542, 405], Auto, []).
% Converted input from Entrezes to Symbols.
% Auto = ['AMBRA1', 'APOL1', 'ARNT'].
%
% ?- gene_family(unknown, Auto, org(hs)).
% ERROR: Unhandled exception: gene_family(cannot_find_input_family_in_the_known_ones(unknown,[autophagy]))
%==
% Family datasets are in pack(bio_analytics/data/families).
%
% @author nicos angelopoulos
% @version  0:1 2019/3/5,  from old code
% @version  0:2 2023/6/7,  added options
% @see go_symbols_reach/4
% @tbd error reporting via print_message/2
%
gene_family( Family, Symbols, Args ) :-
     Self = gene_family,
     options_append( Self, Args, Opts ),
     options( org(OrgPrv), Opts ),
     ( bio_db_organism(OrgPrv,Org) -> true; Org = OrgPrv ),
     gene_family( Family, Org, Symbols, Opts ).

gene_family( Family, Org, Symbols, _Opts ) :-
     gene_family_known( Family, Org, Symbols ),
     !,
     debug( gene_family, 'Located family as the bio_analytics gene collection: ~w', [Family] ).
gene_family( GoTerm, _Org, Symbols, Opts ) :-
     atomic( GoTerm ),
     ( integer(GoTerm) -> GoTerm = UpGoTerm; upcase_atom(GoTerm,UpGoTerm) ),
     go_id( UpGoTerm, _GoAtom, GoInt ),
     % atom_concat( 'GO:', GoAtmInt, UpGOTerm ),
     % atom_number( GoAtmInt, GoInt ),
     !,
     % gene_family_org_go_term( Org, GoInt, Symbols ),
     go_symbols_reach( GoInt, Symbols, Opts ),
     debug( gene_family, 'Located GO term as the family identifier.', [] ).
gene_family( Family, Org, Symbols, _Opts ) :-
     is_list( Family ),
     maplist( is_symbol(Org), Family ),
     % maplist( hgnc_homs_hgnc_symb, _, Family ),
     !,
     Symbols = Family,
     debug( gene_family, 'Ascertained input as list of symbols- passing it through.', [] ).
gene_family( Family, Org, Symbols, _Opts ) :-
     is_list( Family ),
     maplist( entz_symb(Org), Family, Symbols ),
     !,
     debug( gene_family, 'Converted input from Entrezes to Symbols.', [] ).
gene_family( Family, _Org, _Symbols, _Opts ) :-
     is_list( Family ),
     !,
     length( Family, Len ),
     throw( gene_family(family_inpput_not_symbols_or_entrez_ids(length(Len))) ).
gene_family( Family, Org, _Symbols, _Opts ) :-
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
