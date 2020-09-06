
:- lib(suggests(bioc("GOstats"))).
:- lib(suggests(bioc("GSEABase"))).  % installed with GOstats
:- lib(suggests(bioc("GO.db"))).     % installed with GOstats
:- lib(suggests(bioc("Category"))).  % installed with GOstats

:- lib(stoics_lib:kv_decompose/3).

exp_go_over_defaults( Defs ) :-
    Defs = [
                go('BP'),
                go_over_pv_cut(0.05),
                org(hs),
                stem(go_over),
                to_file(false),
                universe(experiment)
    ].

/** exp_go_over( +CsvF, -GoOver, +Opts ).

Perform gene ontology over-representation analysis.

For experimental data in CsvF select de-regulated genes and on<br>
those perform over representation analysis in gene ontology.<br>
Results in GoOver are either as a values list of a csv file.

Opts
  * go('BP')
    gene ontology section (BP,MF,CC)
  * go_over_pv_cut(0.05)
    p value filter for the results
  * org(hs)
    one of bio_db_organism/2 first argument values (hs or mouse for now)
  * stem(Stem=false)
    stem for output csv file. when false use basename of CsvF 
  * to_file(ToF=false)
    when GoOver is unbound, this controls whether the output
    goes to a file or a values list 
  * universe(experiment)
    or genome

Options are also passed to exp_diffex/4.

==
?- absolute_file_name( pack('bio_analytics/data/silac/bt.csv'), CsvF ),
   exp_go_over( CsvF, GoOver, [] ).

CsvF = '.../swipl/pack/bio_analytics/data/silac/bt.csv',
GoOver = [row('GOBPID', 'Pvalue', adj.pvalue, ...), row('GO:0061387', 0.000187, ...)|...].

?- absolute_file_name( pack('bio_analytics/data/silac/bt.csv'), CsvF ),
   exp_go_over( CsvF, GoOver, [stem(here),to_file(true)] ).

CsvF = '.../swipl/pack/bio_analytics/data/silac/bt.csv',
GoOver = here_gontBP_p0.05_univExp.csv.

?- absolute_file_name( pack('bio_analytics/data/silac/bt.csv'), CsvF ),
   exp_go_over( CsvF, 'a_file.csv', [stem(here),to_file(true)] ).

CsvF = '.../swipl/pack/bio_analytics/data/silac/bt.csv',

?- lib(by_unix).
?- @wc(-l,'a_file.csv').
183 a_file.csv

?- @wc(-l,'here_gontBP_p0.05_univExp.csv').
183 here_gontBP_p0.05_univExp.csv
true.

?- absolute_file_name( pack('bio_analytics/data/silac/bt.csv'), CsvF ),                                                  exp_go_over( CsvF, OverF, [to_file(true)] ).
CsvF = '.../swipl/pack/bio_analytics/data/silac/bt.csv',
OverF = go_over_gontBP_p0.05_univExp.csv.

?- absolute_file_name( pack('bio_analytics/data/silac/bt.csv'), CsvF ),                                                  exp_go_over( CsvF, OverF, [to_file(true),stem(false)] ).
CsvF = '.../swipl/pack/bio_analytics/data/silac/bt.csv',
OverF = '.../swipl/pack/bio_analytics/data/silac/bt_gontBP_p0.05_univExp.csv'.
==

@author nicos angelopoulos
@version  0.1 2019/5/2
*/
exp_go_over( CsvF, GoOver, Args ) :-
    Self = exp_go_over,
    options_append( Self, Args, Opts ),
    exp_diffex( CsvF, DEPrs, NDEPrs, Opts ),
    options( org(OrgPrv), Opts ),
    bio_db_organism( OrgPrv, Org ),
    kv_decompose( DEPrs, DEGenes, _ ),
    debug_call( exp_go_over, length, de_pairs/DEPrs ),
    sort( DEGenes, DEGenesSet ),
    debug_call( exp_go_over, length, de_genes_set/DEGenesSet ),
    org_symb_go_over_gene_ids( Org, DEGenesSet, Gids ),
    debug_call( exp_go_over, length, gids/Gids ),
    go_over_frame( Org, goFrameData, GofOrg ),
    goFrame <- 'GOFrame'(goFrameData, organism= +GofOrg),
    goAllFrame <- 'GOAllFrame'(goFrame),
    gsc <- 'GeneSetCollection'(goAllFrame, setType = 'GOCollection'() ),
    genes <- Gids,
    options( universe(UnivOpt), Opts ),
    go_over_universe( UnivOpt, Org, Gids, NDEPrs, Univ ),
    debug_call( exp_go_over, length, universe/Univ ),
    options( go(GoAspect), Opts ),
    universe <- Univ,
    universe <- 'as.character'(universe),
    genes <- 'as.character'(genes),
    options( go_over_pv_cut(PvCut), Opts ),
    gGparams <- 'GSEAGOHyperGParams'(
                    name="bio_analytics_gont",
                    geneSetCollection=gsc,
                    geneIds = genes,
                    universeGeneIds = universe,
                    ontology = + GoAspect,
                    pvalueCutoff = PvCut,
                    conditional = 'FALSE',
                    testDirection = "over"
                ),
    over <- hyperGTest(gGparams),
    dfOver <- 'data.frame'( summary(over) ),
    dfOver$'adj.pvalue' <- 'p.adjust'( dfOver$'Pvalue', method=+'BH' ),
    dfOver <- dfOver[*,c(1,2,8,3,4,5,6,7)],
    Use = use(Self,GoAspect,UnivOpt,PvCut),
    exp_go_over_return( GoOver, dfOver, CsvF, Use, Opts ).

exp_go_over_return( GoOver, DfOveR, _CsvF, _Use, Opts ) :-
    var( GoOver ),
    options( to_file(false), Opts ),
    !,
    Rlist <- 'as.list'(DfOveR),
    findall( List, (member(Head=Tail,Rlist),List=[Head|Tail]), Lists),
	mtx_lists( GoOver, Lists ).

exp_go_over_return( GoOver, DfOveR, CsvF, Use, Opts ) :-
    Use = use(Self,GoAspect,UnivOpt,PvCut),
    SubSep = '',
    options( stem(Stem), Opts ),
    atomic_list_concat( [gont,GoAspect], SubSep, GontAspTkn ),
    capit( UnivOpt, ShUniv ),
    atomic_list_concat( [univ,ShUniv], SubSep, UnivTkn ),
    (Stem == false -> TempF = CsvF; os_ext(csv,Stem,TempF)),
    atom_concat( p, PvCut, PvTkn ),
    Postfixes = [GontAspTkn,PvTkn,UnivTkn],
    (var(GoOver) -> os_postfix(Postfixes, TempF, GoOver) ; true),
    <- 'write.csv'(DfOveR, file=+GoOver, 'row.names'='FALSE'),
    % <- print( warnings() ),
    debug( Self, 'Wrote: ~p', GoOver ).

% go_over_universe( experiment, Org, DEGenes, NDEPrs, Univ ) :-
% we now are passing DeGids, so no need to re-find those
go_over_universe( experiment, Org, DeGids, NDEPrs, Univ ) :-
    go_over_universe_exp( Org, DeGids, NDEPrs, Univ ).
go_over_universe( genome, Org, _DEGids, _NDEPrs, Univ ) :-
    go_over_universe_genome( Org, Univ ).

go_over_universe_genome( hs, Univ ) :-
    findall( Entz, map_hgnc_symb_entz(_Symb,Entz), Entzs ),
    sort( Entzs, Univ ).
go_over_universe_genome( mouse, Univ ) :-
    findall( Mgim, map_mgim_mouse_mgim_symb(Mgim,_Symb), Mgims ),
    sort( Mgims, Univ ).

go_over_universe_exp( hs, DeGids, NDEPrs, Univ ) :-
    findall( Entz, (member(Symb-_,NDEPrs),map_hgnc_symb_entz(Symb,Entz)), NDEEntzs ),
    % findall( Entz1, (member(Symb,DEGenes),map_hgnc_symb_entz(Symb,Entz1)), DEEntzs ),
    % append( DEEntzs, NDEEntzs, Entzs ),
    append( DeGids, NDEEntzs, Entzs ),
    sort( Entzs, Univ ).
go_over_universe_exp( mouse, DeGids, NDEPrs, Univ ) :-
    findall( Mgim, (member(Symb-_,NDEPrs),map_mgim_mouse_mgim_symb(Mgim,Symb)), NDEMgims ),
    append( DeGids, NDEMgims, Mgims ),
    sort( Mgims, Univ ).

org_symb_go_over_gene_ids( hs, Set, Gids ) :-
    findall( Entz,  (member(Symb,Set),map_hgnc_symb_entz(Symb,Entz)), Entzs ),
    sort( Entzs, Gids ).
org_symb_go_over_gene_ids( mouse, Set, Gids ) :-
    findall( Mgim,  (member(Symb,Set),map_mgim_mouse_mgim_symb(Mgim,Symb)), Mgims ),
    sort( Mgims, Gids ).

go_over_frame( mouse, GoFra, GofOrg ) :-
    !,
    findall( row(Gid,E,M), 
            ( map_gont_mouse_mgim_gont(M,E,G),
              go_id(Gid,G)
            ),
                Rows ),
    go_mtx_df( [row(go_id,evidence,gene_id)|Rows], GoFra, [] ),
    GofOrg = "Mus musculus".
go_over_frame( hs, GoFra, GofOrg ) :-
    !,
    findall( row(Gid,E,Entz), 
            ( map_gont_gont_symb(G,E,S),
              go_id(Gid,G),
              map_hgnc_symb_entz(S,Entz)
            ),
                Rows ),
    go_mtx_df( [row(go_id,'Evidence',gene_id)|Rows], GoFra, [] ),
    GofOrg = "Homo sapiens".

capit( Atom, Capit ) :-
    atom_codes( Atom, [A,B,C|_] ),
    atom_codes( Aat, [A] ),
    atom_codes( BCat, [B,C] ),
    upcase_atom( Aat, CapA ),
    downcase_atom( BCat, DwnBC ),
    atom_concat( CapA, DwnBC, Capit ).

go_mtx_df_defaults([check_names(false)]).

go_mtx_df( Csv, Df, Args ) :-
	options_append( go_mtx_df, Args, Opts ),
	mtx_lists( Csv, Lists ),
	findall( Head=Tail, ( member(List,Lists), List=[Head|Tail] ), Pairs ),
	Df <- Pairs,
	( options(check_names(true),Opts) -> Check = 'T'; Check = 'F' ),
	Df <- 'data.frame'(Df,'check.names'=Check).
