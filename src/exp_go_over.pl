
:- lib(promise(exp_go_over_bioc_deps/0,call(exp_go_over_bioc_deps_load))).
:- lib(stoics_lib:kv_decompose/3).

exp_go_over_bioc_deps_load :-
     lib(suggests(bioc("GOstats"))),
     lib(suggests(bioc("GSEABase"))),  % installed with GOstats
     lib(suggests(bioc("GO.db"))),     % installed with GOstats
     lib(suggests(bioc("Category"))),  % installed with GOstats
     assert(exp_go_over_bioc_deps).

exp_go_exp_id_default( hs, symb ).
exp_go_exp_id_default( gallus, symb).
exp_go_exp_id_default( mouse, symb ).

exp_go_over_defaults( Args, Defs ) :-
    Defs = [
                go('BP'),
                go_over_pv_cut(0.05),
                org(hs),
                org_exp_id(ExpId),
                stem(go_over),
                to_file(false),
                universe(go_exp)
    ],
    ( (memberchk(org(InOrg),Args),ground(InOrg)) -> true; InOrg=hs ),
    bio_db_organism( InOrg, BdOrg ),
    exp_go_exp_id_default( BdOrg, ExpId ).

/** exp_go_over( +CsvF, -GoOver, +Opts ).

Perform gene ontology over-representation analysis.

For experimental data in CsvF select de-regulated genes and on<br>
those perform over representation analysis in gene ontology.<br>
Results in GoOver are either as a values list or a csv file, the name of which is 
assumed to be the ground value of GoOVer.


Opts
  * go(GoSec='BP')
    gene ontology section, in: =|[BP,MF,CC]|=
  * go_over_pv_cut(PvCut=0.05)
    p value filter for the results
  * org(Org=hs)
    one of bio_db_organism/2 first argument values (hs, gallus and mouse for now)
  * org_exp_id(OrgExpId)
    the type of the experimental gene ids for the organism. default depends on Org, but currently all map to symb
  * stem(Stem=false)
    stem for output csv file. when false use basename of CsvF 
  * to_file(ToF=false)
    when GoOver is unbound, this controls whether the output
    goes to a file or a values list 
  * universe(Univ=go_exp)
    Univ in : =|[genome,go,go_exp,experiment]|=
    * genome(Gen)
      all gene identifiers in relevant map predicates
    * go(Go)
      all gene identifies appearing in any ontology
    * go_exp(GoExp)
      gene ontology genes that appear in experiment
    * experiment(Exp)
      all identifiers in the experiment
    * ontology(Onto)
      not implemented yet. like OG above, but only include this Ontology branch

Options are also passed to bio_diffex/4.

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
@version  0.2 2022/12/20,   =|Univ=go|= and =|Org=gallus|=
@see go_over_universe/5

*/
exp_go_over( CsvF, GoOver, Args ) :-
    Self = exp_go_over,
    options_append( Self, Args, Opts ),
    exp_go_over_bioc_deps,  % make sure deps are loaded
    debug_call( exp_go_over, options, Self/Opts ),
    bio_diffex( CsvF, DEPrs, NDEPrs, Opts ),
    options( org(OrgPrv), Opts ),
    bio_db_organism( OrgPrv, Org ),
    kv_decompose( DEPrs, DEGenes, _ ),
    debug_call( exp_go_over, length, de_pairs/DEPrs ),
    bio_list_sort_ne( DEGenes, DEGenesSet ),
    debug_call( exp_go_over, length, de_genes_set/DEGenesSet ),
    bio_list_sort_ne( DEGenes, DEGenesSet ),
    options( org_exp_id(ExpId), Opts ),
    org_go_over_std_gene_ids( Org, ExpId, DEGenesSet, Gids ),
    debug_call( exp_go_over, length, gids/Gids ),
    go_over_frame( Org, goFrameData, GofOrg ),
    goFrame <- 'GOFrame'(goFrameData, organism= +GofOrg),
    goAllFrame <- 'GOAllFrame'(goFrame),
    gsc <- 'GeneSetCollection'(goAllFrame, setType = 'GOCollection'() ),
    genes <- Gids,
    options( universe(UnivOpt), Opts ),
    go_over_universe( UnivOpt, Org, ExpId, Gids, NDEPrs, Univ ),
    debug_call( exp_go_over, length, universe/Univ ),
    options( go(GoAspect), Opts ),
    universe <- Univ,
    universe <- 'as.character'(universe),
    genes <- 'as.character'(genes),
    options( go_over_pv_cut(PvCut), Opts ),
    ( debugging(Self) ->
        <- print(paste("length of universe: ", length(universe)))
        ;
        true
    ),
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

bio_list_sort_ne( List, SetNe ) :-
    sort( List, Set ),
    ( select('',Set,SetNe) -> true; Set = SetNe ).

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

%% go_over_universe( +Token, +Org, +ExpIdTkn, +DEGenes, +NDEPrs, -Universe )
%
% Universe is the list of gene identifiers to be used as universe/background for GOstats.
%
% In human (=|Org=hs|=), this is a list of Entrez ids, and in =|Org=mouse|=, a list of Mgim identifiers.<br>
% DEGenes is a list of deferentially expressed gene identifiers, and NDEPrs is a list of non-differential <br>
% expressed Symbol-Pvalue pairs.
% 
% ExpIdTkn identifies the type of ids coming in in NDEPrs, DEGens, are already in standard form for Org.
go_over_universe( experiment, Org, ExpId, DeGids, NDEPrs, Univ ) :-
    go_over_universe_exp( Org, ExpId, DeGids, NDEPrs, Univ ).
go_over_universe( genome, Org, _ExpId, _DEGids, _NDEPrs, Univ ) :-
    go_over_universe_genome( Org, Univ ).
go_over_universe( go, Org, _ExpId, _DEGids, _NDEPrs, Univ ) :-
    go_over_universe_go( Org, Univ ).
go_over_universe( go_exp, Org, ExpId, DEGids, NDEPrs, Univ ) :-
    go_over_universe_go_exp( Org, ExpId, DEGids, NDEPrs, Univ ).

go_over_universe_genome( gallus, Univ ) :-
    findall( Entz, map_cgnc_gallus_cgnc_entz(_Cgnc,Entz), Entzs ),
    sort( Entzs, Univ ).
go_over_universe_genome( hs, Univ ) :-
    findall( Entz, map_hgnc_symb_entz(_Symb,Entz), Entzs ),
    sort( Entzs, Univ ).
go_over_universe_genome( mouse, Univ ) :-
    findall( Mgim, map_mgim_mouse_mgim_symb(Mgim,_Symb), Mgims ),
    sort( Mgims, Univ ).

go_over_universe_go( gallus, Univ ) :-
    findall( Entz, ( map_gont_gallus_symb_gont(Symb,_Rl,_Ev,_Go),
                     map_cgnc_gallus_cgnc_symb(Cgnc,Symb),
                     map_cgnc_gallus_cgnc_entz(Cgnc,Entz)
                   ),
                    Entzs
           ),
    sort( Entzs, Univ ).
go_over_universe_go( hs, Univ ) :-
    findall( Entz,   (   map_gont_gont_symb(_Go,_En,Symb),
                         map_hgnc_symb_entz(Symb,Entz)
                     ),
                         Entzs
           ),
    sort( Entzs, Univ ).
go_over_universe_go( mouse, Univ ) :-
    findall( Mgim, map_gont_mouse_mgim_gont(Mgim,_E,_G), Mgims ),
    sort( Mgims, Univ ).

% fixme: give doc here, what is this for ?
go_over_universe_exp( gallus, DeGids, NDEPrs, Univ ) :-
    findall( Entz, ( member(Symb-NDEPrs),
                     map_cgnc_gallus_cgnc_symb(Cgnc,Symb),
                     map_cgnc_gallus_cgnc_entz(Cgnc,Entz)
                   ), NDEEntzs ),
    append( DeGids, NDEEntzs, Entzs ),
    sort( Entzs, Univ ).
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

go_over_universe_exp( gallus, ExpIdTkn, DEGids, NDEPrs, Univ ) :-
    findall( Id, member(Id-_,NDEPrs), Ids ),
    org_go_over_std_gene_ids_gallus( ExpIdTkn, Ids, NDGids ),
    append( DEGids, NDGids, ExpGids ),
    % fixme: we can add tests here, that the ids exist in some table ? 
    bio_list_sort_ne( ExpGids, Univ ).
go_over_universe_exp( hs, ExpIdTkn, DEGids, NDEPrs, Univ ) :-
    findall( Id, member(Id-_,NDEPrs), Ids ),
    org_go_over_std_gene_ids_hs( ExpIdTkn, Ids, NDGids ),
    append( DEGids, NDGids, ExpGids ),
    % fixme: we can add tests here, that the ids exist in some table ? 
    bio_list_sort_ne( ExpGids, Univ ).
go_over_universe_exp( mouse, ExpIdTkn, DEGids, NDEPrs, Univ ) :-
    findall( Id, member(Id-_,NDEPrs), Ids ),
    org_go_over_std_gene_ids_hs( ExpIdTkn, Ids, NDGids ),
    append( DEGids, NDGids, ExpGids ),
    % fixme: we can add tests here, that the ids exist in some table ? 
    bio_list_sort_ne( ExpGids, Univ ).

go_over_universe_go_exp( gallus, ExpId, DEGids, NDEPrs, Univ ) :-
    findall( Id, member(Id-_,NDEPrs), Ids ),
    org_go_over_std_gene_ids_gallus( ExpId, Ids, NDGids ),
    append( DEGids, NDGids, ExpGids ),
    findall( ExpGid,( member(ExpGid,ExpGids),   % these are std form, here NCBI Gene ids
                      map_cgnc_gallus_cgnc_entz(Cgnc,ExpGid),
                      map_cgnc_gallus_cgnc_symb(Cgnc,Symb),
                      once(map_gont_gallus_symb_gont(Symb,_,_,_))
                    ),
                         List ),
    bio_list_sort_ne( List, Univ ).
go_over_universe_go_exp( hs, ExpIdTkn, DEGids, NDEPrs, Univ ) :-
    findall( Id, member(Id-_,NDEPrs), Ids ),
    org_go_over_std_gene_ids_hs( ExpIdTkn, Ids, NDGids ),
    append( DEGids, NDGids, ExpGids ),
    findall( Entz,  ( member(Entz,ExpGids),
                      map_hgnc_symb_entz(Symb,Entz),
                      once(map_gont_symb(_,_,Symb))
                    ),
                         List ),
    bio_list_sort_ne( List, Univ ).
go_over_universe_go_exp( mouse, ExpIdTkn, DEGids, NDEPrs, Univ ) :-
    findall( Id, member(Id-_,NDEPrs), Ids ),
    org_go_over_std_gene_ids_hs( ExpIdTkn, Ids, NDGids ),
    append( DEGids, NDGids, ExpGids ),
    findall( Mgim,  ( member(Mgim,ExpGids),
                      map_mgim_mouse_mgim_symb(Mgim,Symb),
                      once(map_gont_mouse_gont_symb(_,_,Symb))
                    ), 
                         List ),
    bio_list_sort_ne( List, Univ ).

org_go_over_std_gene_ids( Org, Gtyp, Set, Gids ) :-
     at_con( [org_go_over_std_gene_ids,Org], '_', Pname ),
     ( call(Pname,Gtyp,Set,Gids) ->
          true
          ;
          throw( go_over(could_not_convert(org(Org),gene_id_type(Gtyp))) )
     ).

% it should for ensg and symb
org_go_over_std_gene_ids_gallus( entz, Set, Gids ) :-
     !,  % not really needes as it called from an if() above
     Set = Gids.
org_go_over_std_gene_ids_gallus( cgnc, Set, Gids ) :-
     !,
    findall( Entz,  (member(Cgnc,Set),map_cgnc_gallus_cgnc_entz(Cgnc,Entz)), Entzs ),
    sort( Entzs, Gids ).
org_go_over_std_gene_ids_gallus( SrcT, Set, Gids ) :-
    atom_concat( map_cgnc_gallus_cgnc_, SrcT, SrcNm ),
    findall( Entz,  (member(SrcG,Set),call(SrcNm,Cgnc,SrcG),map_cgnc_gallus_cgnc_entz(Cgnc,Entz)), Entzs ),
    sort( Entzs, Gids ).
% fixme: add more rules... for hs and mouse
org_go_over_std_gene_ids_hs( entz, Set, Gids ) :-
    sort( Set, Gids ).
org_go_over_std_gene_ids_hs( symb, Set, Gids ) :-
    findall( Entz,  (member(Symb,Set),map_hgnc_symb_entz(Symb,Entz)), Entzs ),
    sort( Entzs, Gids ).
org_go_over_std_gene_ids_mouse( entz, Set, Gids ) :-
    findall( Mgim,  (member(Entz,Set),map_mgim_mouse_mgim_entz(Mgim,Entz)), Mgims ),
    sort( Mgims, Gids ).
org_go_over_std_gene_ids_mouse( symb, Set, Gids ) :-
    findall( Mgim,  (member(Symb,Set),map_mgim_mouse_mgim_symb(Mgim,Symb)), Mgims ),
    sort( Mgims, Gids ).

go_over_frame( gallus, GoFra, GofOrg ) :-
     !, 
     findall( row(Gid,E,Entz), (
                                     map_gont_gallus_symb_gont(Symb,_,E,G),
                                     go_id(Gid,G),
                                     map_cgnc_gallus_cgnc_symb(Cgnc,Symb),
                                     map_cgnc_gallus_cgnc_entz(Cgnc,Entz)
                               ),
                              Rows
            ),
    go_mtx_df( [row(go_id,evidence,gene_id)|Rows], GoFra, [] ),
    GofOrg = "Gallus gallus".
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
go_over_frame( mouse, GoFra, GofOrg ) :-
    !,
    findall( row(Gid,E,M), 
            ( map_gont_mouse_mgim_gont(M,E,G),
              go_id(Gid,G)
            ),
                Rows ),
    go_mtx_df( [row(go_id,evidence,gene_id)|Rows], GoFra, [] ),
    GofOrg = "Mus musculus".

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
