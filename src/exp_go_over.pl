
:- lib(promise(exp_go_over_bioc_deps/0,call(exp_go_over_bioc_deps_load))).
:- lib(stoics_lib:break_nth/4).
:- lib(stoics_lib:kv_decompose/3).

:- lib(bio_gid_std/3).
:- lib(org_gid_map/3).
:- lib(bio_list_sort_ne/2).


exp_go_over_bioc_deps_load :-
     lib(suggests(bioc("GOstats"))),
     lib(suggests(bioc("GSEABase"))),  % installed with GOstats
     lib(suggests(bioc("GO.db"))),     % installed with GOstats
     lib(suggests(bioc("Category"))),  % installed with GOstats
     assert(exp_go_over_bioc_deps).

exp_go_gid_default(chicken, symb, ncbi).
exp_go_gid_default(  human, symb, ncbi).
exp_go_gid_default(  mouse, symb, mgim).
exp_go_gid_default(    pig, ensg, ncbi).

exp_go_over_defaults( Args, Defs ) :-
    Defs = [
                go('BP'),
                go_frame(bioc_ann_dbi),
                go_over_pv_cut(0.05),
                org(hs),
                % org_exp_id(ExpId),
                gid(Gid),
                gid_to(Gto),
                min_count(0),
                stem(go_over),
                symbs(false),
                symbs_hdr("symbols"),
                symbs_max(inf),
                to_file(false),
                universe(go_exp)
    ],
    ( (memberchk(org(InOrg),Args),ground(InOrg)) -> true; InOrg=hs ),
    bio_db_organism( InOrg, BdOrg ),
    exp_go_gid_default( BdOrg, Gid, Gto ),
    options_return( gid_to(Gto), Args, [pack(bio_analytics),pred(ex_go_over/3),option(gid_to(Gto))] ).

/** exp_go_over( +Mtx, -GoOver, +Opts ).

Perform gene ontology over-representation analysis.

For experimental data in Mtx select de-regulated genes and on those perform over representation analysis in gene ontology.
Results in GoOver are either as a values list or a csv file, the name of which is assumed to be the ground value of GoOVer.


Opts
  * go(GoSec='BP')
    gene ontology section, in: =|[BP,MF,CC]|=
  * go_frame(GoFrame=bioc_ann_dbi)
    which go frame of GoTerms, Evidence and Gid to use. 
    Default, bioc_ann_dbi, uses r_lib(org.<Org>.eg.db) frames, alternative bio_db uses bio_db tables, 
    before the intro of the option this was the way, so use this for backward compatibility
  * go_over_pv_cut(PvCut=0.05)
    p value filter for the results
  * org(Org=hs)
    one of bio_db_organism/2 first argument values
  * gid(Gid)
    the type of the experimental gene ids for the organism. The default depends on Org, most take symb, apart for pig which takes ensg
    (Caution, this used to be org_exp_id(OrgExpId).)
  * gid_to(GidTo)
    returns the db type for gene ids used in underlying call, mostly ncbi, apart for mouse where it is mgim
  * min_count(MinCount=0)
    min DE symbols in GO term required for a term to be included
  * stem(Stem=false)
    stem for output csv file. When false, use basename of CsvF.
  * symbs(Symbs=false)
    record the symbols for each term (; separated single entry)
  * symbs_hdr(Hymbs="symbols")
    header for the Symbs column
  * symbs_max(SyMax=inf)
    how many symbols to include in the symbols entry (as defined above)
  * to_file(ToF=false)
    when GoOver is unbound, this controls whether the output goes to a file or a values list 
  * universe(Univ=go_exp)
    Univ in : =|[genome,go,go_exp,experiment]|=
    * genome(Gen)
      all gene identifiers in relevant map predicates
    * go(Go)
      all gene identifiers appearing in any ontology
    * go_exp(GoExp)
      gene ontology genes that appear in experiment
    * experiment(Exp)
      all identifiers in the experiment
    * ontology(Onto)
      not implemented yet. like Go above, but only include this Ontology branch

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

?- absolute_file_name( pack('bio_analytics/data/silac/bt.csv'), CsvF ),
   exp_go_over( CsvF, OverF, [to_file(true)] ).

CsvF = '.../swipl/pack/bio_analytics/data/silac/bt.csv',
OverF = go_over_gontBP_p0.05_univExp.csv.

?- absolute_file_name( pack('bio_analytics/data/silac/bt.csv'), CsvF ),
   exp_go_over( CsvF, OverF, [to_file(true),stem(false)] ).

CsvF = '.../swipl/pack/bio_analytics/data/silac/bt.csv',
OverF = '.../swipl/pack/bio_analytics/data/silac/bt_gontBP_p0.05_univExp.csv'.
==

@author nicos angelopoulos
@version  0.1 2019/5/2
@version  0.2 2022/12/20,   =|Univ=go|= and =|Org=gallus|=
@version  0.3 2023/6/5,     option go_frame(GoFrame)
@version  0.4 2024/10/20,   options symbs(Cymbs), symbs_hdr(Hymbs) and symbs_max(SyMax)
@version  0.5 2025/09/24,   option min_count(MinCnt)
@see go_over_universe/6
@tbd 24.10.20, investigate whether the symbs(Cymbs) option can be generalised to save any gene_id
               particularly the NCBI ones that are used by R libs 

*/
exp_go_over( CsvF, GoOver, Args ) :-
    Self = exp_go_over,
    options_append( Self, Args, Opts ),
    exp_go_over_bioc_deps,  % make sure deps are loaded
    debug_call( exp_go_over, options, Self/Opts ),
    bio_diffex( CsvF, DEPrs, NDEPrs, Opts ),
    bio_db_org_in_opts( Org, Opts ),
    kv_decompose( DEPrs, DEGenes, _ ),
    debug_call( exp_go_over, length, de_pairs/DEPrs ),
    bio_list_sort_ne( DEGenes, DEGenesSet ),
    debug_call( exp_go_over, length, de_genes_set/DEGenesSet ),
    bio_list_sort_ne( DEGenes, DEGenesSet ),
    org_gid_map( DEGenesSet, Gids, Opts ),
    debug_call( exp_go_over, length, gids/Gids ),
    options( go_frame(GfOpt), Opts ), 
    go_over_frame( GfOpt, Org, goFrameData, GofOrg ),
    goFrame <- 'GOFrame'(goFrameData, organism= +GofOrg),
    goAllFrame <- 'GOAllFrame'(goFrame),
    gsc <- 'GeneSetCollection'(goAllFrame, setType = 'GOCollection'() ),
    genes <- Gids,
    options( universe(UnivOpt), Opts ),
    go_over_universe( UnivOpt, Org, Gids, NDEPrs, Univ, Opts ),
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
    % this is in r_lib('Category') 
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
    options( symbs(SymbsEnt), Opts ),
    exp_go_over_symbs_entry( SymbsEnt, Org, over, dfOver, Opts ),
    Use = use(Self,GoAspect,UnivOpt,PvCut),
    options( min_count(MinCnt), Opts ),
    exp_go_over_return( GoOver, dfOver, CsvF, Use, MinCnt, Self, Opts ).

exp_go_over_symbs_entry( true, Org, Rov, Rdf, Opts ) :-
    !,
    Signs <- Rdf[*,1],
    Gcats <- geneIdsByCategory(Rov),
    options( gid_to(Gto), Opts ),
    options( symbs_max(Syx), Opts ),
    exp_go_over_symbs( Gcats, Signs, Gto, Syx, GcSymbs, [gid(Gto),gid_to(symb),org(Org)] ),
    Rdf[*,9] <- GcSymbs,
    Names <- names(Rdf),
    once( append(Fames,[_],Names) ),
    options( symbs_hdr(SyHdr), Opts ),
    append( Fames, [SyHdr], NewNames ),
    names(Rdf) <- NewNames.
    % names(Rdf[9]) <- "symbols". -> this or colnames() should work, but it doesnt...
% defaulty:
exp_go_over_symbs_entry( _SymbsEnt, _Org, _Rov, _Rdf, _Opts ).

exp_go_over_symbs( [], _Signs, _Gto, _Syx, [], _Opts ).
exp_go_over_symbs( [NthGc=NthGnsPrv|Prs], Signs, Gto, Syx, SyAtms, Opts ) :-
     % ( NthGc=='GO:0048534' -> trace; true),
     % ( NthGc=='GO:0050789' -> trace; true),
     ( memberchk(NthGc,Signs) ->
          ( is_list(NthGnsPrv) -> NthGnsPrv = NthGns; NthGns = [NthGnsPrv] ),
          bio_gid_std( Gto, NthGns, NthIns ),
          org_gid_map( NthIns, AllSymbs, Opts ),
          ( Syx =:= inf ->
               AllSymbs = Symbs
               ;
               break_nth( Syx, AllSymbs, Symbs, at_short(true) )
          ),
          atomic_list_concat( Symbs, ';', SyAtm ), 
          SyAtms = [SyAtm|TyAtms]
          ;
          SyAtms = TyAtms 
     ),
     exp_go_over_symbs( Prs, Signs, Gto, Syx, TyAtms, Opts ).

exp_go_over_return( GoOver, DfOveR, _CsvF, _Use, MinCnt, Self, Opts ) :-
    var( GoOver ),
    options( to_file(false), Opts ),
    !,
    Rlist <- 'as.list'(DfOveR),
    findall( List, (member(Head=Tail,Rlist),List=[Head|Tail]), Lists),
    mtx_lists( Mtx, Lists ),
    ( MinCnt > 0 ->
          mtx_column_include_rows( Mtx, 'Count', >=(MinCnt), GoOver ),
          length( Mtx, NrBef ),
          length( GoOver, NrAft ),
          debuc( Self, 'Min count of: ~w, reduced matrix rows from: ~d to ~d', [MinCnt,NrBef,NrAft] )
          ;
          debuc( Self, 'Not applying min count constraint, value given: ~w', [MinCnt] ),
          Mtx = GoOver
    ).
exp_go_over_return( GoOver, DfOveR, CsvF, Use, MinCnt, Self, Opts ) :-
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
    ( MinCnt > 0 -> 
          NrBef  <- nrow(DfOveR),
          DfOveR <- DfOveR[DfOveR$'Count'>=MinCnt,*],
          NrAft  <- nrow(DfOveR),
          debuc( Self, 'Min count of: ~w, reduced matrix rows from: ~d to ~d', [MinCnt,NrBef,NrAft] )
          ;
          debuc( Self, 'Not applying min count contraint, value given: ~w', [MinCnt] )
    ),
    <- 'write.csv'(DfOveR, file=+GoOver, 'row.names'='FALSE'),
    % <- print( warnings() ),
    debuc( Self, 'Wrote: ~p', GoOver ).

%% go_over_universe( +Token, +Org, +DEGenes, +NDEPrs, -Universe, +Opts )
%
% Universe is the list of gene identifiers to be used as universe/background for GOstats.
%
% In human (=|Org=hs|=), this is a list of Entrez ids, and in =|Org=mouse|=, a list of Mgim identifiers.<br>
% DEGenes is a list of deferentially expressed gene identifiers, and NDEPrs is a list of non-differential <br>
% expressed Symbol-Pvalue pairs.
% 
% ExpIdTkn identifies the type of ids coming in in NDEPrs, DEGens, are already in standard form for Org.
go_over_universe( experiment, Org, DeGids, NDEPrs, Univ, Opts ) :-
    go_over_universe_exp( Org, DeGids, NDEPrs, Univ, Opts ).
go_over_universe( genome, Org, _DEGids, _NDEPrs, Univ, _Opts ) :-
    go_over_universe_genome( Org, Univ ).
go_over_universe( go, Org, _DEGids, _NDEPrs, Univ, _Opts ) :-
    go_over_universe_go( Org, Univ ).
go_over_universe( go_exp, Org, DEGids, NDEPrs, Univ, Opts ) :-
    go_over_universe_go_exp( Org, DEGids, NDEPrs, Univ, Opts ).

go_over_universe_genome( chicken, Univ ) :-
    findall( Ncbi, cgnc_galg_cgnc_ncbi(_Cgnc,Ncbi), Ncbis ),
    sort( Ncbis, Univ ).
go_over_universe_genome( human, Univ ) :-
    findall( Ncbi, hgnc_homs_symb_ncbi(_Symb,Ncbi), Ncbis ),
    sort( Ncbis, Univ ).
go_over_universe_genome( mouse, Univ ) :-
    findall( Mgim, mgim_musm_mgim_symb(Mgim,_Symb), Mgims ),
    sort( Mgims, Univ ).
go_over_universe_genome( pig, Univ ) :-
     findall( Ncbi, ncbi_suss_ncbi_ensg(Ncbi,_), Ncbis ),
     sort( Ncbis, Univ ).

go_over_universe_go( chicken, Univ ) :-
    findall( Ncbi, ( gont_galg_symb_gont(Symb,_Rl,_Ev,_Go),
                     cgnc_galg_cgnc_symb(Cgnc,Symb),
                     cgnc_galg_cgnc_ncbi(Cgnc,Ncbi)
                   ),
                    Ncbis
           ),
    sort( Ncbis, Univ ).
go_over_universe_go( human, Univ ) :-
    findall( Ncbi,   (   gont_homs_edge_symb(_Go,_En,Symb),
                         hgnc_homs_symb_ncbi(Symb,Ncbi)
                     ),
                         Ncbis
           ),
    sort( Ncbis, Univ ).
go_over_universe_go( mouse, Univ ) :-
    findall( Mgim, gont_musm_mgim_gont(Mgim,_E,_G), Mgims ),
    sort( Mgims, Univ ).
go_over_universe_go( pig, Univ ) :-
     findall( Ncbi, (gont_suss_symb_gont(Symb,_,_,_),ense_suss_ensg_symb(EnsG,Symb),ncbi_suss_ncbi_ensg(Ncbi,EnsG)), Ncbis ),
     sort( Ncbis, Univ ).

% fixme: give doc here, what is this for ?
go_over_universe_exp( chicken, DeGids, NDEPrs, Univ ) :-
    findall( Ncbi, ( member(Symb-NDEPrs),
                     cgnc_galg_cgnc_symb(Cgnc,Symb),
                     cgnc_galg_cgnc_ncbi(Cgnc,Ncbi)
                   ), NDENcbis ),
    append( DeGids, NDENcbis, Ncbis ),
    sort( Ncbis, Univ ).
go_over_universe_exp( human, DeGids, NDEPrs, Univ ) :-
    findall( Ncbi, (member(Symb-_,NDEPrs),hgnc_homs_symb_ncbi(Symb,Ncbi)), NDENcbis ),
    % findall( Entz1, (member(Symb,DEGenes),hgnc_homs_symb_entz(Symb,Entz1)), DEEntzs ),
    % append( DEEntzs, NDEEntzs, Entzs ),
    append( DeGids, NDENcbis, Ncbis ),
    sort( Ncbis, Univ ).
go_over_universe_exp( mouse, DeGids, NDEPrs, Univ ) :-
    findall( Mgim, (member(Symb-_,NDEPrs),mgim_musm_mgim_symb(Mgim,Symb)), NDEMgims ),
    append( DeGids, NDEMgims, Mgims ),
    sort( Mgims, Univ ).

go_over_universe_exp( chicken, DEGids, NDEPrs, Univ, Opts ) :-
     go_over_universe_exp( gallus, DEGids, NDEPrs, Univ, Opts ).
go_over_universe_exp( gallus, DEGids, NDEPrs, Univ, Opts ) :-
    findall( Id, member(Id-_,NDEPrs), Ids ),
    org_gid_map( Ids, NDGids, Opts ), 
    append( DEGids, NDGids, ExpGids ),
    % fixme: we can add tests here, that the ids exist in some table ? 
    bio_list_sort_ne( ExpGids, Univ ).
go_over_universe_exp( human, DEGids, NDEPrs, Univ, Opts ) :-
    findall( Id, member(Id-_,NDEPrs), Ids ),
    org_gid_map( Ids, NDGids, Opts ), 
    append( DEGids, NDGids, ExpGids ),
    % fixme: we can add tests here, that the ids exist in some table ? 
    bio_list_sort_ne( ExpGids, Univ ).
go_over_universe_exp( mouse, DEGids, NDEPrs, Univ, Opts ) :-
    findall( Id, member(Id-_,NDEPrs), Ids ),
    org_gid_map( Ids, NDGids, Opts ), 
    append( DEGids, NDGids, ExpGids ),
    % fixme: we can add tests here, that the ids exist in some table ? 
    bio_list_sort_ne( ExpGids, Univ ).
go_over_universe_exp( pig, DEGids, NDEPrs, Univ, Opts ) :-
    findall( Id, member(Id-_,NDEPrs), Ids ),
    org_gid_map( Ids, NDGids, Opts ), 
    append( DEGids, NDGids, ExpGids ),
    % fixme: we can add tests here, that the ids exist in some table ? 
    bio_list_sort_ne( ExpGids, Univ ).

go_over_universe_go_exp( chicken, DEGids, NDEPrs, Univ, Opts ) :-
    findall( Id, member(Id-_,NDEPrs), Ids ),
    org_gid_map( Ids, NDGids, Opts ),
    append( DEGids, NDGids, ExpGids ),
    findall( ExpGid,( member(ExpGid,ExpGids),   % these are std form, here NCBI Gene ids
                      cgnc_galg_cgnc_ncbi(Cgnc,ExpGid),
                      cgnc_galg_cgnc_symb(Cgnc,Symb),
                      once(gont_galg_symb_gont(Symb,_,_,_))
                    ),
                         List ),
    bio_list_sort_ne( List, Univ ).
go_over_universe_go_exp( human, DEGids, NDEPrs, Univ, Opts ) :-
    findall( Id, member(Id-_,NDEPrs), Ids ),
    org_gid_map( Ids, NDGids, Opts ),
    append( DEGids, NDGids, ExpGids ),
    findall( Ncbi,  ( member(Ncbi,ExpGids),
                      hgnc_homs_symb_ncbi(Symb,Ncbi),
                      once(gont_homs_edge_symb(_,_,Symb))
                    ),
                         List ),
    bio_list_sort_ne( List, Univ ).
go_over_universe_go_exp( mouse, DEGids, NDEPrs, Univ, Opts ) :-
    findall( Id, member(Id-_,NDEPrs), Ids ),
    org_gid_map( Ids, NDGids, Opts ),
    append( DEGids, NDGids, ExpGids ),
    findall( Mgim,  ( member(Mgim,ExpGids),
                      mgim_musm_mgim_symb(Mgim,Symb),
                      once(gont_musm_gont_symb(_,_,Symb))
                    ), 
                         List ),
    bio_list_sort_ne( List, Univ ).
go_over_universe_go_exp( pig, DEGids, NDEPrs, Univ, Opts ) :-
    % fixme: we should be using go_exp probably, because it only comes gont only comes with symbs and unlikely to have a good coverage.
    findall( Id, member(Id-_,NDEPrs), Ids ),
    org_gid_map( Ids, NDGids, Opts ),
    append( DEGids, NDGids, ExpGids ),
    findall( Ncbi, ( member(Ncbi,ExpGids), 
                     ncbi_suss_ncbi_ensg(Ncbi,EnsG),
                     ense_suss_ensg_symb(EnsG,Symb),
                     gont_suss_symb_gont(Symb,_,_,_)
                    ),
                        Ncbis 
           ),
    bio_list_sort_ne( Ncbis, Univ ).

go_over_frame( bioc_ann_dbi, Org, GoFra, DbiOrg ) :-
     bio_conductor_annot_dbi_org( Org, DbiTkn, DbiOrg ),
     bio_conductor_annot_dbi_org_lib( DbiTkn, true, _DbiLib ),
     atom_string( DbiTknAtm, DbiTkn ),
     atomic_list_concat( [org,DbiTknAtm,'egGO'], '.', EgGO ),
     ggframe <- toTable(EgGO),
     GoFra <- 'data.frame'(ggframe$go_id, ggframe$'Evidence', ggframe$gene_id),
     !.
go_over_frame( bio_db, Org, GoFra, DbiOrg ) :-
     go_over_frame_bio_db( Org, GoFra, DbiOrg ),
     !.
go_over_frame( GFopt, OrgOpt, _GoFra, _DbiOrg ) :-
     throw( cannot_coordinate_exp_go_over_options(go_frame(GFopt),org(OrgOpt)) ).

% fixme: these should go to src/bio_org.pl 
%
go_over_frame_bio_db( chicken, GoFra, DbiOrg ) :-
     findall( row(Gid,E,Ncbi), (
                                     gont_galg_symb_gont(Symb,_,E,G),
                                     go_id(Gid,G),
                                     cgnc_galg_cgnc_symb(Cgnc,Symb),
                                     cgnc_galg_cgnc_ncbi(Cgnc,Ncbi)
                               ),
                              Rows
            ),
    go_mtx_df( [row(go_id,evidence,gene_id)|Rows], GoFra, [] ),
    DbiOrg = "Gallus gallus".
go_over_frame_bio_db( human, GoFra, DbiOrg ) :-
    findall( row(Gid,E,Ncbi), 
            ( gont_homs_edge_symb(G,E,S),
              go_id(Gid,G),
              hgnc_homs_symb_ncbi(S,Ncbi)
            ),
                Rows ),
    go_mtx_df( [row(go_id,'Evidence',gene_id)|Rows], GoFra, [] ),
    DbiOrg = "Homo sapiens".
go_over_frame_bio_db( mouse, GoFra, DbiOrg ) :-
    findall( row(Gid,E,M), 
            ( gont_musm_mgim_gont(M,E,G),
              go_id(Gid,G)
            ),
                Rows ),
    go_mtx_df( [row(go_id,evidence,gene_id)|Rows], GoFra, [] ),
    DbiOrg = "Mus musculus".
go_over_frame_bio_db( pig, GoFra, DbiOrg ) :-
    findall( row(Gid,Evi,Ncbi), (  gont_suss_symb_gont(Symb,_,Evi,G), 
                                    go_id(Gid,G),
                                    ense_suss_ensg_symb(EnsG,Symb),
                                    ncbi_suss_ensg_ncbi(EnsG,Ncbi)
                                 ), 
                                       Rows ),
    go_mtx_df( [row(go_id,evidence,gene_id)|Rows], GoFra, [] ),
    DbiOrg = "Sus scrofa".

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
