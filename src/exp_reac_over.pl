
:- lib(stoics_lib:kv_decompose/3).

exp_reac_over_defaults( Args, Defs ) :-
                    Defs = [ debug(true),
                             org(hs),
                             gid(ncbi),
                             gid_to(Gto),
                             ids_de(_IdsDE),
                             ids_de_reac(_ReacIdsDE),
                             min_count(0),
                             mtx_cutoff(_,_,false),
                             pways(univ),
                             rec_clm(symbols),
                             rec_ids(false),
                             rec_sep(';'),
                             % mtx_cutoff(0.05,'adj.pvalue',<),
                             universe(experiment)
                           ],
     % exp_reac_gid_default( Org, Gid ),
     options_return( gid_to(Gto), Args, [pack(bio_analytics),pred(ex_reac_over/3),option(gid_to(Gto))] ),
     Gto = ncbi.

/** exp_reac_over(+Etx, ?ReOver, +Opts).

Perform reactome pathway over-representation analysis.

Etx is the matrix of experimental values (should pass as 1st arg in mtx/2).<br>
If ReOver is ground at call, it is taken to be a filename on which the matrix results are written to using mtx/2.<br>
If ReOver is not ground at call, it is unified to the resulting matrix.<br>
Note, that the representation of ReOver is different to that of the 2nd argument of exp_go_over/3.

Opts
  * debug(Dbg=false)
    informational, progress messages
  * org(Org=hs)
    organism id for experimental data; bio_db_organism/2 first argument values
  * gid(Gid)
    the gene id db token identifier for the genes in Etx. Default depends on Org.
  * gid_to(Gto)
    returns the gene id db token identifier for interrogating the reactome db
    (currently returns ncbi)
  * ids_de(IdsDE)
    returns the de ids
  * ids_de_reac(IdsDE)
    returns the de ids that appear in Reactome pathways
  * min_count(MinCount=0)
    min DE symbols in a Reactome pathway required for it to be included
  * mtx_cutoff(Cnm=_,Val=_,Dir=false)
    filter the output matrix- meant for FDR cutoffs (see mtx_column_threshold/3)
  * pways(Pways=univ)
    pathways to consider, the value only affects the corrected p.values as the longer the 
    list of pathways the stronger the correction. In order of tightness: 
    * dx(Dx) 
      pathways that contain at least one Gid picked by bio_diffex/4
    * exp(Exp) 
      pathways that contain at least one Gid in Etx
    * univ(PwUniv)
      pathways that contain at least one Gid in Univ (below option universe())
    * reac(Reac)
      all reactome pathways for Org
  * rec_clm(Rlm=symbols)
    name for extra ids column, iff _Rids\==false_ 
  * rec_ids(Rids=false)
    whether to record memeber gids of each pathway (give a -(Gid,Rec) pairlist if you want to map)
    (atom recorded is of the form: Rec1<Rep>Rec2<Rep>Rec3)
  * rec_sep(Rep=';')
    separator for Rids concatenation
  * universe(Univ=experiment)
    the genes universe, or background for genes in the statistical test (also: =|reac(tome)|=)

Examples
==
?- exp_reac_over([]).
?- exp_reac_over([mtx_cutoff(0.05,'adj.pvalue',<)]).
==

@author nicos angelopoulos
@version  0.1 2024/03/16
@version  0.2 2025/09/28,  option min_count(MinCount)
@see bio_diffex/4
@see mtx/2
@see bio_db_organism/2
@see org_gid_map/3
@tbd add oddsRatios to make the output equivalent to output from exp_go_over/3.

*/

exp_reac_over( Etx, ReOver, Args ) :-
     Self = exp_reac_over,
     options_append( Self, Args, Opts ),
     bio_diffex( Etx, DEPrs, NDPrs, Opts ),
     debuc( Self, length, [dx,non_dx]/[DEPrs,NDPrs] ),
     kv_decompose( DEPrs, DEs, _ ),
     kv_decompose( NDPrs, NDs, _ ),
     org_gid_map( DEs, IdsDE, Opts ),
     org_gid_map( NDs, IdsND, Opts ),
     % find all background genes in any reactome pathway = Pop
     bio_db_org_in_opts( _Org, [org_tkn(Okn)|Opts] ),
     options( universe(Univ), Opts ),
     at_con( [reac,Okn,ncbi,reap], '_', Func ),
     debuc( Self, length, de_ids/IdsDE ),
     debuc( Self, length, non_de_ids/IdsND ),
     exp_reac_over_universe_ids( Univ, Self, Func, IdsDE, IdsND, IdsUniV ),
     debuc( Self, length, univ_ids/IdsUniV ),
     length( IdsUniV, UniVNof ),
     % find all reactome pathways that have at least 1 DE product
     options( pways(WhcPways), Opts ),
     exp_reac_over_pways( WhcPways, Func, IdsDE, IdsND, IdsUniV, Pways ),
     debuc( Self, length, pways/Pways ),
     exp_reac_over_ncbi_reactome( IdsDE, Func, ReacIdsDE ),
     length( ReacIdsDE, ReacIdsDENof ),
     debuc( Self, length, reac_ids_de/ReacIdsDE ),
     options( [rec_clm(Rlm),rec_ids(Rids),rec_sep(Rep)], Opts ),
     options( ids_de(IdsDE), Opts ),
     options( ids_de_reac(ReacIdsDE), Opts ),
     % fixme: should the following be using ReacIdsDE instead of IdsDE
     maplist( exp_reac_hygeom(Self,IdsDE,IdsUniV,Func,Okn,UniVNof,ReacIdsDENof,Rids,Rep), Pways, ReOverPrs ),
      /* 
      for each pathway 
          find DE genes in pathway (PathDENof)
          find background genes in pathway (PathBkNof)
          <- phyper(PathDENof,PathBkNof,Pop - PathBkNof,DENof,lower.tail='FALSE')
     */
     keysort( ReOverPrs, OrdROPrs ),
     kv_decompose( OrdROPrs, OrdPvs, ReOverRows ),
     OrdQvs <- 'p.adjust'( OrdPvs, method =+'BH' ),
     ( Rids == false ->
          Hdr = row(reactome,'p.value',expected,count,size,pathway)
          ;
          Hdr = row(reactome,'p.value',expected,count,size,pathway,Rlm)
     ),
     mtx_column_add( [Hdr|ReOverRows], 3, ['adj.pvalue'|OrdQvs], AdjMtx ),
     % "GOMFID","Pvalue","adj.pvalue","OddsRatio","ExpCount","Count","Size","Term"
     CutOpt = mtx_cutoff(_CutClm,_CutVal,CutDir),
     options( CutOpt, Opts ),
     ( CutDir == false -> 
          debuc( Self, 'Not curtailing matrix on p values. Given: ~w', [CutOpt] )
          ;
          mtx_column_threshold( AdjMtx, AdjThreshMtx, Opts ),
          length( AdjMtx, NrBefCut ),
          length( AdjThreshMtx, NrAftCut ),
          debuc( Self, 'P values cutoff term: ~w, reduced matrix length from: ~d to ~d', [CutOpt,NrBefCut,NrAftCut] )
     ),
     options( min_count(MinCnt), Opts ),
     ( MinCnt > 0 -> 
          mtx_column_include_rows( AdjThreshMtx, count, >=(MinCnt), ThreshMtx ),
          length( AdjThreshMtx, NrBefMc ),
          length( ThreshMtx, NrAftMc ),
          debuc( Self, 'Min count of: ~w, reduced matrix length from: ~d to ~d', [MinCnt,NrBefMc,NrAftMc] )
          ;
          debuc( Self, 'Not applying min count constraint, value given: ~w', [MinCnt] ),
          AdjThreshMtx = ThreshMtx
     ),
     exp_reac_over_return( ThreshMtx, Etx, ReOver ),
     debuc( Self, end, true ).

exp_reac_over_pways( dx, Func, IdsDE, _IdsND, _IdsUniV, Pways ) :-
     Goal =.. [Func,ADEId,_,APway],
     findall( APway, (member(ADEId,IdsDE),Goal), PwaysL ),
     sort( PwaysL, Pways ).
exp_reac_over_pways( exp, Func, IdsDE, IdsND, _IdsUniV, Pways ) :-
     append( IdsDE, IdsND, AllIds ),
     sort( AllIds, Ids ),
     Goal =.. [Func,ADEId,_,APway],
     findall( APway, (member(ADEId,Ids),Goal), PwaysL ),
     sort( PwaysL, Pways ).
exp_reac_over_pways( univ, Func, _IdsDE, _IdsND, IdsUniV, Pways ) :-
     Goal =.. [Func,ADEId,_,APway],
     findall( APway, (member(ADEId,IdsUniV),Goal), PwaysL ),
     sort( PwaysL, Pways ).
exp_reac_over_pways( reac, Func, _IdsDE, _IdsND, _IdsUniV, Pways ) :-
     Goal =.. [Func,_,_,APway],
     findall( APway, Goal, PwaysL ),
     sort( PwaysL, Pways ).

exp_reac_over_return( Rtx, _Etx, ReOver ) :-
     ground( ReOver ),
     !,
     mtx( ReOver, Rtx ).
exp_reac_over_return( Rtx, _Etx, Rtx ).

/*   ?-
        reac_galg_reap_repn(Reap,Repn),
        findall(Gene,reac_galg_ncbi_reap(Gene,_,Reap),Genes), 
        sort(Genes,Unes), length(Unes,Len), 
        write( Repn:Len), nl, 
        fail.

        https://stackoverflow.com/questions/8382806/hypergeometric-test-phyper
        Pop = number of genes = 5260
        PopSucc = number of background genes in the pathway
        PopFail = number of background genes not in pathway
        Sam = number of DE genes
        SamSucc= sampled successes = how many DE genes are in pathway
        phyper(SamSucc,PopSucc,PopFail,Sam,lower.tail=FALSE)
        phyper(62,1998,5260-1998,131-62,lower.tail=FALSE)
       ?- Ph <- phyper(62,1998,5260-1998,131-62,lower.tail='FALSE').
       Ph = 1.3718295741097734e-20.
        
"GOMFID","Pvalue","adj.pvalue","OddsRatio","ExpCount","Count","Size","Term"

*/

exp_reac_hygeom( _Self, IdsDE, IdsUniv, Func, Okn, UniVNof, ReacDENof, Rid, Rep, Pway, Row ) :-
     GoalDE =.. [Func,InPwayDE,_,Pway],
     findall( InPwayDE, (member(InPwayDE,IdsDE),GoalDE), InPwayDEs ),
     GoalUniv =.. [Func,InPwayUniv,_,Pway],
     findall( InPwayUniv, (member(InPwayUniv,IdsUniv),GoalUniv), InPwayUniVs ),
     maplist( length, [InPwayDEs,InPwayUniVs], [InPwayDEsNof,InPwayUniVsNof] ),
     NotInPwayUniVsNof is UniVNof - InPwayUniVsNof,
     % fixme: use fisher.test() to also test the OddsRatio ? and double check
     % also check dhyper
     Pv <- phyper(InPwayDEsNof,InPwayUniVsNof,NotInPwayUniVsNof,ReacDENof,'lower.tail'='FALSE'),
     at_con( [reac,Okn,reap,repn], '_', RecnFnc ),
     RecnG =.. [RecnFnc,Pway,Pwnm],
     call( RecnG ),
     Exp is (InPwayUniVsNof * ReacDENof) / UniVNof,
     % Row = Pv-row(Pway,Pv,ReacDENof,Exp,InPwayDEsNof,InPwayUniVsNof,Pwnm).
     ( Rid == false -> 
          Row = Pv-row(Pway,Pv,Exp,InPwayDEsNof,InPwayUniVsNof,Pwnm)
          ;
          ( Rid == true -> 
               atomic_list_concat(InPwayDEs,Rep,Ral)
               ;
               ( Rid = [_-_|_] ->
                    findall( Val, (member(InP,InPwayDEs),memberchk(InP-Val,Rid)), Vals ),
                    atomic_list_concat( Vals, Rep, Ral )
                    ;
                    throw( cannot_decipher_rid(Rid) )
               )
          ),
          Row = Pv-row(Pway,Pv,Exp,InPwayDEsNof,InPwayUniVsNof,Pwnm,Ral)
     ).

exp_reac_over_universe_ids( experiment, _Self, Func, IdsDE, IdsND, IdsUniv ) :-
     % findall( Ncbi, ((member(Id,IdsDE);member(Id,IdsND)),reac_homs_ncbi_reap(Ncbi,_,_Reap)), NcbisL ),
     % all the experimental ids that participate in at least 1 pathway
     append( IdsDE, IdsND, Ids ),
     exp_reac_over_ncbi_reactome( Ids, Func, IdsUniv ).
exp_reac_over_universe_ids( reac, Self, Func, IdsDE, IdsND, IdsUniv ) :-
     exp_reac_over_universe_ids( reactome, Self, Func, IdsDE, IdsND, IdsUniv ).
exp_reac_over_universe_ids( reactome, _Self, Func, _IdsDE, _IdsND, IdsUniv ) :-
     Goal =.. [Func,Ncbi,_,_],
     % findall( Ncbi, reac_homs_ncbi_reap(Ncbi,_,Reap), NcbisL ),
     findall( Ncbi, Goal, NcbisL ),
     sort( NcbisL, IdsUniv ).

exp_reac_over_ncbi_reactome( Ncbis, Func, IdsUniv ) :-
     Goal =.. [Func,Ncbi,_,_],
     findall( Ncbi, (member(Ncbi,Ncbis),Goal), NcbisL ),
     sort( NcbisL, IdsUniv ).
