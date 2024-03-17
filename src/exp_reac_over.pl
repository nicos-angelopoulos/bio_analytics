
:- lib(stoics_lib:kv_decompose/3).

exp_reac_over_defaults( Args, Defs ) :-
                    Defs = [ debug(true),
                             org(hs),
                             gid(ncbi),
                             gid_to(Gto),
                             universe(experiment)
                           ],
     % exp_reac_gid_default( Org, Gid ),
     options_return( gid_to(Gto), Args, [pack(bio_analytics),pred(ex_reac_over/3),option(gid_to(Gto))] ),
     Gto = ncbi.

/** exp_reac_over(+Etx, +Opts).

Perform reactome pathway over-representation analysis.

Etx is the matrix of experimental values (should pass as 1st arg in mtx/2).
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
  * universe(Univ=experiment)
    the universe, or background for genes in the statistical test (also: =|reac|=)

Examples
==
?- exp_reac_over([]).
==

@author nicos angelopoulos
@version  0.1 2024/03/16
@see bio_diffex/4
@see mtx/2
@see bio_db_organism/2
@see org_gid_map/3

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
     length( IdsDE, DENof ),
     % length( IdsND, NDNof ),
     % find all background genes in any reactome pathway = Pop
     bio_db_org_in_opts( _Org, [org_tkn(Okn)|Opts] ),
     options( universe(Univ), Opts ),
     at_con( [reac,Okn,ncbi,reap], '_', Func ),
     debuc( Self, length, de_ids/IdsDE ),
     exp_reac_over_universe_ids( Univ, Self, Func, IdsDE, IdsND, IdsUniV ),
     debuc( Self, length, univ_ids/IdsUniV ),
     length( IdsUniV, UniVNof ),
     % find all reactome pathways that have at least 1 DE product
     Goal =.. [Func,ADEId,_,APway],
     findall( APway, (member(ADEId,IdsDE),Goal), PwaysL ),
     sort( PwaysL, Pways ),
     debuc( Self, length, pways/Pways ),
     exp_reac_over_ncbi_reactome( IdsDE, Func, ReacIdsDE ),
     length( ReacIdsDE, ReacIdsDENof ),
     debuc( Self, length, reac_ids_de/ReacIdsDE ),
     maplist( exp_reac_hygeom(Self,IdsDE,IdsUniV,Func,Okn,UniVNof,ReacIdsDENof), Pways, ReOverPrs ),
      /* 
      for each pathway 
          find DE genes in pathway (PathDENof)
          find background genes in pathway (PathBkNof)
          <- phyper(PathDENof,PathBkNof,Pop - PathBkNof,DENof,lower.tail='FALSE')
     */
     keysort( ReOverPrs, OrdROPrs ),
     kv_decompose( OrdROPrs, OrdPvs, ReOverRows ),
     OrdQvs <- 'p.adjust'( OrdPvs, method =+'BH' ),
     Hdr = row(reactome,'p.value',expected,count,size,pathway),
     mtx_column_add( [Hdr|ReOverRows], 3, ['adj.pvalue'|OrdQvs], AdjMtx ),
     % "GOMFID","Pvalue","adj.pvalue","OddsRatio","ExpCount","Count","Size","Term"
     exp_reac_over_return( AdjMtx, ReOver, Etx ),
     debuc( Self, end, true ).

exp_reac_over_return( Rtx, ReOver, _Etx ) :-
     ground( ReOver ),
     !,
     mtx( ReOver, Rtx ).

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

exp_reac_hygeom( Self, IdsDE, IdsUniv, Func, Okn, UniVNof, ReacDENof, Pway, Row ) :-
     GoalDE =.. [Func,InPwayDE,_,Pway],
     findall( InPwayDE, (member(InPwayDE,IdsDE),GoalDE), InPwayDEs ),
     GoalUniv =.. [Func,InPwayUniv,_,Pway],
     findall( InPwayUniv, (member(InPwayUniv,IdsUniv),GoalUniv), InPwayUniVs ),
     maplist( length, [InPwayDEs,InPwayUniVs], [InPwayDEsNof,InPwayUniVsNof] ),
     NotInPwayUniVsNof is UniVNof - InPwayUniVsNof,
     % fixme: use fisher.test() to also test the OddsRatio ? and double check
     % also check dhyper
     Pv <- phyper(InPwayDEsNof,InPwayUniVsNof,NotInPwayUniVsNof,ReacDENof,'lower.tail'='FALSE'),
     debuc( Self, 'Got: ~w', [Pv <- phyper(InPwayDEsNof,InPwayUniVsNof,NotInPwayUniVsNof,ReacDENof,'lower.tail'='FALSE')] ),
     at_con( [reac,Okn,reap,repn], '_', RecnFnc ),
     RecnG =.. [RecnFnc,Pway,Pwnm],
     call( RecnG ),
     Exp is (InPwayUniVsNof * ReacDENof) / UniVNof,
     Row = Pv-row(Pway,Pv,ReacDENof,Exp,InPwayDEsNof,InPwayUniVsNof,Pwnm).

exp_reac_over_universe_ids( experiment, _Self, Func, IdsDE, IdsND, IdsUniv ) :-
     % findall( Ncbi, ((member(Id,IdsDE);member(Id,IdsND)),reac_homs_ncbi_reap(Ncbi,_,_Reap)), NcbisL ),
     % all the experimental ids that participate in at least 1 pathway
     append( IdsDE, IdsND, Ids ),
     exp_reac_over_ncbi_reactome( Ids, Func, IdsUniv ).
exp_reac_over_universe_ids( reac, _Self, Func, _IdsDE, _IdsND, IdsUniv ) :-
     Goal =.. [Func,Ncbi,_,_],
     % findall( Ncbi, reac_homs_ncbi_reap(Ncbi,_,Reap), NcbisL ),
     findall( Ncbi, Goal, NcbisL ),
     sort( NcbisL, IdsUniv ).

exp_reac_over_ncbi_reactome( Ncbis, Func, NcbisSubset ) :-
     Goal =.. [Fun,Ncbi,_,_],
     findall( Ncbi, (member(Ncbi,Ncbis),Goal), NcbisL ),
     sort( NcbisL, IdsUniv ).
