
:- lib(stoics_lib:kv_decompose/3).

exp_reac_gid_default

exp_reac_over_defaults( Args, Defs ) :-
                    Defs = [ debug(true),
                             org(hs),
                             gid(Gid),
                             gid_to(Gto)
                           ],
          ( memberchk(org(Org),Args)->true; Org=hs ),
          exp_reac_gid_default( Org, Gid ),
          options_return( gid_to(Gto), Args, [pack(bio_analytics),pred(ex_reac_over/3),option(gid_to(Gto))] ).

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
@see org_id_map/3

*/

exp_reac_over( Etx, Args ) :-
     Self = exp_reac_over,
     options_append( Self, Args, Opts ),
     bio_diffex( Etx, DEPrs, NDPrs, Opts ),
     debuc( Self, length, [dx,non_dx]/[DEPrs,NDPrs] ),
     kv_decompose( DEPrs, DEs, _ ),
     kv_decompose( NDPrs, NDs, _ ),
     org_id_map( DEs, IdsDE, Opts ),
     org_id_map( NDs, IdsND, Opts ),
     length( IdsDE, DENof ),
     % length( IdsND, NDNof ),
     % find all background genes in any reactome pathway = Pop
     options( universe(Univ), Opts ),
     known( exp_reac_over_universe_ids(Univ,Self,IdsDE,IdsND,IdsUniv) ),

     find all reactome pathways
      for each pathway 
          find DE genes in pathway (PathDENof)
          find background genes in pathway (PathBkNof)
          <- phyper(PathDENof,PathBkNof,Pop - PathBkNof,DENof,lower.tail='FALSE')

     debuc( Self, end, true ).

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
        
*/

exp_reac_over_universe_ids( experiment, Self, IdsDE, IdsND, IdsUniv ) :-
     findall( Ncbi, ((member(Id,IdsDE);member(Id,IdsND)),reac_homs_ncbi_reap(Ncbi,_,_Reap)), NcbisL ),
     sort( Ncbi, NcbisL ),
exp_reac_over_universe_ids( reac, IdsDE, IdsND, IdsUniv ) :-
     findall( Ncbi, reac_homs_ncbi_reap(Ncbi,_,Reap), NcbisL ),
     sort( Ncbi, NcbisL ),
