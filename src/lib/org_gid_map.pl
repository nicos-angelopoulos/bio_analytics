
:- use_module(library(apply)).  % maplist/3.
:- use_module(library(lib)).

:- lib(options).
:- lib(debug_call).

:- lib(stoics_lib:map_succ_list/3).

org_gid_map_defaults( Defs ) :-
     Defs = [  debug(true),
               gid(hgnc),
               org(hs),
               % gid_to(_),
               rmv_e(true),
               sort(true)
            ].

/** org_gid_map(+Ids,-ToIds,+Opts).

Map gene Ids from one type to another type guided by selections in Opts.

Opts
  * debug(Dbg=false)
    informational, progress messages
  * org(Org=hs)
    organism for the Ids and ToIds
  * gid(IdDb=hgnc)
    gene id source database
  * gid_to(ToDb)
    gene id target database
  * rmv_e(RmvE=true)
    whether to remove '' entries
  * sort(Sort=true)
    whether to sort the results

Examples
==
?-  
    hgnc_homs_symb_hgnc('LMTK3', Hgnc), 
    org_gid_map(Hgnc,Ncbi)
    hgnc_homs_symb_ncbi(Symb,Ncbi).

Hgnc = 19295,
Ncbi = 114783,
Symb = 'LMTK3'.

==

@author nicos angelopoulos
@version  0.1 2024/03/16
@version  0.2 2024/10/18,  recipe for chicken: ncbi -> symb
@tbd needs some resilience options, what happens if an id cannot be mapped ? (call pred, or fail, or now=skip)

*/

org_gid_map( Ids, ToIds, Args ) :-
     Self = org_gid_map,
     options_append( Self, Args, Opts ),
     bio_db_org_in_opts( Org, Opts ),
     options( gid(Gid), Opts ),
     options( gid_to(Gto), Opts ),
     org_gid_map_to( Org, Gid, Gto, Ids, ToIdsPrv ),
     options( sort(Sort), Opts ),
     options( rmv_e(RmvE), Opts ),
     once( org_gid_map_clean( Sort, RmvE, ToIdsPrv, ToIds ) ).

org_gid_map_clean( true, true, ToIdsPrv, ToIds ) :-
     sort( ToIdsPrv, Ord ),
     ( select('',Ord,ToIds) -> true; ToIds = Ord ).
org_gid_map_clean( true, false, ToIdsPrv, ToIds ) :-
     sort( ToIdsPrv, ToIds ).
org_gid_map_clean( false, true, ToIdsPrv, ToIds ) :-
     findall( Id, (member(Id,ToIdsPrv),Id \== ''), ToIds ).
org_gid_map_clean( false, false, ToIdsPrv, ToIds ) :-
     ToIdsPrv = ToIds.

org_gid_map_to( _Org, Gid, Gto, Ids, ToIds ) :-
     Gto == Gid,
     !,
     ToIds = Ids.
org_gid_map_to( Org, Gid, Gto, Ids, ToIds ) :-
     is_list( Ids ),
     org_gid_map_1( Org, Gid, Gto, Ids, ToIds ),
     !.
org_gid_map_to( Org, Gid, Gto, Ids, ToIds ) :-
     \+ is_list( Ids ),
     org_gid_map_1( Org, Gid, Gto, [Ids], [ToIds] ),
     !.
org_gid_map_to( Org, Gid, Gto, _Ids, _ToIds ) :-
     throw( cannot_map_gids(Org,Gid,Gto), [pack(bio_analytics),pred(org_gid_map/3)] ).

org_gid_map_1( chicken, Gid, Gto, Ids, ToIds ) :-
     org_gid_map_galg( Gid, Gto, Ids, ToIds ).
org_gid_map_1( human, Gid, Gto, Ids, ToIds ) :-
     org_gid_map_homs( Gid, Gto, Ids, ToIds ).
org_gid_map_1( mouse, Gid, Gto, Ids, ToIds ) :-
     org_gid_map_musm( Gid, Gto, Ids, ToIds ).
org_gid_map_1( pig, Gid, Gto, Ids, ToIds ) :-
     org_gid_map_suss( Gid, Gto, Ids, ToIds ).

% fixme: add ensg -> ncbi, from ifitm project
org_gid_map_galg( cgnc, ncbi, Ids, ToIds ) :-
    findall( Ncbi,  (member(Cgnc,Ids),cgnc_galg_cgnc_ncbi(Cgnc,Ncbi)), ToIds ).
org_gid_map_galg( ncbi, symb, Ids, ToIds ) :-
     map_succ_list( bio_analytics:org_gid_map_galg_ncbi_symb, Ids, ToIds ).

org_gid_map_galg_ncbi_symb( Ncbi, Symb ) :-
     ncbi_galg_ncbi_symb( Ncbi, Symb ),
     !.
org_gid_map_galg_ncbi_symb( Ncbi, Symb ) :-
     cgnc_galg_cgnc_ncbi( Cgnc, Ncbi ),
     cgnc_galg_cgnc_symb( Cgnc, Symb ),
     !.
% fixme: i doubt this gives anything new, we can test in bio_db
org_gid_map_galg_ncbi_symb( Ncbi, Symb ) :-
     unip_galg_unip_ncbi( Unip, Ncbi ),
     unip_galg_unip_symb( Unip, Symb ),
     !.

% org_go_over_std_gene_ids_chicken( SrcT, Set, Gids ) :-
%    atom_concat( cgnc_galg_cgnc_, SrcT, SrcNm ),
%    findall( Ncbi,  (member(SrcG,Set),call(SrcNm,Cgnc,SrcG),cgnc_galg_cgnc_ncbi(Cgnc,Ncbi)), Ncbis ),
%    sort( Ncbis, Gids ).

org_gid_map_homs( symb, ncbi, Ids, ToIds ) :-
    findall( Ncbi,  (member(Symb,Ids),hgnc_homs_symb_ncbi(Symb,Ncbi)), ToIds ).

org_gid_map_musm( ncbi, mgim, Ids, ToIds ) :-
    findall( Mgim,  (member(Ncbi,Ids),mgim_musm_mgim_ncbi(Mgim,Ncbi)), ToIds ).
org_gid_map_musm( symb, mgim, Ids, ToIds ) :-
    findall( Mgim,  (member(Symb,Ids),mgim_musm_mgim_symb(Mgim,Symb)), ToIds ).

org_gid_map_suss( ensg, ncbi, Ids, ToIds ) :-
     findall( Ncbi, (member(EnsG,Ids), org_gid_map_suss_ensg_ncbi(EnsG,Ncbi)), ToIds ).
org_gid_map_suss( symb, ncbi, Ids, ToIds ) :-
     findall( Ncbi, ( member(Symb,Ids), 
                      ense_suss_ensg_symb(EnsG,Symb),
                      org_gid_map_suss_ensg_ncbi(EnsG,Ncbi)
                    ), 
                                             ToIds ).
org_gid_map_suss_ensg_ncbi( EnsG, Ncbi ) :-
     ncbi_suss_ncbi_ensg( Ncbi, EnsG ),
     !.
