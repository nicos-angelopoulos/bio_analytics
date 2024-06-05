
bio_symbols_defaults( Defs ) :-
     Defs = [ 
               debug(false),
               org(hs),
               gid(symb)
     ].

/** bio_symbols( +Vect, -Symbs, +Opts ).

Convert gene identifiers pointed by Vect to Symbols.

Vect can be a list or pairlist (K-V) of which the K element is taken to be the identifier.
Vect can also be an atomic representing  the Nth column on column for Mtx. 

Opts
  * debug(Dbg=false)
    progress, informational messages
  * mtx(Mtx)
    to be used if Vect is atomic
  * org(Org=hs)
    organism the data come from, via bio_db_organism/2
  * gid(Gid=symb)
    the type of the experimental gene ids for the organism. default depends on Org, but currently all map to symb
    (caution, this used to be org_exp_id(), changed in v0.2

==
?- hgnc_homs_hgnc_ncbi( 19295, Ncbi ), bio_symbols( [Ncbi], Symbs, gid(ncbi) ).
Ncbi = 114783,
Symbs = ['LMTK3'].
==

@author nicos angelopoulos
@version  0:1 2022/12/21
@version  0:2 2024/06/05,  changed option org_exp_id() to gid()
@see bio_db_organism/2

*/
bio_symbols( Vect, Symbs, Args ) :-
     Self = bio_symbols,
     options_append( Self, Args, Opts ),
     ground( Vect ),
     bio_symbols_values( Vect, Vals, Opts ),
     options( org(OrgIn), Opts ),
     options( gid(Gid), Opts ),
     bio_db_organism( OrgIn, Org ),
     ( bio_symbols_map( Org, Gid, Vals, Symbs ) ->
          true
          ;
          throw( bio_symbols(3,exp_id_missing_for_org(Gid,Org)) )
     ).

bio_symbols_values( Ids, Ids, _Opts ) :-
     is_list( Ids ),
     !.
bio_symbols_values( Vect, Ids, Opts ) :-
     atomic( Vect ),
     !,
     options( mtx(Mtx), Opts ),
     mtx_column( Mtx, Vect, Ids ).
bio_symbols_values( Oth, _Ids, Opts ) :-
     thrown( bio_symbols(cannot_id_vect(Oth,Opts)) ).
