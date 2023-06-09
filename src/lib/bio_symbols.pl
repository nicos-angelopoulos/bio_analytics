
bio_symbols_defaults( Defs ) :-
     Defs = [ 
               debug(false),
               org(hs),
               org_exp_id(symb)
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
  * org_exp_id(ExpId=symb)
    the type of the experimental gene ids for the organism. default depends on Org, but currently all map to symb

==
?- hgnc_homs_hgnc_ncbi( 19295, Ncbi ), bio_symbols( [Ncbi], Symbs, org_exp_id(ncbi) ).
Ncbi = 114783,
Symbs = ['LMTK3'].
==

@author nicos angelopoulos
@version  0:1 2022/12/21
@see bio_db_organism/2

*/
bio_symbols( Vect, Symbs, Args ) :-
     Self = bio_symbols,
     options_append( Self, Args, Opts ),
     ground( Vect ),
     bio_symbols_values( Vect, Vals, Opts ),
     options( org(OrgIn), Opts ),
     options( org_exp_id(ExpId), Opts ),
     bio_db_organism( OrgIn, Org ),
     ( bio_symbols_map( Org, ExpId, Vals, Symbs ) ->
          true
          ;
          throw( bio_symbols(3,exp_id_missing_for_org(ExpId,Org)) )
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
