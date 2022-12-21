
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

@author nicos angelopoulos
@version  0:1 2022/12/21
@see bio_db_organism/2

*/
bio_symbols( Vect, Symbs, Opts ) :-
     ground( Vect ),
     bio_symbols_values( Vect, GIds, Opts ),
     options( org(OrgIn), Opts ),
     options( org_exp_id(ExpId), Opts ),
     bio_db_organism( OrgIn, Org ),
     ( bio_symbols_map( Org, ExpId, GIds, Symbs ) ->
          true
          ;
          throw( bio_symbols(3,possible_exp_id_missing(Org,ExpId)) )
     ).

bio_symbols_map( _Org, symb, Values, Symbs ) :-
     !,
     Values = Symbs.
bio_symbols_map( gallus, ExpId, Values, Symbs ) :-
     bio_symbols_map_gallus( ExpId, Values, Symbs ).
bio_symbols_map( hs, ExpId, Values, Symbs ) :-
     bio_symbols_map_hs( ExpId, Values, Symbs ).
bio_symbols_map( mouse, ExpId, Values, Symbs ) :-
     bio_symbols_map_mouse( ExpId, Values, Symbs ).

bio_symbols_map_gallus( entz, Values, Symbs ) :-
     % fixme: additionals ?
     findall( Symb, ( member(Valu,Values),
                      map_cgnc_gallus_cgnc_entz(Cgnc,Valu),
                      map_cgnc_gallus_cgnc_symb(Cgnc,Symb)
                    ),
                         Symbs ).

bio_symbols_map_hs( entz, Values, Symbs ) :-
     % fixme: additionals ?
     findall( Symb, ( member(Valu,Values),
                      map_hgnc_hgnc_entz(Cgnc,Valu),
                      map_hgnc_hgnc_symb(Cgnc,Symb)
                    ),
                         Symbs ).
bio_symbols_map_mouse( entz, Values, Symbs ) :-
     % fixme: additionals ?
     findall( Symb, ( member(Valu,Values),
                      map_mgim_mouse_mgim_entz(Cgnc,Valu),
                      map_mgim_mouse_mgim_symb(Cgnc,Symb)
                    ),
                         Symbs ).

bio_symbols_values( Vect, Ids, Opts ) :-
     atomic( Vect ),
     !,
     options( mtx(Mtx), Opts ),
     mtx_column( Mtx, Vect, Ids ).
bio_symbols_values( [H|T], Ids, _Opts ) :-
     ( (compound(H),functor(_,_2)) -> 
          findall( Fst, (member(Elem,[H|T]),arg(1,Elem,Fst)), Ids )
          ;
          Ids = [H|T]
     ).
bio_symbols_values( Oth, _Ids, Opts ) :-
     thrown( bio_symbols(cannot_id_vect(Oth,Opts)) ).
