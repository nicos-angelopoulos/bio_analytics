
:- lib(promise(pl_vector/3,b_real:pl_vector/3)).

bio_p_adjust_defaults( Defs ) :-
               Defs = [
                         debug(true),
                         hdr('q.value'),
                         rel(1)
                      ].

/** bio_p_adjust( +Obj, +Adj, +Opts ).

Adj is the p.adjusted version of Obj. 

Obj can be a list, an R vector, or Cid:Mtx term (see pl_vector/3).
In the later case RelPos (see options) is the relative position for the new column (see mtx_relative_position/4),
if none is given Adj is returned as a list.

The adjustment happens via =|p.adjust()|= function in R.

Opts
  * debug(Dbg=true)
    informational, progress messages
  * hdr(Hdr='q.value')
    the column header if we adding Adj 
  * rel(RelPos=1)
    when Obj is a matrix and this is given, it is taken to be the relative position for a new
    column of the new values to be added to the matrix (see mtx_relative_position/4).
    Def is to append the column one after the pvalue column.
    New matrix returned in Adj.

Opts are passed to b_real:pl_vector/3 (predicate only loaded at run time, so this is a 
promised dependency to pack(b_real)).

Examples
==
?-   Mtx = [row(exp,pv),row(a,0.1),row(b,0.02),row(c,0.001)],   
     bio_p_adjust( pv:Mtx, Adj, cps(Cps) ).
==

@author nicos angelopoulos
@version  0.1 2023/09/02
@see b_real:pl_vector/3

*/

bio_p_adjust( Obj, Adj, Args ) :-
     Self = bio_p_adjust,
     options_append( Self, Args, Opts ),
     pl_vector( Obj, List, [mtx(Mtx),cps(Cps)|Opts] ),
     debuc( Self, length, vector_list/List ),
     AdjV <- p.adjust(List),
     ( var(Mtx) ->
          AdjV = Adj
          ;
          debuc( Self, dims, in_mtx/Mtx ),
          Mtx = [Hdr|_],
          options( rel(Rel), Opts ),
          mtx_relative_pos( Rel, Cps, Hdr, Pos ),
          options( hdr(AdjCnm), Opts ),
          mtx_column_add( Mtx, Pos, [AdjCnm|AdjV], Adj )
     ),
     debuc( Self, end, true ).
