
:- lib(stoics_lib:n_digits_min/3).

go_over_plot_defaults( Defs ) :-
             Defs = [
                         debug(false),
                         pfx_len(50),
                         pos_go_term(8),
                         pos_p_val(2),
                         top_n(30)
                    ].

/** go_over_plot( +GovF, +Opts ).

Plots a lollipop plot for the GO over represented (csv) file or matrix (mtx/2), GovF.

By default when GovF is a file, an .svg is generated at the same name but postfixed by top_TopN (see options for TopN value).

Opts
  * debug(Dbg=false)
    informational, progress messages
  * go(GoSec)
    gene ontology section, in: =|[BP,MF,CC]|=, used for the label. Default is looked for in the filename.
  * pfx_len(PfxLen=50)
    how many initial characters of the GO term name should be included ?
  * pos_go_term(GoTPos=8)
    argument/column position for the GO term name
  * pos_p_val(PvalPos=2)
    argument/column position for the p value
  * stem(Stem)
  * top_n(TopN=30)
    how many of the top terms to include

==
?- go_over_plot( 'res_resist-24.06.05/genes_15kb_resistant_enh_bp_go_over.csv', [] ).
==

@author  nicos angelopoulos
@version 0:1 2024/06/05
@see b_real: gg_lollipop/2
@tbd allow prepending order such as 01- as titles have to be unique
@tbd allow adding the Count/Size on tick label names

*/
go_over_plot( GovF, Args ) :-
     Self = go_over_plot,
     options_append( Self, Args, Opts ),
     % CsvF = 'res_resist-24.06.05/genes_15kb_resistant_enh_bp_go_over.csv',
     mtx( GovF, Mtx ),
     debuc( Self, dims, mtx/Mtx ),
     Mtx = [_|Rows],
     options( pfx_len(PfxLen), Opts ),
     length( PfxCodes, PfxLen ),
     options( top_n(TopN), Opts ),
     options( pos_go_term(GoTPos), Opts ),
     options( pos_p_val(PvalPos), Opts ),
     ( memberchk(go(GoSct),Opts) -> downcase_atom(GoSct,GoLwr); go_over_plot_go_def(GovF,GoLwr) ),
     go_over_plot_got_abv( GoLwr, GoAbv ),
     findall( Derm-Nog,    ( nth1(N,Rows,R),
                             N =< TopN,
                             arg(GoTPos,R,Term),
                             atom_codes(Term, Codes ),
                             (append(PfxCodes,_,Codes) -> PfxCodes = UseCodes; Codes = UseCodes),
                             % number_codes( N, NCs ),
                             % append(NCs,[0'-|UseCodes], Dodes),
                             % atom_codes(Derm, Dodes),
                             atom_codes(Derm, UseCodes),
                             arg(PvalPos,R,Pv),
                             Nog is - log(Pv)
                           ),
                              Pops ),
     reverse( Spop, Pops ),
     ( memberchk(stem(Stem),Opts) -> 
               true
               ;
               ( atomic(GovF) -> os_ext(_,StemPrv,GovF); StemPrv = GoLwr ),
               n_digits_min( 2, TopN, TopTkn ),

               atom_concat( top, TopTkn, TopTknStem ),
               atomic_list_concat( [StemPrv,TopTknStem], '_', Stem )
     ),
     gg_lollipop( Spop, [flip(true),order(true),labels(GoAbv,'-log(p.val)',''),stem(Stem),outputs(svg)] ),
     debuc( Self, length, pops/Pops ).

go_over_plot_go_def( GovF, Sect ) :-
     ( atomic(GovF) -> true; throw(provide_go_opt_if_govf_is_mtx(go_over_plot/2)) ),
     member( Sect, [bp,cc,mf] ),
     atomic_list_concat( ['_',Sect,'_'], '', Patt ),
     sub_atom( GovF, _, _, _, Patt ),
     !.
go_over_plot_go_def( GovF, _Sect) :-
     throw( couldnot_find_one_of_sections_in_filename([bp,cc,mf],GovF,go_over_plot/2) ).

go_over_plot_got_abv(bp, 'GO/Biol.Process').
go_over_plot_got_abv(cc, 'GO/Cell.Component').
go_over_plot_got_abv(mf, 'GO/Mol.Function').
