
:- lib(stoics_lib:en_list/2).
:- lib(stoics_lib:map_succ_list/4).
:- lib(stoics_lib:date_two_digit_dotted/1).

exp_go_over_string_graphs_multi_defaults( Defs ) :-
    Defs = [
        data(data),
        data_ext(csv),
        dated(true),
        keep_mtx_sel(true),
        keep_mtx_viz(true),
        mtx_opts([]),
        viz_diffex_only_constraint(true),
        viz_pv_constraint(true)
    ].

/** exp_go_over_string_graphs_multi( +Opts ).

Run multiple exp_go_over_string_graphs/4 runs.

Can ran across a number of values for parameters and number of datasets.
A directory is created for each combination of input parameters.
Subdirectories are created for each input dataset.

Can be used with single parameter values as well, in order to add that
parameter to the output directory's name.


Opts
  * data(Data=data)
    which datasets to use, can be a directory
  * data_ext(DtExt=csv)
    extension selector for when Data is a directory
  * dated(Dated=true)
    adds date (ate_two_digit_dotted/1) stamp on output directories
  * dfx_opts(DfxOpts)
    invariant options for exp_go_over_string_graphs/4
  * dfx_viz_opts(VizOpts)
    invariant options for viz_de_opts() of exp_go_over_string_graphs/4
  * keep_mtx_sel(KpSel=true)
    record, to file, the mtx of values significant values selected
  * keep_mtx_viz(KpViz=true)
    record, to file, the mtx of values added for vizualisation (non significants)
  * mtx_opts(MtxOpts=[])
    options to pass when reading the input matrices
  * multi(Keys,Vals)
    parameters and associated possible values to pass to exp_go_over_string_graphs/4 (Vals, can be a single value).
    If multiple Keys are given, then they synchronise- always taking the 
  * viz_diffex_only_constraint(DoCon=true)
    enforce that when de_max() and exp_pv_cut() are the same, use viz_de_opts(diffex_only)
  * viz_pv_constraint(PvCon=true)
    enforce that PvCut for non-significants must be equal or larger than that for significants

multi() options are processed in the order they are given, with later ones varied first in the runs.

==
?- Opts = [],
   exp_go_over_string_graphs_multi( Opts ).
==

@author nicos angelopoulos
@version  0:1 2020/09/17
@tbd shorten to exp_multi/1 ?

*/
exp_go_over_string_graphs_multi( Args ) :-
    Self = exp_go_over_string_graphs_multi,
    options_append( Self, Args, Opts ),
    options( data(DataIn), Opts ),
    os_exists( DataIn, type(base(Type)) ),
    exp_multi_data_files( Type, DataIn, DataFs, Opts ),
    debuc( Self, enum, data_files/DataFs ),
    options( dfx_opts(DfxOpts), Opts ),
    options( dfx_viz_opts(VizOpts), Opts ),
    findall( Pairs, ( member(multi(KeyPrv,ValsPrv),Opts),
                      findall( Pair,    (
                                            % en_list(KeysPrv,Keys), 
                                            en_list(ValsPrv,Vals),
                                            % member(Key,Keys), 
                                            member(Val,Vals),
                                            ( is_list(KeyPrv) -> findall(Key-Val,member(Key,KeyPrv),Pair) ; Pair = KeyPrv-Val )
                                        ),
                                            Pairs )
                    ),
                        Nest ),
    debuc( Self, 'Nest: ~w', [Nest] ),
    options( viz_pv_constraint(PvCon), Opts ),
    % findall( Comb, (exp_multi_combine(Nest,Comb),viz_con(PvCon,Comb)), Combs ),
    findall( AComb, combination(Nest,AComb/[]), AllCombs ),
    include( viz_con(PvCon), AllCombs, Combs ),

    debuc( Self, 'Combs: ~w', [Combs] ),
    options( viz_diffex_only_constraint(DoCon), Opts ),
    exp_multi_date( Date, Opts ),
    options( mtx_opts(MtxOpts), Opts ),
    exp_multi( Combs, DataIn, DataFs, DfxOpts, DoCon, VizOpts, MtxOpts, Date ).

   %                  exp_set_constraints( PvCon, Set )
                  % ),
                    % Sets ),

exp_multi_date( Date, Opts ) :-
    options( dated(Dated), Opts ),
    ( Dated == true ->
                    date_two_digit_dotted( Date )
                    ; % defaulty
                    Date = false
    ).  

exp_multi( [], _DataIn, _DataFs, _DfxOpts, _DoCon, _VizOpts, _MtxOpts, _Date ).
exp_multi( [Set|Sets], DataIn, DataFs, DfxOpts, DoCon, VizOpts, MtxOpts, Date ) :-
    findall( Tkn, (member(K-V,Set),atomic_list_concat([K,V],'=',Tkn)), Tkns ),
    atomic_list_concat( Tkns, '_', Stem ),
    exp_multi_dir_name( Date, Stem, Milled ),
    map_succ_list( bio_analytics:viz_option_pair, Set, VizPrs, DfxPrs ),
    debuc( exp_go_over_string_graphs_multi, 'dname: ~w, set:~w', [Milled,Set] ),
    debuc( exp_go_over_string_graphs_multi, 'VizPrs: ~w, DfxPrs:~w', [VizPrs,DfxPrs] ),
    Milts = [   call_options(false),
                outputs_to('debug_messages.txt'),
                type(dir) ],
    length( DataFs, DataFsLen ),
    os_mill( DataIn, bio_analytics:exp_multi_run(DataFs,1,DataFsLen,VizPrs,DfxPrs,DfxOpts,DoCon,VizOpts,MtxOpts), Milled, Milts ),
    exp_multi( Sets, DataIn, DataFs, DfxOpts, DoCon, VizOpts, MtxOpts, Date ).

exp_multi_run( [], _I,_Ln,_VizPrs, _DfxPrs, _DfxOptsIn, _DoCon, _VizOptsIn, _MtxOpts, _Data, _Milled ).
exp_multi_run( [DataF|DataFs], I, Ln, VizPrs, DfxPrs, DfxOptsIn, DoCon, VizOptsIn, MtxOpts, Data, Milled ) :-
    Self = exp_go_over_string_graphs_multi,
    ( Ln > 1 -> debug_consec(Self,magenta,'Doing (~d/~d): ~p',[I,Ln,DataF]) ; true ),
    debuc( Self, task(Milled), true ),
    exp_multi_run_viz_options( DoCon, DfxPrs, VizPrs, DfxOptsIn, VizOpts ),
    findall( Opt, (member(Key-Val,DfxPrs),Opt =.. [Key,Val]), DfxOptsMulti ),
    flatten( [debug(true),viz_de_opts(VizOpts),DfxOptsMulti,DfxOptsIn], DfxOptsPrv ),
    os_path( Data, DataF, RelDataF ),
    os_ext( _Ext, DataStem, DataF ),
    mtx( RelDataF, Mtx, MtxOpts ),
    os_path( Milled, DataStem, Target ),
    ( select(de_mtx(DeMtxIn),DfxOptsPrv,DfxOptsShort) -> os_path(Target,DeMtxIn,DeMtx), DfxOpts=[de_mtx(DeMtx)|DfxOptsShort] ; DfxOptsPrv=DfxOpts ),
    debuc( Self, 'Run options: ~w', [DfxOpts] ),
    Catcher = debuc(ba(info),'Exp_go_over_string_graphs/5 call finished with caught error: ~w', Ball),
    catch( exp_go_over_string_graphs(Mtx,_GoOver,Target,DfxOpts), Ball, Catcher ),
    J is I + 1,
    exp_multi_run( DataFs, J, Ln, VizPrs, DfxPrs, DfxOptsIn, DoCon, VizOptsIn, MtxOpts, Data, Milled ).

exp_multi_run_viz_options( true, DfxPrs, VizPrs, _DfxOptsIn, VizOpts ) :-
% defaulty
    memberchk( exp_pv_cut-PvCut, DfxPrs ),
    memberchk( exp_pv_cut-PvCut, VizPrs ),
    !,
    diffex_only = VizOpts.
exp_multi_run_viz_options( _, _DfxPrs, VizPrs, DfxOptsIn, VizOpts ) :-
    findall( Opt, (member(Key-Val,VizPrs),Opt =.. [Key,Val]), RunPfxOpts ),
    append( RunPfxOpts, DfxOptsIn, VizOpts ).

viz_option_pair( K-V, NkK-V ) :-
    atom_concat( 'viz_', NkK, K ).

exp_multi_dir_name( false, Stem, Dname ) :-
    !,
    atomic_list_concat( [res,Stem], '-', Dname ).
exp_multi_dir_name( Date, Stem, Dname ) :-
    atomic_list_concat( [res,Stem,Date], '-', Dname ).

viz_con( true, Comb ) :-
    !,
    memberchk( exp_pv_cut-PvCut, Comb ),
    memberchk( viz_exp_pv_cut-VzCut, Comb ),
    PvCut =< VzCut.
% defaulty
viz_con( _, _ ).

exp_multi_combine( [], [] ).
exp_multi_combine( [H|T], [Elem|Rest] ) :-
    member( Elem, H ),
    exp_multi_combine( T, Rest ).

exp_multi_data_files( file, Data, DataFs, _Opts ) :-
    DataFs = [Data].
exp_multi_data_files( dir, Data, DataFs, Opts ) :-
    options( data_ext(Ext), Opts ),
    os_sel( os_files, ext(Ext), DataFs, dir(Data) ).

/** combination( +Nest, -List ).


==
?- combination( [[1,2,3],[a,b,c],[[x,y],[s,t]]], List ).
==
*/

combination_naive( [], [] ).
combination_naive( [H|T], Comb ) :-
    member( Elem, H ),
    combination_naive( T, Tomb ),
    ( is_list(Elem) -> append(Elem,Tomb,Comb); Comb=[Elem|Tomb] ).

combination( [], T/T ).
combination( [H|T], Comb/Tomb ) :-
    member( Elem, H ),
    ( is_list(Elem) -> append(Elem,Bomb,Comb); Comb=[Elem|Bomb] ),
    combination( T, Bomb/Tomb ).
