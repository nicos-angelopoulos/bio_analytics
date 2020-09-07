
exp_go_over_string_graphs_defaults( Defs ) :-
    Wplots = [vjust= -1, node_size(3), format(svg)],
    Defs = [ dir_postfix(go_strings), go_id_clm(1),
             max_overs(false),
             stem_type(go_pair_ord),
             viz_de_opts([]), wgraph_plot_opts(Wplots) ].

/** exp_go_over_string_graphs( +Exp, ?GoOver, ?Dir, -Opts )

Create string graphs for all the over-represented terms in GoOver.

When GoOver is a variable expr_go_over/3 is called to generate it.

Opts
  * dir_postfix(Psfx=go_strings)
    postfix for outputs directory (when Dir is a variable)
  * go_id_clm(GoIdClm=1)
    column id for over represented GO terms
  * max_overs(MaxOvs=false)
    when a number, it is taken as the maximal integer of terms to plot graphs for
  * stem_type(Sty=go_pair_ord)
    similar to go_string_graph/3, but different default, others: =|go_name|=, =|go_id|=.
    Here the length of GO terms (Len) is added to form go_pair_ord(I,Len) for forwarding
  * viz_de_opts(VizOpts = [])
    options for restricting genes to visualise via exp_diff/4.
    Default does not restrict what genes are visualised.
  * wgraph_plot_opts(WgOpts=WgOpts)
    defaults to =|[vjust = -1, node_size(3), format(svg)]|=.


Options are passed to exp_gene_family_string_graph/4.

==
?- debug(real).
?- debug(exp_go_over_string_graphs).
?- absolute_file_name( pack('bio_analytics/data/silac/bt.csv'), Exp ),
   exp_go_over_string_graphs( Exp, Gov, Dir, [] ).

% Sending to R: pltv <- ggnet2(lp_adj,vjust = -1,size = 3,label = pl_v_1,color = pl_v_2,edge.size = pl_v_3,edge.color = "#BEAED4")

==

@author nicos angelopoulos
@version  0.1 2019/5/5
@version  0.2 2020/9/6,   option max_overs()
@see exp_go_over/3
@see exp_gene_family_string_graph/4

*/
exp_go_over_string_graphs( Exp, GoOverIn, Dir, Args ) :-
    options_append( exp_go_over_string_graphs, Args, Opts ),
    os_make_path( Dir ),
    exp_go_over_mtx( GoOverIn, Exp, GoOver, Dir, Opts ),
    options( go_id_clm(Cid), Opts ),
    mtx_column( GoOver, Cid, GOs ),
    options( dir_postfix(Psfx), Opts ),
    exp_go_over_string_graphs_to_dir( Exp, Psfx, Dir ),
    debug( exp_go_over_string_graphs, 'Dir: ~p', [Dir] ),
    options( viz_de_opts(VdfOpts), Opts ),
    options( wgraph_plot_opts(WgOpts), Opts ),
    ( VdfOpts == [] -> Exp = RedExp; exp_diffex(Exp,RedExp,_,[as_pairs(false)|VdfOpts]) ),
    mtx( red_exp.csv, RedExp ),
    options( max_overs(MaxOvsPrv), Opts ),
    (number(MaxOvsPrv) -> MaxOvs is integer(MaxOvsPrv),findall(Go,(between(1,MaxOvs,I),nth1(I,GOs,Go)),MaxGOs); GOs = MaxGOs), 
    options( stem_type(Sty), Opts ),
    ( Sty == go_pair_ord ->
        length( MaxGOs, MaxGOsLen ),
        number_codes( MaxGOsLen, MaxGOsLenCs ),
        length( MaxGOsLenCs, PadLen )
        ;
        PadLen = 0
    ),
    exp_diffex( RedExp, DEPrs, NonDEPrs, [as_pairs(true)|Opts] ),
    go_over_string_graphs_dir( MaxGOs, 1, RedExp, Dir, WgOpts, Sty, PadLen, DEPrs, NonDEPrs, Opts ).
    % maplist( go_over_string_graphs_dir(Exp,Dir,WgOpts,Opts), GOs ).

go_over_string_graphs_dir( [], _I, _Exp, _Dir, _WgOpts, _Sty, _Pad, _DEPrs, _NonDEPrs, _Opts ).
go_over_string_graphs_dir( [Go|Gos], I, Exp, Dir, WgOpts, Sty, Pad, DEPrs, NonDEPrs, Opts ) :-
    % ( atom_concat('GO:',Stem,Go) -> atom_concat(go,Stem,GoTkn); GoTkn=Go ),
    debug( exp_go_over_string_graphs, 'Doing: ~w', [Go] ),
    ( Sty == go_pair_ord -> Tty = go_pair_ord(I,Pad) ; Tty = Sty ),
    go_string_graph_stem( Tty, Go, GoTkn ),
    % ( Go == 'GO:0034447' -> trace; true ),
    directory_file_path( Dir, GoTkn, DirGo ),
    debug( exp_go_over_string_graphs, 'File: ~p', [DirGo] ),
    GoWgOpts = [stem(DirGo)|WgOpts],
    % fixme: you can pass DEPrs and NonDEPrs below, to avoid re-finding them...
    exp_gene_family_string_graph( Exp, Go, DEPrs, NonDEPrs, _, [wgraph_plot_opts(GoWgOpts)|Opts] ),
    J is I + 1,
    go_over_string_graphs_dir( Gos, J, Exp, Dir, WgOpts, Sty, Pad, DEPrs, NonDEPrs, Opts ).

exp_go_over_mtx( GoOverIn, Exp, GoOver, Dir, Opts ) :-
    var( GoOverIn ),
    !,
    os_path( Dir, 'go_over', + Stem ),
    append( Opts, [to_file(true),stem(Stem)], GoOpts ),
    exp_go_over( Exp, GoOver, GoOpts ),
    GoOverIn = GoOver.
exp_go_over_mtx( GoOverIn, _Exp, GoOver, _Dir, _Opts ) :-
    debug( exp_go_over_string_graphs, 'Using stored go_over results from: ~p', [GoOverIn] ),
    mtx( GoOverIn, GoOver ). % fixme: pass Opts ?

exp_go_over_string_graphs_to_dir( _Exp, _Psfx, Dir ) :-
    ground( Dir ),
    !.
exp_go_over_string_graphs_to_dir( Exp, Psfx, Dir ) :-
    file_name_extension( Stem, _Ext, Exp ),  % fixme: Ext = csv
    atom_concat( '_', Psfx, PsfxTkn ),
    atom_concat( Stem, PsfxTkn, Dir ).
