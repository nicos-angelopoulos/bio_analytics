
exp_go_over_string_graphs_defaults( Defs ) :-
    Wplots = [vjust= -1, node_size(3), format(svg)],
    Defs = [ dir_postfix(go_strings), go_id_clm(1),
             stem_type(go_pair),
             viz_de_opts([]), wgraph_plot_opts(Wplots) ].

/** exp_go_over_string_graphs( +Exp, ?GoOver, ?Dir, -Opts )

Create string graphs for all the over-represented terms in GoOver.

Opts
  * dir_postfix(Psfx=go_strings)
    postfix for outputs directory (when Dir is a variable)
  * go_id_clm(GoIdClm=1)
    column id for over represented GO terms
  * stem_type(Sty=go_pair)
    as in go_string_graph/3, (possibly different default), others: =|go_name|=, =|go_id|=
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
    options( stem_type(Sty), Opts ),
    go_over_string_graphs_dir( GOs, RedExp, Dir, WgOpts, Sty, Opts ).
    % maplist( go_over_string_graphs_dir(Exp,Dir,WgOpts,Opts), GOs ).

go_over_string_graphs_dir( [], _Exp, _Dir, _WgOpts, _Sty, _Opts ).
go_over_string_graphs_dir( [Go|Gos], Exp, Dir, WgOpts, Sty, Opts ) :-
    % ( atom_concat('GO:',Stem,Go) -> atom_concat(go,Stem,GoTkn); GoTkn=Go ),
    debug( exp_go_over_string_graphs, 'Doing: ~w', [Go] ),
    go_string_graph_stem( Sty, Go, GoTkn ),
    directory_file_path( Dir, GoTkn, DirGo ),
    debug( exp_go_over_string_graphs, 'File: ~p', [DirGo] ),
    GoWgOpts = [stem(DirGo)|WgOpts],
    exp_gene_family_string_graph( Exp, Go, _, [wgraph_plot_opts(GoWgOpts)|Opts] ),
    go_over_string_graphs_dir( Gos, Exp, Dir, WgOpts, Sty, Opts ).

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
