
:- lib(stoics_lib:prefix_atom/3).
:- lib(stoics_lib:kv_compose/3).
:- lib(stoics_lib:kv_decompose/3).

:- lib(col_faint/3).

exp_gene_family_string_graph_defaults( Defs ) :-
    Defs = [ 
                exp_ev_log(true),
                extra_symbol_prefix('p.'),
                extra_symbols([]),
                faint_factor(1.2),
                include_non_significant(true),
                include_non_present(false),
                % node_colours(clrs(red,green,yellow,khaki4,grey)),
                node_colours(clrs(tomato,steelblue,yellow,khaki4,grey)),
                org(hs),
                plot(true),
                wgraph_plot_opts([]),
                wgraph_plot_defs([])
    ].

/** exp_gene_family_string_graph( +Exp, +Family, -Graph, +Opts ).
    exp_gene_family_string_graph( +Exp, +Family, ?DEPrs, ?NonDEPrs, -Graph, +Opts ).

Generate and possibly plot the STRING graph of a known gene family as affected in a biological Exp_eriment.
Each family gene is placed of the following states according to wether it was identified in Exp and 
how it was modified (in relation to background conditions): 

  * Red = upregulated
  * SteelBlue = downregulated
  * Orange = either both directions or unknown direction
  * Khaki4 = for extra nodes
  * Grey = not identified

Note that de-reguation trumps identification, that is, currently there is no distiction
between genes that seen to be both significantly de-regulated and identified versus
just those that are simply significantly de-regulated.

In the exp_gene_family_string_graph/6 version DEPrs and NonDEPrs can be provided, which optimises use in
loops over families (see exp_go_over_string_graphs/4).

Opts
  * exp_ev_log(EvLog=true)
    are the expression values (ev) log values ? (also passed to bio_diffex/4)
  * extra_symbol_prefix(Xrfs='p.')
    for when a symbol is already in graph
  * extra_symbols(Xymbs=[])
    extra symbols for the graph
  * faint_factor(Fctr=1.2)
    faint factor for non-significant nodes
  * include_non_present(IncP=true)
    false excludes non-indentified family genes 
  * include_non_significant(IncS=true)
    false excluded non-significant family genes 
  * node_colours(Clrs=clrs(red,green,yellow,khaki4,grey))
    colours for the different types of nodes
  * org(Org=hs)
    organism in which Family is looked in (and experiment was performed)
  * org_exp_id(ExpID)
    the type of the experimental gene ids for the organism. default depends on Org, but currently all map to symb
  * plot(Plot=true)
    whether to plot the graph via wgraph_plot/2
  * wgraph_plot_opts(WgOpts=[])
    options for wgraph_plot/2 (in preference to any defaults from Self)
  * wgraph_plot_defs(WgDefs=[])
    options for wgraph_plot/2 (added to the end, after any defaults from Self)

These Options are passed to a number of other pack predicates: bio_diffex/4, bio_symbols/3, symbols_string_graph/3 and, selectively, wgraph_plot/2.

See [pack('bio_analytics/examples/bt.pl')].

==
?- lib(real), lib(mtx),
   absolute_file_name( pack('bio_analytics/data/silac/bt.csv'), CsvF ),
   Ppts = [ vjust= -1, node_size(3), mode="fruchtermanreingold", format(svg), stem(bt)],
   Opts = [ exp_ev_cut_let(inf), exp_ev_cut_get(-inf), 
            include_non_present(false), include_non_significant(false), 
            minw(200), wgraph_plot_opts(Ppts) ],
   exp_gene_family_string_graph( CsvF, autophagy, G, Opts ).
==
Produces file: bt.svg

[[../doc/images/bt.svg]]

@author nicos angelopoulos
@version  0.1 2019/4/15
@version  0.2 2024/10/20,    option extra_symbol_prefix(Xfx)
@see bio_diffex/4, bio_symbols/3, symbols_string_graph/3
@see gene_family/3
@see wgraph_plot/2

*/
exp_gene_family_string_graph( Exp, Fam, StGraph, Args ) :-
    exp_gene_family_string_graph( Exp, Fam, _DEPrs, _NonDEPrs, StGraph, Args ).

exp_gene_family_string_graph( Exp, Fam, DEPrs, NonDEPrs, StGraph, Args ) :-
    options_append( exp_gene_family_string_graph, Args, Opts ),
    ( var(DEPrs) ->
        bio_diffex( Exp, DEPrs, NonDEPrs, Opts )
        ;
        true
    ),
    options( org(Org), Opts ),
    gene_family( Fam, Fymbs, org(Org) ),
    sort( Fymbs, Oymbs ),
    options( exp_ev_log(EvLog), Opts ),
    ( EvLog == true -> Zero is 1; Zero is 0 ),
    exp_gene_family_string_graph_node_colours( Nlrs, Flrs, Opts ),
    options( include_non_present(IncP), Opts ),
    options( include_non_significant(IncS), Opts ),
    bio_symbols( DEPrs, DESymbPrs, Opts ),
    bio_symbols( NonDEPrs, NonDESymbPrs, Opts ),
    options( extra_symbols(Xymbs), Opts ),
    sort( Xymbs, Xombs ),
    exp_gene_family_string_graph_nodes( Oymbs, DESymbPrs, NonDESymbPrs, Zero, Nlrs, Flrs, IncP, IncS, Nodes, NdClrs ),
    sort( Nodes, OrdNodes ),
    ord_intersection( Xombs, OrdNodes, ComSymbs ),
    ord_subtract( Xombs, OrdNodes, XonlySymbs ),
    ord_union( Xombs, OrdNodes, SymbNodes ),
    options( extra_symbol_prefix(XybPfx), Opts ),
    maplist( prefix_atom(XybPfx), ArtNodes, ComSymbs ),
    append( XonlySymbs, ArtNodes, AddNodes ),
    arg( 4, Nlrs, XClr ),
    findall( XClr, member(_,AddNodes), AddNdClrs ),
    append( Nodes, AddNodes, GraphNodes ),
    append( NdClrs, AddNdClrs, GraphNdClrs ),
    kv_compose( GraphNodes, GraphNdClrs, GraphPrs ),
    sort( GraphPrs, OrdGrPrs ),
    % fixme: test this is correct and update wgraph documentantion
    kv_decompose( OrdGrPrs, _, OrdGraphNdClrs ),
    symbols_string_graph( SymbNodes, StGraphNat, Opts ),
    exp_gene_family_string_graph_add_edges( ComSymbs, ArtNodes, 750, StGraphNat, StGraphPrv ),
    options( wgraph_plot_opts(DomWGOpts), Opts ),
    options( wgraph_plot_defs(DefWGOpts), Opts ),
    EgfWGOpts = [colours(OrdGraphNdClrs),plotter(ggnet2)],
    flatten( [DomWGOpts,EgfWGOpts,DefWGOpts], WGOpts ),
    string_graph_single_edge_weight_adjusted( StGraphPrv, StGraph ),
    debug( exp_gene_family_string_graph(graph), 'Family: ~w, string graph: ~w', [Fam,StGraph] ),
    exp_wgraph_plot( StGraph, Fam, WGOpts ).

exp_wgraph_plot( [], Fam, _WGOpts ) :-
    !,
    debug( ba(info), 'Empty string graph removed, for family: ~w', [Fam] ).
exp_wgraph_plot( Graph, _Fam, Opts ) :-
     wgraph_plot( Graph, Opts ).

string_graph_single_edge_weight_adjusted( Graph, Adjusted ) :-
    findall( W, member(_X-_Y:W,Graph), AllWs ),
    sort( AllWs, Ws ),
    % when there is a single weight it should be in range 0-10 
    % as ggplot2 (and others i think) have difficulty drawing reasonably widthed edges
    ( Ws = [_Single] -> 
        maplist( string_graph_single_edge_weight_adjust_entry, Graph, Adjusted )
        ;
        Adjusted = Graph
    ).

string_graph_single_edge_weight_adjust_entry( X-Y:W, Edge ) :-
        Wadj is max(1, W // 500),
        Edge = X-Y:Wadj,
        !.
string_graph_single_edge_weight_adjust_entry( Edge, Edge ).
        
exp_gene_family_string_graph_add_edges( [], [], _W, GraphIn, GraphNw ) :-
     GraphIn = GraphNw.
exp_gene_family_string_graph_add_edges( [S|Ss], [A|As], W, GraphIn, GraphNw ) :-
     ( S @< A ->
          wgraph_add_edges(GraphIn, [S-A:W], GraphNx )
          ;
          wgraph_add_edges(GraphIn, [A-S:W], GraphNx )
     ),
     % 0114
     exp_gene_family_string_graph_add_edges( Ss, As, W, GraphNx, GraphNw ).

exp_gene_family_string_graph_nodes( [], _DEPrs, _NonDEPrs, _Zero, _Nlrs, _Flrs, _IncGrey, _IncS, [], [] ).
exp_gene_family_string_graph_nodes( [Gn|Gns], DEPrs, NonDEPrs, Zero, Nlrs, Flrs, IncP, IncS, Nodes, NdClrs ) :-
    debug(exp_gene_family_string_graph(node), 'Gene node: ~w', [Gn] ),
    findall( Ev, member(Gn-Ev,DEPrs), Evs ),
    exp_gene_family_string_graph_nodes_de( Evs, Gn, Zero, NonDEPrs, Nlrs, Flrs, IncP, IncS, Nodes, NdClrs, TNodes, TNdClrs ),
    exp_gene_family_string_graph_nodes( Gns, DEPrs, NonDEPrs, Zero, Nlrs, Flrs, IncP, IncS, TNodes, TNdClrs ).

exp_gene_family_string_graph_nodes_de( [], Gn, Zero, NonDEPrs, Nlrs, Flrs, IncP, IncS, Nodes, NdClrs, TNodes, TNdClrs ) :-
    findall( Ev, member(Gn-Ev,NonDEPrs), Evs ),
    ( Evs == [] ->
        ( IncP == false ->
            TNodes = Nodes,
            TNdClrs = NdClrs
            ;
            arg(5,Nlrs,Grey),
            Nodes = [Gn|TNodes],
            NdClrs = [Grey|TNdClrs]
        )
        ;
        ( IncS == false ->
            TNodes = Nodes,
            TNdClrs = NdClrs
            ; 
            Nodes = [Gn|TNodes],
            exp_gene_family_string_graph_nodes_de_direction( Evs, Zero, Dir ),
            arg( Dir, Flrs, NdClr ),
            NdClrs = [NdClr|TNdClrs]
        )
    ).
exp_gene_family_string_graph_nodes_de( [Ev|Evs], Gn, Zero, _NonDEPrs, Nlrs, _Flrs, _IncP, _IncS, Nodes, NdClrs, TNodes, TNdClrs ) :-
    exp_gene_family_string_graph_nodes_de_direction( [Ev|Evs], Zero, Dir ),
    arg( Dir, Nlrs, NdClr ),
    Nodes = [Gn|TNodes],
    NdClrs = [NdClr|TNdClrs].

exp_gene_family_string_graph_nodes_de_direction( Evs, Zero, Dir ) :-
    partition( <(Zero), Evs, More, LessOrEqual ),
    ( More == [] ->
        ( LessOrEqual == [] ->
            throw( directionless(Evs) )
            ;
            Dir is 2
        )
        ;
        ( LessOrEqual == [] ->
            Dir is 1
            ;
            Dir is 3
        )
    ).

exp_gene_family_string_graph_node_colours( NdClrs, FtClrs, Opts ) :-
    options( node_colours(NdClrs), Opts ),
    options( faint_factor(Fctr), Opts ),
    NdClrs =.. [NcName|ClrsL],
    col_faint( ClrsL, Fctr, FaintsL ),
    FtClrs =.. [NcName|FaintsL].
