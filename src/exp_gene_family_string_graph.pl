
:- lib(col_faint/3).

exp_gene_family_string_graph_defaults( Defs ) :-
    Defs = [ 
                exp_ev_log(true),
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
  * Green = downregulated
  * Orange = either both directions or unknown direction
  * Grey = not identified

Note that de-reguation trumps identification, that is, currently there is no distiction
between genes that seen to be both significantly de-regulated and identified versus
just those that are simply significantly de-regulated.

In the exp_gene_family_string_graph/6 version DEPrs and NonDEPrs can be provided, which optimises use in
loops over families (see exp_go_over_string_graphs/4).

Opts
  * exp_ev_log(EvLog=true)
    are the expression values (ev) log values ? (also passed to bio_diffex/4)
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

These Options are passed to a number of other pack predicates.

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
@see bio_diffex/4
@see gene_family/3
@see wgraph_plot/2
@see bio_symbols/3

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
    %
    sort( Fymbs, Oymbs ),
    options( exp_ev_log(EvLog), Opts ),
    ( EvLog == true -> Zero is 1; Zero is 0 ),
    exp_gene_family_string_graph_node_colours( Nlrs, Flrs, Opts ),
    options( include_non_present(IncP), Opts ),
    options( include_non_significant(IncS), Opts ),
    bio_symbols( DEPrs, DESymbPrs, Opts ),
    bio_symbols( NonDEPrs, NonDESymbPrs, Opts ),
    exp_gene_family_string_graph_nodes( Oymbs, DESymbPrs, NonDESymbPrs, Zero, Nlrs, Flrs, IncP, IncS, Nodes, NdClrs ),
    % 
    % construct_graph( Oymbs, DEPrs, NonDEPrs, Graph, GraphWGOpts ),
    symbols_string_graph( Nodes, StGraphPrv, Opts ),
    options( wgraph_plot_opts(DomWGOpts), Opts ),
    options( wgraph_plot_defs(DefWGOpts), Opts ),
    EgfWGOpts = [colours(NdClrs),plotter(ggnet2)],
    flatten( [DomWGOpts,EgfWGOpts,DefWGOpts], WGOpts ),
    string_graph_single_edge_weight_adjusted( StGraphPrv, StGraph ),
    debug( exp_gene_family_string_graph(graph), 'Family: ~w, string graph: ~w', [Fam,StGraph] ),
    exp_wgraph_plot( StGraph, Fam, WGOpts ).

exp_wgraph_plot( [], Fam, _WGOpts ) :-
    !,
    debug( ba(info), 'Empty string graph removed, for family: ~w', [Fam] ).
% exp_wgraph_plot( [Node], Fam, _WGOpts ) :-
    % atomic(Node),
    % !,
    % debug( exp_gene_family_string_graph, 'Single node graph removed, for family: ~w', [Fam] ).
exp_wgraph_plot( Graph, _Fam, Opts ) :-
    wgraph_plot( Graph, Opts ).

string_graph_single_edge_weight_adjusted( Graph, Adjusted ) :-
    findall( W, member(_X-_Y:W,Graph), AllWs ),
    sort( AllWs, Ws ),
    % when therre is a single weight it should be in range 0-10 
    % as ggplot2 (and others i think) have difficult drawing reasonably widthed edges
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

/*
construct_graph( [], _DEPrs, _NonDEPrs, [], GraphWGOpts ) :-
construct_graph( [Gn|Gns], DEPrs, NonDEPrs, Graph, GraphWGOpts ) :-
    findall( Ev, member(Gn-Ev,DEPrs), Evs ),
    */

exp_gene_family_string_graph_node_colours( NdClrs, FtClrs, Opts ) :-
    options( node_colours(NdClrs), Opts ),
    options( faint_factor(Fctr), Opts ),
    NdClrs =.. [NcName|ClrsL],
    col_faint( ClrsL, Fctr, FaintsL ),
    FtClrs =.. [NcName|FaintsL].
