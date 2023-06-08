%
:- lib(ordsets).

symbols_string_graph_defaults( Defs ) :-
     Defs = [
                cohese(max),
                    include_orphans(true),
                org(hs),
                    minw(500),
                    sort_pairs(true),
                    sort_graph(true)
     ].

/** symbols_string_graph( +Symbols, -Graph, +Opts ).

Create the string database Graph between Symbols.

Opts 
  * cohese(Coh=max)
    method for cohesing multiple edges between two nodes =|[false,max,min,umax,umin]|=
    the =|max|= and =|min|= versions are more efficient but do sort the edges,
    whereas =|umax|= and =|umin|= take much longer but leave edges in found order

  * include_orphans(Orph=true)
    set to false to exclude orphans from Graph

  * org(Org=hs)
    which organism do the gene symbols come from

  * minw(500)
    minimum weight (0 =< x =< 999) - not checked

  * sort_pairs(Spairs=true)
    set to false to leave order of edges dependant on order of Symbols

  * sort_graph(Sort=true)
    set to false for not sorting the results

==
?- Got = 'GO:0043552', gene_family( Got, hs, Symbs ), length( Symbs, SymbsLen ),
   symbols_string_graph(Symbs, Graph, true ), length( Graph, GraphLen ).

Got = 'GO:0043552',
Symbs = ['AMBRA1', 'ATG14', 'CCL19', 'CCL21', 'CCR7'|...],
SymbsLen = 31,
Graph = ['TNFAIP8L3', 'AMBRA1'-'ATG14':981, 'AMBRA1'-'PIK3R4':974|...],
GraphLen = 235.

% the following is flawed as it uses Got from mouse and tries to build the graph in human STRING...
?- Got = 'GO:0043552', gene_family( Got, mouse, Symbs ), length( Symbs, SymbsLen ),
   symbols_string_graph(Symbs, Graph, true ), length( Graph, GraphLen ).

Got = 'GO:0043552',
Symbs = Graph, Graph = ['Ambra1', 'Atg14', 'Ccl19', 'Cd19', 'Cdc42', 'Epha8', 'Fgf2', 'Fgfr3', 'Fgr'|...],
SymbsLen = GraphLen, GraphLen = 28.

% this the correct way for running the first query for mouse:
?- Got = 'GO:0043552', gene_family( Got, mouse, Symbs ), length( Symbs, SymbsLen ),
   symbols_string_graph(Symbs, Graph, org(mouse) ), length( Graph, GraphLen ).

Got = 'GO:0043552',
Symbs = ['Ambra1', 'Atg14', 'Ccl19', 'Cd19', 'Cdc42', 'Epha8'|...],
SymbsLen = 28,
Graph = ['Ambra1', 'Cdc42', 'Lyn', 'Nod2', 'Pdgfra', 'Sh3glb1', 'Tnfaip8l3', 'Vav3', ... : ...|...],
GraphLen = 117.
==

@author nicos angelopoulos
@version  0.1 2016/1/18
@version  0.2 2019/4/8,    added organism to incorporate mouse
@version  0.3 2020/9/5,    option cohese(), and debugs
@tbd  implement cohese() values: min, umax and umin

*/
symbols_string_graph( Symbols, Graph, Args ) :-
     Self = symbols_string_graph,
     options_append( symbols_string_graph, Args, Opts ),
     options( org(OrgIn), Opts ),
     bio_db_organism( OrgIn, Org ),
     options( minw(MinW), Opts ),
     options( sort_pairs(Sprs), Opts ),
     findall( SymbA-SymbB:W, ( member(Symb1,Symbols),
                               member(Symb2,Symbols),
                               Symb1 @< Symb2,
                               % Symb1 \== Symb2,
                               symbols_string_graph_pair(Sprs,Symb1,Symb2,SymbA,SymbB),
                               bio_db:org_edge_strg_symb_ord(Org,Symb1,Symb2,W),
                               MinW =< W
                             ),
                                   Pgraph ),
    options( cohese(Coh), Opts ),
    symbols_string_graph_cohese( Coh, Pgraph, Wgraph ),
    options( include_orphans(IncO), Opts ),
    symbols_string_graph_orphans( IncO, Wgraph, Symbols, Ograph ),
    debuc( Self, length, [string_edges_found,edges_cohesed,w_orphans]/[Pgraph,Wgraph,Ograph] ),
    options( sort_graph(Sgra), Opts ),
    symbols_string_graph_sort( Sgra, Ograph, Graph ).

symbols_string_graph_cohese( false, Vgraph, Wgraph ) :-
    Vgraph = Wgraph.
symbols_string_graph_cohese( max, Vgraph, Wgraph ) :-
    ( sort(Vgraph,[H|Graph]) ->
        H = HSA-HSB:Wei,
        symbols_string_graph_cohese_max( Graph, HSA, HSB, Wei, Wgraph )
        ;
        % only possibility is Vgraph is emtpy list...
        Vgraph = Wgraph
    ).
symbols_string_graph_cohese( min, Vgraph, Wgraph ) :-
     throw( unimplemented_cohesion_method(min) ),
     sort( Vgraph, Ograph ),
     symbols_string_graph_cohese_min( Ograph, Wgraph ).
symbols_string_graph_cohese( umax, Vgraph, Wgraph ) :-
     throw( unimplemented_cohesion_method(umax) ),
     symbols_string_graph_cohese_umax( Vgraph, Wgraph ).
symbols_string_graph_cohese( umin, Vgraph, Wgraph ) :-
     throw( unimplemented_cohesion_method(umin) ),
     symbols_string_graph_cohese_umin( Vgraph, Wgraph ).

symbols_string_graph_cohese_max( [], Csa, Csb, Cwe, [Csa-Csb:Cwe] ).
symbols_string_graph_cohese_max( [Hsa-Hsb:Hwe|T], Csa, Csb, Cwe, Wgraph ) :-
    ( Hsa-Hsb = Csa-Csb ->
        Nsa = Csa,
        Nsb = Csb,
        Nwe is max( Cwe, Hwe ),
        Wgraph = Tgraph,
        % debugging only
        Rwe is min( Cwe, Hwe ),
        debuc( symbols_string_graph(details), info, 'Removing duplicate (pair only) edge: ~w'/[Hsa-Hsb:Rwe] )
        ;
        Nsa = Hsa,
        Nsb = Hsb,
        Nwe = Hwe,
        Wgraph = [Csa-Csb:Cwe|Tgraph]
    ),
    symbols_string_graph_cohese_max( T, Nsa, Nsb, Nwe, Tgraph ).

symbols_string_graph_pair( true, Symb1, Symb2, SymbA, SymbB ) :-
     sort( [Symb1,Symb2], [SymbA,SymbB] ).
symbols_string_graph_pair( false, Symb1, Symb2, Symb1, Symb2 ).

symbols_string_graph_sort( true, Ograph, Graph ) :-
     sort( Ograph, Graph ).
symbols_string_graph_sort( false, Graph, Graph ).

symbols_string_graph_orphans( true, Wgraph, Symbols, Ograph ) :-
     string_add_vertices_1( Symbols, Wgraph, Ograph ).
symbols_string_graph_orphans( false, Wgraph, _Symbols, Wgraph ).

string_add_vertices_1( [], New, New ).
string_add_vertices_1( [V|Vs], G, New ) :-
     atomic( V ),
     string_has_vertex( G, V ),
     !,
     string_add_vertices_1( Vs, G, New ).
string_add_vertices_1( [V|Vs], G, New ) :-
     atomic( V ),
     ord_add_element( G, V, G1 ),
     string_add_vertices_1( Vs, G1, New ).

string_has_vertex( [E|_Es], V ) :-
     string_edge_has_vertex( E, V ),
     !.
string_has_vertex( [_E|Es], V ) :-
     string_has_vertex( Es, V ).

string_edge_has_vertex( V, V ).
string_edge_has_vertex( V-_:_, V ).
string_edge_has_vertex( _-V:_, V ).

