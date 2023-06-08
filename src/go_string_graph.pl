
go_string_graph_defaults( Defs ) :-
     Defs = [    
              org(hs),
              plot(false),           % 
              save(false),           % see wgraph_plot/2.
              stem_type(go_name)
            ].

/** go_string_graph( +GoTerm, -Graph, +Opts ).

Plot the STRING graph of all the Symbols in a GO term.

The Graph can be saved via options to wgraph_plot/2.

Opts are passed to symbols_string_graph/3 and wgraph_plot/2.

Opts 

  * org(Org=hs)
    organism for the operation
  * plot(Plot=false)
    whether to plot the graph
  * save(Save=false)
    whether to save the graph (passed to wgraph_plot/2
  * stem_type(Sty=go_name)
    constructs stem for the output filenames: =|go_id|=, =|go_pair|=, =|go_pair_ord(I,Len)|=.

==
?- go_string_graph( 'GO:0016601', G, true ).
?- go_string_graph( 'GO:0016601', G, plot(true) ).

?- go_string_graph( 'GO:0016601', G, [plot(true),save(true)] ).
G = ['ABI2', 'AIF1', 'CDH13', 'HACD3', 'NISCH', 'ALS2'-'RAC1':892, 'BRK1'-'CYFIP1':999, ... - ... : 904, ... : ...|...].
?- ls.
% Rac_protein_signal_transduction.csv                Rac_protein_signal_transduction_graph.csv          Rac_protein_signal_transduction_layout.csv         
true.

?- go_string_graph( 'GO:0016601', G, [plot(true),save(true),stem_type(go_pair)] ).

==

@author nicos angelopoulos
@version  0.1 2016/4/12
@version  0.2 2020/9/3,  expanded go_string_graph_stem/3
@tbd make sure underlying options are compatible.

*/

go_string_graph( Got, Graph, Args ) :-
     options_append( go_string_graph, Args, Opts ),
     % go_term_symbols( GoTerm, Symbs, Opts ),
     options( org(Org), Opts ),
     go_org_symbols( Got, Org, Symbs ),
     symbols_string_graph( Symbs, Graph, Opts ),
     options( plot(Plot), Opts ),
     go_string_graph_plot( Plot, Graph, Got, Opts ).

go_string_graph_plot( true, Graph, GoTerm, Opts ) :-
     options( stem_type(Stype), Opts ),
     go_string_graph_stem( Stype, GoTerm, Stem ),
     wgraph_plot( Graph, [stem(Stem)|Opts] ).
go_string_graph_plot( false, _Graph, _GoTerm, _Opts ).

% go_string_graph_stem( null, _GoTerm, ). % fixme: breaks backward compatibility
go_string_graph_stem( go_pair_ord(I,Len), GoTerm, Stem ) :-
    go_string_graph_stem( go_id, GoTerm, IdStem ),
    go_string_graph_stem( go_name, GoTerm, NmStem ),
    number_codes( I, ICs ),
    length( ICs, ICsLen ),
    From is ICsLen + 1,
    findall( 0, between(From,Len,_), Zeros ),
    append( Zeros, [I], Numbs ),
    atomic_list_concat( Numbs, '', Pfx ),
    atomic_list_concat( [Pfx,IdStem,NmStem], '-', Stem ).

go_string_graph_stem( go_pair, GoTerm, Stem ) :-
    go_string_graph_stem( go_id, GoTerm, IdStem ),
    go_string_graph_stem( go_name, GoTerm, NmStem ),
    atomic_list_concat( [IdStem,NmStem], '-', Stem ).
go_string_graph_stem( go_id, GoTermIn, Stem ) :-
    ( atom_concat('GO:',GoBasename,GoTermIn) ->
        true
        ;
        GoBasename = GoTermIn
    ),
    atomic_list_concat( [go,GoBasename], '', Stem ).
go_string_graph_stem( go_name, GoTermIn, Stem ) :-
    ( atom_concat('GO:',GoBasename,GoTermIn) ->
        % GoTerm = GoTermIn
        true
        ;
        GoBasename=GoTermIn
        % atom_concat('GO:',GoTermIn,GoTerm)
    ),
    ( atom(GoBasename) -> atom_number(GoBasename,GoNum); GoNum = GoBasename ),
     ( gont_homs_gont_gonm( GoNum, GoName ) -> true ; atom_concat(go,GoBasename,GoName) ),
     % atom_codes( GoName, GoNameCs ),
     % codes_replace( GoNameCs, 0' , 0'_, GoUnderCs ),
     % atom_codes( GoUnder, GoUnderCs ),
     % atomic_list_concat( [go,GoUnder], '_', Stem ).
    atomic_list_concat( PartsA, '/', GoName ),
    atomic_list_concat( PartsA, '_OR_', OredGoName ),
    atomic_list_concat( PartsB, ' ', OredGoName ),
    atomic_list_concat( PartsB, '_', Stem ).  % removing go_ prefix

     
codes_replace( [], _C, _W, [] ).
codes_replace( [C|T], C, W, [W|R] ) :-
     !,
     codes_replace( T, C, W, R ).
codes_replace( [H|T], C, W, [H|R] ) :-
     codes_replace( T, C, W, R ).
