
:- use_module(library(lib)).

% :- set_prolog_flag(lib_suggests_warns,debug).
% use the following to install missing R libs automatically
:- initialization( set_prolog_flag(lib_suggests_warns,install), now ).

:- lib(real).
:- lib(mtx).

:- lib(bio_analytics).

:- debug(bt).

/** bt.

Display the string graph for autophagy from the data in data/silac/bt.csv.

==
?- [pack('bio_analytics/examples/bt')].
?- bt.
==

@author nicos angelopoulos
@version  0:1 2019/4/16
@add ggnet2 save to svg
@see https://briatte.github.io/ggnetwork/ which works with ggrepel

*/
bt :-
    Self = bt,
    debug( Self, 'Starting example ~w ...', [Self] ),
    absolute_file_name( pack('bio_analytics/data/silac/bt.csv'), CsvF ),
    % PvCnm = '', % EvCnm = '',
    Ppts =  [vjust= -1, node_size(3), mode="fruchtermanreingold",
             save(false), format(svg), stem(bt)],
                    % mode="adj"; mode="kamadakawai"
    Opts = [ % exp_pv_cnm(PvCnm), % exp_ev_cnm(EvCnm),
             exp_ev_cut_let(inf), exp_ev_cut_get(-inf),
             exp_pv_cut(0.05),
             include_non_present(false),
             include_non_significant(false),
             minw(200),
             wgraph_plot_opts(Ppts),
             gene_id_cnm('Symbols')
    ],
    exp_gene_family_string_graph( CsvF, autophagy, Graph, Opts ),
    debug( Self, 'Graph: ~w', [Graph] ),
    debug( Self, 'See file: ~w ', ['bt.svg'] ),
    % <- ggsave(file="test.svg", plot=pltv, width=10, height=8),
    debug( Self, '...end of example: ~w', [Self] ).
