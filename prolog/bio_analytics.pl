:- module( bio_analytics, [
                bio_analytics_version/2,          % -Vers, -Date
                bio_diffex/4,                     % +Exp, -DEs -NDEs, +Opts
                bio_p_adjust/3,                   % +Obj, +AdjObj, +Opts
                bio_symbols/3,                    % +Vect, -Symbs, +Opts
                bio_volcano_plot/2,               % +Mtx, Opts
                exp_gene_family_string_graph/4,   % +Exp, +Fam, -Graph, +Opts
                exp_go_over/3,                    % +Exp, -GoOver, +Opts
                exp_go_over_string_graphs/4,      % +Exp, ?GoOver, ?Dir, +Opts
                exp_go_over_string_graphs_multi/1,% +Opts
                exp_reac_over/3,                  % +Exp, -ReOver, +Opts
                gene_family/3,                    % +Alias, -Symbols, +Opts
                go_org_symbol/3,                  % +Org, +GoT, -Symbol
                go_org_symbols/3,                 % +Org, +GoT, -Symbols
                go_over_plot/2,                   % +GovF, +Opts
                go_string_graph/3,                % +GoT, -Graph, +Opts
                go_symbols_reach/3,               % +GoT, -Symbs, +Opts
                symbols_string_graph/3            % +Symbs, -Graph, +Opts
        ] ).

/** <module> Computational biology data analytics.

Collects a number of biological data analytics tasks.
This library provides tools for the bio_db served data,
empowering downstream analyses of user's experimental data.

Installation and loading: 
==
?- pack_install(bio_analytics).  
        % also installs pack(lib). other dependencies are installed at load time via lib
?- use_module(library(lib)).
?- lib(bio_analytics).
==

The library comes with one experimental dataset: data/silac/bt.csv which
is used in the example files in directory examples/ .

@author nicos angelopoulos
@version  0.1 2019/4/22
@version  0.2 2019/5/08
@version  0.3 2019/5/11
@version  0.4 2020/9/18
@version  0.5 2023/1/2
@version  0.6 2023/6/6
@version  0.7 2024/10/7
@license  MIT

*/

/** bio_analytics_version( -Vers, -Date ).

Version and release date.

==
?- bio_analytics_version(V,D).
V = 0:7:0,
D = date(2024,10,17).
==

*/
bio_analytics_version( 0:7:0, date(2024,10,17) ).

:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(debug)).
:- use_module(library(filesex)).
:- use_module(library(lib)).

:- lib(mtx).
:- lib(real).
:- lib(b_real).
:- lib(bio_db).
:- lib(os_lib).
:- lib(options).
:- lib(debug_call).
:- lib(pack_errors).
:- lib(promise(wgraph_plot/2,wgraph)).

:- debug_call:debug(ba(info)).

:- lib(source(bio_analytics), homonyms(true)).

:- lib(gene_family/3).
:- lib(go_org_symbols/3).
:- lib(go_symbols_reach/3).
:- lib(go_string_graph/3).
:- lib(symbols_string_graph/3).
:- lib(bio_diffex/4).
:- lib(exp_gene_family_string_graph/4).
:- lib(exp_go_over/3).
:- lib(exp_go_over_string_graphs/4).
:- lib(exp_go_over_string_graphs_multi/1).
:- lib(bio_volcano_plot/2).
:- lib(bio_symbols/3).
:- lib(bio_org/0).
:- lib(bio_p_adjust/3).
:- lib(exp_reac_over/3).
:- lib(go_over_plot/2).
:- lib(end(bio_analytics), homonyms(true)).

pack_errors:message( cannot_map_gids(Org,Gid,Gto) ) -->
    ['Predicate org_gid_map/3 could not convert gids from: ~w to: ~w in organism: ~w.'-[Gid,Gto,Org]].
