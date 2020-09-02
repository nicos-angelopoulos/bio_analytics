:- module( bio_analytics, [
                bio_analytics_version/2,        % -Vers, -Date
                exp_diffex/4,                   % +Exp, -DEs -NDEs, +Opts
                exp_gene_family_string_graph/4, % +Exp, +Fam, -Graph, +Opts
                exp_go_over/3,                  % +Exp, NDEs, +Opts
                exp_go_over_string_graphs/4,    % +Exp, ?GoOver, ?Dir, -Opts
                gene_family/3,                  % +Alias, +Org, +Symbols
                go_org_symbol/3,                % +Org, +GoT, -Symbol
                go_org_symbols/3,               % +Org, +GoT, -Symbols
                go_string_graph/3,              % +GoT, -Graph, +Opts
                go_symbols_reach/3,             % +GoT, -Symbs, +Opts
                symbols_string_graph/3          % +Symbs, -Graph, +Opts
        ] ).

/** <module> Computational biology data analytics.

Collects a number of biological data analytics tasks.
This library provides tools for the bio_db served data,
empowering downstream analyses of user's experimental data.

Installation and loading: 
==
?- pack_install(bio_analytics).  % also installs pack(lib). other dependencies are at load time via lib
?- use_module(library(lib)).
?- lib(bio_analytics).
==

The library comes with one experimental dataset: data/silac/bt.csv which
is used in the example files in directory examples/ .

@author nicos angelopoulos
@version  0.1 2019/4/22
@version  0.2 2019/5/08
@version  0.3 2019/5/11
@license  MIT

*/

/** bio_analytics_version( +Vers, +Date ).

Version and release date.

==
?- bio_analytics_version(V,D).
V = 0:3:0,
D = date(2019,5,12).
==

*/
bio_analytics_version( 0:3:0, date(2019,5,12) ).

:- ensure_loaded(library(lib)).

:- lib(bio_db).
:- lib(mtx).
:- lib(real).
:- lib(os_lib).
:- lib(wgraph).
:- lib(options).

:- lib(source(bio_analytics), homonyms(true)).

:- lib(gene_family/3).
:- lib(go_org_symbol/3).
:- lib(go_org_symbols/3).
:- lib(go_symbols_reach/3).
:- lib(go_string_graph/3).
:- lib(symbols_string_graph/3).
:- lib(exp_diffex/4).
:- lib(exp_gene_family_string_graph/4).
:- lib(exp_go_over/3).
:- lib(exp_go_over_string_graphs/4).

:- lib(end(bio_analytics), homonyms(true)).
