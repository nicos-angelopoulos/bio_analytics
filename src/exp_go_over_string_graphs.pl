
exp_go_over_string_graphs_defaults( Defs ) :-
    Wplots = [vjust= -1, node_size(3), format(svg)],
    Defs = [ dir_postfix(go_strings), extra_symbols_column(0), go_id_clm(1),
             ov_max(false),
             stem_type(go_pair_ord),
             symbols_column(0),
             viz_de_opts([]), wgraph_plot_opts(Wplots) ].

/** exp_go_over_string_graphs( +Exp, ?GoOver, ?Dir, -Opts ).

Create string graphs for all the over-represented Exp-eriment terms in GoOver.

When GoOver is a variable exp_go_over/3 is called to generate it.

Opts
  * dir_postfix(Psfx=go_strings)
    postfix for outputs directory (when Dir is a variable)
  * extra_symbols_column(Xsymc=0)
    additional symbols for the graph
  * go_id_clm(GoIdClm=1)
    column id for over represented GO terms
  * ov_max(MaxOvs=false)
    when a number, it is taken as the maximal integer of terms to plot graphs for
  * symbols_column(Cymbs=0)
    any non-zero value is the column containing Symb1;Symb2... values which are pushed down 
    instead of the GOid means of identifying the family members (GOid is still used for filename)
  * stem_type(Sty=go_pair_ord)
    similar to go_string_graph/3, but different default, others: =|go_name, go_id|=.
    Here the length of GO terms (Len) is added to form go_pair_ord(I,Len) for forwarding
  * viz_de_opts(VizOpts = [])
    options for restricting genes to visualise via bio_diffex/4.
    Default does not restrict what genes are visualised. =|diffex_only = VizOpts|= restricts genes to significants only.
  * wgraph_plot_opts(WgOpts=WgOpts)
    defaults to =|[vjust = -1, node_size(3), format(svg)]|=.


Options are passed to exp_gene_family_string_graph/4.

==
?- debug(real).
?- debug(exp_go_over_string_graphs).
?- absolute_file_name( pack('bio_analytics/data/silac/bt.csv'), Exp ),
   exp_go_over_string_graphs( Exp, Gov, Dir, [] ).

% Sending to R: pltv <- ggnet2(lp_adj,vjexp_go_over_mtxust = -1,size = 3,
...    label = pl_v_1,color = pl_v_2,edge.size = pl_v_3,edge.color = "#BEAED4")

==

@author nicos angelopoulos
@version  0.1 2019/5/5
@version  0.2 2020/9/6,     option ov_max(OvMax)
@version  0.3 2024/10/20,   options extra_symbols_column(Xsymc), symbols_column(Cymbs)
@see exp_go_over/3
@see exp_gene_family_string_graph/4

*/
exp_go_over_string_graphs( Exp, GoOverIn, Dir, Args ) :-
    Self = exp_go_over_string_graphs,
    options_append( exp_go_over_string_graphs, Args, Opts ),
    options( dir_postfix(Psfx), Opts ),
    exp_go_over_string_graphs_to_dir( Exp, Psfx, Dir ),
    os_make_path( Dir ),
    exp_go_over_mtx( GoOverIn, Exp, GoOver, Dir, Self, Opts ),
    options( go_id_clm(Cid), Opts ),
    mtx_column( GoOver, Cid, GOs ),
    options( symbols_column(Cymb), Opts ),
    ( (number(Cymb),Cymb =:= 0) -> 
          GOPrs = GOs
          ;
          mtx_column( GoOver, Cymb, SymbsClm ),
          maplist( atomic_list_deconcat(';'), SymbsClm, SymbsLs ),
          kv_compose( GOs, SymbsLs, GOPrs )
    ),
    options( extra_symbols_column(Xymb), Opts ),
    ( (number(Xymb),Xymb =:= 0) ->
          findall( [], member(_Ago,GOs), XymbsLs )
          ;
          mtx_column( GoOver, Xymb, XymbsClm ),
          maplist( atomic_list_deconcat(';'), XymbsClm, XymbsLs )
    ),
    debuc( Self, 'Dir: ~p', [Dir] ),
    options( viz_de_opts(VdfOpts), Opts ),
    options( wgraph_plot_opts(WgOpts), Opts ),
    % ( VdfOpts == [] -> Exp = RedExp; bio_diffex(Exp,RedExp,_,[as_pairs(false)|VdfOpts]) ),
    % VdfOpts -> [] or diffex_only  succeed on the following
    ( atomic(VdfOpts) -> 
          debuc( Self, 'Not reducing the input matrix.', [] ),
          Exp = RedExp
          ;
          bio_diffex(Exp,RedExp,_,[as_pairs(false)|VdfOpts]),
          debuc( Self, length, [exp_in,exp_red]/[Exp,RedExp] )
    ),
    % mtx( red_exp.csv, RedExp ), % fixme: do it properly (in subdirectory)
    options( ov_max(MaxOvsPrv), Opts ),
    (number(MaxOvsPrv) -> MaxOvs is integer(MaxOvsPrv),findall(Go,(between(1,MaxOvs,I),nth1(I,GOPrs,Go)),MaxGOs); GOPrs = MaxGOs), 
    options( stem_type(Sty), Opts ),
    ( Sty == go_pair_ord ->
        length( MaxGOs, MaxGOsLen ),
        number_codes( MaxGOsLen, MaxGOsLenCs ),
        length( MaxGOsLenCs, PadLen )
        ;
        PadLen = 0
    ),
    bio_diffex( RedExp, DEPrs, NonDEPrs, [diffex_mtx(DiffMtx),as_pairs(true)|Opts] ),
    debuc( Self, length, [de_prs,non_de_prs]/[DEPrs,NonDEPrs] ),
    ( memberchk(diffex_mtx(PpgUpDiffMtx),Opts ) -> mtx( PpgUpDiffMtx, DiffMtx ) ; true ),
    ( VdfOpts == diffex_only -> RedNonDEPrs = [] ; RedNonDEPrs = NonDEPrs ),
    debuc( Self, length, go_over_non_de_prs/RedNonDEPrs ),
    go_over_string_graphs_dir( MaxGOs, XymbsLs, 1, RedExp, Dir, WgOpts, Sty, PadLen, DEPrs, RedNonDEPrs, Self, Opts ).

go_over_string_graphs_dir( [], [], _I, _Exp, _Dir, _WgOpts, _Sty, _Pad, _DEPrs, _NonDEPrs, _Self, _Opts ).
go_over_string_graphs_dir( [GoPrv|Gos], [Xymb|Xymbs], I, Exp, Dir, WgOpts, Sty, Pad, DEPrs, NonDEPrs, Self, Opts ) :-
    % ( atom_concat('GO:',Stem,Go) -> atom_concat(go,Stem,GoTkn); GoTkn=Go ),
    ( GoPrv = Go-GoSymbs -> Fam = GoSymbs; Fam = GoPrv, Go=GoPrv ),
     
    ( Sty == go_pair_ord -> Tty = go_pair_ord(I,Pad) ; Tty = Sty ),
    go_string_graph_stem( Tty, Go, GoTkn ),
    directory_file_path( Dir, GoTkn, DirGo ),
    debuc( Self, 'GOid: ~w, File: ~p', [Go,GoTkn] ),
    GoWgOpts = [stem(DirGo)|WgOpts],
    ( (is_list(Xymb), Xymb \== []) -> GnOpts=[extra_symbols(Xymb)|Opts] ; GnOpts = Opts ),
    % here check if we are carrying the symbols and extras and plot directly
    exp_gene_family_string_graph( Exp, Fam, DEPrs, NonDEPrs, _, [wgraph_plot_opts(GoWgOpts)|GnOpts] ),
    J is I + 1,
    go_over_string_graphs_dir( Gos, Xymbs, J, Exp, Dir, WgOpts, Sty, Pad, DEPrs, NonDEPrs, Self, Opts ).

exp_go_over_mtx( GoOverIn, Exp, GoOver, Dir, _Self, Opts ) :-
    var( GoOverIn ),
    !,
    os_path( Dir, 'go_over', + Stem ),
    append( Opts, [to_file(true),stem(Stem)], GoOpts ),
    exp_go_over( Exp, GoOver, GoOpts ),
    GoOverIn = GoOver.
exp_go_over_mtx( GoOverIn, _Exp, GoOver, _Dir, Self, _Opts ) :-
    debuc( Self, 'Using stored go_over results from: ~p', [GoOverIn] ),
    mtx( GoOverIn, GoOver ). % fixme: pass Opts ?

exp_go_over_string_graphs_to_dir( _Exp, _Psfx, Dir ) :-
    ground( Dir ),
    !.
exp_go_over_string_graphs_to_dir( Exp, Psfx, Dir ) :-
    file_name_extension( Stem, _Ext, Exp ),  % fixme: Ext = csv
    atom_concat( '_', Psfx, PsfxTkn ),
    atom_concat( Stem, PsfxTkn, Dir ).

atomic_list_deconcat( Sep, Atom, List ) :-
     at_con( ListPrv, Sep, Atom ),
     ( ListPrv = [''] -> List = []; ListPrv = List ).
