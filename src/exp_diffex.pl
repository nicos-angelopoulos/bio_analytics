
:- lib(stoics_lib:kv_decompose/3).

exp_diffex_defaults( Args, Defs ) :-
    % this makes EvLog default: true
    ( (memberchk(exp_ev_log(EvLog),Args),EvLog==false) ->
        EvDefs = [  as_non(pvalue),
                    exp_ev_cnm('experiment'),
                    exp_ev_cut_let(-1),
                    exp_ev_cut_get(2),
                    exp_ev_include_inf(false),
                    exp_ev_log(false)
            ]
        ; 
        EvDefs = [  as_non(all),
                    exp_ev_cnm('log2FC'),
                    exp_ev_cut_let(-1),
                    exp_ev_cut_get(1),
                    exp_ev_include_inf(true),
                    exp_ev_log(true)
            ]
    ),
    Defs = [   
                as_pairs(true),
                diffex_max(false),
                exp_pv_cnm('adj.pvalue'),
                exp_pv_cut(0.05),
                gene_id_cnm('Symbols')
                | EvDefs
    ].

/** exp_diffex( +Exp, -DEs, -NonDEs, +Opts ).

Select a sublist of significantly, differentially expressed genes in an experiment.

Exp is in the form of a mtx/1 matrix. Selected (DEs) and rejected (NonDEs) are returned in common format, 
of either (default) Gene-Ev, where Ev is the expression value for Gene in Exp, or sub-matrices of Exp.
The format is controlled by option =|as_pairs|=. 

We use ev as short for expression value, to avoid confusion with _exp_ which short for experiment.

Opts
  * as_non(AsNon=pvalue)
    what to return as non differentially expressed: either everything in Gcnm (AsNon=all, default when EvLog is true), or only those with numeric pvalue (AsNon=pvalue, default when EvLog is false)
  * as_pairs(AsPairs=true)
    whether to return pairs or matrices
  * diffex_max(DexMax=false)
    puts a cap on the number of differentially expressed genes returned
  * exp_pv_cnm(ExpPcnm='adj.pvalue')
    the experimental column (found in MsF) on which Pcut is applied as a filter
  * exp_pv_cut(Pcut=0.05)
    p.value cut off for experimental input from MSstats
  * exp_ev_log(EvLog=true)
    are the expression values log values 
  * exp_ev_cnm(EvCnm)
    expession value column name (_log2FC_ if EvLog=true, _expression_ otherwise)
  * exp_ev_cut_let(EvLet= -1)
    expression values below (less or equal) which values are selected
  * exp_ev_cut_get(EvGet)
    expression values above (greater or equal) which values are selected (defaults is 2 if EvLog is false, and 1 otherwise)
  * exp_ev_include_inf(IncInf=false)
    include infinity values as diffexs ? (default is _false_ if EvLog is false, and _true_ otherwise)
  * gene_id_cnm(Gcnm='Symbols')
    column name for the key value in the pair lists: DEs and NonDEs.

==
?- lib(mtx),
   absolute_file_name(pack('bio_analytics/data/silac/bt.csv'), CsvF),
   mtx(CsvF, Mtx, convert(true)),
   assert(mtx_data(Mtx)).

CsvF = '/usr/local/users/nicos/local/git/lib/swipl/pack/bio_analytics/data/silac/bt.csv',
Mtx = [row('Protein IDs', 'Symbols', log2FC, adj.pvalue), row('Q9P126', ...), row(..., ..., ..., ...)|...].

?- lib(debug_call),
   debug(testo),
   mtx_data(Mtx),
   debug_call(testo, dims, mtx/Mtx),
   exp_diffex(Mtx, DEPrs, NonDEPrs, []),
   debug_call(testo, length, de/DEPrs),
   debug_call(testo, length, nde/NonDEPrs).

% Dimensions for matrix,  (mtx) nR: 1245, nC: 4.
% Length for list, de: 426.
% Length for list, nde: 818.

DEPrs = ['CLEC1B'- -4.80614893171469, 'LGALS1'- -2.065877096524, ... - ...|...],
NonDEPrs = ['CNN2'- -0.69348026828097, 'CXCR4'-0.73039667221395, ... - ...|...].

?- 
    mtx_data(Mtx),
    exp_diffex(Mtx, DEPrs, NonDEPrs, exp_pv_cut(0.01)),
    debug_call(testo, length, de/DEPrs),
    debug_call(testo, length, nde/NonDEPrs).

% Length for list, de: 286.
% Length for list, nde: 958.

?-
    mtx_data(Mtx),
    Opts = [exp_ev_cut_let(inf),exp_ev_cut_get(-inf)],
    exp_diffex(Mtx, DEPrs, NonDEPrs, Opts),
    debug_call(testo, length, de/DEPrs),
    debug_call(testo, length, nde/NonDEPrs).

% Length for list, de: 581.
% Length for list, nde: 663.

?- mtx_data(Mtx),
    Opts = [exp_ev_cut_let(inf),exp_ev_cut_get(-inf),as_pairs(false)],
    exp_diffex(Mtx, DEs, NonDEs, Opts),
    debug_call(testo, length, de/DEPrs),
    debug_call(testo, length, nde/NonDEPrs).

% Length for list, de: 582.
% Length for list, nde: 664.

Mtx = [row('Protein IDs', 'Symbols', log2FC, adj.pvalue), row(..., ..., ..., ...)|...],
Opts = [exp_ev_cut_let(inf), exp_ev_cut_get(-inf), as_pairs(false)],
DEs = [row('Protein IDs', 'Symbols', log2FC, adj.pvalue), row('Q9P126', 'CLEC1B', -4.80614893171469, 2.057e-37), row(..., ..., ..., ...)|...],
NonDEs = [row('Protein IDs', 'Symbols', log2FC, adj.pvalue), row('B4DUT8;Q6FHC3;Q6FHE4;Q99439;B4DDF4;Q53GK7;B4DN57;B4DHU5;K7ES69;H3BVI6;H3BQH0', 'CNN2', -0.69348026828097, 0.0814675), row(..., ..., ..., ...)|...].

==

@author nicos angelopoulos
@version  0.1 2019/5/2
@version  0.2 2020/9/3,  ability to return sub-matrices

*/

exp_diffex( MtxIn, DEs, NonDEs, Args ) :-
    Self = exp_diffex,
    options_append( Self, Args, Opts ),
    options( as_non(AsNon), Opts ),
    options( as_pairs(AsPrs), Opts ),
    options( exp_pv_cnm(ExpPCnm), Opts ),
    options( exp_pv_cut(Pcof), Opts ),
    options( exp_ev_cnm(ExpEvCnm), Opts ),
    options( exp_ev_cut_let(EvLet), Opts ),
    options( exp_ev_cut_get(EvGet), Opts ),
    options( exp_ev_include_inf(InfInc), Opts ),
    options( gene_id_cnm(Gcnm), Opts ),
    options( diffex_max(DEMaxPrv), Opts ),
    mtx( MtxIn, Mtx, convert(true) ),
    Mtx = [Hdr|Rows],
    Cids = [ExpPCnm,ExpEvCnm,Gcnm],
    CPos = [PvPos,EvPos,GnPos],
    maplist( mtx_header_column_name_pos(Hdr), Cids, _Cnms, CPos ),
    options( diffex_max(DEMaxPrv), Opts ),
    ( number(DEMaxPrv) -> DEMaxPrv = DEMax; length(Rows,RsLen), DEMax is (RsLen + 2) * 2 ),  % ideally we want infinity ...
    exp_diffex_separate( Rows, DEMax, PvPos, EvPos, GnPos, Pcof, EvLet, EvGet, InfInc, AsNon, AsPrs, TDEs, TNonDEs ),
    exp_diff_add_header( AsPrs, Hdr, TDEs, TNonDEs, DEs, NonDEs ).

exp_diffex_separate( [], _I, _PvPos, _EvPos, _GnPos, _Pcof, _EvLet, _EvGet, _InfInc, _AsNon, _AsPrs, [], [] ).
exp_diffex_separate( [Row|Rows], I, PvPos, EvPos, GnPos, Pcof, EvLet, EvGet, InfInc, AsNon, AsPrs, DEs, NonDEs ) :-
    arg( PvPos, Row, Pv ), 
    arg( EvPos, Row, Ev ),
    arg( GnPos, Row, Gn ),
    ( (I>0,exp_diffex_select(Pv,Ev,Pcof,EvLet,EvGet,InfInc)) ->
        exp_diffex_ret_elem( AsPrs, Row, Gn, Ev, DEs, TDEs ),
        J is I - 1,
        NonDEs = TNonDEs
        ;
        J is I,
        DEs = TDEs,
        ( (AsNon == pvalue, \+ number(Pv)) -> NonDEs = TNonDEs; exp_diffex_ret_elem(AsPrs,Row,Gn,Ev,NonDEs,TNonDEs)
                % NonDEs = [Gn-Ev|TNonDEs] 
        )
    ),
    exp_diffex_separate( Rows, J, PvPos, EvPos, GnPos, Pcof, EvLet, EvGet, InfInc, AsNon, AsPrs, TDEs, TNonDEs ).

exp_diff_add_header( false, Hdr, TDEs, TNonDEs, DEs, NonDEs ) :-
    !,
    DEs = [Hdr|TDEs],
    NonDEs = [Hdr|TNonDEs].
exp_diff_add_header( _Defaulty, _Hdr, TDEs, TNonDEs, DEs, NonDEs ) :-
    DEs = TDEs,
    NonDEs = TNonDEs.

exp_diffex_ret_elem( false, Row, _Gn, _Ev, List, Tail ) :-
    !,
    List = [Row|Tail].
exp_diffex_ret_elem( _Defaulty, _Row, Gn, Ev, List, Tail ) :-
    % DEs = [Gn-Ev|TDEs],
    List = [Gn-Ev|Tail].

exp_diffex_select( Pv, Ev, Pcof, EvLet, EvGet, InfInc ) :-
    number( Pv ), 
    Pv < Pcof,
    number(Ev),
    ( Ev =< EvLet ; EvGet =< Ev ),
    exp_diffex_if_inf_include( Ev, InfInc ).

exp_diffex_if_inf_include( Val, Incl ) :-
    ( Val is inf; Val is - inf ),
    !,
    Incl == true.
exp_diffex_if_inf_include( _Val, _Incl ).
