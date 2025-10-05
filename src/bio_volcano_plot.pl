
:- lib(b_real:options_rvar_rmv/2).
:- lib(stoics_lib:at_con/3).
:- lib(stoics_lib:lexi/2).
:- lib(stoics_lib:colour_hex/2).

bio_volcano_plot_defaults( Args, Defs ) :-
     bio_analytics:bio_diffex_defaults( Args, DfxDefs ),
     Defs = [
               clr_down(brandeisblue),
               clr_inv(cadetgrey),
               clr_line(bazaar),
               clr_up(cadmiumred),
               ext(pdf),
               debug(true),
               labels(false),
               legend_cnm(regulation),
               legend_down("down"),
               legend_inv("invariant"),
               legend_title(_),
               legend_up("up"),
               % lim_x_min(_),
               % lim_x_max(_),
               lim_y_min(0),
               % lim_y_max(_),
               plot_file(_),
               repel_box_pad(0.5),
               repel_max_overlaps('Inf'),
               rvar_cv(bvp_cv),
               rvar_df(bvp_df),
               rvar_gt(bvp_gt),
               rvar_rmv(true),
               stem(volc_plot),
               theme(classic)
               | DfxDefs
     ].

/** bio_volcano_plot( +Mtx, +Opts ).

Create a volcano plot of differential expressed values (x) against p-values (y).

Uses ggplot2 and bio_analytics:bio_diffex/4. 
The predicate creates an R data.frame and creates a ggplot2 object from which
typically output files are produced.

Opts
  * clr_down(ClrDw=brandeisblue)
    bio_colour_hex/2 compatible colour for down regulated entries
  * clr_inv(ClrInv=cadetgrey)
    colour for invariant entries
  * clr_line(ClrLn=bazaar)
    colour for cut-off lines
  * clr_up(ClrUp=cadmiumred)
    colour for up regulated entries
  * dir(Dir)
    as understood by os_dir_stem_ext/2- see Odir, below too
  * debug(Dbg=true)
    informational, progress messages
  * ext(Ext=pdf)
    extension that defines the type of output (passed to os_dir_stem_ext/2)
  * labels(Lbls=false)
    labels column id, or =|false|= or free variable for no labels
  * legend_cnm(LgCnm='regulation')
    the legend frame column name (also appears as legend title, if none is given)
  * legend_down(LegUp="down")
    token for up-regulation
  * legend_inv(LegUp="invariant")
    token for invariants
  * legend_title(LgTitle)
    by default LgCnm is used.
  * legend_up(LegUp="up")
    token for up-regulation
  * lim_x_max(LimXmax)
    if LimXmin is given this defaults to -LimXmin. By default ggplot2 decides this.
  * lim_x_min(LimXmin)
    if LimXmax is given this defaults to -LimXmax. By default ggplot2 decides this.
  * lim_y_max(LimXmax)
    By default ggplot2 decides this.
  * lim_y_min(LimY=0)
    this is not used unless LimYmax is given
  * odir(Odir)
    as understood by os_dir_stem_ext/2- preferred to Dir above
  * plot_file(File)
    returns the output file
  * repel_box_pad(Rbp=0.5)
    box.padding argument for geom_text_repel() R function from library(ggrepel)
  * repel_max_overlaps(Rmo='Inf')
    max.overlaps argument for geom_text_repel() R function from library(ggrepel)
  * rvar_cv(Rvar=bvp_cv)
    R var for colors vector
  * rvar_df(Rvar=bvp_df)
    R variable for the 3 column data frame
  * rvar_gt(Rvar=bvp_gt)
    R variable holding the plot term
  * rvar_rmv(RRmv=true)
    remove the three R variables from above
  * stem(Stem=volc_plot)
    stem for output, set to =|false|= if no output is wanted
  * theme(Theme=classic)
    theme_Theme should be a ggplot() theme

Opts are passed to bio_diffex/4. This predicate also pinches the defaults from there.

Examples, fixme: use iris or other R dataset.
==
?- bio_volcano_plot([]).
==

As of v0.2 labels can be added by including a label column and passing option 

@author nicos angelopoulos
@version  0.1 2022/12/15
@version  0.2 2025/10/5   added support for labels with options: labels(Lbls), repel_box_padding(Rbp) and repel_max_overlaps(Rmo)
@see bio_diffex/4
@see os_dir_stem_ext/2
@tbd add args for influencing of the gg term and the output-to-file call

*/
bio_volcano_plot( Mtx, Args ) :-
     Self = bio_volcano_plot,
     <- library("ggplot2"),
     options_append( Self, Args, Opts ),
     % get the up-down regulated row indices
     bio_diffex( Mtx, _, _, [which(dx(Up,Dw))|Opts] ),
     debuc( Self, length, [up,down]/[Up,Dw] ),
     % construct the data.frame
     options( exp_ev_cnm(EvCnm), Opts ),
     options( exp_pv_cnm(PvCnm), Opts ),
     mtx_column( Mtx, EvCnm, Evs ),
     mtx_column( Mtx, PvCnm, Pvs ),
     options( legend_cnm(LgCnm), Opts ),
     options( legend_inv(LgIv), Opts ),
     options( [rvar_cv(Rcv),rvar_df(Rdf),rvar_gt(Rgt)], Opts ),
     options( labels(LblsCnm), Opts ),
     mtx_column( Mtx, LblsCnm, Lbls ),
     ( (var(Lbls);Lbls==false) ->
          HasLbls = false,
          Rdf <- 'data.frame'(EvCnm=Evs,PvCnm=Pvs,LgCnm=LgIv),
          Aes = aes(x=EvCnm, y= -log10(PvCnm), col=LgCnmAtm)
          ;
          <- library("ggrepel"),
          HasLbls = true,
          Rdf <- 'data.frame'(EvCnm=Evs,PvCnm=Pvs,LgCnm=LgIv,labels=Lbls),
          Aes = aes(x=EvCnm, y= -log10(PvCnm), col=LgCnmAtm, label=labels)
     ),
     options( legend_up(LgUp), Opts ),
     options( legend_down(LgDw), Opts ),
     lexi( LgCnm, +(LgCnmAtm) ),
     Rdf$LgCnmAtm <- replace( Rdf$LgCnmAtm, Up, LgUp ),
     Rdf$LgCnmAtm <- replace( Rdf$LgCnmAtm, Dw, LgDw ),
     Rdf$LgCnmAtm <- factor(Rdf$LgCnmAtm, levels=c(LgUp,LgIv,LgDw)),
     % build the plot term
     options( theme(GGthemeTkn), Opts ),
     atom_concat( theme_, GGthemeTkn, GGthemeNm ),
     atom_concat( GGthemeNm, '()', GGthemeTerm ),
     Rgt <- ggplot(data=Rdf, Aes) + geom_point() + GGthemeTerm,
     ( HasLbls == true ->
               % Rgt <- Rgt + geom_text(hjust=0, vjust=0)
               options( repel_box_pad(Rbp), Opts ),
               Rgt <- Rgt + geom_text_repel('box.padding' = Rbp, 'max.overlaps' = 'Inf')
               ;
               true
     ),
     options( [clr_down(ClrDw),clr_inv(ClrIv),clr_line(ClrLn),clr_up(ClrUp)], Opts ),
     maplist( colour_hex, [ClrDw,ClrIv,ClrLn,ClrUp], [HexDw,HexIv,HexLn,HexUp] ),
     options( [exp_ev_cut_let(CutLw),exp_ev_cut_get(CutUp),exp_pv_cut(CutPv)], Opts ),
     Rgt <- Rgt + geom_vline(xintercept=c(CutLw,CutUp), col=HexLn)
                      + geom_hline(yintercept= -log10(CutPv), col=HexLn),
     bio_volcano_gg_limits_term( Rgt, Opts ),
     Rcv <- c(HexUp,HexIv,HexDw),
     names(Rcv) <- c(+LgUp,+LgIv,+LgDw),
     options( legend_title(LgTtl), Opts ),
     ( var(LgTtl) -> LgTtl = LgCnm; true ),
     Rgt <- Rgt + scale_color_manual(+LgTtl,values=Rcv),

     % options( plot_save_width(Width), Opts ),
     bio_volcano_plot_save( Opts ),
     % fixme: keep rvars ? 
     options_rvar_rmv( [Rcv,Rdf,Rgt], Opts ),
     debuc( Self, end, true ).

bio_volcano_plot_save( Opts ) :-
     options( stem(false), Opts ),
     !.
bio_volcano_plot_save( Opts ) :-
     os_dir_stem_ext( File, Opts ),
     options( plot_file(File), Opts ),
     % <- ggsave( +File, width=Width )
     % fixme: allow for arbitary = pairs to be passed ?
     <- ggsave( +File ).

bio_volcano_gg_limits_term( Rv, Opts ) :-
     bio_volcano_gg_limits_term_x( Xt, Opts ),
     bio_volcano_gg_limits_term_y( Yt, Opts ),
     bio_volcano_gg_limits_term_conc( Xt, Yt, Rv ).

bio_volcano_gg_limits_term_conc( null, null, _Rv ) :-  !.
bio_volcano_gg_limits_term_conc( null, Yt, Rv ) :-  !,
     Rv <- Rv + coord_cartesian(ylim=Yt).
bio_volcano_gg_limits_term_conc( Xt, null, Rv ) :-  !,
     Rv <- Rv + coord_cartesian(xlim=Xt).
bio_volcano_gg_limits_term_conc( Xt, Yt, Rv ) :-
     Rv <- Rv + coord_cartesian(xlim=Xt, ylim=Yt).

bio_volcano_gg_limits_term_x( Xt, Opts ) :-
     bio_volcano_gg_limits_term_x_given( Low, High, Opts ),
     !,
     Xt = c(Low,High).
bio_volcano_gg_limits_term_x( null, _Opts ).

bio_volcano_gg_limits_term_x_given( Low, High, Opts ) :-
    memberchk( lim_x_max(High), Opts ),
    memberchk( lim_x_min(Low), Opts ),
    !.
bio_volcano_gg_limits_term_x_given( Low, High, Opts ) :-
    memberchk( lim_x_max(High), Opts ),
    !,
    Low is - High.
bio_volcano_gg_limits_term_x_given( Low, High, Opts ) :-
    memberchk( lim_x_min(Low), Opts ),
    !,
    High is - Low.

bio_volcano_gg_limits_term_y( Yt, Opts ) :-
     bio_volcano_gg_limits_term_y_given( Low, High, Opts ),
     !,
     Yt = c(Low,High).
bio_volcano_gg_limits_term_y( null, _Opts ).

bio_volcano_gg_limits_term_y_given( Low, High, Opts ) :-
    memberchk( lim_y_max(High), Opts ),
    !,
    options( lim_y_min(Low), Opts ).
