
:- lib( un_list_singleton/3 ).
:- lib( stoics_lib:en_list/2 ).

%% faint( +Cols, +Scale, -Faint ).
%% faint( +Cols, -Faint ).
%
% Create Faint versions of Cols RGB palette by a factor of Scale (default = 2). 
% 0 =< Scale with 1 making each Faint colour equal to White and 2 spliting the difference
% between the colour and white in to half in creating each Faint.
% Cols is a single, or a list of, rgb/3 terms (see col_rgb/3).
%
%==
% ?- use_module( library(real) ).
% ?- requires( col_rgb/3 ).
% ?- Pal <- palette(), col_rgb( Pal, hexed, PalHex ), colours_pie( PalHex ),
% 	col_rgb( Pal, term_dec, PalRGBs ), faint( PalRGBs, 1.5, New ),
% 	maplist( type_term_rgb(hexed), New, NewHex ),
% 	<- x11(),
%    colours_pie( NewHex ).
%  ?- requires( interleave/3 ).
%  ?- Pal <- palette(), col_rgb( Pal, hexed, PalHex ), colours_pie( PalHex ), 
%     col_rgb( Pal, term_dec, PalRGBs ), faint( PalRGBs, 2, New ),
%     maplist( type_term_rgb(hexed), New, NewHex ), assert( new_hex(NewHex) ),
%     interleave( PalHex, NewHex, Hexes ), colours_pie( Hexes )
%  ?- <- library( "RColorBrewer" ).
%  ?- Pal <- 'brewer.pal'(9,"Set1"), col_rgb( Pal, term_dec, PalRGBs ), faint( PalRGBs, 1.5, New ), maplist( type_term_rgb(hexed), New, NewHex ), interleave( Pal, NewHex, Inter ), colours_pie( Inter ).
%== 
%
% @author nicos angelopoulos
% @version  0.1 2014/03/17
% @see col_rgb/3
% @tbd used input(Type), output(Type)
% 
faint( Cols, New ) :-
	faint( Cols, 2, New ).

faint( ColS, By, New ) :-
	en_list( ColS, Cols ),
	debug( faint, 'In colours: ~w', [Cols] ),
	findall( rgb(Rf,Gf,Bf), (member(rgb(R,G,B),Cols),maplist(faint_scale(By),[R,G,B],[Rf,Gf,Bf])), NewList ),
	un_list_singleton( ColS, NewList, New ),
	debug( faint, 'new: ~w', [New] ).

faint_scale( By, V, New ) :-
	X is 255 - V,
	Off is integer( X * ( 1 / By) ),
	New is V + Off.
