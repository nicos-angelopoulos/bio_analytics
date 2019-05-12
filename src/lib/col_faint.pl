
:- lib(col_rgb/3).
:- lib(faint/3).

%% col_faint( +Colours, +By, -Faints ).
%
% Faint a Colour or Colours by a factor of By.
%
% Faints is a list if Colours was one.
%
% See col_rgb/3 and faint/3.
%
%==
% ?- 
%==
col_faint( Clrs, By, Faints ) :-
	col_rgb( Clrs, term_dec, Rgbs ),
	faint( Rgbs, By, RgbFaints ),
	col_rgb( Faints, RgbFaints ).
