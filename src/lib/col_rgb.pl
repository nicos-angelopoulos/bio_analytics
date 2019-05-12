:- lib(decimal_hexadecimal/2).
:- lib(stoics_lib:n_digits_min/3).
:- lib(stoics_lib:en_list/2).

%% col_rgb( +Colours, RGBtype, -RGBs ).
%% col_rgb( +Colours, -RGBs ).
%
% Convert a Colour to corresponding RGB term of RGBtype. Colours should be
% digestable by R as a vector. If RGBs is singleton the list is stripped except if
% Colours was a list.
% 
% RGBtype should define the form of RGBs. 
%  * *hexed* return a string or strings of the form "#RRGGBB"
%  * hexed_atom  at hexed but atom(s) rather than strings
%  * hex    return string without leading #
%  * term_hex  return rgb(RR,GG,BB), with each colour as hex atom
%  * term_dec return rgb(RR,GG,BB), with each RGB as decimal number this is the
%     cannonical representation, see rgb_term/2.
%  * dec    to decimal strings
%  * dec_atom  as decimal atoms
% 
%==
% ?- CC <- colours(), length( CC, Len ), col_rgb( CC, RGBs ), length( RGBs, Len ).
% ?- requires( colours_pie/1 ).
% ?- PP <- palette(), length( PP, Len ), col_rgb( PP, RGBs ), length( RGBs, Len ), colours_pie( RGBs ), maplist( rgb_term, RGBs, Terms ).
%
%==
%
% @author nicos angelopoulos
% @version  0.1 2014/03/17
% @tbd      not all RGBType types are implemented...
% @tbd      documentation for other predicates
%
col_rgb( Clrs, RGBs ) :-
	ground( Clrs ),
	!,
	col_rgb( Clrs, hexed, RGBs ).
col_rgb( Clrs, RGBS ) :-
	en_list( RGBS, RGBs ),
	findall( Clr, (member(rgb(R,G,B),RGBs),Clr <- rgb(R/255,G/255,B/255)), ClrsList ),
	( is_list(RGBS) -> Clrs = ClrsList; ClrsList = [Clrs] ).

col_rgb( Clrs, Type, RGBs ) :-
	RGBNest <- col2rgb( Clrs ),
	RGBNest = [Rs,Gs,Bs],
	findall( rgb(R,G,B), (nth1(N,Rs,R),nth1(N,Gs,G),nth1(N,Bs,B)), RGBterms ), 
	maplist( type_term_rgb(Type), RGBterms, RGBsList ),
	( is_list(Clrs) -> RGBs = RGBsList;
		(RGBsList = [RGBs] -> true; 
			RGBs = RGBsList 
		)
	).

type_term_rgb( hexed, rgb(R,G,B), RGB ) :-	
	maplist( decimal_hexadecimal, [R,G,B], Hexas ),
	maplist( n_digits_min(2), Hexas, Nexas ),
	atomic_list_concat( ['#'|Nexas], Atom ),
	atom_string( Atom, RGB ).

type_term_rgb( term_dec, RGB, RGB ).

rgb_term( RGB, Term ) :-
	term_codes( RGB, RGBcodes ),
	% type_rgb_hash_bool( Type, RGB, Hash ),
	codes_rgb_term( RGBcodes, Term ).

codes_rgb_term( [0'#|Codes], RGB ) :-
	!,
	non_hash_codes_rgb_term( Codes, RGB ).
codes_rgb_term( Codes, RGB ) :-
	non_hash_codes_rgb_term( Codes, RGB ).

non_hash_codes_rgb_term( Codes, rgb(R,G,B) ) :-
	Codes = [R1,R2,G1,G2,B1,B2],
	decimal_hexadecimal( R, [R1,R2] ),
	decimal_hexadecimal( G, [G1,G2] ),
	decimal_hexadecimal( B, [B1,B2] ).
	/* number_codes( R, [R1,R2] ),
	   number_codes( G, [G1,G2] ),
	   number_codes( B, [B1,B2] ).
	   */
