%% decimal_hexadecimal( -Dec, -Hex ).
%
% Convert a decimal number or atom to its hexadecimal atom or a hexadecimal
% (may it be anything that term_codes/2 can convert to codes list) to it decimal 
% number counterpart.
% 
% This is based on old code that works for 10 -> all bases
%
%==
%   between( 0, 256, I ), decimal_hexadecimal( I, Hex ), write( I:Hex ), nl, fail.
%
%   between( -32, 0, I ), decimal_hexadecimal( I, Hex ), write( I:Hex ), nl, fail.
%
%   between( 0, 256, I ), decimal_hexadecimal( I, Hex ), 
%   decimal_hexadecimal( New, Hex ), write( I:Hex:New ), nl, fail.
%==
%
% @author nicos angelopoulos
% @version  0.1 2014/03/17
% @tbd implement Hex -> Dec
%
decimal_hexadecimal( DecPrv, Hex ) :-
	ground( DecPrv ),
	!,
	atom_number_number( DecPrv, DecNum ),
	( DecNum < 0 -> Dec is DecNum * (-1); Dec is DecNum ),
	hex_digits( HexDigs ),
	decimal_to_hexadecimal( Dec, HexDigs, RevHexCodes ),
	reverse_decimal_to_hexadecimal( RevHexCodes, HexCodesPsv ),
	( DecNum < 0 -> HexCodes = [0'-|HexCodesPsv]; HexCodes = HexCodesPsv ),
	atom_codes( Hex, HexCodes ).
decimal_hexadecimal( Dec, HexPrv ) :-
	ground( HexPrv ),
    term_to_atom( HexPrv, HexAtm ),
    atom_codes( HexAtm, HexRev ),
	% term_codes( HexPrv, HexRev ),
	reverse( HexRev, Hex ),
	hex_digits( HexDigs ),
	hexadecimal_to_decimal( Hex, 1, 0, HexDigs, Dec ).

decimal_to_hexadecimal( 0, _HexDigs, [] ) :- !.
decimal_to_hexadecimal( Dec, Digs, [Hc|HCs] ) :-
	Based is Dec mod 16,
	Pos   is Based + 1,
	nth1( Pos, Digs, Hc ),
	Nxt is Dec // 16,
	decimal_to_hexadecimal( Nxt, Digs, HCs ).

reverse_decimal_to_hexadecimal( [], [0'0] ) :- !.
reverse_decimal_to_hexadecimal( Reverse, Hex ) :-
	reverse( Reverse, Hex ).

hexadecimal_to_decimal( [], _Mult, Dec, _HexDigs, Dec ).
hexadecimal_to_decimal( [H|T], Mult, Acc, HexDigs, Dec ) :-
	once( nth0(V,HexDigs,H) ),
	Nxt is Acc + (V * Mult),
	Pow is Mult * 16,
	hexadecimal_to_decimal( T, Pow, Nxt, HexDigs, Dec ).
	
hex_digits( `0123456789ABCDEF` ).

atom_number_number( Dec, Num ) :-
    atom( Dec ),
    !,
    atom_number( Dec, Num ).
atom_number_number( Dec, Num ) :-
    number( Dec ),
    !,
    Num = Dec.
atom_number_number( Dec, _Num ) :-
    throw( decimal_hexadecimal(cannot_convert_to_number(Dec)) ).
