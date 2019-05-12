%% un_list_singleton( +Either, -Element ).
%
% True iff Removed is the element of singleton list Either, 
% or is unified with Either if the latter is not a singleton list.
%
% @author nicos angelopoulos
% @version  0.1 2014/4/30  (time before)
%
un_list_singleton( MayBeSingleton, ListRmved ) :-
     ( (\+ var(MayBeSingleton),MayBeSingleton=[ListRmved]) ->
          true
          ;
          ListRmved = MayBeSingleton
     ).
     
%% un_list_singleton( +Guide, +Either, -Element ).
%
% Element is the single element of singleton Either when Guide is not a list.
% Else, Element is unified to Either. This is the permissive version.
% A stricter version could fail, if Guide is not a list and Either is not
% singleton.
%
% This predicate allows polymorphism on predicate's inputs and outputs.
% Matching input Guide to output Element.
%
% @author nicos angelopoulos
% @version  0.1 2014/4/
%
un_list_singleton( Guide, Either, Element ) :-
	\+ var( Either ),
	\+ var( Guide ),
	\+ Guide = [_|_],
	Either = [Element],
	!.
un_list_singleton( _Guide, Either, Either ).
