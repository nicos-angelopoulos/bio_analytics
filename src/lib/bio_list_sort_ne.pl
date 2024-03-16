bio_list_sort_ne( List, SetNe ) :-
    sort( List, Set ),
    ( select('',Set,SetNe) -> true; Set = SetNe ).
