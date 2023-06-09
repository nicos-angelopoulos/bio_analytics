
/** bio_org.

Bio analytics supports a number of organism specific predicates.

This is a doc predicate to provide access to them:
  * bio_conductor_annot_dbi_org/3

Examples
==
?- bio_org.
==

@author nicos angelopoulos
@version  0.1 2023/06/07

*/

bio_org.

/** bio_conductor_annot_dbi_org(+Org, -DbiToken, -DbiOrg).

The Annotation DBI token and organism strings for bio_db Org-anism.

==
?- bio_conductor_annot_dbi_org(human, Tkn, Dorg).
Tkn = "Hs",
Dorg = "Homo sapiens".

?- bio_db_organism(gallus,Ggallus), bio_conductor_annot_dbi_org(Ggallus, Tkn, Dorg).
Ggallus = chicken,
Tkn = "Gg",
Dorg = "Gallus gallus".
==

@author nicos angelopoulos
@version  0:1 2023/06/07

*/
bio_conductor_annot_dbi_org(chicken, "Gg", "Gallus gallus").
bio_conductor_annot_dbi_org(  human, "Hs", "Homo sapiens").
bio_conductor_annot_dbi_org(  mouse, "Mm", "Mus musculus").
bio_conductor_annot_dbi_org(    pig, "Ss", "Sus scrofa").

% maybe we should move the pred below to a bio_conductor_annot_dbi.pl source file ?

/** bio_conductor_annot_dbi_org_lib(+DbiToken, +EnsLoaded, -Lib).

From a annotation dbi, organism token, construct Dbi library and possibly load it.

When EnsLoaded is grounded to =|true\= the library is loaded. 

==
?- bio_conductor_annot_dbi_org(pig, PigTkn, PigOrg), 
   bio_conductor_annot_dbi_org_lib(PigTkn, false, PigDbiLib ).
      
PigTkn = "Ss",
PigOrg = "Sus scrofa",
PigDbiLib = "org.Ss.eg.db".
==

@author nicos angelopoulos
@version  0:1 2023/06/07
@see bio_conductor_annot_dbi_org/3
@tbd ensure loaded rather than load everytime (most likely the message says loading only)

*/
bio_conductor_annot_dbi_org_lib( DbiTkn, Load, DbiLib ) :-
     string_concat( "org.", DbiTkn, DbiLibPfx ),
     string_concat( DbiLibPfx, ".eg.db", DbiLib ),
     ( Load == true -> 
               lib( bioc(DbiLib) )      % fixme: ensure loaded, instead
               ;
               true
     ).


