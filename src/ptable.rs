/// default valence charge, indexed by atomic number. NOTE probably should have
/// de-duplicated this. see atomic_data.org if these need to be updated. I
/// assumed the default valence was the first entry in the allowed valences
pub(crate) const DEFAULT_VALENCE: [isize; 121] = [
    -1, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, 3, 4, 3, 2, 1, 0, 1, 2, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, 3, 2, 3, 2, 1, 0, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, 2, 3, 2, 1,
    0, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
];

/// table reporting whether or not a given atom is left of Carbon, indexed by
/// atomic number
pub(crate) const IS_EARLY_ATOM: [bool; 119] = [
    false, // #0 *
    false, // #1 H
    false, // #2 He
    true,  // #3 Li
    true,  // #4 Be
    true,  // #5 B
    false, // #6 C
    false, // #7 N
    false, // #8 O
    false, // #9 F
    false, // #10 Ne
    true,  // #11 Na
    true,  // #12 Mg
    true,  // #13 Al
    false, // #14 Si
    false, // #15 P
    false, // #16 S
    false, // #17 Cl
    false, // #18 Ar
    true,  // #19 K
    true,  // #20 Ca
    true,  // #21 Sc
    true,  // #22 Ti
    false, // #23 V
    false, // #24 Cr
    false, // #25 Mn
    false, // #26 Fe
    false, // #27 Co
    false, // #28 Ni
    false, // #29 Cu
    true,  // #30 Zn
    true,  // #31 Ga
    true,  // #32 Ge  see github #2606
    false, // #33 As
    false, // #34 Se
    false, // #35 Br
    false, // #36 Kr
    true,  // #37 Rb
    true,  // #38 Sr
    true,  // #39 Y
    true,  // #40 Zr
    true,  // #41 Nb
    false, // #42 Mo
    false, // #43 Tc
    false, // #44 Ru
    false, // #45 Rh
    false, // #46 Pd
    false, // #47 Ag
    true,  // #48 Cd
    true,  // #49 In
    true,  // #50 Sn  see github #2606
    true,  // #51 Sb  see github #2775
    false, // #52 Te
    false, // #53 I
    false, // #54 Xe
    true,  // #55 Cs
    true,  // #56 Ba
    true,  // #57 La
    true,  // #58 Ce
    true,  // #59 Pr
    true,  // #60 Nd
    true,  // #61 Pm
    false, // #62 Sm
    false, // #63 Eu
    false, // #64 Gd
    false, // #65 Tb
    false, // #66 Dy
    false, // #67 Ho
    false, // #68 Er
    false, // #69 Tm
    false, // #70 Yb
    false, // #71 Lu
    true,  // #72 Hf
    true,  // #73 Ta
    false, // #74 W
    false, // #75 Re
    false, // #76 Os
    false, // #77 Ir
    false, // #78 Pt
    false, // #79 Au
    true,  // #80 Hg
    true,  // #81 Tl
    true,  // #82 Pb  see github #2606
    true,  // #83 Bi  see github #2775
    false, // #84 Po
    false, // #85 At
    false, // #86 Rn
    true,  // #87 Fr
    true,  // #88 Ra
    true,  // #89 Ac
    true,  // #90 Th
    true,  // #91 Pa
    true,  // #92 U
    true,  // #93 Np
    false, // #94 Pu
    false, // #95 Am
    false, // #96 Cm
    false, // #97 Bk
    false, // #98 Cf
    false, // #99 Es
    false, // #100 Fm
    false, // #101 Md
    false, // #102 No
    false, // #103 Lr
    true,  // #104 Rf
    true,  // #105 Db
    true,  // #106 Sg
    true,  // #107 Bh
    true,  // #108 Hs
    true,  // #109 Mt
    true,  // #110 Ds
    true,  // #111 Rg
    true,  // #112 Cn
    true,  // #113 Nh
    true,  // #114 Fl
    true,  // #115 Mc
    true,  // #116 Lv
    true,  // #117 Ts
    true,  // #118 Og
];

/// number of outer electrons. TODO also not de-duplicated
pub(crate) const OUTER_ELECS: [usize; 121] = [
    0, 1, 2, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6,
    7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2,
    3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 3, 4, 5, 6, 7,
    8, 9, 10, 11, 12, 13, 14, 15, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2,
];

/// atomic number to symbol
pub(crate) const SYMBOL: [&str; 121] = [
    "*", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Nh", "Uut", "Fl", "Mc", "Uup", "Lv", "Ts", "Og",
];
