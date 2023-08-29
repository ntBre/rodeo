#[derive(Debug, Default, PartialEq)]
pub enum BondDir {
    EitherDouble,
    Wedge,
    Dash,
    None,
    #[default]
    Unknown,
}

#[derive(Debug, PartialEq)]
pub enum BondStereo {
    Any,
    E,
    Z,
    None,
}

#[derive(Default, Clone, Copy)]
pub enum BondType {
    Zero,
    Ionic,
    Single,
    Aromatic,
    Double,
    Triple,
    Quadruple,
    Quintuple,
    Hextuple,
    OneAndAHalf,
    TwoAndAHalf,
    ThreeAndAHalf,
    FourAndAHalf,
    FiveAndAHalf,
    Dative,
    DativeOne,
    Hydrogen,
    #[default]
    Unspecified,
}

impl BondType {
    /// Returns `true` if the bond type is [`Double`].
    ///
    /// [`Double`]: BondType::Double
    #[must_use]
    pub fn is_double(&self) -> bool {
        matches!(self, Self::Double)
    }

    /// Returns `true` if the bond type is [`Triple`].
    ///
    /// [`Triple`]: BondType::Triple
    #[must_use]
    pub fn is_triple(&self) -> bool {
        matches!(self, Self::Triple)
    }

    /// Returns `true` if the bond type is [`Aromatic`].
    ///
    /// [`Aromatic`]: BondType::Aromatic
    #[must_use]
    pub fn is_aromatic(&self) -> bool {
        matches!(self, Self::Aromatic)
    }
}

impl From<BondType> for f64 {
    fn from(val: BondType) -> Self {
        use BondType::*;
        match val {
            Unspecified | Ionic | Zero | Hydrogen => 0.0,
            Single => 1.0,
            Double => 2.0,
            Triple => 3.0,
            Quadruple => 4.0,
            Quintuple => 5.0,
            Hextuple => 6.0,
            OneAndAHalf | Aromatic => 1.5,
            TwoAndAHalf => 2.5,
            ThreeAndAHalf => 3.5,
            FourAndAHalf => 4.5,
            FiveAndAHalf => 5.5,
            // there's a FIXME comment saying these are both probably wrong
            Dative | DativeOne => 1.0,
        }
    }
}

impl From<usize> for BondType {
    fn from(value: usize) -> Self {
        match value {
            1 => Self::Single,
            2 => Self::Double,
            3 => Self::Triple,
            4 => Self::Aromatic,
            0 => Self::Unspecified,
            _ => unimplemented!("MolFileParser.cpp:1733"),
        }
    }
}

#[derive(Default)]
pub struct Bond {
    pub(crate) bond_type: BondType,

    pub(crate) bond_dir: BondDir,

    pub(crate) react_status: usize,

    negation: bool,

    /// index in the owning molecule's bonds vector
    index: usize,

    /// first atom index in the owning molecule's atoms vector
    pub(crate) begin_atom_index: usize,

    /// second atom index in the owning molecule's atoms vector
    pub(crate) end_atom_index: usize,

    is_aromatic: bool,

    // pretty sure this is a bool, but who knows
    pub(crate) stereo_care: usize,
}

impl Bond {
    pub fn get_index(&self) -> usize {
        todo!();
    }

    pub fn get_stereo(&self) -> BondStereo {
        todo!();
    }

    pub fn set_stereo(&mut self, _stereo: BondStereo) {
        todo!();
    }

    pub fn set_index(&mut self, index: usize) {
        self.index = index;
    }

    pub fn set_begin_atom_index(&mut self, begin_atom_index: usize) {
        self.begin_atom_index = begin_atom_index;
    }

    pub fn set_end_atom_index(&mut self, end_atom_index: usize) {
        self.end_atom_index = end_atom_index;
    }

    pub fn set_is_aromatic(&mut self, is_aromatic: bool) {
        self.is_aromatic = is_aromatic;
    }

    pub fn set_bond_type(&mut self, bond_type: BondType) {
        self.bond_type = bond_type;
    }

    pub(crate) fn get_valence_contrib(&self, atom_idx: usize) -> f64 {
        if atom_idx != self.begin_atom_index && atom_idx != self.end_atom_index
        {
            return 0.0;
        }

        use BondType::*;
        if matches!(self.bond_type, Dative | DativeOne)
            && atom_idx != self.end_atom_index
        {
            0.0
        } else {
            self.bond_type.into()
        }
    }

    pub(crate) fn has_complex_bond_type(&self) -> bool {
        if !self.has_query() {
            return false;
        }
        todo!();
    }

    const fn has_query(&self) -> bool {
        false
    }

    pub fn set_bond_dir(&mut self, bond_dir: BondDir) {
        self.bond_dir = bond_dir;
    }

    pub fn set_negation(&mut self, negation: bool) {
        self.negation = negation;
    }
}
