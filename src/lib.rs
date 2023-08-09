use std::collections::HashMap;

pub enum Chi {
    TetrahedralCCW,
    TetrahedralCW,
}

pub struct Atom;

impl Atom {
    pub fn new(_num: usize) -> Self {
        Self
    }

    pub fn get_index(&self) -> usize {
        todo!();
    }

    pub fn get_stereo(&self) -> String {
        todo!();
    }

    pub fn set_formal_charge(&mut self, _charge: isize) {
        todo!()
    }

    pub fn set_is_aromatic(&mut self, _aromatic: bool) {
        todo!()
    }

    pub fn set_name(&mut self, _name: String) {
        todo!();
    }

    pub fn set_chiral_tag(&mut self, _chirality: Chi) {
        todo!();
    }
}

#[derive(Debug, PartialEq)]
pub enum BondStereo {
    E,
    Z,
    None,
}

pub enum BondType {
    Single,
    Aromatic,
    Double,
    Triple,
    Quadruple,
    Quintuple,
    Hextuple,
    OneAndAHalf,
}

pub struct Bond;

impl Bond {
    pub fn set_is_aromatic(&mut self, _aromatic: bool) {
        todo!();
    }

    pub fn set_bond_type(&mut self, _typ: BondType) {
        todo!();
    }

    pub fn get_index(&self) -> usize {
        todo!();
    }

    pub fn get_stereo(&self) -> BondStereo {
        todo!();
    }

    pub fn set_stereo(&mut self, _stereo: BondStereo) {
        todo!();
    }
}

bitflags::bitflags! {
    #[repr(transparent)]
    pub struct SanitizeOptions: u16 {
        const None = 0x0;
        const Cleanup = 0x1;
        const Properties = 0x2;
        const SymmRings = 0x4;
        const Kekulize = 0x8;
        const FindRadicals = 0x10;
        const SetAromaticity = 0x20;
        const SetConjugation = 0x40;
        const SetHybridization = 0x80;
        const CleanupChirality = 0x100;
        const AdjustHS = 0x200;
        const CleanupOrganometallics = 0x400;
        const All = 0b11111111111; // TODO not sure about this
    }
}

pub enum AromaticityModel {
    MDL,
}

pub struct RWMol {
    name: String,
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
}

impl RWMol {
    pub fn new() -> Self {
        Self {
            name: String::new(),
            atoms: Vec::new(),
            bonds: Vec::new(),
        }
    }

    pub fn set_name(&mut self, name: String) {
        self.name = name;
        todo!();
    }

    pub fn add_atom(&mut self, _atom: Atom) -> usize {
        todo!()
    }

    pub fn add_bond(&mut self, _atom1: usize, _atom2: usize) {
        todo!();
    }

    pub fn set_aromaticity(&self, _model: AromaticityModel) {
        todo!()
    }

    pub fn get_bond_between_atoms_mut(
        &mut self,
        _atom1: usize,
        _atom2: usize,
    ) -> &mut Bond {
        todo!()
    }

    pub fn sanitize(&mut self, _options: SanitizeOptions) {
        todo!();
    }

    pub fn add_conformer(&mut self, _coordinates: &[f64]) {
        todo!();
    }

    pub fn assign_stereochemistry_from_3d(&mut self) {
        todo!();
    }

    pub fn get_atoms(&mut self) -> impl Iterator<Item = &Atom> {
        self.atoms.iter()
    }

    pub fn get_bonds(&self) -> impl Iterator<Item = &Bond> {
        self.bonds.iter()
    }

    pub fn get_bonds_mut(&mut self) -> impl Iterator<Item = &mut Bond> {
        self.bonds.iter_mut()
    }

    pub fn find_symmetry_classes(&self) -> HashMap<usize, String> {
        todo!();
    }
}
