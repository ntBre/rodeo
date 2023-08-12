use std::collections::HashMap;

use graph::Graph;

mod graph;

#[derive(Default)]
pub enum Chi {
    TetrahedralCCW,
    TetrahedralCW,
    #[default]
    None,
}

#[derive(Default)]
pub struct Atom {
    name: String,
    formal_charge: isize,
    is_aromatic: bool,
    index: usize,
    chiral_tag: Chi,
}

impl Atom {
    pub fn new(_num: usize) -> Self {
        Self {
            name: String::new(),
            formal_charge: 0,
            is_aromatic: false,
            index: 0,
            ..Default::default()
        }
    }

    pub fn get_index(&self) -> usize {
        todo!();
    }

    pub fn get_stereo(&self) -> String {
        todo!();
    }

    pub fn set_formal_charge(&mut self, charge: isize) {
        self.formal_charge = charge;
    }

    pub fn set_is_aromatic(&mut self, aromatic: bool) {
        self.is_aromatic = aromatic;
    }

    pub fn set_name(&mut self, name: String) {
        self.name = name;
    }

    pub fn set_index(&mut self, index: usize) {
        self.index = index;
    }

    pub fn set_chiral_tag(&mut self, chiral_tag: Chi) {
        self.chiral_tag = chiral_tag;
    }
}

#[derive(Debug, PartialEq)]
pub enum BondStereo {
    E,
    Z,
    None,
}

#[derive(Default)]
pub enum BondType {
    Single,
    Aromatic,
    Double,
    Triple,
    Quadruple,
    Quintuple,
    Hextuple,
    OneAndAHalf,
    #[default]
    Unspecified,
}

#[derive(Default)]
pub struct Bond {
    bond_type: BondType,
}

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

#[derive(Default)]
struct Point3D {
    x: f64,
    y: f64,
    z: f64,
}

impl Point3D {
    fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

/// Representation of a 2D or 3D conformation of a molecule as a vector of 3D
/// points.
struct Conformer {
    positions: Vec<Point3D>,
}

impl Conformer {
    /// panics if `atom_id` is greater than usize::MAX lmao
    fn set_atom_pos(&mut self, atom_id: usize, new: Point3D) {
        if atom_id >= self.positions.len() {
            self.positions.resize_with(atom_id + 1, Point3D::default);
        }
        self.positions[atom_id] = new;
    }
}

#[derive(Default)]
pub struct RWMol {
    name: String,
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    d_graph: Graph<Atom>,

    /// some collection of conformers. looks like this is supposed to be some
    /// kind of dictionary? called `d_conf` in C++
    conformers: Vec<Conformer>,
}

impl RWMol {
    pub fn new() -> Self {
        Self {
            name: String::new(),
            atoms: Vec::new(),
            bonds: Vec::new(),
            d_graph: Graph::default(),
            ..Self::default()
        }
    }

    /// corresponds to setProp("_Name", name)
    pub fn set_name(&mut self, name: String) {
        self.name = name;
    }

    pub fn add_atom(&mut self, mut atom: Atom) -> usize {
        let which = self.d_graph.add_vertex();
        atom.set_index(which);
        self.d_graph[which] = atom;
        for cfi in self.conformers_mut() {
            cfi.set_atom_pos(which, Point3D::new(0.0, 0.0, 0.0));
        }
        which
    }

    pub fn add_bond(&mut self, atom1: usize, atom2: usize) {
        assert!(atom1 != atom2, "attempted to add self-bond");
        let b = Bond::default();
        let which = self.d_graph.add_edge(atom1, atom2);
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

    fn conformers_mut(&mut self) -> impl Iterator<Item = &mut Conformer> {
        self.conformers.iter_mut()
    }
}
