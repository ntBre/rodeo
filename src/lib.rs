use std::collections::HashMap;

use bond::{Bond, BondType};
use graph::Graph;

pub mod bond;
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

    fn calc_explicit_valence(&self, _arg: bool) -> usize {
        todo!()
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
        const AdjustHs = 0x200;
        const CleanupOrganometallics = 0x400;
        const All = 0b11111111111; // TODO not sure about this
    }
}

pub enum AromaticityModel {
    MDL,
    Default,
}

#[derive(Default)]
#[allow(unused)]
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

    /// isn't this just the length of bonds?
    num_bonds: usize,

    d_graph: Graph<Atom, Bond>,

    /// some collection of conformers. looks like this is supposed to be some
    /// kind of dictionary? called `d_conf` in C++
    conformers: Vec<Conformer>,
}

impl RWMol {
    pub fn new() -> Self {
        Self {
            name: String::new(),
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
        atom.set_index(which.0);
        self.d_graph[which] = atom;
        for cfi in self.conformers_mut() {
            cfi.set_atom_pos(which.0, Point3D::new(0.0, 0.0, 0.0));
        }
        which.0
    }

    pub fn add_bond(&mut self, atom1: usize, atom2: usize) {
        assert!(atom1 != atom2, "attempted to add self-bond");
        let which = self.d_graph.add_edge();
        let mut b = Bond::default();
        self.num_bonds += 1;
        b.set_index(self.num_bonds - 1);
        b.set_begin_atom_index(atom1);
        b.set_end_atom_index(atom2);
        self.d_graph[which] = b;
    }

    pub fn set_aromaticity(&self, _model: AromaticityModel) {
        todo!()
    }

    pub fn get_bond_between_atoms_mut(
        &mut self,
        atom1: usize,
        atom2: usize,
    ) -> Option<&mut Bond> {
        self.d_graph.edges.iter_mut().find(|b| {
            (b.begin_atom_index, b.end_atom_index) == (atom1, atom2)
                || (b.end_atom_index, b.begin_atom_index) == (atom1, atom2)
        })
    }

    pub fn sanitize(&mut self, options: SanitizeOptions) {
        use SanitizeOptions as SO;

        // clear out any cached properties
        self.clear_computed_props();

        // clean up things like nitro groups
        if options.contains(SO::Cleanup) {
            self.cleanup();
        }

        // update computed properties on atoms and bonds
        if options.contains(SO::Properties) {
            self.update_property_cache(true);
        } else {
            self.update_property_cache(false);
        }

        if options.contains(SO::SymmRings) {
            self.symmetrize_sssr();
        }

        // kekulizations
        if options.contains(SO::Kekulize) {
            self.kekulize();
        }

        // look for radicals:
        //
        // We do this now because we need to know that the N in [N]1C=CC=C1 has
        // a radical before we move into set_aromaticity. It's important that
        // this happen post-Kekulization because there's no way of telling what
        // to do with the same molecule if it's in the form [n]1cccc1
        if options.contains(SO::FindRadicals) {
            self.assign_radicals();
        }

        // then do aromaticity perception
        if options.contains(SO::SetAromaticity) {
            self.set_aromaticity(AromaticityModel::Default);
        }

        // set conjugation
        if options.contains(SO::SetConjugation) {
            self.set_conjugation()
        }

        // set hybridization
        if options.contains(SO::SetHybridization) {
            self.set_hybridization();
        }

        // remove bogus chirality specs
        if options.contains(SO::CleanupChirality) {
            self.cleanup_chirality();
        }

        // adjust hydrogen counts
        if options.contains(SO::AdjustHs) {
            self.adjust_hs();
        }

        // fix things like non-metal to metal bonds that should be dative.
        if options.contains(SO::CleanupOrganometallics) {
            self.cleanup_organometallics();
        }

        // now that everything has been cleaned up, go through and check/update
        // the computed valences on atoms and bonds one more time
        if options.contains(SO::Properties) {
            self.update_property_cache(true);
        }
    }

    pub fn add_conformer(&mut self, _coordinates: &[f64]) {
        todo!();
    }

    pub fn assign_stereochemistry_from_3d(&mut self) {
        todo!();
    }

    pub fn get_atoms(&mut self) -> impl Iterator<Item = &Atom> {
        self.d_graph.vertices.iter()
    }

    pub fn get_bonds(&self) -> impl Iterator<Item = &Bond> {
        self.d_graph.edges.iter()
    }

    pub fn get_bonds_mut(&mut self) -> impl Iterator<Item = &mut Bond> {
        self.d_graph.edges.iter_mut()
    }

    pub fn find_symmetry_classes(&self) -> HashMap<usize, String> {
        todo!();
    }

    fn conformers_mut(&mut self) -> impl Iterator<Item = &mut Conformer> {
        self.conformers.iter_mut()
    }

    fn clear_computed_props(&self) {
        todo!()
    }

    fn cleanup(&self) {
        todo!()
    }

    fn update_property_cache(&self, _arg: bool) {
        todo!()
    }

    fn symmetrize_sssr(&self) {
        todo!()
    }

    fn kekulize(&self) {
        todo!()
    }

    fn assign_radicals(&self) {
        todo!()
    }

    fn set_conjugation(&self) {
        todo!()
    }

    fn set_hybridization(&self) {
        todo!()
    }

    fn cleanup_chirality(&self) {
        todo!()
    }

    fn adjust_hs(&self) {
        todo!()
    }

    fn cleanup_organometallics(&self) {
        todo!()
    }
}
