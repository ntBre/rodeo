#![feature(lazy_cell)]

use std::collections::HashMap;

use bond::{Bond, BondType};
use graph::Graph;
use ptable::{DEFAULT_VALENCE, IS_EARLY_ATOM, VALENCE_LIST};

use crate::ptable::{OUTER_ELECS, SYMBOL};

pub mod bond;
mod graph;
mod ptable;

#[derive(Default)]
pub enum Chi {
    TetrahedralCCW,
    TetrahedralCW,
    #[default]
    None,
}

#[derive(Default)]
pub struct Atom {
    atomic_number: usize,
    name: String,
    formal_charge: isize,
    is_aromatic: bool,
    index: usize,
    chiral_tag: Chi,
    explicit_valence: isize,
    num_explicit_hs: usize,
}

impl Atom {
    pub fn new(atomic_number: usize) -> Self {
        Self {
            atomic_number,
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

    fn get_bond_between_atoms(
        &self,
        atom1: usize,
        atom2: usize,
    ) -> Option<&Bond> {
        self.d_graph.edges.iter().find(|b| {
            (b.begin_atom_index, b.end_atom_index) == (atom1, atom2)
                || (b.end_atom_index, b.begin_atom_index) == (atom1, atom2)
        })
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

    pub fn get_atoms(&self) -> impl Iterator<Item = &Atom> {
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

    fn clear_computed_props(&mut self) {
        // no op for now, not sure which props are computed
    }

    fn cleanup(&mut self) {
        for atom_idx in 0..self.d_graph.vertices.len() {
            match self.d_graph.vertices[atom_idx].atomic_number {
                7 => self.nitrogen_cleanup(atom_idx),
                15 => self.phosphorus_cleanup(atom_idx),
                17 | 35 | 53 => self.halogen_cleanup(atom_idx),
                _ => {}
            }
        }
    }

    fn update_property_cache(&mut self, strict: bool) {
        for atom_idx in 0..self.d_graph.vertices.len() {
            // in-line atom.updatePropertyCache
            self.calc_explicit_valence(atom_idx, strict);
            self.calc_implicit_valence(atom_idx, strict);
        }
        // in-line bond.updatePropertyCache, which does nothing
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

    #[inline]
    fn atoms(&self) -> &Vec<Atom> {
        &self.d_graph.vertices
    }

    #[inline]
    fn atoms_mut(&mut self) -> &mut Vec<Atom> {
        &mut self.d_graph.vertices
    }

    /// conversions here:
    /// - neutral 5 coordinate N's with double bonds to O's to the zwitterionic
    /// form. eg: CN(=O)=O -> C[N+](=O)[O-]
    /// and:
    /// C1=CC=CN(=O)=C1 -> C1=CC=C[N+]([O-])=C1
    /// - neutral 5 coordinate N's with triple bonds to N's to the zwitterionic
    /// form. eg: C-N=N#N -> C-N=[N+]=[N-]
    fn nitrogen_cleanup(&mut self, atom_idx: usize) {
        // NOTE we have to take this out to get around calling other borrowing
        // methods in the loop. DO NOT forget to put it back at the end of the
        // method
        let mut atom = std::mem::take(&mut self.atoms_mut()[atom_idx]);

        // we only want to do neutrals so that things like this don't get
        // munged:
        //  O=[n+]1occcc1
        // this was sf.net issue 1811276
        if atom.formal_charge != 0 {
            return;
        }

        // we need to play this little aromaticity game because the explicit
        // valence code modifies its results for aromatic atoms
        let arom_holder = atom.is_aromatic;
        atom.is_aromatic = false;

        // NOTE that we are calling calcExplicitValence() here, we do this
        // because we cannot be sure that it has already been called on the atom
        // (cleanUp() gets called pretty early in the sanitization process):
        let valence = self.calc_explicit_valence(atom_idx, false);
        if valence == 5 {
            let aid = atom.get_index();
            for nbr_idx in self.atom_neighbors(atom_idx) {
                let mut is_break = false;
                // NOTE as with atom, we are taking this out of self.atoms to
                // mess with it in the loop. Don't forget to put it back at the
                // bottom of the loop. This is why we also can't use break
                // directly and have to set is_break instead
                let mut nbr = std::mem::take(&mut self.atoms_mut()[nbr_idx]);
                if nbr.atomic_number == 8
                    && nbr.formal_charge == 0
                    && (self
                        .get_bond_between_atoms(aid, nbr.index)
                        .unwrap()
                        .bond_type
                        .is_double())
                {
                    // here's the double-bonded oxygen
                    let b = self
                        .get_bond_between_atoms_mut(aid, nbr.index)
                        .unwrap();
                    b.set_bond_type(BondType::Single);
                    atom.set_formal_charge(1);
                    nbr.set_formal_charge(-1);
                    is_break = true;
                } else if nbr.atomic_number == 7
                    && nbr.formal_charge == 0
                    && self
                        .get_bond_between_atoms(aid, nbr.index)
                        .unwrap()
                        .bond_type
                        .is_triple()
                {
                    // here's the triple-bonded nitrogens
                    let b = self
                        .get_bond_between_atoms_mut(aid, nbr.index)
                        .unwrap();
                    b.set_bond_type(BondType::Double);
                    atom.set_formal_charge(1);
                    nbr.set_formal_charge(-1);
                    is_break = true;
                }
                std::mem::swap(&mut self.atoms_mut()[nbr_idx], &mut nbr);
                if is_break {
                    break;
                }
            }
        }
        // force a recalculation of the explicit valence here
        atom.set_is_aromatic(arom_holder);
        std::mem::swap(&mut self.atoms_mut()[atom_idx], &mut atom);
        self.calc_explicit_valence(atom_idx, false);
    }

    fn phosphorus_cleanup(&mut self, _atom_idx: usize) {
        todo!()
    }

    fn halogen_cleanup(&mut self, _atom_idx: usize) {
        todo!()
    }

    /// return a vector of atom indices that are neighbors to `atom_idx`
    fn atom_neighbors(&self, _atom_idx: usize) -> Vec<usize> {
        todo!();
    }

    /// in rdkit this is a method on `Atom`, but it requires an owning Molecule,
    /// which is a reference nightmare in Rust.
    fn calc_explicit_valence(
        &mut self,
        atom_idx: usize,
        strict: bool,
    ) -> isize {
        let mut acc = 0.0;
        for bond in self.get_bonds() {
            acc += bond.get_valence_contrib(atom_idx);
        }
        acc += self.get_num_explicit_hs(atom_idx);
        // check if acc is greater than the default valence
        let atomic_number = self.atoms()[atom_idx].atomic_number;
        let dv = DEFAULT_VALENCE[atomic_number] as f64;

        let formal_charge = self.atoms()[atom_idx].formal_charge;
        let mut chr = formal_charge as f64;
        if IS_EARLY_ATOM[atomic_number] {
            chr *= -1.0; // the usual correction for early atoms
        }

        // special case for carbon - see github #539
        if atomic_number == 6 && chr > 0.0 {
            chr = -chr; // no multiply this time? lol
        }

        if acc > (dv + chr) && self.atoms()[atom_idx].is_aromatic {
            // this needs some explanation : if the atom is aromatic and accum >
            // (dv + chr) we assume that no hydrogen can be added to this atom.
            // We set x = (v + chr) such that x is the closest possible integer
            // to "accum" but less than "accum".
            //
            // "v" here is one of the allowed valences. For example:
            //    sulfur here : O=c1ccs(=O)cc1
            //    nitrogen here : c1cccn1C
            let mut pval = dv + chr;
            let valens = &VALENCE_LIST[atomic_number];
            for val in valens {
                if *val == -1 {
                    break;
                }
                let mut val = *val as f64;
                val += chr;
                if val > acc {
                    break;
                } else {
                    pval = val;
                }
            }
            // if we're within 1.5 of the allowed valence, go ahead and take it.
            // this reflects things like the N in c1cccn1C, which starts with
            // accum of 4, but which can be kekulized to C1=CC=CN1C, where the
            // valence is 3 or the bridging N in c1ccn2cncc2c1, which starts
            // with a valence of 4.5, but can be happily kekulized down to a
            // valence of 3
            if (acc - pval) as f64 <= 1.5 {
                acc = pval;
            }
        }
        // despite promising to not to blame it on him - this a trick Greg came
        // up with: if we have a bond order sum of x.5 (i.e. 1.5, 2.5 etc) we
        // would like it to round to the higher integer value -- 2.5 to 3
        // instead of 2 -- so we will add 0.1 to accum. this plays a role in the
        // number of hydrogen that are implicitly added. This will only happen
        // when the accum is a non-integer value and less than the default
        // valence (otherwise the above if statement should have caught it). An
        // example of where this can happen is the following smiles:
        //
        //    C1ccccC1
        //
        // Daylight accepts this smiles and we should be able to Kekulize
        // correctly.
        acc += 0.1;

        let res = acc.round() as isize;

        if strict {
            let effective_valence;
            if OUTER_ELECS[atomic_number] >= 4 {
                effective_valence = res - formal_charge;
            } else {
                // for boron and co, we move to the right in PT, so adding extra
                // valence means adding negative charge
                effective_valence = res + formal_charge;
            }
            let valens = &VALENCE_LIST[atomic_number];

            let max_valence = valens.last().unwrap();
            // max_valence == -1 signifies that we'll take anything at the high
            // end
            if *max_valence > 0 && effective_valence > *max_valence {
                panic!(
                    "Explicit valence for atom # {} {atomic_number} {}, \
                    {effective_valence}, is greater than permitted",
                    atom_idx, SYMBOL[atomic_number]
                );
            }
        }

        self.atoms_mut()[atom_idx].explicit_valence = res;
        res
    }

    fn calc_implicit_valence(
        &mut self,
        _atom_idx: usize,
        _strict: bool,
    ) -> isize {
        todo!()
    }

    fn get_num_explicit_hs(&self, atom_idx: usize) -> f64 {
        self.atoms()[atom_idx].num_explicit_hs as f64
    }
}
