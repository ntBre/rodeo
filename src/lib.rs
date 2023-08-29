#![allow(unused)]
#![feature(lazy_cell)]

use std::{
    collections::HashMap,
    fs::read_to_string,
    ops::{Index, IndexMut},
    path::Path,
};

use bond::{Bond, BondType};
use graph::Graph;
use ptable::{DEFAULT_VALENCE, IS_EARLY_ATOM, VALENCE_LIST};

use crate::ptable::{OUTER_ELECS, SYMBOL};

pub mod bond;
mod graph;
mod ptable;

const COMPLEX_QUERIES: [&str; 8] = ["A", "AH", "Q", "QH", "X", "XH", "M", "MH"];

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
    implicit_valence: isize,
    num_explicit_hs: usize,
    num_radical_electrons: usize,
}

impl Atom {
    pub fn new(atomic_number: usize) -> Self {
        Self {
            atomic_number,
            name: String::new(),
            formal_charge: 0,
            is_aromatic: false,
            index: 0,
            // signal for not having been set. might make more sense for this to
            // be optional then
            explicit_valence: -1,
            // guessing the same is true for implicit
            implicit_valence: -1,
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

#[derive(Clone, Debug, Default)]
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

    #[inline]
    const fn zero() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
}

/// Representation of a 2D or 3D conformation of a molecule as a vector of 3D
/// points.
#[derive(Default)]
struct Conformer {
    positions: Vec<Point3D>,
    is_3d: bool,
}

impl Conformer {
    fn new(num_atoms: usize) -> Self {
        Self {
            positions: vec![Point3D::zero(); num_atoms],
            ..Self::default()
        }
    }

    /// panics if `atom_id` is greater than usize::MAX lmao
    fn set_atom_pos(&mut self, atom_id: usize, new: Point3D) {
        if atom_id >= self.positions.len() {
            self.positions.resize_with(atom_id + 1, Point3D::default);
        }
        self.positions[atom_id] = new;
    }

    fn set_3d(&mut self, v: bool) {
        self.is_3d = v;
    }

    fn has_non_zero_z_coords(&self) -> bool {
        todo!()
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

    extra_rings: Vec<Vec<isize>>,

    is_3d: bool,
}

impl RWMol {
    pub fn new() -> Self {
        Self {
            name: String::new(),
            d_graph: Graph::default(),
            ..Self::default()
        }
    }

    /// load an [RWMol] from an SDF file
    pub fn from_sdf(file: impl AsRef<Path>) -> Self {
        let s = read_to_string(dbg!(file.as_ref())).unwrap();
        println!("{s}");
        let mut lines = s.lines();

        let name = lines.next().unwrap();
        let info = lines.next().unwrap();
        let comm = lines.next().unwrap();

        // counts line
        let line = lines.next().unwrap();
        let natoms = line[0..3].trim().parse::<usize>().unwrap();
        let nbonds = line[3..6].trim().parse::<usize>().unwrap();
        let nlists = line[6..9].trim().parse::<usize>().unwrap();
        let chiral_flag = line[12..15].trim().parse::<usize>().unwrap();
        let ns_text = line[15..18].trim().parse::<usize>().unwrap();
        let n_rxn_components = line[18..21].trim().parse::<usize>().unwrap();
        let n_reactants = line[21..24].trim().parse::<usize>().unwrap();
        let n_products = line[24..27].trim().parse::<usize>().unwrap();
        let n_intermediates = line[27..30].trim().parse::<usize>().unwrap();

        let mut ctab_version = 2000;
        if line.len() > 35 {
            if line.len() < 39 || line.chars().nth(34).unwrap() != 'V' {
                panic!("CTAB version string invalid");
            } else if &line[34..39] == "V3000" {
                ctab_version = 3000;
            } else if &line[34..39] != "V2000" {
                panic!("unsupported CTAB version {}", &line[34..39]);
            }
        }

        let mut mol = Self::new();
        if ctab_version == 2000 {
            // in-line ParseV2000TAB
            let mut conf = Conformer::new(natoms);
            if natoms == 0 {
                conf.set_3d(false);
            } else {
                // in-line ParseMolBlockAtoms
                for i in 1..=natoms {
                    let line = lines.next().unwrap();
                    // in-line ParseMolFileAtomLine
                    let (atom, pos) = {
                        if line.len() < 34 {
                            // we're always strict, or we could check < 32
                            // without
                            panic!("atom line too short");
                        }
                        let mut pos = Point3D::zero();
                        pos.x = line[0..10].trim().parse().unwrap();
                        pos.y = line[10..20].trim().parse().unwrap();
                        pos.z = line[20..30].trim().parse().unwrap();
                        let symb = line[31..34].trim();
                        // REVIEW: should we handle missing fields at the end of
                        // the line?
                        let mut mass_diff = 0;
                        if line.len() >= 36 && &line[34..36] != " 0" {
                            mass_diff = line[34..36].trim().parse().unwrap();
                        }
                        let mut chg = 0;
                        if line.len() >= 39 && &line[36..39] != "  0" {
                            chg = line[36..39].trim().parse().unwrap();
                        }
                        let mut hcount = 0;
                        if line.len() >= 45 && &line[42..45] != "  0" {
                            hcount = line[42..45].trim().parse().unwrap();
                        }
                        let mut res = Atom::default();
                        let is_complex_query_name = COMPLEX_QUERIES
                            .iter()
                            .find(|&&s| s == symb)
                            .is_some();
                        if is_complex_query_name
                            || symb == "L"
                            || symb == "LP"
                            || symb == "R"
                            || symb == "R#"
                            || (symb.chars().nth(0).unwrap() == 'R'
                                && symb >= "R0"
                                && symb <= "R99")
                        {
                            if is_complex_query_name
                                || symb == "*"
                                || symb == "R"
                            {
                                let mut query = Atom::new(0);
                                if symb == "*" || symb == "R" {
                                    // according to the MDL spec, these match
                                    // anything
                                    query.set_query(make_atom_null_query());
                                } else if is_complex_query_name {
                                    query.convert_complex_name(symb);
                                }
                                res = query;
                                res.set_no_implicit(true);
                            } else {
                                res.set_atomic_num(0);
                            }
                            if mass_diff == 0
                                && symb.chars().nth(0).unwrap() == 'R'
                            {
                                if symb.len() > 1
                                    && symb >= "R0"
                                    && symb <= "R99"
                                {
                                    let rnumber = &symb[1..symb.len() - 1]
                                        .parse()
                                        .unwrap_or(-1);
                                    if *rnumber >= 0 {
                                        res.set_isotope(*rnumber);
                                    }
                                }
                            }
                            if symb.chars().nth(0).unwrap() == 'R' {
                                todo!("setRGPProps(symb, res)");
                            }
                        } else if symb == "D" {
                            res.set_atomic_num(1);
                            res.set_isotope(2);
                        } else if symb == "T" {
                            res.set_atomic_num(1);
                            res.set_isotope(3);
                        } else if symb == "Pol" || symb == "Mod" {
                            res.set_atomic_num(0);
                            res.dummy_label = symb.to_owned();
                        } else {
                            let s: Vec<_> = symb.chars().collect();
                            let symb = if symb.len() == 2
                                && s[1] >= 'A'
                                && s[1] <= 'Z'
                            {
                                format!("{}{}", s[0], s[1].to_ascii_lowercase())
                            } else {
                                symb.to_string()
                            };
                            lookup_atomic_number(res, symb);
                        }
                        if chg != 0 {
                            res.set_formal_charge(chg);
                        }
                        if (hcount >= 1) {
                            // NOTE: skipping res->hasQuery() check here, 1531
                            // MolFileParser.cpp
                            res.set_no_implicit(true);
                            if (hcount > 1) {
                                todo!()
                            //                       ATOM_EQUALS_QUERY *oq = makeAtomImplicitHCountQuery(hcount - 1);
                            //                       auto nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>(
                            // hcount - 1, oq->getDataFunc(),
                            // std::string("less_") + oq->getDescription());
                            //                       res->expandQuery(nq);
                            //                       delete oq;
                            } else {
                                // res->expandQuery(makeAtomImplicitHCountQuery(0));
                                todo!();
                            }
                        }

                        if (mass_diff != 0) {
                            todo!()
                            // int defIso =
                            //     PeriodicTable::getTable()->getMostCommonIsotope(res->getAtomicNum());
                            // int dIso = defIso + mass_diff;
                            // if (dIso < 0) {
                            //   BOOST_LOG(rdWarningLog)
                            //       << " atom " << res->getIdx()
                            //       << " has a negative isotope offset. line:  " << line << std::endl;
                            // }
                            // res->setIsotope(dIso);
                            // res->setProp(common_properties::_hasMassQuery, true);
                        }

                        if (line.len() >= 42 && &line[39..42] != "  0") {
                            let mut parity = 0;
                            parity = line[39..42].trim().parse().unwrap();
                            res.set_parity(parity);
                        }

                        if (line.len() >= 48 && &line[45..45 + 3] != "  0") {
                            let mut stereoCare = 0;
                            stereoCare =
                                line[45..45 + 3].trim().parse().unwrap();
                            res.set_stereo_care(stereo_care);
                        }
                        if (line.len() >= 51 && &line[48..48 + 3] != "  0") {
                            let mut totValence = 0;
                            totValence =
                                line[48..48 + 3].trim().parse().unwrap();
                            if (totValence != 0) {
                                // only set if it's a non-default value
                                res.set_tot_valence(totValence);
                            }
                        }
                        if (line.len() >= 57 && &line[54..54 + 3] != "  0") {
                            let mut rxnRole = 0;
                            rxnRole = line[54..54 + 3].trim().parse().unwrap();
                            if (rxnRole != 0) {
                                // only set if it's a non-default value
                                res.set_molRxnRole(rxnRole);
                            }
                        }
                        if (line.len() >= 60 && &line[57..57 + 3] != "  0") {
                            let mut rxnComponent = 0;
                            rxnComponent =
                                line[57..57 + 3].trim().parse().unwrap();
                            if (rxnComponent != 0) {
                                // only set if it's a non-default value
                                res.set_molRxnComponent(rxnComponent);
                            }
                        }
                        if (line.len() >= 63 && &line[60..60 + 3] != "  0") {
                            let mut atomMapNumber = 0;
                            atomMapNumber =
                                line[60..60 + 3].trim().parse().unwrap();
                            res.set_molAtomMapNumber(atomMapNumber);
                        }
                        if (line.len() >= 66 && &line[63..63 + 3] != "  0") {
                            let mut inversionFlag = 0;
                            inversionFlag =
                                line[63..63 + 3].trim().parse().unwrap();
                            res.set_molInversionFlag(inversionFlag);
                        }
                        if (line.len() >= 69 && &line[66..66 + 3] != "  0") {
                            let mut exactChangeFlag = 0;
                            exactChangeFlag =
                                line[66..66 + 3].trim().parse().unwrap();
                            res.set_molExactChangeFlag(exactChangeFlag);
                        }
                        return res;

                        todo!();
                    };
                    let aid = mol.add_atom(atom);
                    conf.set_atom_pos(aid, pos);
                    mol.set_atom_bookmark(atom, i);
                }

                let nonzero_z = conf.has_non_zero_z_coords();
                if mol.is_3d {
                    conf.set_3d(true);
                    if !nonzero_z {
                        eprintln!(
                            "warning: molecule is tagged as 3d,\
                                but all Z coords are zero"
                        );
                    }
                } else {
                    conf.set_3d(nonzero_z);
                }
            }
            mol.add_conformer2(conf, true);
            todo!("ParseMolBlockBonds");
            todo!("ParseMolBlockProperties");
        } else {
            todo!("unhandled ctab version: {}", ctab_version);
        }
        println!("{}", s);
        todo!("finishMolProcessing");
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

    fn add_conformer2(&mut self, conf: Conformer, arg: bool) {
        todo!()
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

    fn symmetrize_sssr(&mut self) -> usize {
        let mut res: Vec<Vec<isize>> = Vec::new();
        let mut sssrs: Vec<Vec<isize>> = Vec::new();

        // FIX: need to set flag here the symmetrization has been done in order
        // to avoid repeating this work
        let include_dative_bonds = false; // default and not passed
        self.find_sssr(&mut sssrs, include_dative_bonds);

        res.reserve(sssrs.len());
        for r in &sssrs {
            res.push(r.clone());
        }

        // now check if there are any extra rings on the molecule
        if self.extra_rings.is_empty() {
            // no extra rings nothing to be done
            return res.len();
        }
        let extras = self.extra_rings.clone();

        // convert the rings to bond ids
        let bondsssrs = Vec::new();
        self.convert_to_bonds(&sssrs, &bondsssrs);

        //
        // For each "extra" ring, figure out if it could replace a single ring
        // in the SSSR. A ring could be swapped out if:
        //
        // * They are the same size
        // * The replacement doesn't remove any bonds from the union of the
        //   bonds in the SSSR.
        //
        // The latter can be checked by determining if the SSSR ring is the
        // unique provider of any ring bond. If it is, the replacement ring must
        // also provide that bond.
        //
        // May miss extra rings that would need to swap two (or three...) rings
        // to be included.

        // counts of each bond
        let mut bond_counts = vec![0; self.num_bonds];
        for r in &bondsssrs {
            for b in r {
                bond_counts[*b as usize] += 1;
            }
        }

        let extra_ring: Vec<isize> = Vec::new();
        for extra_atom_ring in &extras {
            self.convert_to_bonds2(extra_atom_ring, &extra_ring);
            for ring in &bondsssrs {
                if ring.len() != extra_ring.len() {
                    continue;
                }

                // If `ring` is the only provider of some bond, extraRing must
                // also provide that bond.
                let mut share_bond = false;
                let mut replaces_all_unique_bonds = true;
                for bond_id in ring {
                    let bond_count = bond_counts[*bond_id as usize];
                    if bond_count == 1 || !share_bond {
                        let position =
                            extra_ring.iter().find(|&b| b == bond_id);
                        // TODO this is going to be some kind o if let on
                        // whether the find passed
                        if position.is_some() {
                            share_bond = true;
                        } else if bond_count == 1 {
                            // 1 means `ring` is the only ring in the SSSR to
                            // provide this bond, and extraRing did not provide
                            // it (so extraRing is not an acceptable
                            // substitution in the SSSR for ring)
                            replaces_all_unique_bonds = false;
                        }
                    }
                }

                if share_bond && replaces_all_unique_bonds {
                    res.push(extra_atom_ring.clone());
                    self.store_ring_info(extra_atom_ring.to_vec());
                    break;
                }
            }
        }

        // we keep the property as just the vector
        if !self.extra_rings.is_empty() {
            self.extra_rings.clear();
        }

        res.len()
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
    const fn atoms(&self) -> &Vec<Atom> {
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
            if (acc - pval) <= 1.5 {
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
            let effective_valence = if OUTER_ELECS[atomic_number] >= 4 {
                res - formal_charge
            } else {
                // for boron and co, we move to the right in PT, so adding extra
                // valence means adding negative charge
                res + formal_charge
            };
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
        atom_idx: usize,
        strict: bool,
    ) -> isize {
        if self[atom_idx].explicit_valence == -1 {
            self.calc_explicit_valence(atom_idx, strict);
        }

        // special cases
        if self[atom_idx].atomic_number == 0 {
            self[atom_idx].implicit_valence = 0;
            return 0;
        }

        // sad but true
        let bonds: Vec<&Bond> = self.get_atom_bonds(atom_idx).collect();
        for bnd in bonds {
            if bnd.has_complex_bond_type() {
                self[atom_idx].implicit_valence = 0;
                return 0;
            }
        }
        if self[atom_idx].explicit_valence == 0
            && self[atom_idx].atomic_number == 1
            && self[atom_idx].num_radical_electrons == 0
        {
            if self[atom_idx].formal_charge == 1
                || self[atom_idx].formal_charge == -1
            {
                self[atom_idx].implicit_valence = 0;
                return 0;
            } else if self[atom_idx].formal_charge == 0 {
                self[atom_idx].implicit_valence = 1;
                return 1;
            } else if strict {
                panic!("Unreasonable formal charge on hydrogen # {atom_idx}");
            } else {
                self[atom_idx].implicit_valence = 0;
                return 0;
            }
        }

        // this is basically the difference between the allowed valence of the
        // atom and the explicit valence already specified - tells how many H's
        // to add
        let mut res;

        // The d-block and f-block of the periodic table (i.e. transition metals,
        // lanthanoids and actinoids) have no default valence.
        let atomic_num = self[atom_idx].atomic_number;
        let dv = DEFAULT_VALENCE[atomic_num];
        if dv == -1 {
            self[atom_idx].implicit_valence = 0;
            return 0;
        }

        // here is how we are going to deal with the possibility of
        // multiple valences
        // - check the explicit valence "ev"
        // - if it is already equal to one of the allowed valences for the
        //    atom return 0
        // - otherwise take return difference between next larger allowed
        //   valence and "ev"
        // if "ev" is greater than all allowed valences for the atom raise an
        // exception finally aromatic cases are dealt with differently - these
        // atoms are allowed only default valences
        let valens = &VALENCE_LIST[atomic_num];

        let atom = &self[atom_idx];
        let explicit_plus_rad_v =
            atom.explicit_valence + atom.num_radical_electrons as isize;
        let mut chg = atom.formal_charge;

        // NOTE: this is here to take care of the difference in element on
        // the right side of the carbon vs left side of carbon
        // For elements on the right side of the periodic table
        // (electronegative elements):
        //     NHYD = V - SBO + CHG
        // For elements on the left side of the periodic table
        // (electropositive elements):
        //      NHYD = V - SBO - CHG
        // This reflects that hydrogen adds to, for example, O as H+ while
        // it adds to Na as H-.

        // V = valence
        // SBO = Sum of bond orders
        // CHG = Formal charge

        //  It seems reasonable that the line is drawn at Carbon (in Group
        //  IV), but we must assume on which side of the line C
        //  falls... an assumption which will not always be correct.  For
        //  example:
        //  - Electropositive Carbon: a C with three singly-bonded neighbors (DV
        //    = 4, SBO = 3, CHG = 1) and a positive charge (a 'stable'
        //    carbocation) should not have any hydrogens added.
        //  - Electronegative Carbon: C in isonitrile, R[N+]#[C-] (DV = 4, SBO =
        //    3, CHG = -1), also should not have any hydrogens added.
        //  Because isonitrile seems more relevant to pharma problems, we'll be
        //  making the second assumption: *Carbon is electronegative*.
        //
        // So assuming you read all the above stuff - you know why we are
        // changing signs for "chg" here
        if IS_EARLY_ATOM[atomic_num] {
            chg *= -1;
        }
        // special case for carbon - see GitHub #539
        if atomic_num == 6 && chg > 0 {
            chg = -chg;
        }

        // if we have an aromatic case treat it differently
        if atom.is_aromatic {
            if explicit_plus_rad_v <= dv + chg {
                res = dv + chg - explicit_plus_rad_v;
            } else {
                // As we assume when finding the explicitPlusRadValence if we
                // are aromatic we should not be adding any hydrogen and already
                // be at an accepted valence state,

                // FIX: this is just ERROR checking and probably moot - the
                // explicitPlusRadValence function called above should assure us
                // that we satisfy one of the accepted valence states for the
                // atom. The only diff I can think of is in the way we handle
                // formal charge here vs the explicit valence function.
                let mut satis = false;
                for vi in valens {
                    if *vi < 0 {
                        break;
                    }
                    if explicit_plus_rad_v == ((*vi) + chg) {
                        satis = true;
                        break;
                    }
                }
                if strict && !satis {
                    panic!(
                        "Explicit valence for aromatic atom # {atom_idx} \
                        not equal to any accepted valence"
                    );
                }
                res = 0;
            }
        } else {
            // non-aromatic case we are allowed to have non default valences
            // and be able to add hydrogens
            res = -1;
            for vi in valens {
                if *vi < 0 {
                    break;
                }
                let tot = (*vi) + chg;
                if explicit_plus_rad_v <= tot {
                    res = tot - explicit_plus_rad_v;
                    break;
                }
            }
            if res < 0 {
                if strict && *valens.last().unwrap() != -1 {
                    // this means that the explicit valence is greater than any
                    // allowed valence for the atoms - raise an error
                    panic!(
                        "Explicit valence for atom # {atom_idx} {} \
                        greater than permitted",
                        SYMBOL[atomic_num]
                    );
                } else {
                    res = 0;
                }
            }
        }

        self[atom_idx].implicit_valence = res;
        res
    }

    fn get_num_explicit_hs(&self, atom_idx: usize) -> f64 {
        self[atom_idx].num_explicit_hs as f64
    }

    /// assuming that this gives the bonds involving `atom_idx`
    fn get_atom_bonds(&self, atom_idx: usize) -> impl Iterator<Item = &Bond> {
        self.get_bonds().filter(move |bond| {
            bond.begin_atom_index == atom_idx || bond.end_atom_index == atom_idx
        })
    }

    fn find_sssr(
        &self,
        _sssrs: &mut [Vec<isize>],
        _include_dative_bonds: bool,
    ) {
        todo!()
    }

    fn convert_to_bonds(
        &self,
        _sssrs: &[Vec<isize>],
        _bondsssrs: &[Vec<isize>],
    ) {
        todo!()
    }

    fn store_ring_info(&self, _extra_atom_ring: Vec<isize>) {
        todo!()
    }

    fn convert_to_bonds2(
        &self,
        _extra_atom_ring: &[isize],
        _extra_ring: &[isize],
    ) {
        todo!()
    }

    fn set_atom_bookmark(&mut self, atom: Atom, i: usize) {
        todo!()
    }
}

impl Index<usize> for RWMol {
    type Output = Atom;

    fn index(&self, index: usize) -> &Self::Output {
        &self.atoms()[index]
    }
}

impl IndexMut<usize> for RWMol {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.atoms_mut()[index]
    }
}
