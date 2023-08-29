//! Utilities for parsing SDF files

use crate::{ptable::SYMBOL, Atom, Conformer, Point3D, RWMol, COMPLEX_QUERIES};

pub(crate) fn parse_mol_block_atoms(
    natoms: usize,
    mut lines: std::str::Lines<'_>,
    mol: &mut RWMol,
    conf: &mut Conformer,
) {
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
            let is_complex_query_name =
                COMPLEX_QUERIES.iter().find(|&&s| s == symb).is_some();
            if is_complex_query_name
                || symb == "L"
                || symb == "LP"
                || symb == "R"
                || symb == "R#"
                || (symb.chars().nth(0).unwrap() == 'R'
                    && symb >= "R0"
                    && symb <= "R99")
            {
                if is_complex_query_name || symb == "*" || symb == "R" {
                    let mut query = Atom::new(0);
                    if symb == "*" || symb == "R" {
                        // according to the MDL spec, these match
                        // anything
                        todo!("query.set_query(make_atom_null_query())");
                    } else if is_complex_query_name {
                        todo!("query.convert_complex_name(symb);");
                    }
                    res = query;
                    res.set_no_implicit(true);
                } else {
                    res.set_atomic_number(0);
                }
                if mass_diff == 0 && symb.chars().nth(0).unwrap() == 'R' {
                    if symb.len() > 1 && symb >= "R0" && symb <= "R99" {
                        let rnumber =
                            &symb[1..symb.len() - 1].parse().unwrap_or(-1);
                        if *rnumber >= 0 {
                            res.set_isotope(*rnumber);
                        }
                    }
                }
                if symb.chars().nth(0).unwrap() == 'R' {
                    todo!("setRGPProps(symb, res)");
                }
            } else if symb == "D" {
                res.set_atomic_number(1);
                res.set_isotope(2);
            } else if symb == "T" {
                res.set_atomic_number(1);
                res.set_isotope(3);
            } else if symb == "Pol" || symb == "Mod" {
                res.set_atomic_number(0);
                res.dummy_label = symb.to_owned();
            } else {
                let s: Vec<_> = symb.chars().collect();
                let symb = if symb.len() == 2 && s[1] >= 'A' && s[1] <= 'Z' {
                    format!("{}{}", s[0], s[1].to_ascii_lowercase())
                } else {
                    symb.to_string()
                };
                let num = SYMBOL.iter().position(|&s| s == symb).unwrap();
                res.set_atomic_number(num);
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
                let mut stereo_care = 0;
                stereo_care = line[45..45 + 3].trim().parse().unwrap();
                res.set_stereo_care(stereo_care);
            }
            if (line.len() >= 51 && &line[48..48 + 3] != "  0") {
                let mut tot_valence = 0;
                tot_valence = line[48..48 + 3].trim().parse().unwrap();
                if (tot_valence != 0) {
                    // only set if it's a non-default value
                    res.set_tot_valence(tot_valence);
                }
            }
            if (line.len() >= 57 && &line[54..54 + 3] != "  0") {
                let mut rxn_role = 0;
                rxn_role = line[54..54 + 3].trim().parse().unwrap();
                if (rxn_role != 0) {
                    // only set if it's a non-default value
                    res.set_rxn_role(rxn_role);
                }
            }
            if (line.len() >= 60 && &line[57..57 + 3] != "  0") {
                let mut rxn_component = 0;
                rxn_component = line[57..57 + 3].trim().parse().unwrap();
                if (rxn_component != 0) {
                    // only set if it's a non-default value
                    res.set_rxn_component(rxn_component);
                }
            }
            if (line.len() >= 63 && &line[60..60 + 3] != "  0") {
                let mut atom_map_number = 0;
                atom_map_number = line[60..60 + 3].trim().parse().unwrap();
                res.set_map_number(atom_map_number);
            }
            if (line.len() >= 66 && &line[63..63 + 3] != "  0") {
                let mut inversion_flag = 0;
                inversion_flag = line[63..63 + 3].trim().parse().unwrap();
                res.set_inversion_flag(inversion_flag);
            }
            if (line.len() >= 69 && &line[66..66 + 3] != "  0") {
                let mut exact_change_flag = 0;
                exact_change_flag = line[66..66 + 3].trim().parse().unwrap();
                res.set_exact_change_flag(exact_change_flag);
            }
            (res, pos)
        };
        let aid = mol.add_atom(atom.clone());
        conf.set_atom_pos(aid, pos);
        mol.set_atom_bookmark(aid, i);
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
