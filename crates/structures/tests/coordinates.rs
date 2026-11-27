//! Authored coordinate fixtures and geometric expectations, not golden output
//! copied from another viewer or parser.
use collagene_structures::{
    scene::{self, Alternates, Chains, Selection},
    *,
};
const A: &str = "ATOM      1  CA AALA A  -2B     11.104  13.207   2.100  0.60 20.00           C  ";
const B: &str = "HETATM    2 CA    CA Z  17      12.000  13.000   2.000  1.00 10.00          CA2+";
const C: &str = "ATOM      3  CA BALA A  -2B     11.204  13.307   2.200  0.40 21.00           C  ";
fn selection(model: u32) -> Selection {
    Selection {
        model,
        chains: Chains::All,
        alternates: Alternates::All,
        hydrogens: true,
    }
}
fn pdb(lines: &[&str]) -> Structure {
    import(lines.join("\n").as_bytes(), Format::Pdb).unwrap()
}
const CIF_HEADER: &str = "data_authored\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\n";
const ROW_A: &str = "ATOM 1 C CA A ALA X 1 5 B 1.0 2.0 3.0 0.6 12.0 -2 ALA Author CA 1\n";
const ROW_B: &str = "HETATM 2 CA CA . CA Y 2 . ? -1.0 -2.0 -3.0 1.0 10.0 17 CA Author CA 1\n";
fn cif(rows: &str) -> Structure {
    import(format!("{CIF_HEADER}{rows}").as_bytes(), Format::Mmcif).unwrap()
}

#[test]
fn pdb_columns_preserve_models_alternates_author_identity_and_declared_elements() {
    let source = format!(
        "HEADER    AUTHORED FIXTURE\r\nMODEL        1\r\n{A}\r\n{B}\r\n{C}\r\nENDMDL\r\nMODEL        2\r\n{A}\r\nENDMDL\r\nEND\r\n"
    );
    let s = import(source.as_bytes(), Format::Pdb).unwrap();
    assert_eq!(s.original(), source.as_bytes());
    assert_eq!(
        s.source_hash(),
        blake3::hash(source.as_bytes()).to_hex().as_str()
    );
    s.validate().unwrap();
    assert_eq!(s.models().into_iter().collect::<Vec<_>>(), vec![1, 2]);
    assert_eq!(s.atoms().len(), 4);
    let a = &s.atoms()[0];
    assert_eq!(a.position, [11.104, 13.207, 2.1]);
    assert_eq!(a.occupancy, Some(0.6));
    assert_eq!(a.b_factor, Some(20.0));
    assert_eq!(a.author.sequence.as_deref(), Some("-2"));
    assert_eq!(a.insertion.as_deref(), Some("B"));
    assert_eq!(a.alternate.as_deref(), Some("A"));
    assert_eq!(a.element.as_deref(), Some("C"));
    assert_eq!(a.label.chain, None);
    assert_eq!(s.atoms()[1].element.as_deref(), Some("CA"));
    assert_eq!(s.atoms()[1].formal_charge, Some(2));
    assert_eq!(s.notes().ignored_pdb_records.get("HEADER"), Some(&1));
    let scene = scene::build(
        &s,
        Selection {
            alternates: Alternates::Named("A".into()),
            ..selection(1)
        },
    )
    .unwrap();
    assert_eq!(
        scene
            .points
            .iter()
            .map(|p| p.atom_index)
            .collect::<Vec<_>>(),
        vec![0, 1]
    );
    assert!(scene::build(&s, selection(3)).is_err());
    assert!(
        scene::build(
            &s,
            Selection {
                alternates: Alternates::Named("missing".into()),
                ..selection(1)
            }
        )
        .is_err()
    );
}
#[test]
fn mmcif_label_and_author_namespaces_nonpolymer_missing_ids_and_source_remain_separate() {
    let s = cif(&format!("{ROW_A}{ROW_B}"));
    let a = &s.atoms()[0];
    assert_eq!(a.position, [1.0, 2.0, 3.0]);
    assert_eq!(a.label.chain.as_deref(), Some("X"));
    assert_eq!(a.author.chain.as_deref(), Some("Author"));
    assert_eq!(a.label.sequence.as_deref(), Some("5"));
    assert_eq!(a.author.sequence.as_deref(), Some("-2"));
    assert_eq!(a.insertion.as_deref(), Some("B"));
    assert_eq!(s.atoms()[1].label.sequence, None);
    assert_eq!(s.atoms()[1].insertion, None);
    assert_eq!(s.atoms()[1].element.as_deref(), Some("CA"));
    assert_eq!(
        scene::build(
            &s,
            Selection {
                chains: Chains::Label("X".into()),
                ..selection(1)
            }
        )
        .unwrap()
        .points
        .len(),
        1
    );
    assert_eq!(
        scene::build(
            &s,
            Selection {
                chains: Chains::Author(Some("Author".into())),
                ..selection(1)
            }
        )
        .unwrap()
        .points
        .len(),
        2
    );
    assert!(
        scene::build(
            &s,
            Selection {
                chains: Chains::Label("Author".into()),
                ..selection(1)
            }
        )
        .is_err()
    );
    s.validate().unwrap();
}
#[test]
fn cif_lexical_boundaries_quotes_comments_multiline_case_and_column_order() {
    let source = format!(
        "{}{}\n_struct.title\n;first line\n  indented 'quoted' text\n;\n# a comment\n_struct_keywords.text 'a dog\'s life'\n",
        CIF_HEADER.to_ascii_uppercase(),
        ROW_A.replace("Author", "'loop_'")
    );
    let s = import(source.as_bytes(), Format::Mmcif).unwrap();
    assert_eq!(s.atoms()[0].author.chain.as_deref(), Some("loop_"));
    assert_eq!(s.block(), Some("AUTHORED"));
    assert!(s.notes().unprojected_cif_categories.contains("_struct"));
    let tags = CIF_HEADER.lines().skip(2).collect::<Vec<_>>();
    let values = ROW_A.split_whitespace().collect::<Vec<_>>();
    let source = format!(
        "data_reordered\nloop_\n{}\n{}\n",
        tags.iter().rev().copied().collect::<Vec<_>>().join("\n"),
        values.iter().rev().copied().collect::<Vec<_>>().join(" ")
    );
    assert_eq!(
        import(source.as_bytes(), Format::Mmcif).unwrap().atoms(),
        cif(ROW_A).atoms()
    );
    for newline in ["\r\n", "\r"] {
        let text = format!("{CIF_HEADER}{ROW_A}").replace('\n', newline);
        assert_eq!(
            import(text.as_bytes(), Format::Mmcif).unwrap().atoms(),
            cif(ROW_A).atoms()
        );
    }
}
#[test]
fn malformed_and_unsupported_coordinate_files_refuse_without_partial_success() {
    for text in [
        format!("MODEL        1\n{A}"),
        format!("{A}\n{A}"),
        format!("ENDMDL\n{A}"),
        format!("{A}\nMODEL        2\n{B}\nENDMDL"),
        format!("END\n{A}"),
        format!("{A}\nATOM      9"),
    ] {
        assert!(import(text.as_bytes(), Format::Pdb).is_err(), "{text}");
    }
    for text in [
        format!(
            "{CIF_HEADER}{}",
            ROW_A.replace("1.0 2.0 3.0", "NaN 2.0 3.0")
        ),
        format!(
            "{CIF_HEADER}{}",
            ROW_A.replace("1.0 2.0 3.0", "'1.0' 2.0 3.0")
        ),
        format!("{CIF_HEADER}{}", ROW_A.replace("1.0 2.0 3.0", "? 2.0 3.0")),
        format!("{CIF_HEADER}ATOM 1 C"),
        format!("{CIF_HEADER}{ROW_A}\ndata_second\n"),
        format!("{CIF_HEADER}{ROW_A}\n_atom_site.id 2"),
        format!("{CIF_HEADER}{ROW_A}\n_struct.title 'unclosed"),
        format!("{CIF_HEADER}{ROW_A}\n_struct.title\n;unclosed"),
        format!("{CIF_HEADER}{ROW_A}\nloop_\n_bad.one\n_bad.two\nx\n"),
        format!("{CIF_HEADER}{ROW_A}\nglobal_"),
        format!("{CIF_HEADER}{ROW_A}\n_ value"),
        format!("{CIF_HEADER}{ROW_A}\nloop_\n_\nvalue"),
        format!("{CIF_HEADER}{ROW_A}\n_struct.title [unsupported]"),
        format!("#\\#CIF_2.0\n{CIF_HEADER}{ROW_A}"),
    ] {
        assert!(import(text.as_bytes(), Format::Mmcif).is_err(), "{text}");
    }
    assert!(matches!(
        import(&vec![b' '; MAX_SOURCE_BYTES + 1], Format::Pdb),
        Err(Error::Limit(_))
    ));
    assert!(import(b"data_unicode\n_struct.title '\xc3\xa9'", Format::Mmcif).is_err());
    assert!(import(b"data_empty\nloop_\n_atom_site.id\n", Format::Mmcif).is_err());
}
#[test]
fn scene_uses_explicit_connections_only_and_never_merges_models_or_alternates() {
    let no_bonds = pdb(&[A, B, C]);
    assert!(
        scene::build(&no_bonds, selection(1))
            .unwrap()
            .edges
            .is_empty()
    );
    let s = pdb(&[
        A,
        B,
        C,
        "CONECT    1    2    2    3    9",
        "CONECT    2    1",
    ]);
    let scene = scene::build(&s, selection(1)).unwrap();
    assert_eq!(
        scene.edges.iter().map(|e| e.atoms).collect::<Vec<_>>(),
        vec![[0, 1]]
    );
    assert_eq!(scene.unresolved_connections, 2);
    assert_eq!(s.connections().len(), 5); // multiplicity retained, no bond order guessed
    let s = pdb(&[
        "MODEL        1",
        A,
        B,
        "ENDMDL",
        "MODEL        2",
        A,
        "ENDMDL",
        "CONECT    1    2",
    ]);
    assert_eq!(scene::build(&s, selection(1)).unwrap().edges.len(), 1);
    let second = scene::build(&s, selection(2)).unwrap();
    assert!(second.edges.is_empty());
    assert_eq!(second.unresolved_connections, 1);
}
#[test]
fn cif_connection_selectors_and_symmetry_are_not_flattened_to_guessed_bonds() {
    let connections = "loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_symmetry\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_symmetry\na covale X CA 1_555 Y CA 1_555\nb covale X CA 1_555 Y CA 2_555\nc hydrog X CA 1_555 Y CA 1_555\nd covale Missing CA 1_555 Y CA 1_555\n";
    let source = format!("{CIF_HEADER}{ROW_A}{ROW_B}{connections}");
    let s = import(source.as_bytes(), Format::Mmcif).unwrap();
    let scene = scene::build(&s, selection(1)).unwrap();
    assert_eq!(scene.edges.len(), 1);
    assert_eq!(scene.edges[0].atoms, [0, 1]);
    assert_eq!(scene.withheld_symmetry_or_relation, 2);
    assert_eq!(scene.unresolved_connections, 1);
    // A second candidate with no deposited alternate selector makes the endpoint ambiguous.
    let row = ROW_A
        .replace("ATOM 1 ", "ATOM 3 ")
        .replace("CA A ALA", "CA B ALA");
    let s = import(
        format!("{CIF_HEADER}{ROW_A}{row}{ROW_B}{connections}").as_bytes(),
        Format::Mmcif,
    )
    .unwrap();
    assert!(scene::build(&s, selection(1)).unwrap().edges.is_empty());
}
#[test]
fn cif_connection_missing_unknown_and_quoted_identifiers_have_distinct_meanings() {
    let header = "loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_symmetry\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_symmetry\n_struct_conn.pdbx_ptnr1_label_alt_id\n";
    for (value, edges) in [(".", 0), ("?", 0), ("A", 1), ("'A'", 1), ("'.'", 0)] {
        let source =
            format!("{CIF_HEADER}{ROW_A}{ROW_B}{header}a covale X CA 1_555 Y CA 1_555 {value}\n");
        let s = import(source.as_bytes(), Format::Mmcif).unwrap();
        assert_eq!(
            scene::build(&s, selection(1)).unwrap().edges.len(),
            edges,
            "alternate selector {value}"
        );
        if value == "." {
            assert!(
                matches!(&s.connections()[0],Connection::Mmcif{fields} if fields["_struct_conn.pdbx_ptnr1_label_alt_id"]==CifValue::NotApplicable)
            );
        }
        if value == "?" {
            assert!(
                matches!(&s.connections()[0],Connection::Mmcif{fields} if fields["_struct_conn.pdbx_ptnr1_label_alt_id"]==CifValue::Unknown)
            );
        }
        if value == "'.'" {
            assert!(
                matches!(&s.connections()[0],Connection::Mmcif{fields} if fields["_struct_conn.pdbx_ptnr1_label_alt_id"]==CifValue::Text(".".into()))
            );
        }
    }
    let conflicting = format!(
        "{CIF_HEADER}{ROW_A}{ROW_B}{header}_struct_conn.ptnr1_auth_asym_id\na covale X CA 1_555 Y CA 1_555 A WrongAuthor\n"
    );
    let s = import(conflicting.as_bytes(), Format::Mmcif).unwrap();
    assert!(scene::build(&s, selection(1)).unwrap().edges.is_empty());
}
#[test]
fn geometry_is_finite_source_indexed_and_preserves_physical_coordinates() {
    let s = cif(&format!("{ROW_A}{ROW_B}"));
    let before = s.atoms().to_vec();
    let scene = scene::build(&s, selection(1)).unwrap();
    assert_eq!(scene.scaled_center, [0.0; 3]);
    assert_eq!(scene.coordinate_magnitude, 3.0);
    let expected = [1.0_f64, 2.0, 3.0];
    let length = 14.0_f64.sqrt();
    for axis in 0..3 {
        assert!((scene.points[0].position[axis] as f64 - expected[axis] / length).abs() < 1e-7);
        assert_eq!(
            scene.points[1].position[axis],
            -scene.points[0].position[axis]
        );
    }
    assert_eq!(s.atoms(), before);
    assert_eq!(scene.points[1].atom_index, 1);
    let extreme = ROW_A.replace("1.0 2.0 3.0", "1e308 -1e308 0");
    let opposite = ROW_B.replace("-1.0 -2.0 -3.0", "-1e308 1e308 0");
    let s = cif(&format!("{extreme}{opposite}"));
    let scene = scene::build(&s, selection(1)).unwrap();
    assert!(
        scene
            .points
            .iter()
            .flat_map(|p| p.position)
            .all(f32::is_finite)
    );
}
#[test]
fn category_budget_accepts_last_distinct_category_and_existing_rows_then_refuses_new_category() {
    let mut source = format!("{CIF_HEADER}{ROW_A}");
    for n in 0..4096 {
        source.push_str(&format!("_aux{n}.note retained\n"));
    }
    source.push_str("loop_\n_aux4095.extra\none\ntwo\n");
    let s = import(source.as_bytes(), Format::Mmcif).unwrap();
    assert_eq!(s.notes().unprojected_cif_categories.len(), 4096);
    source.push_str("_beyond.note refused\n");
    assert!(matches!(
        import(source.as_bytes(), Format::Mmcif),
        Err(Error::Limit("CIF categories"))
    ));
}
#[test]
fn serialized_projection_must_reconstruct_from_exact_current_source() {
    let s = pdb(&[A, B]);
    let bytes = serde_json::to_vec(&s).unwrap();
    let loaded: Structure = serde_json::from_slice(&bytes).unwrap();
    loaded.validate().unwrap();
    for field in ["engine", "source_hash"] {
        let mut v = serde_json::to_value(&s).unwrap();
        v[field] = "forged".into();
        assert_eq!(
            serde_json::from_value::<Structure>(v).unwrap().validate(),
            Err(Error::Binding)
        );
    }
    let mut v = serde_json::to_value(&s).unwrap();
    v["atoms"][0]["position"][0] = 999.0.into();
    let forged: Structure = serde_json::from_value(v).unwrap();
    assert_eq!(forged.validate(), Err(Error::Binding));
    assert!(scene::build(&forged, selection(1)).is_err());
}
