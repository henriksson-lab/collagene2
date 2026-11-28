#!/usr/bin/env python3
"""Maintainer-only deterministic compact Corpus emitter. Never called by the app.

Python mirrors the Rust Prepared factual projection. Rust crosschecking and the
independent, previously captured translation oracle remain separate validation.
"""
import hashlib
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
JSON_SHA = "a557dc47e4946b015fcb5414c29c275783c464848fd13b8f4a0386c4d6e686e8"
FASTA_SHA = "27403baa82141e88a065ddedbc9fc6b4bc4594fc9ecf64d8c4056e98faf4d160"
URL = "https://ogrdb.airr-community.org/api_v2/germline/set/10090.C57BL.sl.6%20IGH.C57BL.sl.6/5"
V_POLICY = "Ordinals 1..104 of the paired provider-declared IMGT-gapped triplet columns, crosschecked against exact source coding DNA and CDR nucleotide intervals; JSON labels are absent."
J_POLICY = "Pinned-source observation: j_cdr3_end-1 locates W and has phase frame-1. The AIRR textual frame interpretation (4-frame)%3 differs for these records; both are retained. No generic frame reinterpretation or J numbering labels are supplied."
CODE = dict(zip(map("".join, itertools.product("TCAG", repeat=3)),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))

def require(condition, reason):
    if not condition:
        raise ValueError(reason)

def pinned(name, sha, limit):
    data = (ROOT / name).read_bytes()
    require(len(data) <= limit and hashlib.sha256(data).hexdigest() == sha, name + " byte/hash pin")
    return data

def translate(dna):
    require(len(dna) % 3 == 0, "partial translated codon")
    return "".join("-" if dna[i:i+3] == "..." else CODE[dna[i:i+3]] for i in range(0, len(dna), 3))

def main():
    raw = json.loads(pinned("source.json", JSON_SHA, 512 * 1024))["GermlineSet"]
    require(len(raw) == 1, "one source set")
    raw = raw[0]
    require((raw["germline_set_id"], raw["germline_set_ref"], raw["release_version"], raw["species"]["id"], raw["species_subgroup"], raw["locus"]) ==
            ("OGRDB:G00004", "OGRDB:C57BL/6 IGH.5", 5, "NCBITAXON:10090", "C57BL/6", "IGH"), "source identity")
    records = []
    for line in pinned("gapped.fasta", FASTA_SHA, 64 * 1024).decode("ascii").splitlines():
        if line.startswith(">"):
            records.append([line[1:], ""])
        else:
            require(records and line and set(line) <= set("ACGT."), "FASTA shape/alphabet")
            records[-1][1] += line
    alleles = raw["allele_descriptions"]
    require(len(alleles) == len(records) == 114, "paired inventory")
    by_dna = {}
    for record in records:
        by_dna.setdefault(record[1].replace(".", ""), []).append(record)
    require(len(by_dna) == 114 and all(len(value) == 1 for value in by_dna.values()), "ambiguous paired DNA")
    v_entries, j_entries, excluded, atypical = [], [], [], []
    joined = set()
    for allele in alleles:
        dna = allele["coding_sequence"]
        require(dna in by_dna and dna not in joined, "nonunique/missing source join")
        joined.add(dna)
        require(allele["species"]["id"] == "NCBITAXON:10090" and allele["species_subgroup"] == "C57BL/6" and allele["locus"] == "IGH", "allele species/strain/locus")
        if not allele["functional"] or allele["sequence_type"] == "D":
            excluded.append(dict(allele_ref=allele["allele_description_ref"], reason="provider-not-functional" if not allele["functional"] else "D-not-a-V-or-J-reference-model", sequence_type=allele["sequence_type"]))
            continue
        common = dict(source=0, allele_ref=allele["allele_description_ref"], allele_id=allele["allele_description_id"], label=allele["label"], release_version=int(allele["release_version"]))
        require(common["release_version"] == allele["release_version"], "nonintegral allele release")
        if allele["sequence_type"] == "V":
            aligned = by_dna[dna][0][1][:312]
            require(len(aligned) == 312, "missing V prefix columns")
            protein = translate(aligned)
            require("*" not in protein and dna.startswith(aligned.replace(".", "")), "V prefix coding correspondence")
            delineations = [d for d in allele["v_gene_delineations"] if d["delineation_scheme"] == "IMGT"]
            require(len(delineations) == 1, "single IMGT delineation")
            delineation = delineations[0]
            require(delineation["alignment_labels"] is None and delineation["aligned_sequence"] is None, "unexpected upstream labels")
            labels = [i+1 for i, aa in enumerate(protein) if aa != "-"]
            cdrs = []
            for kind, low, high in [(1, 27, 38), (2, 56, 65)]:
                indices = [i for i, column in enumerate(labels) if low <= column <= high]
                start, end = indices[0], indices[-1] + 1
                require(delineation[f"cdr{kind}_start"] == start * 3 + 1 and delineation[f"cdr{kind}_end"] == end * 3, "CDR source coordinates")
                cdrs.append(dict(kind=kind, start=start, end=end))
            require(delineation["fwr3_end"] == len(aligned.replace(".", "")), "V endpoint")
            common.update(aligned_nucleotides=aligned, aligned_protein=protein, alignment_labels=list(map(str, range(1, 105))), cdrs=cdrs,
                          non_cysteine_anchors=[n for n in (23, 104) if protein[n-1] != "C"])
            v_entries.append(common)
            for column, canonical in [(23, "C"), (41, "W"), (104, "C")]:
                if protein[column-1] != canonical:
                    atypical.append(dict(allele_ref=common["allele_ref"], imgt_column=column, codon=aligned[(column-1)*3:column*3], residue=protein[column-1], canonical_residue=canonical))
        else:
            require(allele["sequence_type"] == "J", "unsupported record type")
            frame, anchor = allele["j_codon_frame"], allele["j_cdr3_end"] - 1
            require(1 <= frame <= 3 and anchor % 3 == frame-1 and anchor % 3 != (4-frame) % 3, "pinned J phase facts")
            end = anchor + (len(dna)-anchor)//3*3
            fr4 = translate(dna[anchor:end])
            require(len(fr4) == 11 and fr4[0] == "W" and fr4[1] == fr4[3] == "G" and dna[end:] == "G", "FR4 observed anchor/extent")
            require(allele["sequence"][allele["gene_start"]-1:allele["gene_end"]] == dna, "J coding source interval")
            common.update(original_sequence=allele["sequence"], coding_nucleotides=dna, gene_start=allele["gene_start"], gene_end=allele["gene_end"],
                          provider_codon_frame=frame, provider_j_cdr3_end=allele["j_cdr3_end"], anchor_offset=anchor,
                          fr4_nucleotides=dna[anchor:end], fr4_protein=fr4,
                          fr4_positions=list(map(str, range(118, 129))), trailing_nucleotides=dna[end:])
            j_entries.append(common)
    for entries in [v_entries, j_entries, excluded]:
        entries.sort(key=lambda entry: entry["allele_ref"])
    atypical.sort(key=lambda a: (a["allele_ref"], a["imgt_column"]))
    require((len(v_entries), len(j_entries), len(excluded), len(atypical)) == (95, 4, 15, 4), "derived inventory")
    oracle = json.loads((ROOT / "independent-oracle.json").read_bytes())
    require({e["label"]:e["aligned_protein"] for e in v_entries} == {e["label"]:e["aligned104_protein"] for e in oracle["v_entries"]}, "independent V translation comparison")
    require({e["label"]:e["fr4_protein"] for e in j_entries} == {e["label"]:e["fr4_protein"] for e in oracle["j_entries"]}, "independent J translation comparison")
    source = {key: raw[key] for key in ["release_date", "author", "locus"]}
    source.update(set_id=raw["germline_set_id"], set_ref=raw["germline_set_ref"], release_version=5, species="Mus musculus", taxonomy_id="NCBITAXON:10090", url=URL, sha256=JSON_SHA, original_records=114, v_records=95, j_records=4, excluded=excluded)
    corpus = dict(version=3, species="MouseC57Bl6", source_policy={"MousePinnedColumns":dict(gapped_sha256=FASTA_SHA, v_label_evidence=V_POLICY, atypical_anchors=atypical, j_phase_policy=J_POLICY)},
                  license="CC0-1.0", license_url="https://ogrdb.airr-community.org/render_page/licensing_statement.html", scope="mouse-C57BL6-functional-IGH-V-prefix-and-observed-J-FR4",
                  j_label_rule="IMGT-118-sequential-complete-J-codons", sources=[source], v_entries=v_entries, j_entries=j_entries)
    output = (json.dumps(corpus, separators=(",", ":"), ensure_ascii=False) + "\n").encode()
    require(len(output) < 1024*1024, "compact runtime bound")
    (ROOT / "references.json").write_bytes(output)
    provenance = dict(emitter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), runtime_sha256=hashlib.sha256(output).hexdigest(), source_sha256=JSON_SHA, gapped_sha256=FASTA_SHA,
                      emitter="Maintainer-only Python stdlib; Rust Prepared projection and independent translation oracle crosscheck separately", bytes=len(output), v_records=95, j_records=4,
                      model_admission="92 canonical C23/W41/C104 V references; three atypical functional alleles retained in source corpus and excluded from this mouse model", j_rule_origin="Authored conditional IMGT 118..128 complete-codon FR4 labels; provider does not supply J numbering labels")
    (ROOT / "runtime-provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    print(json.dumps(provenance))

if __name__ == "__main__":
    main()
