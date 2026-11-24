#!/usr/bin/env python3
"""Development-only check of the real MCP executable over its standard stdio transport."""
import argparse
import base64
import json
import os
import queue
import re
import subprocess
import sys
import threading
import xml.etree.ElementTree as ET
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser(description=__doc__)
    target_dir=Path(os.environ.get("CARGO_TARGET_DIR", "target"))
    if not target_dir.is_absolute():
        target_dir=root / target_dir
    parser.add_argument("binary",nargs="?",type=Path,default=target_dir / "debug" / ("collagene-mcp.exe" if sys.platform == "win32" else "collagene-mcp"))
    parser.add_argument("--storage-root",type=Path,help="Existing managed storage parent passed to the native server")
    parser.add_argument("--workspace",type=Path,help="Browser-generated synthetic QA snapshot")
    parser.add_argument("--candidate",type=Path,help="Write the binary candidate for the browser QA helper")
    parser.add_argument("--native-workspace",type=Path,help="Synthetic source/target native workspace for retained-run checks")
    args=parser.parse_args()
    binary=args.binary
    server_command=[str(binary.resolve())]
    if args.storage_root is not None:
        server_command.extend(["--storage-root",str(args.storage_root)])
    process = subprocess.Popen(server_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8")
    lines = queue.Queue()

    def read():
        for line in process.stdout:
            lines.put(json.loads(line))
        lines.put(None)

    threading.Thread(target=read, daemon=True).start()
    request_id = 0

    def send(method, params, notification=False, allow_error=False):
        nonlocal request_id
        request_id += 1
        request = {"jsonrpc": "2.0", "method": method, "params": params}
        if not notification:
            request["id"] = request_id
        process.stdin.write(json.dumps(request) + "\n")
        process.stdin.flush()
        if notification:
            return
        while True:
            response = lines.get(timeout=20)
            assert response is not None, "server closed stdout unexpectedly"
            if response.get("id") == request_id:
                if allow_error and "error" in response:
                    return {"rpc_error":response["error"]}
                assert "error" not in response, response
                return response["result"]

    def tool(name, arguments):
        result = send("tools/call", {"name": name, "arguments": arguments})
        assert not result.get("isError"), result
        return result["structuredContent"]

    def reject_tool(name,arguments):
        result=send("tools/call",{"name":name,"arguments":arguments},allow_error=True)
        assert "rpc_error" in result or result.get("isError"),result

    def full_json_artifact(artifact_id):
        parts=[];offset=0
        while True:
            page=tool("read_artifact",{"artifact_id":artifact_id,"start":offset,"max_bytes":16384})
            parts.append(page["text"])
            if page["done"]:break
            offset=page["next_start"]
        return json.loads("".join(parts))

    try:
        initialized = send("initialize", {"protocolVersion": "2025-11-25", "capabilities": {}, "clientInfo": {"name": "collagene-stdio-smoke", "version": "1"}})
        assert initialized["protocolVersion"] == "2025-11-25", initialized
        send("notifications/initialized", {}, notification=True)
        tools = send("tools/list", {})["tools"]
        count_source = (root / "crates/mcp/src/lib.rs").read_text()
        current_tool_count = int(re.search(r"pub const CURRENT_TOOL_COUNT: usize = (\d+);", count_source).group(1))
        assert len(tools) == current_tool_count, len(tools)
        assert {"propose_overlap_assembly", "replay_variant_preview", "import_variants", "inspect_variants", "check_variant_reference", "preview_variant_allele", "export_variants", "dispose_variants", "simulate_topo", "inspect_topo_orientations", "project_topo_product", "check_topo_target", "inspect_recombination", "check_recombination_target", "select_recombination_component", "lookup_enzyme_catalog", "resolve_enzyme_evidence", "calculate_dilution", "balance_mixture", "calculate_dna_molar_ratios", "inspect_restriction", "propose_plan", "verify_plan", "read_sequence_region", "validate_reads", "import_fastq", "export_oligos", "import_workspace_archive", "export_workspace_archive", "begin_workspace_upload", "append_workspace_upload", "finalize_workspace_upload", "abort_workspace_upload", "create_feature_library", "import_feature_library", "export_feature_library", "scan_feature_library", "apply_feature_library", "import_annotations", "apply_annotations", "export_annotations", "export_annotation_dataset", "simulate_end_repair", "simulate_oligo_annealing", "search_document", "view_tracks", "import_alignment", "import_alignment_document", "list_alignments", "read_alignment", "export_alignment", "extract_alignment_rows", "analyze_protein", "design_codons", "check_codon_design", "analyze_dotplot", "assemble_contigs", "assess_methylation", "design_crispr_guides", "search_crispr_off_targets", "design_edit_template", "check_edit_template", "begin_genome_import", "append_genome_import", "finish_genome_import", "inspect_genome", "read_genome_export", "read_genome_interval", "design_genome_interval_guides", "start_genome_off_target_search", "advance_genome_off_target_search", "read_genome_off_target_hits", "cancel_genome_job", "simulate_tailing", "check_overlap_profile", "record_workspace_run", "list_workspace_runs", "replay_workspace_run"} <= {t["name"] for t in tools}
        profile_resource=send("resources/read",{"uri":"collagene://recombination/profile/manufacturer-context-v1"})
        recombination_profile=json.loads(profile_resource["contents"][0]["text"])["snapshot"]
        recombination_docs=tool("import_document",{"text":(root/"fixtures/recombination-bp-demo.gb").read_text()})["documents"]
        def immutable_ref(d): return {"id":d["id"],"revision":d["revision"]}
        def recombination_named(name): return next(d for d in recombination_docs if d["name"].split()[0]==name)
        assert recombination_named("BP_insert")["length"]==59 and recombination_named("BP_donor")["length"]==4761
        recombination_sources=[{"role":"insert","document":immutable_ref(recombination_named("BP_insert"))},{"role":"donor","document":immutable_ref(recombination_named("BP_donor"))}]
        recombination_request={"profile":recombination_profile,"reaction":"BP","limits":{"max_input_bases":1000000,"max_comparisons":20000000,"max_sites":256,"max_combinations":4096,"max_candidates":64,"max_components_per_candidate":4,"max_product_bases":4000000,"max_report_bytes":16777216},"annotation_limits":{"max_work":4000000,"max_copy_bytes":16777216,"max_features":16384,"max_report_bytes":33554432}}
        recombination=tool("inspect_recombination",{"sources":recombination_sources,"request":recombination_request,"limit":1})
        assert recombination["complete"] and recombination["candidate_count"]==1
        assert sorted(p["length"] for p in recombination["candidates"][0]["components"])==[2267,2553]
        recombination_target=immutable_ref(recombination_named("BP_target"))
        recombination_policy={"allow_reverse_complement":False}
        checked=tool("check_recombination_target",{"sources":recombination_sources,"request":recombination_request,"target":recombination_target,"policy":recombination_policy})
        matches=checked["assessment"]["reference"]["matches"]
        assert len(matches)==1 and matches[0]["fully_specified_dna"]
        selection={"batch_revision":recombination["batch_revision"],"candidate_id":matches[0]["candidate"],"component_id":matches[0]["product"]}
        selected=tool("select_recombination_component",{"sources":recombination_sources,"request":recombination_request,"selection":selection,"target":recombination_target,"policy":recombination_policy})
        assert selected["document"]["length"]==2553 and len(selected["reciprocal_components"])==2
        selected_evidence=full_json_artifact(selected["artifact_id"])
        assert selected_evidence["product"]["molecule"]["metadata"]["recombination:reference_only"]=="true"
        assert selected_evidence["product"]["molecule"]["metadata"]["source:legacy:chemistry_unmapped"]
        plan_sources=[{"Input":"recombination_insert"},{"Input":"recombination_donor"}]
        plan_inspection=tool("inspect_recombination",{"sources":recombination_sources,"request":recombination_request,"plan_sources":plan_sources})
        plan_check=tool("check_recombination_target",{"sources":recombination_sources,"request":recombination_request,"plan_sources":plan_sources,"target":recombination_target,"policy":recombination_policy})
        plan_match=plan_check["assessment"]["reference"]["matches"][0]
        plan_selection={"batch_revision":plan_inspection["batch_revision"],"candidate_id":plan_match["candidate"],"component_id":plan_match["product"]}
        reference_plan=tool("propose_plan",{"name":"Stdio reference-only BP","target":recombination_target,"inputs":[{"id":"recombination_insert","document":recombination_sources[0]["document"]},{"id":"recombination_donor","document":recombination_sources[1]["document"]}],"steps":[{"RecombineReference":{"sources":plan_sources,"request":recombination_request,"selection":plan_selection}}]})
        reference_report=tool("verify_plan",{"plan_revision":reference_plan["plan_revision"]})
        assert reference_report["status"]=="Incomplete" and reference_report["engine_version"].endswith("/science-9"),reference_report
        catalog=tool("lookup_enzyme_catalog",{"query":"HaeIII","mode":"Enzyme"})
        assert catalog["status"]=="Unique" and catalog["enzymes"][0]["id"]=="haeiii" and catalog["enzymes"][0]["geometry"]=={"recognition":"GGCC","top_cut":2,"bottom_cut":2}
        resource=send("resources/read",{"uri":catalog["catalog_resource_uri"]})
        assert json.loads(resource["contents"][0]["text"])["identity"]["catalog_content_hash"]==catalog["catalog_content_hash"]
        fact=tool("resolve_enzyme_evidence",{"formulation":"neb-r3101","condition":{"kind":"buffer_activity","buffer":"NEBuffer r3.1","temperature_milli_c":37000}})
        assert fact["resolution"]["status"]=="Known" and fact["resolution"]["value"]=={"kind":"less_than_percent","value":10} and fact["resolution"]["evidence"][0]["qualifiers"]
        missing_fact=tool("resolve_enzyme_evidence",{"formulation":"neb-r3101","condition":{"kind":"buffer_activity","buffer":"not-curated","temperature_milli_c":37000}})
        assert missing_fact["resolution"]["status"]=="Unknown" and "value" not in missing_fact["resolution"]
        dilution_args={"stock":{"kind":"molar","value":"10","unit":"uM"},"target":{"kind":"molar","value":"100","unit":"nM"},"final_volume":{"value":"20","unit":"uL"}}
        dilution=tool("calculate_dilution",dilution_args)
        assert dilution["exact"]["stock_volume"]=={"unit":"nL","value":{"numerator":"200","denominator":"1"}} and dilution["target_reached"]
        rounded=tool("calculate_dilution",dict(dilution_args,quantization={"increment":{"value":"0.3","unit":"uL"},"mode":"NearestTiesUp"}))
        assert not rounded["target_reached"] and rounded["achieved_concentration"]=={"unit":"pmol/nL","value":{"numerator":"3","denominator":"20000"}}
        mixture=tool("balance_mixture",{"final_volume":{"value":"20","unit":"uL"},"components":[{"role":"buffer","volume":{"value":"2","unit":"uL"}},{"role":"DNA","volume":{"value":"0.3","unit":"uL"}}]})
        assert mixture["plan"]["remaining_diluent"]=={"unit":"nL","value":{"numerator":"17700","denominator":"1"}}
        dna_args={"model":{"model_id":"explicit650","grams_per_mole_per_base_pair":"650","assumptions":"Synthetic length-only duplex molecule arithmetic; no end correction."},"anchor":{"role":"vector","length_bp":"3000","mass":{"value":"50","unit":"ng"}},"targets":[{"role":"insert","length_bp":"1000","molar_ratio":"3"}],"require_stock_volumes":False}
        ratio=tool("calculate_dna_molar_ratios",dna_args)
        assert ratio["plan"]["anchor"]["amount"]=={"unit":"pmol","value":{"numerator":"1","denominator":"39"}} and ratio["plan"]["targets"][0]["mass"]=={"unit":"ng","value":{"numerator":"50","denominator":"1"}}
        reject_tool("calculate_dna_molar_ratios",{k:v for k,v in dna_args.items() if k!="model"})
        reject_tool("calculate_dilution",dict(dilution_args,target={"kind":"molar","value":"11","unit":"uM"}))
        reject_tool("balance_mixture",{"final_volume":{"value":"20","unit":"ng"},"components":[]})
        if args.workspace:
            imported=tool("import_workspace",{"data_base64":base64.b64encode(args.workspace.read_bytes()).decode("ascii")})
            docs=imported["documents"]
        else:
            docs = tool("import_document", {"text": ">left\nAAAACCCC\n>right\nCCCCGGGG\n>target\nAAAACCCCGGGG\n"})["documents"]

        def ref(name):
            m = next(m for m in docs if m["name"] == name)
            return {"id": m["id"], "revision": m["revision"]}

        if args.workspace:
            assert imported["target"]==ref("target"), "The browser must independently select the target before export"
        if not args.workspace:
            inspection=tool("inspect_restriction",{"document":ref("target"),"source_role":"mvp-source","enzymes":["EcoRI","BamHI"],"max_sites":1,"max_fragments":1})
            assert inspection["scan_complete"] and inspection["products_determined"]
            assert inspection["site_count"]==0 and inspection["fragment_count"]==1
            evidence=full_json_artifact(inspection["artifact_id"])
            assert evidence["batch_revision"]==inspection["batch_revision"] and evidence["source_revision"]==ref("target")["revision"]
            assert len(evidence["enzymes"])==2 and evidence["fragments"][0]["duplex"]["top_length"]==12
            reject_tool("inspect_restriction",{"document":ref("target"),"source_role":"mvp-source","enzymes":["invented"]})
        digest_docs=tool("import_document",{"text":">digest-source\nAAAAGAATTCCCCCGGATCCAAAA\n>digest-possible\nAAAAGAANTCAAAA\n"})["documents"]
        digest_ref=lambda i:{"id":digest_docs[i]["id"],"revision":digest_docs[i]["revision"]}
        digest_args={"document":digest_ref(0),"source_role":"mvp-source-event","enzymes":["ecori","BamHI"],"max_sites":1,"max_fragments":1}
        digest=tool("inspect_restriction",digest_args)
        assert digest["scan_complete"] and digest["products_determined"] and digest["site_count"]==2 and digest["fragment_count"]==3
        digest_full=full_json_artifact(digest["artifact_id"])
        assert [f["duplex"]["top_length"] for f in digest_full["fragments"]]==[5,10,9]
        digest_page=tool("inspect_restriction",dict(digest_args,site_start=1,fragment_start=1))
        assert digest_page["artifact_id"]==digest["artifact_id"] and digest_page["sites"][0]["recognition_start"]==14 and digest_page["fragments"][0]["top_start"]==5
        possible=tool("inspect_restriction",{"document":digest_ref(1),"source_role":"mvp-possible","enzymes":["EcoRI"]})
        assert possible["scan_complete"] and not possible["products_determined"] and not possible["sites"][0]["certain"]
        reject_tool("analyze_restriction",{"document":digest_ref(1),"enzyme":"EcoRI"})
        catalog_doc=tool("import_document",{"text":">catalog-digest\nAAGGCCTT\n"})["documents"][0]
        catalog_args={"document":{"id":catalog_doc["id"],"revision":catalog_doc["revision"]},"source_role":"catalog-source","enzymes":["HaeIII"]}
        catalog_digest=tool("inspect_restriction",catalog_args)
        assert catalog_digest["products_determined"] and [f["duplex"]["top_length"] for f in catalog_digest["fragments"]]==[4,4]
        assert catalog_digest["catalog_pin"]["catalog_content_hash"]==catalog["catalog_content_hash"] and catalog_digest["catalog_pin"]["enzymes"][0]["id"]=="haeiii"
        wrapped=full_json_artifact(catalog_digest["catalog_artifact_id"])
        raw=full_json_artifact(catalog_digest["artifact_id"])
        assert wrapped["inspection"]==raw and "catalog_pin" not in raw and wrapped["batch_revision"]==catalog_digest["catalog_batch_revision"] and wrapped["selection"]==catalog_digest["catalog_pin"]
        stale_pin=dict(catalog_digest["catalog_pin"],catalog_content_hash="0"*64)
        reject_tool("inspect_restriction",dict(catalog_args,catalog_pin=stale_pin))
        search_options={"query":"CCCC","mode":"Sequence","both_strands":True,"case_sensitive":False,"whole_word":False,"enzymes":[],"max_hits":100,"max_rows":1000,"max_work":10000}
        searched=tool("search_document",{"document":ref("target"),"options":search_options})
        assert searched["source_revision"]==ref("target")["revision"] and searched["complete"] and searched["hit_count"]>=1, searched
        search_artifact=tool("read_artifact",{"artifact_id":searched["artifact_id"],"start":0,"max_bytes":16384})
        assert json.loads(search_artifact["text"])["hits"][0]["parts"][0]["start"]==4
        track_options={"features":True,"primers":True,"enzyme_recognition":False,"enzyme_cuts":False,"translations":True,"orfs":False,"enzyme_names":[],"genetic_code":1,"frames":[0],"orientation":"Forward","known_cds":"Off","max_items":100,"max_residues":100,"max_window_bases":1000,"max_work_units":10000,"min_orf_amino_acids":20}
        tracks=tool("view_tracks",{"document":ref("target"),"window":{"start":0,"end":12},"options":track_options})
        assert tracks["source_revision"]==ref("target")["revision"] and tracks["complete"] and tracks["item_count"]==1, tracks
        assert tracks["items"][0]["residues"][0]["codon_parts"][0]["start"]==0, tracks
        alignment_bytes=b">row_one first\r\nA-CG\r\n>row_two\r\nAT-G\r\n"
        aligned=tool("import_alignment",{"data_base64":base64.b64encode(alignment_bytes).decode("ascii"),"format":"AlignedFasta","alphabet":"Dna","name":"Stdio alignment"})
        alignment_ref=aligned["alignment"]
        assert aligned["row_count"]==2 and aligned["column_count"]==4 and aligned["source_is_current"]
        page=tool("read_alignment",{"alignment":alignment_ref,"row_start":0,"max_rows":1,"column_start":1,"max_columns":2})
        assert page["rows"][0]["aligned"]=="-C" and page["rows"][0]["ungapped_positions"]==[None,1],page
        source=tool("export_alignment",{"alignment":alignment_ref,"format":"Original"})
        original=tool("read_artifact",{"artifact_id":source["artifact_id"],"start":0,"max_bytes":16384})
        assert original["text"].encode("ascii")==alignment_bytes
        extracted=tool("extract_alignment_rows",{"alignment":alignment_ref,"row_ids":[aligned["rows"][0]["id"]]})
        assert extracted["documents"][0]["length"]==3 and extracted["documents"][0]["feature_count"]==0,extracted
        protein=tool("import_document",{"text":">protein [molecule=protein]\nACW"})["documents"][0]
        protein_ref={"id":protein["id"],"revision":protein["revision"]}
        protein_options={"stop_policy":"AllowSingleTerminal","ambiguity":"Unavailable","disulfides":"Reduced","modification_state":"AssumeUnmodified","charge_ph":7.0,"max_residues":1000000}
        properties=tool("analyze_protein",{"document":protein_ref,"options":protein_options})
        assert properties["source_revision"]==protein_ref["revision"] and properties["residue_count"]==3,properties
        assert properties["mass"]["average_da"]["minimum"]>0 and properties["artifact_id"]
        codon_source=tool("import_document",{"text":">required_codon_protein [molecule=protein]\nMK*\n"})["documents"][0]
        codon_source_ref={"id":codon_source["id"],"revision":codon_source["revision"]}
        codon_options={"genetic_code":1,"codon_weights":[{"codon":"AAG","weight":100}],"forbidden_motifs":[],"left_context":"","right_context":"","gc":{"minimum_bp":0,"maximum_bp":10000},"max_homopolymer":12,"window_gc":None,"beam_width":128,"max_work":2000000}
        codon_design=tool("design_codons",{"document":codon_source_ref,"options":codon_options})
        assert codon_design["source_revision"]==codon_source_ref["revision"] and codon_design["candidate_sequence_preview"]=="ATGAAGTAA",codon_design
        codon_candidate={"id":codon_design["candidate"]["id"],"revision":codon_design["candidate"]["revision"]}
        codon_check=tool("check_codon_design",{"source":codon_source_ref,"candidate":codon_candidate,"options":codon_options})
        assert codon_check["complete"] and codon_check["satisfies"] and codon_check["candidate_revision"]==codon_candidate["revision"],codon_check
        codon_artifact=tool("read_artifact",{"artifact_id":codon_check["artifact_id"],"start":0,"max_bytes":16384})
        assert json.loads(codon_artifact["text"])["translated_protein"]=="MK*"
        dot_docs=tool("import_document",{"text":">dot_reference\nACGTTT\n>dot_query\nACG\n"})["documents"]
        dot_reference=next(d for d in dot_docs if d["name"]=="dot_reference")
        dot_query_document=next(d for d in dot_docs if d["name"]=="dot_query")
        dot_ref={"id":dot_reference["id"],"revision":dot_reference["revision"]}
        dot_query={"id":dot_query_document["id"],"revision":dot_query_document["revision"]}
        dot_options={"window":3,"min_matches":3,"stride":1,"both_strands":True,"wrap_circular":True,"max_work":20000000,"max_hits":100000,"max_index_entries":1000000}
        dot=tool("analyze_dotplot",{"reference":dot_ref,"query":dot_query,"options":dot_options,"max_hits":1})
        assert dot["complete"] and not dot["sampled"] and dot["hit_count"]==2 and len(dot["hits"])==1,dot
        dot_next=tool("analyze_dotplot",{"reference":dot_ref,"query":dot_query,"options":dot_options,"hit_start":1,"max_hits":1})
        assert dot_next["done"] and dot_next["artifact_id"]==dot["artifact_id"]
        dot_full=tool("read_artifact",{"artifact_id":dot["artifact_id"],"start":0,"max_bytes":16384})
        dot_report=json.loads(dot_full["text"])
        assert any(h["reference_start"]==1 and h["query_start"]==2 and h["reverse"] for h in dot_report["hits"]),dot_report
        dot_svg=tool("read_artifact",{"artifact_id":dot["svg_artifact_id"],"start":0,"max_bytes":16384})
        assert "<svg" in dot_svg["text"] and dot_svg["done"]
        dot_metadata_node=ET.fromstring(dot_svg["text"]).find(".//{http://www.w3.org/2000/svg}metadata")
        assert dot_metadata_node is not None
        dot_metadata=json.loads(dot_metadata_node.text)
        assert dot_metadata["reference_revision"]==dot_ref["revision"] and dot_metadata["query_revision"]==dot_query["revision"]
        assert dot_metadata["engine_version"]==dot_report["engine_version"] and dot_metadata["options"]==dot_options
        assert dot_metadata["complete"] and not dot_metadata["sampled"] and "Unverified" in dot_metadata["authority"]
        contig_docs=tool("import_document",{"text":">contig_left\nAGTCGACTA\n>contig_middle\nACTAGGCTT\n>contig_right\nGCTTACG\n"})["documents"]
        contig_reads=[]
        for role in ["left","middle","right"]:
            document=next(d for d in contig_docs if d["name"]=="contig_"+role)
            contig_reads.append({"id":role,"document":{"id":document["id"],"revision":document["revision"]}})
        contig_options={"trim":{"quality_windows":False,"raw_score_threshold":20,"window_size":5,"minimum_passing":4,"min_retained_bases":1,"allow_unscored":True},"scoring":{"match_score":2,"ambiguous_score":0,"mismatch_score":-3,"gap_open":-4,"gap_extend":-1},"min_overlap":4,"max_error_bp":0,"both_strands":False,"min_support":1,"require_both_strands":False,"max_pair_cells":4000000,"max_total_cells":128000000,"max_edges":4096,"max_columns":100000}
        contig=tool("assemble_contigs",{"reads":contig_reads,"options":contig_options,"max_contigs":1,"max_overlaps":1})
        assert contig["contig_count"]==1 and contig["contigs"][0]["sequence_preview"]=="AGTCGACTAGGCTTACG",contig
        contig_cached=tool("assemble_contigs",{"reads":contig_reads,"options":contig_options,"overlap_start":1,"max_overlaps":1})
        assert contig_cached["reused_result"] and contig_cached["artifact_id"]==contig["artifact_id"]
        contig_text=[];contig_at=0
        while True:
            page=tool("read_artifact",{"artifact_id":contig["artifact_id"],"start":contig_at,"max_bytes":16384})
            contig_text.append(page["text"])
            if page["done"]:break
            contig_at=page["next_start"]
        contig_report=json.loads("".join(contig_text))
        assert contig_report["batch_revision"]==contig["batch_revision"] and len(contig_report["contigs"][0]["placements"])==3
        consensus=next(d for d in contig["consensus_documents"] if d["length"]==17)
        consensus_ref={"id":consensus["id"],"revision":consensus["revision"]}
        assert tool("read_sequence_region",{"document":consensus_ref,"start":0,"end":17})["sequence"]=="AGTCGACTAGGCTTACG"
        methylation_document=tool("import_document",{"text":">methylation_source\nTCTAGATCGAATTCGGATCCCTCGAGGGTCTCCGTCTC\n"})["documents"][0]
        methylation_ref={"id":methylation_document["id"],"revision":methylation_document["revision"]}
        methylation_options={"enzymes":["EcoRI","BamHI","HindIII","EcoRV","SmaI","NotI","XbaI","BsaI","PstI","BsmBI","BbsI","SapI"],"dam":"Unknown","dcm":"Unknown","cpg":"Unknown","max_comparisons":20000000,"max_sites":4096}
        methylation=tool("assess_methylation",{"document":methylation_ref,"options":methylation_options,"max_sites":1})
        assert methylation["source_revision"]==methylation_ref["revision"] and len(methylation["profiles"])==12 and methylation["site_count"]>1,methylation
        methylation_next=tool("assess_methylation",{"document":methylation_ref,"options":methylation_options,"site_start":1,"max_sites":1})
        assert methylation_next["artifact_id"]==methylation["artifact_id"] and methylation_next["next_site_start"]==2
        methylation_text=[];methylation_at=0
        while True:
            page=tool("read_artifact",{"artifact_id":methylation["artifact_id"],"start":methylation_at,"max_bytes":16384})
            methylation_text.append(page["text"])
            if page["done"]:break
            methylation_at=page["next_start"]
        methylation_report=json.loads("".join(methylation_text))
        assert methylation_report["batch_revision"]==methylation["batch_revision"] and any(site["outcome"]=="Unknown" for site in methylation_report["sites"])
        crispr_spacer="ACGTACGTACGTACGTACGT"
        crispr_document=tool("import_document",{"text":">crispr_reference\n"+crispr_spacer+"AGGTTT"+crispr_spacer+"AGG\n"})["documents"][0]
        crispr_ref={"id":crispr_document["id"],"revision":crispr_document["revision"]}
        crispr_guides=tool("design_crispr_guides",{"document":crispr_ref,"max_candidates":1})
        assert crispr_guides["source_revision"]==crispr_ref["revision"] and crispr_guides["candidate_count"]>=2,crispr_guides
        crispr_guide=crispr_guides["candidates"][0]
        crispr_refs=[{"id":"intended","document":crispr_ref,"build":"declared-stdio-build","accession":"test.1"},{"id":"duplicate-record","document":crispr_ref,"build":None,"accession":None}]
        crispr_options={"profile":crispr_guides["options"]["profile"],"max_mismatches":0,"intended":{"reference_id":"intended","reference_revision":crispr_ref["revision"],"guide_positions":crispr_guide["guide_positions"],"pam_positions":crispr_guide["pam_positions"],"reverse":crispr_guide["reverse"]},"reference_scope":"Explicit stdio selected fixtures","max_comparisons":100000000,"max_index_entries":1000000,"max_hits":4096}
        crispr_query={"id":crispr_guide["id"],"sequence":crispr_guide["guide_sequence"]}
        crispr_off=tool("search_crispr_off_targets",{"query":crispr_query,"references":crispr_refs,"options":crispr_options,"max_hits":1})
        crispr_full=full_json_artifact(crispr_off["artifact_id"])
        assert crispr_off["hit_count"]==4 and sum(hit["is_intended"] for hit in crispr_full["hits"])==1,crispr_full
        edit_documents=tool("import_document",{"text":">edit_source\nAAAACCCCGGGGTTTT\n>edit_target\nAAAACCTCGGGGTTTT\n"})["documents"]
        edit_refs={d["name"]:{"id":d["id"],"revision":d["revision"]} for d in edit_documents}
        edit_options={"edit_start":6,"delete_length":1,"replacement":"T","left_arm":3,"right_arm":3,"reverse_donor":False,"max_cut_distance":100}
        edit=tool("design_edit_template",{"source":edit_refs["edit_source"],"target":edit_refs["edit_target"],"options":edit_options})
        assert edit["consistent"] and edit["donor_preview"]=="ACCTCGG" and edit["donor_document"],edit
        edit_donor={"id":edit["donor_document"]["id"],"revision":edit["donor_document"]["revision"]}
        edit_check=tool("check_edit_template",{"source":edit_refs["edit_source"],"target":edit_refs["edit_target"],"donor":edit_donor,"options":edit_options})
        assert edit_check["check"]["consistent"] and full_json_artifact(edit_check["artifact_id"])["target_matches"],edit_check
        # Managed-reference lifecycle uses bounded transport and actual private files.
        genome_dna="A"*65526+crispr_spacer+"AGG"+"A"*77
        genome_bytes=(">chr\n"+genome_dna+"\n").encode("ascii")
        genome_build={"name":"Declared stdio reference","assembly_accession":None,"version":None,"provenance":[],"expected_components":[{"id":"chr","source_id":"chr","accession":None,"length":len(genome_dna),"sequence_hash":None,"topology":"Linear"}]}
        genome_import=tool("begin_genome_import",{"format":"Fasta","expected_bytes":len(genome_bytes),"build":genome_build})
        for offset in range(0,len(genome_bytes),7001):
            chunk=genome_bytes[offset:offset+7001]
            progress=tool("append_genome_import",{"import_id":genome_import["import_id"],"offset":offset,"data_base64":base64.b64encode(chunk).decode("ascii")})
            assert progress["received"]==offset+len(chunk)
        while True:
            genome_ready=tool("finish_genome_import",{"import_id":genome_import["import_id"],"max_blocks":1})
            if genome_ready["done"]:break
            assert genome_ready["dataset"] is None
        genome_dataset=genome_ready["dataset"];genome_ref=genome_dataset["genome"]
        assert genome_dataset["semantic_verified_at_import"] and genome_dataset["expected_components_complete"]
        reject_tool("read_genome_export",{"genome":genome_ref,"filesystem_path":"/etc/passwd"})
        genome_inspect=tool("inspect_genome",{"genome":genome_ref,"max_contigs":1})
        assert genome_inspect["contigs"][0]["block_count"]==2
        genome_interval={"contig_id":"chr","start":65526,"end":65549}
        genome_read=tool("read_genome_interval",{"genome":genome_ref,"options":{"source_role":"reference","interval":genome_interval}})
        assert genome_read["read"]["sequence"]==crispr_spacer+"AGG"
        genome_guide_options={"source_role":"reference","interval":{"contig_id":"chr","start":65526,"end":65527},"guide":crispr_guides["options"]}
        genome_guides=tool("design_genome_interval_guides",{"genome":genome_ref,"options":genome_guide_options,"max_candidates":1})
        genome_guide=genome_guides["candidates"][0]
        assert 65536 in genome_guide["guide_positions"]
        genome_search_options={"reference_role":"reference","profile":crispr_guides["options"]["profile"],"max_mismatches":0,"intended":{"reference_role":"reference","manifest_hash":genome_ref["manifest_hash"],"content_revision":genome_ref["content_revision"],"contig_id":"chr","guide_positions":genome_guide["guide_positions"],"pam_positions":genome_guide["pam_positions"],"reverse":genome_guide["reverse"]},"max_total_comparisons":1000000000000,"max_index_entries_per_block":1000000,"max_hits":4096}
        genome_search=tool("start_genome_off_target_search",{"genome":genome_ref,"query":{"id":"stdio-global","sequence":crispr_spacer},"options":genome_search_options})
        assert not genome_search["progress"]["source_integrity_verified"]
        reject_tool("advance_genome_off_target_search",{"job_id":genome_search["job_id"],"checkpoint":{"finished":True,"source_integrity_verified":True}})
        while not genome_search["progress"]["finished"]:
            genome_search=tool("advance_genome_off_target_search",{"job_id":genome_search["job_id"],"max_blocks":1})
        assert genome_search["progress"]["complete"] and genome_search["progress"]["source_integrity_verified"]
        genome_hits=tool("read_genome_off_target_hits",{"job_id":genome_search["job_id"],"limit":1})
        assert len(genome_hits["hits"])==1 and genome_hits["hits"][0]["is_intended"]
        genome_native=[];genome_at=0
        while True:
            page=tool("read_genome_export",{"genome":genome_ref,"offset":genome_at,"max_bytes":4096})
            genome_native.append(base64.b64decode(page["data_base64"],validate=True))
            if page["done"]:break
            genome_at=page["next_offset"]
        genome_native=b"".join(genome_native)
        assert genome_native.startswith(b"CGREF\0\x01\0") and len(genome_native)==genome_dataset["file_length"]
        native_import=tool("begin_genome_import",{"format":"Cgref","expected_bytes":len(genome_native)})
        for offset in range(0,len(genome_native),7001):
            tool("append_genome_import",{"import_id":native_import["import_id"],"offset":offset,"data_base64":base64.b64encode(genome_native[offset:offset+7001]).decode("ascii")})
        while True:
            native_ready=tool("finish_genome_import",{"import_id":native_import["import_id"],"max_blocks":1})
            if native_ready["done"]:break
        assert native_ready["dataset"]["genome"]["manifest_hash"]==genome_ref["manifest_hash"]
        tool("cancel_genome_job",{"id":genome_search["job_id"]})
        tool("cancel_genome_job",{"id":genome_ref["id"],"genome":genome_ref})
        reject_tool("inspect_genome",{"genome":genome_ref})
        native_ref=native_ready["dataset"]["genome"]
        tool("cancel_genome_job",{"id":native_ref["id"],"genome":native_ref})
        assert not tool("inspect_genome",{})["datasets"]
        tracked_docs=tool("import_document",{"text":">tracked_source\nAAAAGAATTCCCC\n>tracked_target\nAAAAG\n"})["documents"]
        tracked_refs={d["name"]:{"id":d["id"],"revision":d["revision"]} for d in tracked_docs}
        tracked_profile=tool("list_enzymes",{})["checked_digest_profile_revision"]
        tracked_steps=[{"DeclareSourceModifications":{"template":{"Input":"s"},"assumptions":{"dam":"Unmodified","dcm":"Unmodified","cpg":"Unmodified"}}},{"CheckedDigest":{"template":{"Step":0},"enzymes":["EcoRI"],"profile_revision":tracked_profile,"fragment_index":0}}]
        tracked=tool("propose_plan",{"name":"Stdio checked digest","target":tracked_refs["tracked_target"],"inputs":[{"id":"s","document":tracked_refs["tracked_source"]}],"steps":tracked_steps})
        tracked_report=tool("verify_plan",{"plan_revision":tracked["plan_revision"]})
        assert tracked_report["status"]=="ConsistentUnderAssumptions" and tracked_report["engine_version"].endswith("/science-9"),tracked_report
        assert tracked_report["products"][1]["length"]==5 and tracked_report["products"][1]["tracked_modifications"]["generation"]
        candidate = tool("propose_plan", {"name": "Stdio assembly", "target": ref("target"), "inputs": [{"id": name, "document": ref(name)} for name in ["left", "right"]], "steps": [{"Overlap": {"fragments": [{"Input": "left"}, {"Input": "right"}], "overlaps": [4], "circular": False}}]})
        assert candidate["verified"] is False
        recipe_candidate=tool("propose_overlap_assembly",{
            "inputs":[{"id":name,"document":ref(name)} for name in ["left","right"]],
            "target":ref("target"),
            "recipe":{"name":"Stdio prepared overlap recipe","fragments":[
                {"role":"left-fragment","input_id":"left","preparation":"WholeLinear","reverse":False},
                {"role":"right-fragment","input_id":"right","preparation":"WholeLinear","reverse":False}],
                "overlaps":[4],"circular":False,
                "options":{"method":"HiFi","scan_repeated_overlaps":False,"max_comparisons":10000,"max_occurrences":100}}})
        assert recipe_candidate["verified"] is False and recipe_candidate["preparation_complete"],recipe_candidate
        assert recipe_candidate["verification"]["status"]=="ConsistentUnderAssumptions",recipe_candidate
        recipe_evidence=full_json_artifact(recipe_candidate["artifact_id"])
        assert len(recipe_evidence["compiled"]["rows"])==2 and len(recipe_evidence["prepared"])==2
        assert recipe_evidence["profile"]["geometry_valid"] and not recipe_evidence["profile"]["exhaustive_alternative_search"]
        assert tool("verify_plan",{"plan_revision":recipe_candidate["plan_revision"]})["status"]=="ConsistentUnderAssumptions"
        plan = {"plan_revision": candidate["plan_revision"]}
        report = tool("verify_plan", plan)
        assert report["status"] == "ConsistentUnderAssumptions", report
        assert report["products"][0]["length"] == 12
        resource = send("resources/read", {"uri": "collagene://report/" + plan["plan_revision"]})
        assert json.loads(resource["contents"][0]["text"])["plan_revision"] == plan["plan_revision"]
        artifact = tool("export_protocol", plan)
        html = tool("read_artifact", {"artifact_id": artifact["artifact_id"], "start": 0, "max_bytes": 16384})
        assert html["done"] and "ConsistentUnderAssumptions" in html["text"]
        candidate_artifact=tool("export_candidate",plan)
        encoded=""
        offset=0
        while True:
            chunk=tool("read_artifact",{"artifact_id":candidate_artifact["artifact_id"],"start":offset,"max_bytes":16384})
            assert chunk["encoding"]=="base64"
            encoded+=chunk["text"]
            if chunk["done"]:
                break
            offset=chunk["next_start"]
        candidate_bytes=base64.b64decode(encoded,validate=True)
        assert candidate_bytes.startswith(b"CGP\0") and candidate_bytes[4:6]==bytes([7,0])
        if args.candidate:
            args.candidate.write_bytes(candidate_bytes)
        # Route source-bound annotation/library tools over the actual stdio executable.
        # Every application returns a new session document; the original plan remains bound.
        annotations=tool("import_annotations",{"document":ref("target"),"text":"stdio_target\t0\t4\tStdio local part\t0\t+\n","format":"Bed","seqid":"stdio_target"})
        assert annotations["proposal_count"]==1 and annotations["applied"] is False
        applied=tool("apply_annotations",{"document":ref("target"),"dataset":annotations["dataset"],"selected":[0]})["documents"][0]
        applied_ref={"id":applied["id"],"revision":applied["revision"]}
        library=tool("create_feature_library",{"name":"Stdio own features","version":"1","features":[{"document":applied_ref,"feature_index":applied["feature_count"]-1}]})["library"]
        scan=tool("scan_feature_library",{"document":ref("target"),"library":library,"options":{"both_strands":True,"max_mismatches":0,"max_comparisons":20000,"max_hits":100}})
        assert scan["complete"] and scan["hit_count"]>=1 and scan["applied"] is False
        first=scan["hits"][0]
        assert first["start"]==0 and first["length"]==4 and first["needs_review"] is True
        new_document=tool("apply_feature_library",{"document":ref("target"),"library":library,"scan":scan["scan"],"selected":[first["id"]]})["documents"][0]
        assert new_document["revision"]!=ref("target")["revision"]
        repair_text="LOCUS       repair_source             10 bp    DNA     circular SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 gaattcaaaa\n//\n"
        repair_source=tool("import_document",{"text":repair_text})["documents"][0]
        repair_ref={"id":repair_source["id"],"revision":repair_source["revision"]}
        repair=tool("simulate_end_repair",{"document":repair_ref,"enzymes":["EcoRI"],"fragment_index":0,"method":"Blunt"})
        assert repair["length"]==14 and repair["top_added_count"]==4 and repair["bottom_added_count"]==4
        tailed=tool("simulate_tailing",{"document":repair_ref,"recipe":{"Digest":{"enzymes":["EcoRI"],"fragment_index":0,"blunt_before":True}},"base":"T"})
        assert tailed["top_length"]==15 and tailed["bottom_start"]==-1 and len(tailed["additions"])==2,tailed
        profile_options={"method":"HiFi","scan_repeated_overlaps":False,"max_comparisons":10000,"max_occurrences":100}
        profiled=tool("check_overlap_profile",{"documents":[ref("left"),ref("right")],"reverse":[False,False],"overlaps":[4],"circular":False,"options":profile_options})
        assert profiled["geometry_valid"] and not profiled["recommendation_compliant"] and not profiled["exhaustive_alternative_search"],profiled
        oligos=tool("import_document",{"text":">anneal_top\nACGTCG\n>anneal_bottom\nCGACGT\n"})["documents"]
        top={"id":oligos[0]["id"],"revision":oligos[0]["revision"]};bottom={"id":oligos[1]["id"],"revision":oligos[1]["revision"]}
        anneal=tool("simulate_oligo_annealing",{"top":top,"bottom":bottom,"options":{"min_paired_bases":4,"max_comparisons":20000,"max_candidates":100,"top_chemistry":{"five_prime":"Hydroxyl","three_prime":"Hydroxyl"},"bottom_chemistry":{"five_prime":"Hydroxyl","three_prime":"Hydroxyl"}}})
        assert anneal["unique"] and anneal["complete"] and anneal["candidate_count"]==1
        if args.native_workspace:
            native=tool("import_workspace_archive",{"data_base64":base64.b64encode(args.native_workspace.read_bytes()).decode("ascii")})["workspace"]
            workspace_ref={"id":native["id"],"revision":native["revision"]}
            initial=tool("list_workspace_runs",workspace_ref)
            assert initial["run_count"]==0 and initial["event_count"]==2
            events=initial["events"]
            source_event=next(e for e in events if e["document_id"]=="document-1")
            target_event=next(e for e in events if e["document_id"]=="document-2")
            event_ref=lambda e:{k:e[k] for k in ["document_id","event_id","revision"]}
            document_ref=lambda e:{"id":e["revision"],"revision":e["revision"]}
            run_plan=tool("propose_plan",{"name":"Stdio retained run","target":document_ref(target_event),"inputs":[{"id":"s","document":document_ref(source_event)}],"steps":[{"Overlap":{"fragments":[{"Input":"s"}],"overlaps":[],"circular":False}}]})
            recorded=tool("record_workspace_run",{"workspace":workspace_ref,"plan_revision":run_plan["plan_revision"],"target_binding":event_ref(target_event),"input_bindings":[{"input_id":"s","event":event_ref(source_event)}]})
            assert recorded["workspace"]["id"]!=native["id"] and len(recorded["receipt"]["products"])==1
            new_ref={k:recorded["workspace"][k] for k in ["id","revision"]}
            assert tool("list_workspace_runs",workspace_ref)["run_count"]==0
            assert tool("list_workspace_runs",new_ref)["run_count"]==1
            replay=tool("replay_workspace_run",{"workspace":new_ref,"run_id":recorded["receipt"]["run_id"]})
            assert replay["products_match"] and replay["engine_matches"] and replay["status_matches"]
        print(f"MCP stdio smoke passed: {len(tools)} tools, 2025-11-25 negotiation, tools/resources, candidate, deterministic verification, HTML artifact, revision-bound unified search/browser tracks, retained alignment exchange, protein properties, codon design and independent candidate checks, annotation evidence, feature-library scan/application, reconstructed end repair/tailing, explicit overlap profiles and declared oligo annealing; optional immutable retained-run record/list/replay, source-bound dot plots/SVG, and reference-free source-bound read contigs and explicitly declared methylation-context assessment and reconstructed tracked-generation checked digest plans, source-bound CRISPR guides/selected-reference off-target evidence and independently checked literal editing donors and file-backed streamed references with global seam guides/private cooperative full-manifest search; exact unit-aware dilution/quantization, named volume balance and caller-model DNA molecule ratios; versioned catalog lookup/contextual evidence and compatible catalog-pinned digest artifacts.")
    finally:
        process.stdin.close()
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.terminate()
            process.wait(timeout=5)
        error = process.stderr.read()
        assert process.returncode == 0, (process.returncode, error)


if __name__ == "__main__":
    main()
