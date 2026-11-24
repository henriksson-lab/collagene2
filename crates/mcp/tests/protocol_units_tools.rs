//! Actual MCP SDK; literal expected arithmetic is independent of library result builders.
use collagene_mcp::CollageneServer;
use rmcp::{
    RoleClient, ServiceExt,
    model::CallToolRequestParams,
    service::{RunningService, ServiceError},
};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
async fn setup() -> (Client, tokio::task::JoinHandle<()>) {
    let (a, b) = tokio::io::duplex(65536);
    let server = tokio::spawn(async move {
        CollageneServer::default()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap();
    });
    (().serve(b).await.unwrap(), server)
}
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(r.is_error, Some(true), "{name}: {:?}", r.content);
    r.structured_content.unwrap()
}
async fn reject(c: &Client, name: &str, args: Value, code: Option<&str>) -> Value {
    let result = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    let data = match result {
        Err(ServiceError::McpError(e)) => e.data.unwrap_or(Value::Null),
        Ok(r) => {
            assert_eq!(r.is_error, Some(true), "expected rejection: {r:?}");
            assert!(
                code.is_none(),
                "library error must retain structured variant data: {r:?}"
            );
            Value::Null
        }
        Err(e) => panic!("expected semantic or parameter rejection, not transport failure: {e}"),
    };
    if let Some(code) = code {
        assert_eq!(data["domain"], "protocol_units");
        assert_eq!(data["code"], code);
    }
    data
}
fn q(unit: &str, n: &str, d: &str) -> Value {
    json!({"unit":unit,"value":{"numerator":n,"denominator":d}})
}
fn dilution() -> Value {
    json!({"stock":{"kind":"molar","value":"10","unit":"uM"},"target":{"kind":"molar","value":"100","unit":"nM"},"final_volume":{"value":"20","unit":"uL"}})
}
#[tokio::test]
async fn dilution_exact_units_quantization_achieved_concentration_and_typed_errors() {
    let (c, s) = setup().await;
    let tools = c.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), 105);
    let tool = tools
        .iter()
        .find(|t| t.name == "calculate_dilution")
        .unwrap();
    assert_eq!(tool.input_schema["additionalProperties"], false);
    let before = call(&c, "list_documents", json!({})).await;
    let args = dilution();
    let exact = call(&c, "calculate_dilution", args.clone()).await;
    assert_eq!(exact["engine_version"], "protocol-units-1");
    assert_eq!(exact["exact"]["stock_volume"], q("nL", "200", "1"));
    assert_eq!(exact["exact"]["diluent_volume"], q("nL", "19800", "1"));
    assert!(exact["quantized"].is_null());
    assert_eq!(exact["achieved_concentration"], q("pmol/nL", "1", "10000"));
    assert_eq!(exact["target_reached"], true);
    let mut rounded = args.clone();
    rounded["quantization"] =
        json!({"increment":{"value":"0.3","unit":"uL"},"mode":"NearestTiesUp"});
    let rounded = call(&c, "calculate_dilution", rounded).await;
    assert_eq!(
        rounded["quantized"]["stock"]["rounded"],
        q("nL", "300", "1")
    );
    assert_eq!(
        rounded["quantized"]["stock"]["difference"],
        q("nL", "100", "1")
    );
    assert_eq!(rounded["quantized"]["stock"]["direction"], "Up");
    assert_eq!(
        rounded["quantized"]["diluent_volume"],
        q("nL", "19700", "1")
    );
    assert_eq!(
        rounded["achieved_concentration"],
        q("pmol/nL", "3", "20000")
    );
    assert_eq!(rounded["target_reached"], false);
    let floored=call(&c,"calculate_dilution",json!({"stock":{"kind":"molar","value":"1000","unit":"nM"},"target":{"kind":"molar","value":"1","unit":"nM"},"final_volume":{"value":"1","unit":"uL"},"quantization":{"increment":{"value":"10","unit":"nL"},"mode":"Floor"}})).await;
    assert_eq!(floored["quantized"]["stock"]["rounded"], q("nL", "0", "1"));
    assert_eq!(floored["achieved_concentration"], q("pmol/nL", "0", "1"));
    assert_eq!(floored["quantized"]["diluent_volume"], q("nL", "1000", "1"));
    assert_eq!(floored["target_reached"], false);
    let mass=call(&c,"calculate_dilution",json!({"stock":{"kind":"mass","value":"50","mass_unit":"ng","volume_unit":"uL"},"target":{"kind":"mass","value":"10","mass_unit":"ng","volume_unit":"uL"},"final_volume":{"value":"20","unit":"uL"}})).await;
    assert_eq!(mass["exact"]["stock_volume"], q("nL", "4000", "1"));
    assert_eq!(mass["achieved_concentration"], q("ng/nL", "1", "100"));
    for (path, value, code) in [
        (
            "target",
            json!({"kind":"molar","value":"11","unit":"uM"}),
            Some("ImpossibleDilution"),
        ),
        ("target", Value::Null, Some("MissingInput")),
        (
            "target",
            json!({"kind":"mass","value":"1","mass_unit":"ng","volume_unit":"uL"}),
            Some("UnsupportedUnit"),
        ),
        (
            "stock",
            json!({"kind":"molar","value":"1e3","unit":"uM"}),
            Some("InvalidNumber"),
        ),
        (
            "stock",
            json!({"kind":"molar","value":"-1","unit":"uM"}),
            Some("NegativeQuantity"),
        ),
        (
            "final_volume",
            json!({"value":"18446744073709551615","unit":"uL"}),
            Some("Overflow"),
        ),
        ("final_volume", json!({"value":"20","unit":"ng"}), None),
        (
            "quantization",
            json!({"increment":{"value":"0.3","unit":"uL"},"mode":"ExactOnly"}),
            Some("RoundingRequired"),
        ),
        (
            "quantization",
            json!({"increment":{"value":"0","unit":"uL"},"mode":"Floor"}),
            Some("ZeroIncrement"),
        ),
        ("trusted_stock_inventory", json!(true), None),
    ] {
        let mut bad = args.clone();
        bad[path] = value;
        reject(&c, "calculate_dilution", bad, code).await;
    }
    let missing = reject(&c, "calculate_dilution", json!({}), Some("MissingInput")).await;
    assert_eq!(missing["field"], "StockConcentration");
    assert!(missing["component_index"].is_null());
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    c.cancel().await.unwrap();
    s.await.unwrap();
}
#[tokio::test]
async fn named_mixture_exact_balance_missing_component_and_dimension_guard() {
    let (c, s) = setup().await;
    let args = json!({"final_volume":{"value":"20","unit":"uL"},"components":[{"role":"buffer","volume":{"value":"2","unit":"uL"}},{"role":"enzyme","volume":{"value":"0.5","unit":"uL"}},{"role":"DNA","volume":{"value":"300","unit":"nL"}}]});
    let r = call(&c, "balance_mixture", args.clone()).await;
    assert_eq!(r["plan"]["specified_volume"], q("nL", "2800", "1"));
    assert_eq!(r["plan"]["remaining_diluent"], q("nL", "17200", "1"));
    assert_eq!(r["plan"]["recipe"]["components"][2]["role"], "DNA");
    let mut over = args.clone();
    over["final_volume"]["value"] = json!("2");
    reject(&c, "balance_mixture", over, Some("OverfilledMixture")).await;
    let mut missing = args.clone();
    missing["components"][1]["volume"] = Value::Null;
    let e = reject(&c, "balance_mixture", missing, Some("MissingInput")).await;
    assert_eq!(e["field"], "ComponentVolume");
    assert_eq!(e["component_index"], 1);
    let mut dup = args.clone();
    dup["components"][1]["role"] = json!("buffer");
    reject(&c, "balance_mixture", dup, Some("DuplicateRole")).await;
    let mut invalid = args.clone();
    invalid["components"][0]["volume"]["unit"] = json!("pmol");
    reject(&c, "balance_mixture", invalid, None).await;
    let mut many = args;
    many["components"] = json!(
        (0..129)
            .map(|i| json!({"role":format!("c{i}"),"volume":{"value":"0","unit":"nL"}}))
            .collect::<Vec<_>>()
    );
    reject(&c, "balance_mixture", many, Some("TooManyComponents")).await;
    c.cancel().await.unwrap();
    s.await.unwrap();
}
fn dna() -> Value {
    json!({"model":{"model_id":"explicit650","grams_per_mole_per_base_pair":"650","assumptions":"Declared length-only duplex model; no end or modification correction.","source":"Synthetic arithmetic fixture"},"anchor":{"role":"vector","length_bp":"3000","mass":{"value":"50","unit":"ng"},"stock":{"value":"25","mass_unit":"ng","volume_unit":"uL"}},"targets":[{"role":"insert","length_bp":"1000","molar_ratio":"3","stock":{"value":"10","mass_unit":"ng","volume_unit":"uL"}}],"require_stock_volumes":true})
}
#[tokio::test]
async fn dna_ratios_keep_explicit_mass_model_molecule_counts_and_missing_stock_policy() {
    let (c, s) = setup().await;
    let args = dna();
    let r = call(&c, "calculate_dna_molar_ratios", args.clone()).await;
    //50ng*1000/(3000bp*650g/mol/bp)=1/39pmol; threefold1000bp target needs50ng.
    assert_eq!(r["plan"]["anchor"]["amount"], q("pmol", "1", "39"));
    assert_eq!(r["plan"]["anchor"]["volume"], q("nL", "2000", "1"));
    assert_eq!(r["plan"]["targets"][0]["mass"], q("ng", "50", "1"));
    assert_eq!(r["plan"]["targets"][0]["amount"], q("pmol", "1", "13"));
    assert_eq!(r["plan"]["targets"][0]["volume"], q("nL", "5000", "1"));
    assert_eq!(r["plan"]["recipe"]["model"]["model_id"], "explicit650");
    let mut other = args.clone();
    other["model"]["grams_per_mole_per_base_pair"] = json!("500");
    let other = call(&c, "calculate_dna_molar_ratios", other).await;
    assert_eq!(other["plan"]["anchor"]["amount"], q("pmol", "1", "30"));
    assert_eq!(other["plan"]["targets"][0]["mass"], q("ng", "50", "1"));
    let mut no_model = args.clone();
    no_model.as_object_mut().unwrap().remove("model");
    reject(&c, "calculate_dna_molar_ratios", no_model, None).await;
    let mut bad_model = args.clone();
    bad_model["model"]["grams_per_mole_per_base_pair"] = json!("0");
    reject(
        &c,
        "calculate_dna_molar_ratios",
        bad_model,
        Some("InvalidDnaMassModel"),
    )
    .await;
    let mut missing = args.clone();
    missing["targets"][0]["stock"] = Value::Null;
    let e = reject(
        &c,
        "calculate_dna_molar_ratios",
        missing.clone(),
        Some("MissingInput"),
    )
    .await;
    assert_eq!(e["field"], "StockConcentration");
    assert_eq!(e["component_index"], 1);
    missing["require_stock_volumes"] = json!(false);
    let missing = call(&c, "calculate_dna_molar_ratios", missing).await;
    assert!(missing["plan"]["targets"][0]["volume"].is_null());
    for (key, value, code) in [
        ("length_bp", json!("0"), "ZeroDnaLength"),
        ("molar_ratio", json!("0"), "ZeroMolarRatio"),
        ("length_bp", json!("1000.5"), "InvalidNumber"),
    ] {
        let mut bad = args.clone();
        bad["targets"][0][key] = value;
        reject(&c, "calculate_dna_molar_ratios", bad, Some(code)).await;
    }
    let mut dimensional = args;
    dimensional["anchor"]["mass"]["unit"] = json!("uL");
    reject(&c, "calculate_dna_molar_ratios", dimensional, None).await;
    c.cancel().await.unwrap();
    s.await.unwrap();
}
