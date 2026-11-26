# SciLifeLab OpenLLM planning integration assessment

Assessed 2026-11-26. Feasible as an optional Rust HTTP planning provider; account
entitlement, current model capabilities and data-processing terms still need
confirmation before deployment. No credential values were read, no model/API
requests made, no biological/workspace data transmitted, and no external source
files changed. This assessment does not report inference tests or cloning accuracy.

## What the scraper actually uses

Read-only evidence from `/home/mahogny/github/claude/btyper_scraper`:

- `docs/scilifelab_llm_usage.md:3–28`: OpenLLM pilot, default base URL
  `https://open-llm.scilifelab.se/api`; documents
  `https://openllm.scilifelab.se/api` as an alias. Pilot Users membership and a
  personal API key are required. Four-week key expiry and shared-service/no-SLA
  limits are repository statements, not independently verified service promises.
- `crates/llm-extract/src/scilifelab.rs`: native `reqwest` client, Bearer auth,
  `GET /models`, `POST /chat/completions`, 90-second whole-request timeout.
  Config names are `SCILIFELAB_API_KEY`, `SCILIFELAB_BASE_URL`; the documented
  fallback credential filename exists as a convention but was **not opened**.
- `crates/llm-extract/src/command.rs`: explicit model wins, then
  `SCILIFELAB_MODEL`, then `gemma3-27b`. Local usage notes also name
  `Qwen3.6-35B-A3B-FP8`. Neither is a live model-discovery result from this review.
- Request DTO supports only model, text role/content messages, max_tokens and
  optional temperature. It contains no tools, tool_choice, tool_calls,
  response_format or json_schema. Extraction runner requests compact JSON in a
  prompt and parses `choices[0].message.content`. Successful extraction is not
  evidence of server-enforced structured output or MCP support.

No applicable AGENTS.md was found in the inspected external repository or its
ancestor directories.

## Authoritative public evidence and its limits

SciLifeLab's [deployment repository](https://github.com/ScilifelabDataCentre/openllm-setup)
describes an Open WebUI front end with vLLM backends, including Swedish KTH and
SafeSpring deployments. Its [C3SE deployment guide](https://github.com/ScilifelabDataCentre/openllm-setup/blob/main/vllm-03-c3se-podman/README.md)
also describes a LiteLLM proxy, scoped model keys and per-consumer limits.
These are deployment descriptions, not proof of this account's current routing.

Configured examples are [Qwen3-32B](https://github.com/ScilifelabDataCentre/openllm-setup/blob/main/vllm-03-c3se-podman/single-node/models/qwen3-32b.conf)
with 32,768 tokens and automatic tool choice/Hermes parsing, and
[Qwen3-235B-A22B](https://github.com/ScilifelabDataCentre/openllm-setup/blob/main/vllm-03-c3se-podman/multi-node/conf/qwen3-235b.conf)
with 131,072 tokens and tool parsing. Those configured context sizes are **not**
verified limits for the scraper's Gemma/Qwen3.6 names, nor guaranteed usable input
budgets after generation/reasoning/tool schemas. No exact context/output limit
for those names was established.

SciLifeLab's [MCP workshop](https://github.com/ScilifelabDataCentre/scilifelab-ai-agent-mcp-workshop-2026-11-02/blob/main/Section_2_MCP/README.md)
connects OpenLLM-backed agents to a separate MCP literature-search server.
This supports the architecture of client-side tool orchestration; it does not
establish that the inference endpoint itself accepts MCP requests or exposes
hosted MCP tools. Strict JSON-schema enforcement remains unconfirmed.

The general [Serve service](https://serve.scilifelab.se/about/) describes free
access for Swedish-affiliated life-science researchers and collaborators.
However, the [roadmap](https://serve.scilifelab.se/about/roadmap/) says its planned
application-developer LLM API requires apps to be available through Serve.
Do not assume pilot access authorizes arbitrary public-web/Tauri redistribution,
a shared bundled key, or unrestricted multi-user traffic. Ask the service team
which policy applies to this specific pilot and deployment.

## Minimal integration proposal

1. Add an optional Rust provider adapter behind a provider-neutral planning
   interface. Tauri and a trusted server hold credentials; browser/WASM receives
   no shared service secret. Keep credentials outside archives, catalog metadata,
   prompts and logs. Offline editing and deterministic analysis remain independent.
   Preserve local binary IPC; translate to the provider's HTTP/JSON only at this
   boundary. No Python, Docker or subprocess is needed in Collagene.
2. Bind an explicit planning session to selected immutable native input/event
   pins and an independently declared target/accepted intent. Send only the
   necessary selected sequences, feature summaries and constraints. Do not upload
   entire workspaces, private history, inventories, credentials or unrelated files.
3. Ask for the existing typed operation DAG/intent or bounded tool calls referencing
   those pins. A Rust orchestrator may expose a small read/preview/design tool
   allowlist through the existing in-process MCP/session APIs. It controls the
   round/output/time budget; the model cannot choose filesystem paths, arbitrary
   HTTP destinations, database writes, publication or source/material authority.
4. Parse strict bounded DTOs, reject unknown fields and invalid IDs/coordinates,
   and locally resolve all small handles. Run normal native preflight, compiler
   and actual-source verifier against the unchanged independent target. Present
   the proposed DAG and scientific findings; preserve Incomplete/Failed states.
   Model prose/confidence, plausible JSON or a tool-call response never confers
   verified material status. Saving/adopting a design remains the ordinary native
   workflow, with source/epoch checks and real persistence acknowledgements.
5. First support plain JSON proposals with strict local validation. If authorized
   synthetic capability tests later confirm tool calling/JSON schema, enable those
   per model; do not make correctness depend on them. Record provider/model ID,
   settings and verifier version for reproducibility. Temperature zero reduces
   sampling variability but is not a determinism or correctness guarantee.

Implement fresh boundary code rather than copying the scraper client unchanged:
validate configured HTTPS origin, prevent redirect/credential forwarding, redact
and mark sensitive headers, cap response bytes before JSON allocation, handle
truncated output/finish_reason, and bound retries/tool loops. The existing client
alone does not establish these protections.

## Remaining decisions and unknowns

- Current authorized model IDs, per-model tool/structured-output support and exact
  input/output/reasoning token limits; discover only through an authorized later
  capability check, not assumptions from public deployment files.
- Pilot enrollment, key renewal and model-scoped permissions; desktop/public-web
  use, end-user keys versus a separately authorized server credential.
- Rate/concurrency quotas, billing/budgets, latency, availability and supported
  cancellation/streaming semantics. General Serve being free is not a quota grant.
- Prompt/output logging, retention, training use, subprocessors and research/IP
  data policy. Swedish hosting does not by itself resolve confidential sequence
  or sensitive-data permissions; the scraper explicitly advises against sensitive
  and patient data.
- Cloning performance: use synthetic/public benchmark plans with independent
  intended products, negative cases and native verifier outcomes after permission
  to call the service. No such evaluation was conducted here.
