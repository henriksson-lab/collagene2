# Editing and applying reaction quantities

Open a saved quantity report, then choose **Revise saved quantities as new report**.
Its exact declarations return to the form. Review the edited values and choose
**Save as new quantity report**. The original attachment and biological document
heads remain unchanged. Export controls continue to refer to the explicitly
opened or saved report.

The form preserves IDs, labels, optional material references, missing values,
rounding declarations, concentration dimensions, reaction totals and each
component's own dilution final volume. Canonical nL, ng/nL and pmol/nL fields
prevent overflow during display conversion. Numeric fields accept decimals and
exact unsigned fractions such as `1/3`; no floating-point rounding is introduced.
Missing or inconsistent saved declarations remain missing or inconsistent.

Dilution final volume is independently editable. **Copy declared reaction total**
is an explicit action; opening a report never silently supplies it. Native
arithmetic still requires it to match the reaction total for a complete
calculation. Direct component volume and dilution final volume have separate raw
fields, so switching amount modes does not reinterpret unfinished inactive text.

Under **Apply quantities to selected reactions**, choose one source reaction and
explicit destination rows. Applying copies the amount, concentration, rounding
and fill declarations and replaces those destination declarations. The source
and unselected reactions remain unchanged. Destination reaction IDs remain
stable; copied component/fill IDs are newly allocated without colliding with
retained IDs.

Every selected destination must reference a distinct operation of the same native
type with the same operand layout. Planned material references map by operation
operand position. Repeated source operands that imply different destination
references are ambiguous and refuse the entire action. No DNA-content matching,
new source choice or chemistry inference is performed. Any incompatible or
missing selection leaves the complete form unchanged.

Copying is only form editing. Exact plan/source/target bindings and workspace
epoch guards remain in force; prior quantity review is cleared. Saving still
requires fresh native quantity review and the existing atomic installation
acknowledgment. These numeric declarations do not establish measured
concentrations, inventory, molecular consistency or experimental readiness.
Molecular verification remains separate, including in combined HTML/PDF exports.

This changes only the web form and review controls. Native quantity records,
worker messages and MCP requests retain their single existing representation;
no compatibility reader or extra installation mechanism is added.

Six focused model cases pass, covering exact restoration, missing values,
immutable capture/open, selected operand remapping and atomic refusal. The
integrated web/worker WASM check passes. Full graphical acceptance remains
deferred under the implementation workflow.
