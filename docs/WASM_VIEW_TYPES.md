# Bounding dialog view types

The desktop App builds many dialog branches inside one `Show`. Leptos0.8.21's
`Show` takes `TypedChildrenFn<C>` and returns a branch containing `C`; the concrete
child view type therefore propagates into the parent render/state types. One
`into_any()` conversion now sits around the existing modal backdrop **inside**
that `Show`. The dialog remains lazy, with its existing markup, callbacks, source
bindings, transaction guards and cleanup ownership. This is a framework type
boundary, not a new HTML wrapper or another dialog workflow.

Only this boundary changes in application code. An exact source comparison
against6d78c69 confirmed that the dialog body and all other App code are identical.
Scientific operations, worker/protocol/archive schemas and other view types are
unchanged. The type-erasure cost is limited to this dialog boundary; inner canvas
painting and scientific algorithms remain concrete. Framework type erasure adds
its usual stored value and build/rebuild dispatch, so it is not applied broadly.

## Measured output

The full existing offline development build passed, using the same lean target,
`CARGO_INCREMENTAL=0`, debug settings and two-job build limit. No extra compiler
profile, browser check or scientific-test pass was introduced for this reversible
rendering refactor.

| Artifact | Before (6d78c69) | After | Reduction |
| --- | ---: | ---: | ---: |
| Generated frontend WASM | 662,440,140 bytes | 651,694,648 bytes | 10,745,492 bytes |
| Optional WASM `name` section | 609,017,735 bytes | 598,281,602 bytes | 10,736,133 bytes |
| Shipped frontend after existing name stripping | 53,422,405 bytes | 53,413,046 bytes | 9,359 bytes |

The existing web-assets hook removes only the optional non-executable `name`
section. This confirms a modest reduction in generated type-name data; it does
not establish a meaningful runtime loading/rendering improvement. The scientific
worker hash remains69df84345300ee73de69f9d261d90866f49008ab10e2d7394dedd1df68359aa8.

The candidate Cargo build reported6m26s. A two-second process sampler observed
about411.8 CPU seconds and a maximum observed5,114.49MiB RSS for the web compiler.
The previous9m32s build also rebuilt browser bindings, and its complete CPU/peak
memory baseline was not captured. No causal compile-time or memory speedup is
claimed, and no repeated build was run solely to manufacture a favorable result.

The small boundary is retained because it demonstrably limits parent type-name
propagation while preserving the existing source/view architecture. The wider
application and its remaining cloning/desktop milestones are still incomplete.
Visual acceptance remains deferred under the user's instruction.
