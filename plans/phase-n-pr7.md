# Phase N PR-7 — `KroneckerComponent` and SPDE-friendly coordinate forms

## Context

Carried forward from [`plans/phase-n.md`](phase-n.md) — PR-7 is the
Phase N stretch item: extend `@lgm` from index-column random-effect
syntax (`f(idx, IID(n))`) to coordinate-column syntax for SPDE and
separable space-time models. Target: `v0.2.2`.

Phase N PRs 1–6 closed in three weeks (2026-04-29 → 2026-05-06). PR-6
landed the Scotland and Tokyo `@lgm` migration story but had to leave
Meuse as "explicit constructor only" because SPDE projectors are
mesh-barycentric, not column-indexing — the very thing PR-7 fixes.

PR-7 is non-trivial enough to warrant a subplan because it forces
four ADRs that were out of scope for PR-1..PR-6: a parser change to
accept tuple-typed first arguments to `f(...)`, a new schema-side
runtime helper for mesh-barycentric projectors, a question about who
owns the mesh reference, and a dependency policy call about whether
LGMFormula needs `INLASPDE` knowledge.

## State of the world (read-only audit)

What ships today, against PR-7's target surface:

| Piece | Today | Gap |
|---|---|---|
| `KroneckerComponent(space, time)` | Lives at [`packages/LatentGaussianModels.jl/src/components/kronecker.jl`](../packages/LatentGaussianModels.jl/src/components/kronecker.jl); concrete struct, ships in v0.2.0 (Phase M PR-5). | `@lgm` cannot construct one — components tuple slots only emit literal `Component(args)` calls or runtime `_wrap_term(...)` helpers (PR-5). |
| `KroneckerMapping(A_space, A_time)` | Lives at [`packages/LatentGaussianModels.jl/src/observation_mapping.jl`](../packages/LatentGaussianModels.jl/src/observation_mapping.jl) (Phase M PR-1). | `@lgm` always builds a single `LinearProjector(A)` via `hcat(blocks…)`; the macro has no path to emit a `KroneckerMapping`. |
| `SPDE2(points, triangles; …)` | [`packages/INLASPDE.jl/src/components/spde2.jl`](../packages/INLASPDE.jl/src/components/spde2.jl). Stores `FEMMatrices` and a `GMRFGraph`; **does not** retain an `INLAMesh` reference. | The macro needs the mesh to build a `MeshProjector` at runtime — currently the mesh is "lost" once `SPDE2` is constructed. |
| `MeshProjector(mesh, loc)` | [`packages/INLASPDE.jl/src/projector.jl:30`](../packages/INLASPDE.jl/src/projector.jl). Takes `INLAMesh` + `n_obs × 2` matrix. | Has no awareness of `SPDE2` — needs a way to bridge `(component, coordinate columns) → projector block`. |
| `_parse_f_term` | [`packages/LGMFormula.jl/src/parse.jl:142`](../packages/LGMFormula.jl/src/parse.jl). Hard-coded `col isa Symbol` check on the first positional argument. | Cannot accept `(east, north)` tuple-coordinate first arg. |
| `_build_term_block` (schema) | [`packages/LGMFormula.jl/src/schema.jl`](../packages/LGMFormula.jl/src/schema.jl). Emits `sparse(1:n_obs, idx_col, 1.0, n_obs, length(comp))` per term. | No path for mesh-barycentric or Kronecker projector blocks. |

Two cross-cutting realities the audit surfaces:

1. **`SPDE2` does not currently retain `INLAMesh`.** It stores
   `points, triangles` only as the input to FEM assembly; the
   `INLAMesh` object (which `MeshProjector` requires) is never
   captured. The current public API for projector construction is
   "user calls `inla_mesh_2d(...)` once, threads the mesh into both
   `SPDE2(mesh, …)` and `MeshProjector(mesh, locations)` separately."

2. **The macro has no precedent for shape-changing the projector.**
   Every PR-1..PR-6 expansion produces a single
   `LinearProjector(A::SparseMatrixCSC)`. A `KroneckerMapping` block
   inside an otherwise-`hcat`-stacked design matrix is not
   representable in today's expansion — it would require either
   a `StackedMapping` row-partition (which is for multi-likelihood,
   not column-blocks) or a redesign of how the macro emits the
   final mapping.

## Design calls (ADRs)

PR-7 requires four new ADRs in [`plans/decisions.md`](decisions.md).
Order by dependency — ADR-036 unlocks the rest.

### ADR-036 — How does the macro acquire the mesh?

The macro needs to call `MeshProjector(mesh, locations)` at runtime,
but `SPDE2` does not retain `INLAMesh`. Three candidates:

(a) **`SPDE2` retains its `INLAMesh`.** Add `mesh::INLAMesh` field
    (or accept `INLAMesh` in the constructor and store both
    `(mesh, fem)`). Backwards-compatible additive change in
    `INLASPDE.jl`. Once `SPDE2` knows its mesh, a runtime helper
    `_build_spde_block(spde, data, coord_cols)` is straightforward.

(b) **The macro extracts mesh from the `f(...)` expression at parse
    time.** Pattern-match `SPDE2(mesh, ...)` in the AST and emit
    `MeshProjector($mesh, ...)` directly. Brittle — only works for
    the literal `SPDE2(...)` constructor call shape; breaks under
    aliasing (`spde = SPDE2(...); @lgm ... f(coords, spde)`).

(c) **A separate `LGMFormula._SpatialTerm(mesh, component)` wrapper.**
    User writes `@lgm ... f(coords, _SpatialTerm(mesh, SPDE2(mesh, …)))`.
    Verbose; mesh redundancy invites bugs.

**Recommendation: (a).** Adding `mesh` to `SPDE2` is the cleanest
seam, and the v0.1.x API already takes `(points, triangles)` which
are the inputs `inla_mesh_2d` produces — the upgrade path is "users
call `inla_mesh_2d(loc, ...)` to get an `INLAMesh`, then pass that to
both `SPDE2(mesh; ...)` and any projector calls." This is the same
pattern R-INLA uses.

This is a contract change in `INLASPDE.jl`, **not** a contract change
in `LatentGaussianModels.jl`'s `AbstractLatentComponent` — `SPDE2` is
the only component that carries a mesh, and the cross-package coupling
stays inside `INLASPDE.jl`.

### ADR-037 — Tuple-coordinate parser semantics

`f((east, north), SPDE2(mesh; ...))` — the first argument is now a
2-tuple of column symbols. Two questions:

1. **What does the tuple mean for non-SPDE components?** The 2-tuple
   has no meaning for `IID(n)` or `AR1(n)`; it's specifically a
   coordinate-pair for mesh-barycentric projection. Reject at parse
   time when the second arg is not "spatial-shaped" — but the parser
   doesn't know component types. Compromise: reject only the obvious
   shape errors at parse time (`(s, t)` with `length > 3` etc.) and
   defer the type-fit check to schema-time when `data` is bound.

2. **2-tuple vs N-tuple.** R-INLA's SPDE always uses 2D coordinates;
   3D SPDE is out of scope (deferred to v0.3+ per
   [`packages/INLASPDE.jl/plans/plan.md:213-219`](../packages/INLASPDE.jl/plans/plan.md)).
   1D SPDE uses a single coordinate column — `f(t, SPDE1D(mesh))`
   already works under the existing parser since `t` is a `Symbol`,
   not a tuple. So PR-7 only needs to accept 2-tuples.

**Recommendation:** Accept `(s_col::Symbol, t_col::Symbol)` 2-tuples
as the first positional argument of `f(...)`; reject other arities
(`(s,)`, `(s, t, u)`) at parse time with a user-readable error. The
type-fit check ("component must accept a coordinate matrix") is a
runtime error from the schema-side helper.

### ADR-038 — `KroneckerComponent`: explicit-only or `group =` overload

R-INLA writes separable space-time as
`f(s, model = "spde", group = t, control.group = list(model = "ar1"))`.
PR-5 already routes `f(t, AR1; group = grp_col)` to a runtime
`Group(AR1, data.grp_col)` — which is **not** the same as Kronecker
(Group is one inner component per group label, with the latent
flattening being block-diagonal; Kronecker is `Q_s ⊗ Q_t`).

Three candidates:

(a) **Explicit only.** User writes
    `f((east, north), KroneckerComponent(SPDE2(mesh), AR1(T))) data=df`.
    Tuple-coordinate first arg, `KroneckerComponent` second arg, time
    column is implicit (it's whatever `AR1(T)`'s indexing column is —
    but wait, the coordinate tuple is only `(east, north)`, so where
    does the time column go?). Doesn't fit cleanly — KroneckerComponent
    needs *three* columns: two spatial coords + one time index.

(b) **Three-tuple coord.** User writes
    `f((east, north, time), KroneckerComponent(SPDE2(mesh), AR1(T)))`.
    Parser accepts 3-tuple when the component is `KroneckerComponent`.
    First two are mesh coordinates; third is time index. Mirrors the
    Cameletti space-time fixture's natural data layout.

(c) **R-INLA-style `group = time` + auto-Kronecker.** User writes
    `f((east, north), SPDE2(mesh); group = time)` and the macro
    auto-wraps as `KroneckerComponent(SPDE2(...), AR1(T))` —
    but then the macro picks the `time` model and inherits R-INLA's
    `control.group = list(model = "ar1")` ambiguity. Hidden behavior.

**Recommendation: (b) for PR-7.** The 3-tuple form is the most
explicit and least surprising, and the parse-time arity check from
ADR-037 generalises cleanly to "≤ 3 columns, with the trailing
column being a time index when the component is
`KroneckerComponent`." Defer (c) — `group = time` as syntax sugar
over `KroneckerComponent` — to a follow-up if user demand surfaces.

### ADR-039 — LGMFormula ↔ INLASPDE dependency

The mesh-barycentric runtime helper lives where? Three placements:

(a) **`LGMFormula` hard-deps `INLASPDE`.** Simplest, but pulls the
    SPDE FEM stack into every `using LGMFormula` import — wrong
    layering.

(b) **`LGMFormula` weakdeps `INLASPDE` via a Julia 1.9 extension.**
    `LGMFormulaINLASPDEExt` adds the spatial-projector helper when
    the user has `using INLASPDE`. Macros from `LGMFormula` cannot
    export new symbols from extensions, but the extension can attach
    methods to existing `LGMFormula` functions — ADR-008 / Julia
    1.9 extension contract.

(c) **`INLASPDE` adds a hook the macro calls via duck typing.**
    Define a method
    `LGMFormula._build_spatial_block(c::AbstractLatentComponent, data, cols)`
    that throws by default, and `INLASPDE` overloads it for `SPDE2`.
    Same effect as (b) but the bookkeeping is on the `INLASPDE` side.

**Recommendation: (b).** Julia extensions are the standard
Anthropic-ecosystem pattern; ADR-008 already establishes that
`LGMFormula` lives in its own package precisely so optional
integrations layer cleanly. The extension `LGMFormulaINLASPDEExt`
adds methods to a `LGMFormula._build_spatial_block` function defined
(but not exported) in `LGMFormula`'s core, and the schema-side
runtime helper dispatches into it.

## PR sequence

PR-7 ships as a single PR if the audit holds, but the audit surfaces
enough cross-package work that a 7a/7b/7c split is the realistic
shape. Default to a split; collapse back to one PR only if 7a turns
out trivial.

### PR-7a — `SPDE2` retains `INLAMesh` (ADR-036 implementation)

Contract change inside `INLASPDE.jl`:
- New `SPDE2(mesh::INLAMesh; α, pc)` constructor that stores
  `mesh::INLAMesh` alongside `fem`, `graph`, `pc`. Keep the existing
  `SPDE2(points, triangles; …)` constructor as a thin wrapper that
  builds an `INLAMesh` internally via `inla_mesh_2d` or accepts the
  raw points + triangles and constructs a minimal `INLAMesh` that
  preserves the original mesh data.
- Bump `INLASPDE.jl` 0.2.0 → 0.3.0 (additive but the field set
  changes; conservative bump).

**Tests**: existing `SPDE2` regression suite continues to pass; new
test `test_spde2_mesh_field.jl` asserts the mesh round-trips through
construction with `MeshProjector(spde.mesh, loc)` agreeing with
`MeshProjector(mesh, loc)` built from the original mesh.

**Risk**: `SPDE2` is the central type for the entire SPDE stack;
adding a field is risky if downstream code (oracle test fixtures,
raster prediction) compares structs by field set. Audit before merging.

### PR-7b — Tuple-coordinate parser + spatial-projector schema (ADR-037, ADR-039)

Inside `LGMFormula.jl`:
- Extend `_parse_f_term` to accept `Expr(:tuple, Symbol[, Symbol[, Symbol]])`
  as the first positional argument. Reject `length > 3`,
  `length == 1` (the bare `Symbol` case is already handled).
- Lower tuple-coordinate f-terms to a NamedTuple with a new
  `coord_cols::Tuple{Symbol, …}` field instead of `col::Symbol`.
- Add `_build_spatial_block(comp, data, coord_cols)` as an
  unexported function; throw `ArgumentError("@lgm: component
  $(typeof(comp)) does not accept coordinate columns; install
  INLASPDE.jl and load it for SPDE support")` from the core.
- Schema-side dispatch: when the random-effect term has
  `coord_cols`, call `_build_spatial_block` instead of the existing
  `sparse(1:n_obs, idx, 1.0, ...)` block.

In a new `ext/LGMFormulaINLASPDEExt.jl`:
- Overload `LGMFormula._build_spatial_block(c::INLASPDE.SPDE2,
  data, (east, north))` to build
  `MeshProjector(c.mesh, hcat(data.east, data.north)).A`.
- Add `LGMFormula = "..."` weakdep + `INLASPDE = "..."` to the
  extension target in `Project.toml`.

**Tests** (`test/regression/test_spatial_term.jl`):
- Roundtrip: Meuse-shape `@lgm` (Gaussian + Intercept + dist + SPDE)
  produces an `_struct_isequal` model to the explicit constructor
  using `MeshProjector(mesh, hcat(east, north))`.
- Macroexpand: tuple-coord f-term emits a `_build_spatial_block(...)`
  call, not a literal `sparse(...)` block.
- Errors: tuple of length 1 / 4+ rejected at parse; non-Symbol tuple
  entries rejected; calling `_build_spatial_block` without
  `INLASPDE` loaded raises a user-readable error pointing at
  installation.

**Compat bumps**: `LGMFormula.jl` 0.2.0 → 0.3.0 (extension target
addition; new public schema function `_build_spatial_block`).

### PR-7c — `KroneckerComponent` 3-tuple + space-time roundtrip (ADR-038)

Extend the parser further to accept 3-tuples when the second arg is
`KroneckerComponent(...)`. The schema-side runtime helper builds a
`KroneckerMapping(A_space, A_time)` from
`MeshProjector(spatial.mesh, hcat(data.s_col, data.t_col))` and
`sparse(1:n_obs, data.time_col, 1.0, n_obs, length(temporal))`.

This is also the first PR where the macro emits a non-`LinearProjector`
mapping — until now every expansion produced
`LatentGaussianModel(family, components, A::SparseMatrixCSC)`
relying on the implicit `LinearProjector(A)` wrap. PR-7c emits
`StackedMapping` when the design matrix mixes a `KroneckerMapping`
block with non-Kronecker blocks, or `KroneckerMapping` directly when
the entire RHS is one Kronecker term.

**Tests**:
- Cameletti-shape roundtrip: `@lgm` produces an `_struct_isequal`
  model to the explicit `LatentGaussianModel(GaussianLikelihood(),
  (Intercept(), KroneckerComponent(SPDE2(mesh), AR1(T))),
  StackedMapping(...))` form.
- Macroexpand structural assertion on the `KroneckerComponent` slot
  + `KroneckerMapping` mapping.
- Pure `KroneckerComponent` with no other terms: mapping is
  `KroneckerMapping(A_space, A_time)` directly, no `StackedMapping`
  wrap.

**Compat bumps**: `LGMFormula.jl` 0.3.0 → 0.4.0 (semantically
visible change — mapping shape can now be Kronecker).

## Test surface (PR-7 close gate)

- 1 new oracle-shape roundtrip per PR-7b (Meuse) and PR-7c
  (Cameletti). No new oracle fixtures — the LGM-core / INLASPDE
  oracles already exist; the macro just needs to produce
  `_struct_isequal` models.
- Macroexpand structural tests for all new term shapes.
- Error-message tests for: tuple arity bounds, non-Symbol tuple
  entries, missing INLASPDE extension, type-mismatch (`f((s, t), IID(n))`).
- Aqua + JET clean across the new code paths.
- Migration guide updated (`docs/src/lgmformula-tutorial.md`): the
  PR-6 "PR-7" stub paragraphs convert into live examples; the Meuse
  vignette gains a "Same model, written with `@lgm`" `@example`
  block matching the Scotland / Tokyo treatment.

## Out of scope for PR-7

- 3D SPDE coordinates (`(x, y, z)` 3-tuples for spatial-only) —
  defer with the rest of 3D SPDE per
  [`packages/INLASPDE.jl/plans/plan.md:213-219`](../packages/INLASPDE.jl/plans/plan.md).
- Non-separable space-time SPDE (Lindgren et al. 2024) — Phase M+1
  per [`plans/phase-m.md`](phase-m.md "Out of scope").
- R-INLA-style `f(s, model = "spde", group = t, control.group = ...)`
  syntax sugar (ADR-038 candidate (c)).
- Coordinate columns from `Rasters.RasterStack` extraction — that's
  `INLASPDERasters.jl` territory.
- `lgmformula` function-form parity for tuple-coordinate forms.
  PR-7c can ship with the function-form supporting only the macro
  surface; full parity is a follow-up.

## Critical files

| Concern | Path |
|---|---|
| Parser (target of PR-7b changes) | [`packages/LGMFormula.jl/src/parse.jl:142`](../packages/LGMFormula.jl/src/parse.jl) |
| Schema (target of PR-7b changes) | [`packages/LGMFormula.jl/src/schema.jl`](../packages/LGMFormula.jl/src/schema.jl) |
| Expansion AST emit | [`packages/LGMFormula.jl/src/expand.jl`](../packages/LGMFormula.jl/src/expand.jl) |
| `SPDE2` (target of PR-7a) | [`packages/INLASPDE.jl/src/components/spde2.jl`](../packages/INLASPDE.jl/src/components/spde2.jl) |
| `MeshProjector` | [`packages/INLASPDE.jl/src/projector.jl:30`](../packages/INLASPDE.jl/src/projector.jl) |
| `KroneckerComponent` | [`packages/LatentGaussianModels.jl/src/components/kronecker.jl`](../packages/LatentGaussianModels.jl/src/components/kronecker.jl) |
| `KroneckerMapping` | [`packages/LatentGaussianModels.jl/src/observation_mapping.jl`](../packages/LatentGaussianModels.jl/src/observation_mapping.jl) |
| Migration guide (PR-6 stub → PR-7c live) | [`docs/src/lgmformula-tutorial.md`](../docs/src/lgmformula-tutorial.md) |
| Meuse vignette (PR-6 stub → PR-7c live) | [`docs/src/vignettes/meuse-spde.md`](../docs/src/vignettes/meuse-spde.md) |
| ADR registry | [`plans/decisions.md`](decisions.md) — append ADR-036, ADR-037, ADR-038, ADR-039 |
| Macro policy | [`plans/macro-policy.md`](macro-policy.md) |
| Phase N PR ledger | [`plans/phase-n.md`](phase-n.md) |

## Estimated cadence

- PR-7a: ~3 days. `SPDE2` field addition + audit of downstream
  consumers (oracle fixtures). Risk lives here, not in the macro work.
- PR-7b: ~4 days. Parser change, extension scaffolding, new
  regression tests, migration-guide live example for Meuse.
- PR-7c: ~3 days. KroneckerComponent path + Cameletti roundtrip;
  hardest part is matching `StackedMapping` vs `KroneckerMapping`
  emit logic to existing PR-1..PR-6 invariants.

Total: ~2 weeks. Phase N's stretch budget per
[`plans/phase-n.md:339-346`](phase-n.md#estimated-cadence) was
5 weeks for PR-1..PR-6 with PR-7 as slack; PR-1..PR-6 closed in
~3 weeks, leaving a 2-week buffer that aligns with this estimate.

If PR-7a's audit surfaces a downstream-consumer issue that needs its
own PR, the cadence stretches; in that case PR-7 defers to Phase N+1
and the stub paragraphs in PR-6 stay as-is.

## Phase close + release target

PR-7 closes Phase N with `v0.2.2`:
- `INLASPDE.jl` 0.3.0 (PR-7a contract change).
- `LGMFormula.jl` 0.4.0 (PR-7b extension target + PR-7c mapping
  shape change).
- Umbrella `INLA.jl` v0.2.2 with a CHANGELOG entry covering all
  three sub-PRs.
- Tag at PR-7 close: `v0.2.2` on `main`, release commit
  `chore(release): v0.2.2 — Phase N close (LGMFormula PR-7)`.

If PR-7 defers to Phase N+1, the v0.2.1 close tag (which lands at
PR-6 plus ADR housekeeping) stands as Phase N's actual close, and
PR-7 graduates to a Phase O scope item alongside any post-Phase-M
follow-ups (fractional-α SPDE per
[ADR-030](decisions.md#adr-030-fractional-α-spde-deferred-to-v021)).
