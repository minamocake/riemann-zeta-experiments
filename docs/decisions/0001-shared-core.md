# 0001. One shared core for all experiments

**Status:** accepted (v0.1)

**Context.** The audited code duplicated kernel and validation logic in two modules, and both
defects that produced invalid bounds were missing input validation (partition coverage).

**Decision.** Kernels, the Montgomery functional, partition validation, certificates and status
logic live in `src/riemann_research/`. Experiments import it and may not reimplement validation.

**Consequences.** A validation fix reaches every experiment at once. New test-function families
should be added as new kernel modules with exact counterparts so that certificates stay possible.
