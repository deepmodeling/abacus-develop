# ESolver Module Rules

`source_esolver` is the cockpit dashboard of ABACUS: it presents the high-level
algorithm workflow, while the low-level implementation details stay "under the
hood" in their owning modules.

1. `esolver` functions only orchestrate the algorithm flow within one ionic
   step (e.g., the SCF iteration); the outer ionic loop lives in the calling
   driver (`source_relax`, `source_md`). Any computation of a concrete
   physical quantity must be implemented in the owning module (e.g.,
   `source_estate`, `source_hamilt`) and called from `esolver` as a named
   step.
2. Do not add ad-hoc function names to `esolver`. New steps must reuse the
   naming conventions already shared across solvers, so that all `esolver`
   implementations read uniformly. If no suitable name exists, propose the
   convention change first.
