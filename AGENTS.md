# Agent Instructions

## First Steps (Read on Every Session)
1. **Read this file completely** before making any changes.
2. **Read DEVELOPMENT.md sections** "Current Status" and "Known Issues" for latest package state.
3. **Only read README.md** if you need to understand architecture/API details.

## Project Quick Reference
- **Julia 1.10 ONLY** - Do not attempt to use Julia 1.11+. APIs are incompatible.
- **Threading**: Always ensure `JULIA_NUM_THREADS` is set (default: 6) or pass `--threads=6` to julia commands.
- **BLAS**: Run `BLAS.set_num_threads(1)` before any optimization work to avoid thread contention.
- **Package environment**: Fixed via locked Manifest.toml. Do NOT modify dependency versions.

## Commands
- **Install deps**: `julia --project -e 'using Pkg; Pkg.instantiate()'`
- **Run tests**: `julia --project -e 'using Pkg; Pkg.test()'`
- **Run example**: `julia --project examples/main.jl`
- **Build docs**: `julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'`

## Code Conventions
- **Indentation**: 4 spaces, no tabs. Keep lines reasonably short.
- **Naming**: `snake_case` for functions/variables, `CamelCase` for types/modules.
- **Mutating functions**: End with `!` (e.g., `calculate_fitness!`).
- **Exports**: Keep explicit in `src/OscillatorOptimization.jl` with docstrings for public APIs.
- **RNG**: Use deterministic seeds (`MersenneTwister`/`Xoshiro`) in examples and tests.

## Testing
- **Framework**: stdlib `Test` with `@testset`/`@test`.
- **Location**: Add tests under `test/` and include from `test/runtests.jl`.
- **Pass threshold**: All **29/29 tests must pass** before committing changes.
- **Determinism**: Set seeds; avoid long-running ODE grids.

## Critical Warnings
- **DO NOT delete Manifest.toml** - It locks ALL dependencies including transitive ones (e.g., NonlinearSolve v3.15.2).
- **DO NOT run Pkg.resolve()** - It regenerates Manifest and breaks transitive dependency locks.
- **DO NOT relax [compat] constraints** in Project.toml - Exact versions required for SciML ecosystem stability.

## Commit Style
- **Imperative, concise** commit messages. Scope if useful.
- Examples: `add solver from DataFrame row`, `fix CI: pin SciML deps`, `update row solvers`

## After Making Changes
1. **Run full test suite** and confirm all 29/29 tests pass.
2. **Update DEVELOPMENT.md** "Current Status" date and any affected sections.
3. **Move resolved items** from "Known Issues" to "Resolved Issues Log" with date.
