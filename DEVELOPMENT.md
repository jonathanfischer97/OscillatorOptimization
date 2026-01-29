# Development Status & Progress Tracking

## 📋 Current Status (January 28, 2026)

### Package State: ✅ FULLY FUNCTIONAL
- Package environment is fixed and all core functionality works
- **All tests pass**: 29/29 (100%)
- **Julia version**: 1.10 ONLY (1.11+ require source code changes for API compatibility)
- **Multi-threading**: Configured (6 threads via `~/.claude/env-setup.sh`)

### Julia Version Compatibility
**Currently supports Julia 1.10 ONLY.** Julia 1.11+ are not yet supported due to breaking API changes in:
- `SymbolicIndexingInterface.jl` (parameter/species access APIs changed)
- `DiffEqCallbacks.jl` (SavingCallback API changes)

See the Julia 1.11+ Migration roadmap section below for details.

---

## ⚠️ Known Issues

### Deprecation Warnings (Cosmetic - Do Not Affect Functionality)
- `Product(v)` → `product_distribution(v)` in population_generation.jl:70
- `sys.name` → `get_name()` and `sys.var_to_name` → `get_var_to_name()`

These warnings are from ModelingToolkit/Catalyst deprecations and do not affect package functionality. They will be addressed during the Julia 1.11+ migration.

---

## 📊 Test Results Summary

### Overall: 29/29 tests passing (100%) ✅

```
Test Summary:
✅ Trait mapping & constructor:             8/8  passed  (100%)
✅ ODE solver test:                         1/1  passed  (100%)
✅ Fitness function test:                   6/6  passed  (100%)
✅ Optimization test:                       3/3  passed  (100%)
✅ Version comparison tests:                9/9  passed  (100%)
   ✅ Individual ODE solution:              3/3  passed
   ✅ Fitness calculation reproducibility:  3/3  passed
   ✅ Oscillation detection:                4/4  passed
   ✅ Medium population test:               1/1  passed
```

### Detailed Test Coverage

**✅ ALL Tests Passing:**
- **Trait system**: Observable selection works for both `fullrn` and `trimer_rn` models
- **ODE solving**: Solves correctly with Rodas5P (verified peak position)
- **Fitness calculation**: Produces EXACT expected values:
  - Fitness = 0.520 (within 1% tolerance)
  - Period = 49.2 (within 1% tolerance)
  - Amplitude = 0.630 (within 1% tolerance)
- **Optimization execution**: Runs successfully and finds 643 solutions as expected
- **Version comparison**: All ODE solution, fitness, and oscillation detection tests pass
- **Medium population test**: Runs successfully without errors

---

## ✅ Resolved Issues Log

### January 28, 2026 - Test Fixes
1. **Exported missing functions** - Added `check_oscillatory` and `check_oscillation_regularity` to exports
2. **Fixed ReturnCode comparison** - Changed `sol.retcode == :Success` to `successful_retcode(sol)` to handle modern OrdinaryDiffEq enum return codes
3. **Fixed FFTW plan stride issue** - Removed pre-computed plan from `getFrequencies()` call in test to avoid stride mismatches

### December 19, 2025 - Environment & Threading
1. **Fixed package environment** (THE BIG ONE) - Copied working Manifest.toml from GeometricallyTunableOscillator to lock all dependencies including transitive ones (NonlinearSolve v3.15.2)
2. **Added exact version constraints** - Set exact compat constraints in Project.toml `[compat]` for all direct SciML dependencies
3. **Configured multi-threading** - Created `~/.claude/env-setup.sh` with `JULIA_NUM_THREADS=6`

### December 17-19, 2025 - Code Fixes
1. **Fixed duplicate `DF` parameter** - Removed duplicate declaration in trimer_model.jl line 38
2. **Added missing Symbolics import** - Added `using Symbolics` to main module to fix `UndefVarError`

---

## 📦 Dependency Management

### The Transitive Dependency Problem

**Key Discovery**: Julia's `[compat]` system can only constrain packages in your `[deps]` section. You **cannot** constrain indirect (transitive) dependencies.

**Example:**
- We constrained `ModelingToolkit = "=9.41.0"` in `[compat]`
- But ModelingToolkit depends on `NonlinearSolve` (transitive dependency)
- NonlinearSolve resolved to v3.14.0 (incompatible) because we can't add it to `[compat]`
- No way to fix this with `[compat]` alone!

**Solution**: Use **both**:
1. **Exact compat constraints** on direct dependencies (Project.toml `[compat]`) - guide the resolver
2. **Manifest.toml** to lock ALL dependencies including transitive ones - **this is where the actual version locks live**

This two-part approach is acceptable and recommended for:
- Development packages
- Packages not yet published to General registry
- Research code requiring exact reproducibility

### Critical Package Versions

| Package | Working (GTO) | Previous (Broken) | Current (Fixed) | Status |
|---------|---------------|-------------------|-----------------|--------|
| Catalyst | 14.4.1 | 15.0.11 | **14.4.1** | ✅ |
| ModelingToolkit | 9.41.0 | 9.59.0 | **9.41.0** | ✅ |
| Symbolics | 6.44.0 | 6.58.0 | **6.44.0** | ✅ |
| OrdinaryDiffEq | 6.92.0 | 6.92.0 | **6.92.0** | ✅ |
| SciMLBase | 2.77.2 | 2.77.2 | **2.77.2** | ✅ |
| DiffEqCallbacks | 3.9.1 | 3.9.1 | **3.9.1** | ✅ |
| ADTypes | 1.15.0 | 1.20.0 | **1.15.0** | ✅ |
| SymbolicIndexingInterface | 0.3.37 | 0.3.46 | **0.3.37** | ✅ |
| **NonlinearSolve** | **3.15.2** | **3.14.0** | **3.15.2** | ✅ |

**Key Insight**: NonlinearSolve (transitive dependency) was the culprit. It resolved to v3.14.0 which was missing `AbstractNonlinearTerminationMode`, causing compilation failure. Manifest.toml now locks it to v3.15.2.

### Publishing Strategy

When ready to publish to General registry, consider:

1. **Option A: Keep exact constraints** (simplest, most reliable, but inflexible)
   - Users must accept the version constraints
   - Best for specialized research code

2. **Option B: Relax to minor version ranges** (more flexible but higher maintenance)
   - Test extensively against multiple SciML package versions
   - More work to maintain compatibility

3. **Option C: Vendor critical functionality** (most flexible but highest maintenance)
   - Copy essential code from SciML packages to avoid dependency
   - Maximum flexibility but significant maintenance burden

For now, **Option A with Manifest.toml** provides the stability needed for development and internal use.

---

## 🎯 Roadmaps

### Julia 1.11+ Migration

**Status**: Not started (planned)

**Goal**: Migrate from Julia 1.10 LTS with pinned SciML packages to Julia 1.11+ with current package versions, while maintaining identical numerical behavior through Test-Driven Development.

**Known Breaking APIs**:
- `SymbolicIndexingInterface.jl`: `getu`, `setp`, `setu` usage throughout codebase
- `DiffEqCallbacks.jl`: SavingCallback API changes in `evaluate_individual.jl`
- `Catalyst.jl`/`ModelingToolkit.jl`: Observable definitions and model syntax changes

**5-Phase TDD Plan**:
1. **Establish Numerical Baseline** - Capture exact numerical outputs on Julia 1.10 + pinned packages as reference
2. **Julia 1.11 Migration** - Install identical pinned package versions on Julia 1.11 and identify numerical differences
3. **SciML Package Migration** - Update DiffEqCallbacks (3.9.1 → current), SymbolicIndexingInterface (0.3.37 → current), ModelingToolkit (9.41.0 → current) one at a time
4. **Catalyst.jl Migration** - Update model definition syntax and observable definitions to current Catalyst version
5. **Integration Testing** - Run full optimization suite and validate numerical equivalence with baseline

**Success Criteria**:
- All test cases pass with `rtol=1e-12` or better (numerical equivalence)
- Runtime within 10% of baseline performance
- No more pinned package versions; compatible with latest Julia LTS

**Estimated Effort**: 6-9 weeks of focused development

---

### Package Generalization

**Status**: Not started (planned)

**Goal**: Transform OscillatorOptimization from a specialized package for specific biochemical models into a generic Quality Diversity optimizer for **any** Catalyst.jl ReactionSystem with oscillatory behavior.

**5 Hardcoded Assumptions to Address**:
1. **Observable names** - Currently hardcoded `:Amem_old`, `:Amem`, `:TrimerYield` in trait system
2. **Model-specific symbol extraction** - Assumes "LpA" and "A" patterns in `get_amem_symbols()`
3. **Fixed parameter ranges** - Hardcoded KF_RANGE, KR_RANGE, BOUNDS dictionaries in models
4. **Fixed time domain settings** - Hardcoded `tspan=(0.0, 1968.2)`, `dt=0.1`
5. **Hardcoded fitness thresholds** - Assumes `dt=0.1` in `check_oscillatory()` with `last_200_points = end-2000:end`

**Target Generic API**:
```julia
OptimizationReactionSystem(
    rx_sys::ReactionSystem;
    fft_observable::Symbol,              # User specifies
    time_domain_observable::Symbol,      # User specifies
    time_config::TimeConfig = TimeConfig(),
    oscillation_config::OscillationConfig = OscillationConfig(),
    fitness_function::AbstractFitnessFunction = FFTAmplitudeFitness(),
    fixed_params::Dict{Symbol, Float64} = Dict{Symbol, Float64}(),
    constraint_set::ConstraintSet = ConstraintSet([null_constraint()]),
    alg = Rodas5P(autodiff = AutoForwardDiff(chunksize=length(species(rx_sys)))),
    remove_conserved::Bool = false
)
```

**7 Implementation Phases**:
1. **Test Infrastructure** - Set up comprehensive test suite for existing behavior as regression baseline
2. **Observable Parameterization** - Add observable parameters to constructor, remove hardcoded traits
3. **Time Configuration** - Parameterize tspan, dt, tolerances with configurable defaults
4. **Bounds Generalization** - Extract bounds from Catalyst metadata instead of hardcoded dictionaries
5. **Fitness Function Modularity** - Create AbstractFitnessFunction interface, refactor existing logic into FFTAmplitudeFitness
6. **User Experience** - Add convenience constructors, comprehensive examples, tutorial documentation
7. **Advanced Features** - Support custom constraints, alternative fitness functions, multi-objective optimization

**Example Target Usage for Custom Models**:
```julia
# User defines their own oscillatory model
@reaction_network my_oscillator begin
    @parameters k1 = 1.0, [bounds=(0.1, 10.0)]
                k2 = 2.0, [bounds=(1.0, 20.0)]
    @species A(t) = 5.0, [bounds=(1.0, 100.0)]
             B(t) = 2.0, [bounds=(0.1, 50.0)]

    @observables begin
        A_ratio ~ A / (A + B)
        total ~ A + B
    end

    k1, A + B --> 2A
    k2, 2A --> A + B
end

# Configure optimizer for their specific observables
opt_sys = OptimizationReactionSystem(my_oscillator;
    fft_observable=:A_ratio,
    time_domain_observable=:total,
    tspan=(0.0, 100.0),
    dt=0.1
)

results = run_optimization(500, opt_sys)
```

**Success Criteria**:
- Zero breaking changes for existing users (backward compatibility)
- Simple API for new users with custom models
- Comprehensive test coverage (>90%)

---

## 💻 Environment Reference

### System
- **Platform**: macOS (darwin 25.2.0)
- **Julia Version**: 1.10.10 (aarch64)
- **Working Directory**: `/Users/jonathanfischer/.julia/dev/OscillatorOptimization`

### Threading
- **Threads available**: 6 (multi-threaded)
- **Configuration**: `~/.claude/env-setup.sh` with `JULIA_NUM_THREADS=6`
- **Environment variable**: Automatically loaded via `CLAUDE_ENV_FILE` in `~/.zshrc`
- **Impact**: Optimization uses proper parallelization with `:threadprogress`

---

## 📝 Notes for Maintainers and AI Agents

### Critical Rules
- **DO NOT delete Manifest.toml** - It contains the actual version locks for ALL dependencies including transitive ones (e.g., NonlinearSolve v3.15.2)
- **DO NOT run Pkg.resolve()** - It would regenerate Manifest and lose transitive dependency locks
- **DO NOT relax [compat] constraints** in Project.toml - Exact versions are required for SciML ecosystem stability

### Key Concepts
- **`[compat]` constraints** (in Project.toml) - Guide the resolver on what versions to consider
- **Manifest.toml locks** - The actual resolved versions that are used (includes transitive dependencies)
- **Transitive dependencies** - Cannot be constrained via `[compat]`, only via Manifest.toml

### Development Workflow
- Package environment is **FIXED** - package loads and works correctly
- Manifest.toml approach **works** - don't try to remove or regenerate it
- Multi-threading **configured** - 6 threads by default
- All tests **pass** - maintain 29/29 passing before committing changes

### Version Locking Philosophy
Given this package's research use case and SciML ecosystem fragility:
- ✅ Exact compat constraints are correct - guides resolver to compatible versions
- ✅ Manifest.toml is essential - contains the actual locked versions (including transitive deps)
- ⚠️ Trade inflexibility for stability - the right choice for research code

---

*Last Updated: January 28, 2026*
