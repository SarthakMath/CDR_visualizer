# CDR FEM FCT lab

An interactive finite element laboratory on Ω = (0,1)². It solves real P1 triangular Galerkin systems in a browser worker, displays a rotatable solution surface or a 2D field, and lets you change the mesh, source, boundary data, reduced matrix A, and reduced load vector f.

Choose Steady or Transient, then C, D, R, CD, CR, DR, or CDR, and a matching example. C means convection/advection, D means physical diffusion, and R means reaction. Transient problems include playback and time scrubbing; Galerkin/isotropic/SUPG use the selected implicit integrator, while FCT uses its documented SSPRK2 substeps.

No package installation or build step is needed. The only external resource is an optional Google Fonts stylesheet; system fonts work offline. All numerical computation runs locally in the browser. Session edits live in memory and reset on reload. **Export data** downloads the complete experiment, mesh, sparse matrices, and solutions as JSON.


## Organized experiment setup

The sidebar follows the order **Problem → Mesh → Coefficients/forcing → Time/initial data (transient only) → Boundary data → Stabilization**. First select the time setting and physical terms, then choose a matching example. The large workspace shows only the current experiment summary, equation, and solution; the former rows of overlapping example buttons have been removed. The system inspector, source/boundary expressions, matrix edits, exact reference, μ control, mass choices, and exports remain available.

| Terms | Steady starter | Transient starter |
| --- | --- | --- |
| C: convection/advection | Compatible affine transport | Moving front · FCT |
| D: diffusion | Smooth sine field | Diffusing pulse |
| R: reaction | Reaction equilibrium | Pure decay / growth |
| CD | Outflow layer · BJK AFC | Moving and diffusing pulse |
| CR | Affine transport + reaction | Moving and decaying pulse |
| DR | Diffusion–reaction sine field | Spreading and decay |
| CDR | CDR sine field | Moving, spreading, and decay |

Only examples matching the time setting and equation terms appear. Inactive physical coefficients are disabled: choose a family containing R to edit reaction, and one containing D to edit physical diffusion. SUPG is unavailable without a convection term. The method selector identifies **AFC · BJK-type limiter** for steady problems and **FCT · Zalesak-style limiter (SSPRK2)** for transient problems. Existing numerical methods are retained; no additional limiter has been added in this reorganization.

The four added steady starter cases cover C, R, CR, and DR. The C/R/CR examples use `u = 1+x+y`, matching source `vx+vy+c*(1+x+y)`, and matching traces. The DR example uses `u = sin(pi*x)*sin(pi*y)` with source `(2*eps*pi^2+c)*u`. Pure steady convection usually takes inflow boundary data: this laboratory retains the full Dirichlet machinery and supplies compatible traces in its affine starter. Its default n = 9 gives a solvable Galerkin system, while some other meshes (including n = 16) produce a singular baseline. That condition is reported explicitly; the app does not insert physical diffusion silently. Use CD with small positive ε for a robust boundary-layer experiment. Pure steady reaction requires nonzero c.

`dist/experiments.js` contains the example catalogue, family mapping, defaults, and method availability. Five checks in `tests/experiments.test.mjs` cover all 14 starter combinations, physical-term activation, new manufactured equilibria, the singular pure-convection case, and method labels/availability. The full suite has 64 checks. Direct links such as `?example=afc-affine` and `?example=afc-layer` still load their matching time setting and equation family.

## Run and test

From this directory:

```bash
python3 -m http.server 4173 --bind 127.0.0.1 --directory dist
# Open http://127.0.0.1:4173
```

Use a web server: ES modules and Web Workers do not work reliably from a file:// URL. With npm installed, `npm start` runs the same command. Tests need Node.js 18 or newer, with no third-party packages:

```bash
node --test tests/*.test.mjs
# Or: npm test
```

## Time-dependent experiments: moving layers and mass lumping

In **Problem**, choose **Transient**, then the desired equation terms, and an example from the filtered list. Every C/D/R combination has a starter example. Changing the time setting or equation terms loads a matching default; choosing another example loads that example’s documented parameters and comparison method.

Each time-dependent experiment solves

```text
uₜ − εΔu + β·∇u + cu = f,        u(x,y,0) = u₀(x,y),
MII ḃ + AII b = FI − AIB g − MIB ġ.
```

The **consistent baseline** always uses consistent mass and δ = 0. The **comparison** uses the mass matrix selected in the sidebar and the current added diffusion δ. Both use the same initial interior field, boundary values, source, mesh, time step, and time integrator. This lets you isolate mass lumping by keeping δ = 0, or artificial diffusion by choosing consistent mass for the comparison.

### Moving front: convection and dispersion

1. Choose **Moving front**. The defaults are n = 24, β = (0.55,0), physical ε = 0, Δt = 0.005, T = 0.8, Crank–Nicolson, and lumped mass for the comparison.
2. Press **Play**, or drag the time slider. A smooth but steep-sided pulse moves to the right; the exact reference translates its initial profile without spreading. At t = 0.8, its center has moved from x ≈ 0.25 to x ≈ 0.69.
3. Inspect the orange consistent baseline and teal comparison in the cross-section. Both can oscillate. Mass lumping changes dispersion and is **not** sufficient stabilization for central Galerkin convection.
4. At the same time, click **Apply mesh-based value**. This gives δ = 0.55/(2×24) ≈ 0.0114583. The comparison becomes much smoother and also substantially wider and lower. Small undershoots can remain; the 1D-inspired formula is not a general 2D positivity guarantee.
5. Set comparison mass to **Consistent M** to isolate the effect of δK. Set δ back to 0 to recover the baseline exactly. Switch to backward Euler to observe its additional temporal damping; both curves then use backward Euler.

The initial profile is

```text
u₀(x,y) = ½[tanh((x−0.15)/0.015) − tanh((x−0.35)/0.015)]
          × exp(−((y−0.5)/0.18)²).
```

The reference is u₀(x−βx t, y−βy t). Its trace is prescribed on all four edges for this compatible analytic example. This first version does not implement periodic boundaries or a separate inflow-only boundary formulation for arbitrary hyperbolic data.

### Diffusing pulse: see where lumping helps positivity

1. Choose **Diffusing pulse**. Keep n = 16, ε = 0.01, Δt = 0.001, T = 0.1, **backward Euler**, δ = 0, and comparison mass **Row-sum lumped Mᴸ**.
2. Click **Next time step**, then compare **Consistent baseline** and **Comparison**. The former develops negative values outside the narrow positive pulse, while the latter stays nonnegative in this experiment.
3. Scrub near t = 0.005, or play the whole evolution. Across the default run the consistent baseline reaches a minimum of approximately −0.01611; the lumped backward-Euler comparison has no negative nodes below the diagnostic tolerance of −10⁻¹⁰.
4. In the matrix inspector, choose **Consistent mass M**, then **Lumped mass Mᴸ**. The first has off-diagonal entries; the second is diagonal. These assembled matrices are read only; **Spatial operator A** remains editable.
5. Compare the L² errors too. This deliberately narrow initial Gaussian is under-resolved on n = 16. Lumping removes negative values here but does not necessarily decrease error, conserve the exact initial Gaussian's mass, or reproduce its peak better. Refine the mesh to improve initial resolution.

This is a particular positivity demonstration on the supplied right-triangle mesh with a diffusion operator and backward Euler. Mass lumping alone is not a guarantee after arbitrary A edits, convection, changed mesh geometry, or a change to Crank–Nicolson. For the mathematical context, see [Thomée and Wahlbin, “On the existence of maximum principles in parabolic finite element equations”](https://research.chalmers.se/en/publication/79818).

### Moving + diffusing: physical versus artificial spreading

Choose **Moving + diffusing**. The default Gaussian has center (0.25,0.4), initial squared width s² = 0.0049, β = (0.5,0.1), and ε = 0.002. The exact whole-plane Gaussian is

```text
w(t) = s² + 4εt,
u_exact(x,y,t) = [s²/w(t)] exp(−[(x−cx−βx t)²+(y−cy−βy t)²]/w(t)).
```

Its trace is imposed on the square's boundary, so it is also the exact reference for this bounded-domain test. The pure diffusion example uses the same formula with β = 0, center (0.5,0.5), and s² = 0.000625. Change ε to alter physical spreading; add δ to alter the comparison's spatial operator; change the mass choice to see its separate effect. References follow changed ε and velocities while the matching initial/source/boundary expressions remain intact.

### Time controls and diagnostics

- **Play/Pause**, single-step buttons, rewind, and the slider operate on the computed trajectory. Playback speed is measured in frames per second; it does not alter the numerical Δt.
- Changing a parameter pauses playback and recomputes from the initial field. The displayed physical time is preserved as closely as the new time grid permits. Loading a preset returns to t = 0.
- Both surface/field plots and cross-sections share scales locked over the full current trajectory. This prevents shrinking peaks from appearing unchanged through automatic time-by-time rescaling.
- Source and edge expressions accept `t`, `vx`, and `vy` as well as `x`, `y`, `eps`, `c`, and `pi`. `eps` is the physical diffusion (zero in pure convection), and `vx`, `vy` are the active convection velocity (both zero in pure diffusion).
- Initial data are interpolated at nodes. At boundary nodes, prescribed Dirichlet data take precedence over u₀; a warning reports disagreements at t = 0.
- The integral ∫uₕ, negative-node count (threshold −10⁻¹⁰), Courant number |β|Δt/h, and physical diffusion number εΔt/h² are shown for time-dependent runs. The integral can change through sources and prescribed boundaries; it is not labeled a conserved quantity. At t = 0 no step has been solved, so backward error is left blank.
- At most 240 steps are computed. A final shorter step lands exactly at T. The final matrix inspector uses that shorter step where appropriate. Reference and difference views do not report a solve residual.
- JSON export includes the whole trajectory, the current snapshot, consistent and lumped masses, spatial operators, and the actual current step matrices and right-hand sides. Algebraic A/f overrides are held throughout the run; f overrides replace the reduced spatial forcing at each time before the changing-boundary mass contribution.

### Time integration and boundary lifting

For each element T, the P1 consistent mass is

```text
Mᵀ = (|T|/12) [[2,1,1],[1,2,1],[1,1,2]].
```

The lumped mass is diagonal with entries ∑ⱼ Mᵢⱼ = ∫φᵢ. **The sum is taken before boundary elimination**; summing only MII would lose the boundary-neighbor contributions and fail the constant-state patch test.

For θ = 1 (backward Euler) or θ = ½ (Crank–Nicolson), let Aδ = A + δK and let fδ(t) include spatial boundary lifting. The consistent-mass step is

```text
(MII + θΔt Aδ) bⁿ⁺¹ = (MII − (1−θ)Δt Aδ) bⁿ
+ Δt[θ fδ(tⁿ⁺¹) + (1−θ) fδ(tⁿ)]
− MIB (gⁿ⁺¹ − gⁿ).
```

For the lumped comparison, replace MII by the full-row-sum diagonal Mᴸ and omit the MIB term, since the lumped mass has no off-diagonal boundary block. Adding δK changes spatial diffusion; it does not change M. Backward Euler and Crank–Nicolson have distinct damping and accuracy properties, so the UI identifies the time integrator alongside the mass treatment.

The worker reuses partial-pivot LU factors (n≤32) or sparse preconditioners (n>32) for each distinct step size and uses sparse row products for time-step loads. This avoids rebuilding or refactoring a dense system at every frame. There is no explicit Euler option in this version; the displayed Courant and diffusion numbers are diagnostics, not blanket stability or positivity guarantees.

## The steady experiment: a thin layer and artificial diffusion

**Important distinction:** the default problem is Poisson. A thin, under-resolved Poisson solution can be inaccurate, but pure Poisson on this uniform right-triangle mesh does not generally exhibit the convection-driven oscillation mechanism. The **Transport layer** preset explicitly adds a convection term to demonstrate that mechanism honestly. The **Sharp Poisson layer** preset remains a pure Poisson example.

1. Click **Transport layer**. This loads:

   ```text
   n = 16, slope +1 diagonals
   ε = 0.005, βx = 1, βy = 0
   source f(x,y) = 1, multiplier = 1
   artificial diffusion δ = 0
   ```

   The PDE is −εΔu + ∂u/∂x = 1. Its exact solution is independent of y:

   ```text
   q(x) = x − [exp((x−1)/ε) − exp(−1/ε)] / [1 − exp(−1/ε)].
   ```

   Set u = 0 on the left/right, and u = q(x) on the bottom/top. The preset fills these expressions automatically. The trace depends on `eps`, so changing ε updates the prescribed data too.

2. Keep **Original** selected and inspect the cross-section at y = 0.5. The physical layer near x = 1 has width O(ε), whereas h = 1/16 = 0.0625. The coordinate mesh Péclet number Peₕ = |β|h/(2ε) is 6.25. The discrete solution oscillates; its maximum is approximately **1.84579**, although the exact solution stays between 0 and 1.

3. Click **Apply mesh-based value**. This sets:

   ```text
   δ = max(0, |β|/(2n) − ε) = 0.02625.
   ```

   The app adds the assembled diffusion matrix δK. The new coordinate mesh Péclet number is 1. **With diffusion** becomes selected; the cross-section simultaneously shows the original orange curve, stabilized teal curve, and exact dashed curve. The stabilized maximum is approximately **0.93750**. The pronounced oscillations are suppressed, but the layer is spread out.

4. Select **Δ Change** to locate where the solution changed. **Original** and **With diffusion** use the same color and height range, so switching views cannot hide a change by rescaling the plot.

5. Enter **δ = 0.1** to see excess diffusion. Then change **n to 32** and click **Apply mesh-based value** again; the suggested δ becomes **0.010625**. Refinement changes the amount of stabilization needed. The suggestion is not recalculated automatically when the mesh changes, so you can also keep δ fixed for controlled comparisons.

The mesh-based formula is a 1D-inspired starting point for this example, not a positivity guarantee for all 2D velocity directions, diagonal patterns, or manually edited matrices. This is isotropic artificial diffusion, **not SUPG**. Stabilization changes the differential equation and can worsen accuracy even when the picture is smoother.

### Pure Poisson layer

Click **Sharp Poisson layer**. Let a = 0.025 and

```text
qₐ(x) = x − [exp((x−1)/a) − exp(−1/a)] / [1 − exp(−1/a)],
u_exact(x,y) = qₐ(x) sin(πy),
f(x,y) = [exp((x−1)/a)/(a²(1−exp(−1/a))) + π²qₐ(x)] sin(πy).
```

This solves −Δu = f with zero Dirichlet data. Compare n = 4, 8, 16, and 32 at δ = 0. A coarse mesh misses the sharp right-hand layer. Refinement improves its representation. Adding δK with fixed f scales the discrete zero-boundary Poisson solution by 1/(1+δ) when ε = 1: it lowers the amplitude and does **not** repair the missing resolution. For arbitrary physical ε, this preset's exact solution is qₐ(x) sin(πy)/ε.

### Smooth Poisson baseline

The default solves −Δu = 2π² sin(πx) sin(πy), with u = 0 on ∂Ω. The exact solution is sin(πx) sin(πy). Approximate L² errors with the default diagonal pattern and δ = 0 are:

| n | Interior unknowns | L² error (quadrature estimate) |
|---:|---:|---:|
| 4 | 9 | 0.0790906 |
| 8 | 49 | 0.0211337 |
| 16 | 225 | 0.00537749 |
| 32 | 961 | 0.00135044 |

## Matrix, forcing, and boundary experiments

- **Edit A:** choose Smooth Poisson, set n = 4, select Matrix A, row 4, column 4. Change its value from 4 to 8 and click **Apply edit**. The center value falls. The mirror checkbox also writes A[j,i]; turn it off to make a nonsymmetric perturbation. Transport matrices are naturally nonsymmetric, so the transport preset disables mirroring by default.
- **Edit f:** select Load f and DOF 4. Increase the displayed value, then apply it. The center rises. This edits the assembled reduced load, including any original boundary lift; it is different from changing the continuous source function.
- **Edit the continuous source:** enter an expression such as `1`, `100*exp(-100*((x-0.5)^2+(y-0.5)^2))`, or `sin(pi*x)*sin(pi*y)`. Source multiplier scales the function before integration.
- **Edit boundary data:** each edge accepts a number or an expression. For example, source `0` with all four boundary expressions set to `x+y` reproduces an affine solution. Conflicting edge values at a corner are averaged. A boundary node override takes precedence.
- **Inspect nodes:** click the mesh preview or the solution plot. The inspector shows the global node ID, coordinates, and free DOF if interior. Boundary node values can be prescribed individually; interior nodal values are solved unknowns, so edit their A or f instead.
- **Mesh rules:** changing n or the diagonal pattern clears all indexed overrides because their locations/connectivity change. PDE/source/boundary changes retain indexed overrides. **Clear edits**, presets, and **Reset** clear them.

Indices are zero based. Nodes are ordered left to right, bottom to top: node(i,j) = j(n+1)+i. Interior DOF(i,j) = (j−1)(n−1)+(i−1), for 1 ≤ i,j ≤ n−1. The default splits each square along its bottom-left to top-right diagonal (slope +1). Slope −1 and alternating patterns are also available. The range is 2 ≤ n ≤ 128.

## Discretization and stabilization

For steady modes, the supported PDE is −εΔu + β·∇u + cu = f with prescribed Dirichlet data. Poisson mode sets β = 0. The time-dependent extension adds uₜ as described above. For local nodal basis functions φᵢ on a triangle T:

```text
Kᵀᵢⱼ = |T| ∇φᵢ · ∇φⱼ
Cᵀᵢⱼ = (|T|/3) β · ∇φⱼ
Fᵀᵢ  = ∫T f φᵢ
```

Source integrals use a seven-point, degree-five triangle quadrature rule. Diffusion and constant-velocity convection matrices are integrated exactly. If I denotes interior nodes and B boundary nodes, the reduced system is

```text
A = (εK + C + cM)II
f = FI − (εK + C + cM)IB g
Ab = f.
```

Manual A and f overrides are applied to this reduced base system. Stabilization solves

```text
(A + δKII) bδ = f − δKIB g.
```

The extra boundary lifting is essential for nonzero Dirichlet data. Both solutions share the same prescribed boundary values. K is reassembled from the current mesh; manual A edits do not change K.

The solver uses double-precision Gaussian elimination with partial pivoting inside a Web Worker. Singular or numerically unsafe systems show an error and dim the last successful plot; they are not silently regularized. An exact reference is shown only when the current data still match a known preset. Adding artificial diffusion retains the exact reference for the **original** PDE, so the error measures stabilization bias too.

The L² value is a seven-point quadrature estimate of the continuous piecewise-linear field's error, not a nodal RMS and not a certified bound. Thin unresolved layers can also challenge the error quadrature. Backward error is ‖Ab−f‖∞ / (‖A‖∞‖b‖∞ + ‖f‖∞), using the appropriate original/stabilized system. A tiny algebraic residual does not imply an accurate PDE solution. The displayed Pe uses h = 1/n (square spacing), not the triangle diameter √2/n.

## Source and verification

- `dist/fem.js`: mesh, expression parser, element assembly, boundary elimination, direct solve, interpolation, and diagnostics.
- `dist/transient.js`: consistent/lumped mass time stepping, time-dependent boundary lifting, and trajectory generation.
- `dist/worker.js`: asynchronous steady solve and trajectory computation with progress updates.
- `dist/plots.js`: dependency-free Canvas surface/heatmap, mesh, cross-section, and sparsity plot.
- `dist/app.js`: controls, matrix edits, inspectors, and JSON export.
- `tests/fem.test.mjs`: 12 steady-solver checks covering stencil/load, affine patches, boundary lifting, convergence, layer stabilization, edits, parsing, and singular systems.
- `tests/transient.test.mjs`: 10 checks covering exact local masses, reused pivoted LU, changing-boundary affine patches, Crank–Nicolson with quadratic time data, constant-state preservation, heat-pulse positivity, isolated mass comparison, zero-diffusion transport, exact Gaussian data, time grids, and resource limits.

The tests validate these finite cases at floating-point tolerances; they are not a proof of stability for arbitrary input. A useful background example is the [FEniCS convection–diffusion tutorial](https://www.karlin.mff.cuni.cz/~hron/fenics-tutorial/convection_diffusion/doc.html). For the distinction between suppressing oscillations and smearing layers, see the [classical artificial-diffusion discussion in this research article](https://link.springer.com/article/10.1007/s00419-026-03047-y).

## Reaction experiments: decay, growth, and transport

The constant reaction coefficient **c** is available in every mode, including steady problems. The sign convention is **+cu on the left**: with no forcing or transport, uₜ + cu = 0 gives u = exp(−ct)u₀. This is linear reaction; nonlinear kinetics and coupled chemical species are not implemented.

Five buttons in **Add reaction** provide starting points:

| Preset | Physical equation |
| --- | --- |
| Pure decay | uₜ + cu = 0 |
| Spreading + decay | uₜ − εΔu + cu = 0 |
| Moving + decay | uₜ + β·∇u + cu = 0 |
| Moving + spreading + decay | uₜ − εΔu + β·∇u + cu = 0 |
| Steady ADR | −εΔu + β·∇u + cu = f |

1. Choose **Pure decay**, with c = 2, consistent comparison mass, δ = 0, Crank–Nicolson. Play to T = 0.8. The exact amplitude falls by exp(−1.6) ≈ 0.202 without moving or spreading. The fixed color scale makes decay visible. Set c = 0 to preserve the initial field, or c = −1 to see growth by exp(0.8) ≈ 2.226.
2. Choose **Moving + spreading + decay**. The defaults are n = 24, ε = 0.002, β = (0.5, 0.1), c = 2, Δt = 0.005, T = 0.8. Compare c = 0, 2, and −1 at the same time. Motion, spreading, and reaction amplitude changes have separate causes. Add δK to see extra spreading, then switch comparison mass to isolate the effect of lumping.
3. To examine a reacting layer, choose **Moving front**, set c = 2, and multiply every boundary expression by `exp(-c*t)`. This adds decay to the original steep front. Its old exact reference is hidden after changing c. Compare δ = 0 with the mesh-based δ; reduced amplitude is not evidence that oscillations are resolved.
4. To isolate stiff reaction, choose **Pure decay**, set all four boundary expressions to `0`, c = 100, Δt = 0.04, T = 0.4, and keep both masses consistent and δ = 0. Here cΔt = 4. Crank–Nicolson's amplification (1−cΔt/2)/(1+cΔt/2) = −1/3 produces alternating signs. Backward Euler's factor 1/(1+cΔt) = 1/5 stays positive, though the exact factor is exp(−4). Reduce Δt for accuracy. Adding diffusion cannot repair an unresolved reaction time scale. Edited boundaries hide the preset reference.

For the Gaussian presets, the exact reference is

```text
u(x,y,t) = exp(−ct) s²/(s² + 4εt)
           × exp(−((x−0.25−βx t)²+(y−0.4−βy t)²)/(s²+4εt)),
s² = 0.0049.
```

Inactive diffusion/velocity terms are zero. Compatible exact traces are prescribed on all four edges. Expressions for source, initial data, and boundaries accept `c`; changes to coefficients update these traces and references. Other presets lose their reference when reaction invalidates it. The steady ADR preset uses u = sin(πx)sin(πy) and the corresponding coefficient-dependent manufactured source.

The reaction element matrix is **Rᵀ = cMᵀ**, assembled exactly including its Dirichlet boundary contribution. **A includes reaction before manual edits and before δK**. Comparison mass lumping changes only the time derivative, not R: both runs retain consistent cM. Therefore the diffusion-only positivity demonstration must not be generalized to stiff reaction. Negative c models physical growth and may also make a steady operator singular or indefinite.

`tests/reaction.test.mjs` checks reaction assembly and boundary lifting, steady/transient manufactured patches, exact BE/CN amplification for decay and growth, PDE residuals of all four transient analytic profiles, invalid inputs/reference suppression, and steady refinement.

## Layers and SUPG: three reproducible comparisons

Use **Steady → CD → Outflow layer · SUPG**, **Transient → C → Moving front · SUPG**, or **Transient → CDR → Developing outflow layer**. The orange solution is always unstabilized Galerkin with consistent time mass. The teal **Comparison** uses the selected method: none, isotropic δK, or SUPG. These methods are mutually exclusive; an old δ value is ignored in SUPG/none mode. Keep comparison mass consistent to isolate stabilization. Switching methods rebuilds the run and preserves the selected time.

### Steady outflow layer

1. Select **Steady outflow layer**: n = 16, ε = 0.005, β = (1,0), c = 0, f = 1, SUPG strength = 1. Its exact q(x) is the same outflow-layer profile documented above; width is O(ε), much smaller than h = 1/16.
2. Compare the cross-section at y = 0.5. The baseline maximum is about 1.846; SUPG lowers it to about 0.9375. The quadrature L² error falls from 0.1648 to 0.1186. The layer is still unresolved: a low residual or a smooth curve does not recover its sub-element shape.
3. Switch **Comparison method** to **Isotropic artificial diffusion**, then apply the mesh-based value. In this essentially one-dimensional example the results are similar. This example demonstrates oscillation control, not SUPG's crosswind advantage.
4. Switch back to SUPG, try strengths 0, 0.5, 1, and 2, then n = 32. Strength 0 recovers the unmodified comparison. τ is recalculated for the current mesh, velocity, diffusion and reaction.

### Moving interior layer

1. Select **Moving interior layer**: n = 24, ε = 0, β = (0.55,0), c = 0, Δt = 0.005, T = 0.8, Crank–Nicolson. It uses the steep-sided pulse from Moving front, with exact translating edge traces; the comparison uses consistent M plus full transient SUPG.
2. Play to t = 0.8. The baseline min/max are about −0.2322 / 1.1771; SUPG gives −0.0665 / 1.0440. L² error falls from about 0.05657 to 0.04803. SUPG helps but does **not** remove all overshoot or undershoot. The count of negative nodes can rise even as undershoot magnitude falls, so compare extrema and error as well as counts.
3. Switch to isotropic diffusion and apply the mesh-based δ ≈ 0.0114583. Use **2D field** to compare spreading perpendicular to the flow as well as the horizontal cross-section. Isotropic δK spreads across y; SUPG's leading artificial diffusion acts only along β. Neither method restores detail already lost by interpolating a thin front on a coarse mesh.
4. Return to SUPG and change c to 1: the exact front also decays by exp(−t). Try backward Euler separately to see additional temporal damping. With lumped mass selected, the comparison uses Mᴸ + H, not a diagonal time operator; the SUPG H term remains consistent.

### Developing boundary layer

1. Select **Developing boundary layer**: n = 16, ε = 0.005, β = (1,0), c = 0.5, Δt = 0.01, T = 1, Crank–Nicolson, consistent mass and SUPG.
2. The exact solution is u(x,y,t) = (1−exp(−t))q(x). The source is exp(−t)q(x) + (1−exp(−t))(1+cq(x)); bottom/top traces match u and left/right values are zero. This manufactured example develops an outflow layer from zero and tests changing-boundary and forcing terms. It is not an unforced physical startup.
3. Scrub from t = 0 to 1. At t = 1, baseline maximum ≈ 1.148 and SUPG maximum ≈ 0.5926; L² error is approximately 0.1008 versus 0.0750. Refine the mesh and compare the exact curve. Changing β invalidates this particular manufactured reference; changing ε or c updates its compatible source/traces.

### What the implementation adds

For constant physical coefficients, SUPG adds the element residual tested along the flow:

```text
ΣT τT ∫T (β·∇v) [uₜ − εΔu + β·∇u + cu − f].
```

P1 functions are affine on each triangle, so Δuₕ = 0 inside each element. Writing di = β·∇φi, the assembled additions are

```text
Hᵀij = τT |T| di / 3
Sᵀij = τT |T| di dj + c Hᵀij
Qᵀi(t) = τT di ∫T f(t).
```

The steady comparison is (A+S)b = F+Q−SIBg. The transient comparison is

```text
(M*II + HII) ḃ + (AII + SII)b
    = FI + QI − (AIB + SIB)g − (M*IB + HIB)ġ,
```

where M* is consistent or full-row-sum lumped mass. The same θ step is applied to every term. Even with lumped M*, HIB contributes to changing-boundary lifting. Source corrections use the same degree-five quadrature as Galerkin. Reaction uses consistent cM; SUPG also includes its cH residual contribution. Matrix inspectors expose S, H, and the actual comparison time mass; JSON v4 exports them.

The default coefficient is a documented heuristic based on the classical 1D coth choice:

```text
hT = 2|β| / Σi |β·∇φi|,    PeT = |β|hT/(2ε),
τadv = hT/(2|β|) [coth(PeT) − 1/PeT],
τT = strength × τadv / sqrt(1 + (c τadv)²).
```

The last factor is this app's reaction time-scale cap, not a universal optimum. The small-Pe limit is evaluated with a series; ε = 0 uses the large-Pe limit; β = 0 gives no SUPG. τ is semidiscrete and independent of Δt; reducing time-step error is a separate task. The displayed scalar mesh Péclet number includes physical/isotropic diffusivity only, not directional SUPG; its τ range is shown separately.

See the [pyoomph SUPG derivation](https://pyoomph.readthedocs.io/en/latest/tutorial/pde/convdiffu/supg.html) for the streamline test function, transient residual, and coth parameter. Our triangle length and reaction cap are explicitly specified above. The [ASPECT SUPG description](https://aspect-documentation.readthedocs.io/en/stable/user/methods/advection-stabilization/supg.html) also explains residual testing along the flow.

SUPG is not shock capturing, flux correction, or a discrete maximum-principle guarantee. It adds no crosswind diffusion; layers perpendicular to its stabilizing direction can remain troublesome. It vanishes when β = 0, so it cannot fix a pure reaction–diffusion layer. Resolve such layers with mesh refinement. Matrix/load overrides retain the physical PDE's S, H, and Q; after arbitrary algebraic edits the combination is no longer residual-consistent with the edited equation.

`tests/supg.test.mjs` checks zero-flow/zero-strength limits, τ limits, steady and time-dependent affine ADR patches across all mesh patterns and both mass/time choices, changing-boundary H lifting, absence of crosswind stiffness, and the measured layer comparisons above.

## Algebraic flux correction (AFC)

Select **Steady AFC layer** or **Nonnegative moving front** in the filtered **Problem → Example** selector. The comparison selector now offers **AFC · low-order graph diffusion** and **AFC · limited flux correction**. These are alternatives to SUPG and isotropic diffusion, not extra terms stacked on top of them.

AFC is nonlinear even for a linear PDE: its limiter depends on the solution. That does not always require a nonlinear system solve. This app uses iterative steady AFC and explicit transient flux-corrected transport (FCT):

| Method | Comparison computation | Cost reported |
| --- | --- | --- |
| Steady low order | One linear solve with graph diffusion | Solve time, zero limiter iterations |
| Steady limited AFC | Damped fixed-point iterations; reused LU or sparse preconditioner | Iterations, scaled nonlinear residual, convergence, AFC time |
| Transient low order | Lumped-mass SSPRK2 with graph diffusion | Internal substeps and stepping time |
| Transient limited FCT | Same SSPRK2 stages plus pairwise flux limiting | Internal substeps, limiter passes, stepping time |

### Steady layer experiment

1. Choose **Steady AFC layer**. Defaults are n = 16, ε = 0.005, β = (1,0), c = 0, source = 1, tolerance = 1e−8, relaxation = 0.5, and at most 500 iterations.
2. Select **With AFC** above the plot. The uncorrected Galerkin maximum is about 1.846. Limited AFC gives about 0.9643; switching to **AFC · low-order graph diffusion** gives about 0.9375. Compare the cross-section with the exact curve: the thin boundary layer remains unresolved.
3. On the tested machine the corrected steady solve took about 55 ms overall, including about 38 ms inside AFC, and converged in 110 iterations. These are illustrative timings, not a speed guarantee. Mesh size, tolerance, and relaxation affect cost; refinement increases linear-solver work (dense LU through n=32, sparse Krylov above).
4. Set the iteration limit to 1 to inspect failure reporting. The app labels the result **AFC iteration limit** and reports the last iterate and its residual; it does not mark it converged. Restore 500. If another case stalls, try smaller relaxation or mesh refinement and inspect the residual; convergence is not guaranteed by increasing the limit.

The stopping residual is the infinity norm of `L b − (F_low + limited_flux(b))`, scaled by the largest rowwise sum of absolute matrix contributions plus the absolute corrected right-hand side. A converged algebraic residual is separate from discretization error.

### Moving-front experiment

1. Choose **Nonnegative moving front**: n = 24, β = (0.55,0), ε = c = 0, Δt = 0.005, T = 0.8. Select **Comparison**, then play or scrub to the last frame.
2. The consistent-mass Crank–Nicolson baseline reaches approximately −0.232 / 1.177. Limited FCT stays nonnegative in this experiment, with maximum about 0.934 and zero negative nodes. Values are not clipped after solving.
3. Switch to **AFC · low-order graph diffusion** at the same time: the peak falls to about 0.428. Flux correction restores much of the sharp profile lost to graph diffusion. Switch to SUPG separately to compare its remaining undershoot and smaller error in this example.
4. Check accuracy as well as extrema. At T = 0.8, quadrature L² errors are approximately 0.0566 for the baseline, 0.1273 for low order, and 0.0685 for FCT. Removing negative values does not necessarily minimize L² error.
5. Increase Δt or physical diffusion to see automatic internal substeps. The requested output times remain the same, while the comparison divides each interval to satisfy its explicit low-order positivity restriction. The status shows substeps and limiter passes for the selected interval, cumulative substeps, and cumulative AFC stepping time. At time zero no steps have occurred.

The default moving-front trajectory took about 0.9 s overall in one local run, including about 44 ms of AFC stepping over 160 internal steps. The remaining work includes the consistent baseline, assembly, and diagnostics. The two SSPRK2 stages each evaluate a nonlinear limiter directly; they do not run Newton or fixed-point iterations. Stiff diffusion/reaction or finer meshes can require many more substeps, so explicit FCT may then be slower than an implicit method. The app stops with an explanation if more than 50,000 internal steps would be needed.

**Comparison convention:** the baseline still uses the selected backward-Euler or Crank–Nicolson method with consistent mass. AFC uses lumped mass and SSPRK2 regardless of that selection. Thus comparing those two curves changes both spatial stabilization and time discretization. To isolate flux correction, compare the two AFC options; they share the same mass, time integrator, and internal step sizes. Transient antidiffusion corrects the spatial graph viscosity only; it does not restore consistent mass.

### Operator, limiter, and limits

The full nodal operator is `A = εK + C + cM`, assembled before Dirichlet elimination. Symmetric graph diffusion has off-diagonals `D_ij = −max(0, A_ij, A_ji)` and diagonal entries chosen for zero row sums. Consequently `L = A + D` has nonpositive off-diagonals. Graph diffusion includes boundary edges and the positive off-diagonal contributions of consistent reaction; it is not just a scalar multiple of K. The reduced D shown in the inspector need not have zero row sums after boundary columns are eliminated.

Steady AFC uses a BJK-type symmetric local-extremum limiter with an adjustable μ (default 1). It restores a limited portion of each antisymmetric edge flux and solves repeatedly with the same L. This is an educational variant; it makes no claim of linearity preservation on every mesh patch. See [this primary article's BJK limiter discussion](https://link.springer.com/article/10.1007/s11075-023-01511-2) and [Jha and John's study of nonlinear AFC iteration schemes](https://www.wias-berlin.de/publications/wias-publ/run.jsp?number=2533&template=abstract&type=Preprint&year=).

Each transient Euler stage starts with a lumped-mass low-order update. Zalesak-style nodal flux budgets constrain pairwise corrections using neighboring old and low-order values. Equal and opposite corrections preserve the interior pair contribution to lumped mass; Dirichlet boundaries, forcing, and reaction can change the total integral. SSPRK2 combines two such stages. Internal steps satisfy `h ≤ 0.8 min_i(m_i / L_ii)` over positive interior diagonals.

This implementation accepts c ≥ 0 and assembled PDE coefficients. Clear manual A/f overrides before selecting AFC; changing the physical source expression or boundary data remains supported. Signed initial data, sources, and boundary values are allowed and may legitimately produce negative solutions. Transient nonnegativity relies on nonnegative discrete input/load/boundary data and the enforced step restriction, up to floating-point tolerance. The status flags signed data rather than claiming positivity. AFC is neither a universal accuracy improvement nor a cure for insufficient mesh resolution.

JSON export v4 includes AFC convergence/substep metadata and graph D. For steady AFC, `stabilizedA` contains the retained graph diffusion at the returned iterate. For transient AFC, `stabilizedA` and `stabilizedF` describe the **low-order backbone**; nonlinear corrections are applied separately at each explicit stage. There is no single implicit comparison step matrix or right-hand side, so those export fields are null and the step-matrix inspector is disabled. The comparison's implicit backward-error metric is also unavailable. SUPG exports remain included.

`dist/afc.js` implements graph construction, steady limiting, FCT stages, and automatic substeps. `dist/linear.js` holds the shared LU routines. Twelve checks in `tests/afc.test.mjs` cover graph signs and conservation, bounded pair fluxes, nonlinear convergence and iteration-limit reporting, moving-front nonnegativity/sharpness, changing-boundary manufactured solutions on all three meshes, signed data, and unsupported input rejection. Together with the existing checks, the test suite contains 64 tests. These finite tests do not establish accuracy or convergence for arbitrary input.

### Adjustable μ and affine preservation

Choose **Affine preservation · μ** in the filtered **Problem → Example** selector. Its exact solution is **u = 1+x+y**, with ε = 0.005, β = (1,0.3), c = 0, source = βx+βy+c(1+x+y), and exact Dirichlet traces. It starts on the original slope +1 mesh with n = 8. Positive reaction can also be tested without changing its reference because the source expression includes c. The default solve tolerance is 1e−11.

The **Limiter budget μ** control multiplies the steady nodal budgets: `R± = min(1, μQ±/P±)`. Its range is 0–10. μ = 0 is an exploratory low-order endpoint; the published limiter uses positive parameters. Larger μ relaxes limiting. This control has no effect on transient Zalesak limiting and is disabled for the explicit low-order comparison option.

**The original three mesh patterns have centrally symmetric interior vertex patches. On these meshes, μ = 1 is sufficient for linearity preservation for every affine field.** This is a geometric statement, not merely an empirical observation for the selected profile. The key identity is `u_i − u_i^min = u_i^max − u_i` for affine u on a centrally symmetric patch. With nonnegative graph weights, both flux budgets dominate their corresponding flux sums, so the steady limiter leaves affine fluxes unlimited.

This result is proved in **Barrenechea, John & Knobloch (2017), “An algebraic flux correction scheme satisfying the discrete maximum principle and linearity preservation on general meshes,” §6, Lemma 6.1**, in the [accepted author manuscript, printed p. 14](https://pure.strath.ac.uk/ws/portalfiles/portal/90418296/Barrenechea_etal_MMMAS_2017_An_algebraic_flux_correction_scheme_satisfying_the_discrete_maximum_principle_and_linearity_preservation.pdf). Their γ scales the same budget that this app scales with μ; their signed q is the negative of our positive sum of graph weights. **Knobloch (2023), §4.2, following equation (30)** explicitly states the μ = 1 symmetric-patch conclusion in the notation used here: [published article, Numerical Algorithms 94, 547–580](https://link.springer.com/article/10.1007/s11075-023-01511-2). These primary texts and the translation to the implemented formulas were checked on 2026-10-28. Our boundary graph construction remains the documented BJK-type variant; the affine budget argument only needs the complete symmetric neighbor set and nonnegative graph weights.

To explore the distinction between a reproduced affine solution and an inactive limiter:

1. On the default uniform mesh, compare μ = 0.25 and μ = 1. At μ = 1, the affine edge fluxes are fully restored. With stronger limiting, a small nodal error can persist because the retained graph diffusion cancels on affine fields. Therefore nodal error alone does not establish that the limiter is inactive.
2. The status separately shows **maximum nodal error**, **maximum retained affine edge flux**, and the **retained-diffusion residual evaluated on the exact affine field**. The nonlinear solve residual measures convergence of the computed solution, a different question.
3. Select **Shifted interior node · affine test** in Triangulation. This optional mesh moves the node at indices `(floor(n/2), floor(n/2))` by `(0.35h,0.2h)` while preserving connectivity and boundary positions. It breaks local symmetry. At n = 8, μ = 0, 0.5, and 1 give maximum nodal errors about 0.05473, 0.02830, and 5e−10, respectively; all three solves converge. The baseline Galerkin error is about 5e−15. At μ = 1 the exact-field retained flux and residual are zero; the remaining nodal error is iteration error.
4. For the shifted example, the displayed threshold is about **0.7341544**. It is computed from the exact profile's flux-budget ratios and is sufficient to remove limiting **for this field and these coefficients only**. It is not an all-direction geometry bound. Changing mesh size or coefficients changes it. μ = 1 is not asserted to suffice on arbitrary nonsymmetric meshes.

The shifted mesh uses actual triangle areas for source quadrature, L² error, and integrals; interpolation searches adjacent cells to find its physical triangles. The original uniform mesh generation is unchanged. Regression checks verify central symmetry and unlimited affine fluxes on all three original patterns for several affine directions, plus the shifted example, geometric quadrature, interpolation, μ validation, and the unchanged transient limiter.


### Open the convection-dominated examples directly

The app now starts with **Convection-dominated AFC layer**, with ε = 0.005, β = (1,0), n = 16 and mesh Péclet number 6.25. Use **Problem → Example** to select **Affine preservation · μ** without returning to the buttons above the plot. That affine preset solves `−0.005Δu + (1,0.3)·∇u = 1.3` with exact traces `u = 1+x+y`; n = 8 gives mesh Péclet number approximately 13.05. It is a convection-dominated consistency test, whereas the layer preset tests oscillation control.

With the local server running, open `http://127.0.0.1:4173/?example=afc-affine` for the affine test or `http://127.0.0.1:4173/?example=afc-layer` for the layer. Selecting any preset updates this URL parameter, so the chosen experiment survives reloads and can be bookmarked. Smooth Poisson remains available through the experiment selector and Reset. If the new selector is absent, the displayed page is an older version; reload the updated local checkout or deploy the current files before using the hosted site.

## Solid-body rotation (transient C)

Choose **Transient → C → Solid-body rotation · three shapes**, or open
`?example=solid-rotation`. This adds the benchmark in Kuzmin, *Explicit and implicit
FEM-FCT algorithms with flux linearization*, JCP 228 (2009), 2517–2534,
§7.1, equations (64)–(65), printed page 2526.

The velocity is β(x,y) = (0.5−y, x−0.5), with zero diffusion, reaction and source.
It rotates counterclockwise about (0.5,0.5), with period 2π. The initial field
`rotation(x,y)` contains three disjoint radius-0.15 bodies:

- A unit-height cylinder centered at (0.5,0.75), with the slot removed where
  |x−0.5| < 0.025 and y < 0.85.
- A cone centered at (0.5,0.25), height 1−r.
- A hump centered at (0.25,0.5), height (1+cos(πr))/4.

Here r is distance from each body's center divided by 0.15, and each profile is
zero outside its disk. The exact reference evaluates the initial field at the
inverse-rotated point, so it returns to the initial field at T=2π. Editing source,
boundary data or initial data disables this benchmark reference.

The default is n=32, Δt=π/100, T=2π (200 output intervals), a consistent-mass
Crank–Nicolson baseline and the existing lumped SSPRK2 spatial FCT comparison.
Play the trajectory, then scrub to the last frame. Compare FCT with low-order
graph diffusion to see how limiting restores sharpness, and with Galerkin to
see undershoots/overshoots. Try y=0.75 for the slot cross-section at t=0 or 2π.
The coarse mesh resolves the 0.05-wide slot poorly; refinement matters.

This is the paper's **test data**, not an implementation of its CN-FCT or flux
linearization algorithms. The app uses P1 triangles and its existing methods.
The paper prescribes zero inflow values; this app retains zero Dirichlet traces
on **all** edges. The exact support stays inside the domain, but numerical
spreading can reach the boundary, so monitor integral drift; exact discrete mass
conservation is not promised. The mesh now supports n≤128; n=32 remains the rotation preset default.

Convection integrates the affine velocity exactly against the P1 test functions.
The full AFC graph and reduced Galerkin matrix use the same element entries;
SUPG uses local centroid τ and quadrature for the varying streamline factors.
CFL and Péclet summaries use the maximum speed over the square, √(1/2).


### Rotating and diffusing (transient CD)

Choose **Transient → CD → Rotating and diffusing · three shapes**, or use
`?example=rotation-diffusion`. This extension keeps the same rotation field,
initial profiles, zero boundary traces, and T=2π, and adds adjustable physical
ε (default 0.0001) to both the baseline and comparison operators.
Try ε=0, 0.0001, and 0.001 to distinguish physical spreading from numerical
smoothing. Compare FCT, low-order graph diffusion, SUPG, and artificial diffusion.
Physical spreading prevents an exact return after one revolution.

For positive ε, the app deliberately supplies no exact reference or L² error:
the inviscid rotated profile is not the solution of this diffusion problem, and
whole-plane heat convolution would not satisfy the square's zero boundary data.
Setting ε=0 restores the inviscid reference for unedited benchmark data.
Boundary absorption can change the integral as the diffusing field reaches edges.


## Fine meshes: 64 and 128

The mesh slider supports 2–128 subdivisions per direction. At n=128 the mesh has
16,641 nodes, 32,768 P1 triangles and 16,129 interior unknowns. Defaults stay small;
choose 64 or 128 explicitly for refinement experiments.

At n≤32, existing dense pivoted LU remains in use. Above 32, all assembled and
inspection matrices use sparse rows. Linear solves use ILU(0)-preconditioned
BiCGSTAB, with restarted GMRES as a fallback; the original matrix is never
regularized. Guarded pivots affect the preconditioner only. Each solution is
checked against the true residual (relative 2-norm target 1e-11, acceptance up
to 5e-11). Failure to converge raises an error and leaves the last successful
plot marked stale. Difficult pure-convection or edited systems can fail.
Time-step preconditioners are reused, with the preceding solution as a starting
estimate. Steady AFC reuses its low-order preconditioner across limiter iterations.
Its nonlinear convergence tolerance remains a separate check.

All requested numerical time steps are computed (up to 240 as before). To bound
playback memory, n=33–64 stores at most 121 frames and n=65–128 at most 81,
including both endpoints. The timeline shows saved frames and computed steps;
exports include `computedSteps` and `savedFrameStride`. Changing inputs terminates
an obsolete worker calculation before starting the latest one. Sparse matrix
inspection, arbitrary entry edits, and coordinate-triplet JSON exports remain
available. Dense reconstruction is not needed for drawing the sparsity preview.

Measured locally in Node on this workstation, the default 200-step full revolution
at n=128 took approximately 15 s for pure rotation and 9 s for rotation with
ε=0.0001. These are solver-only examples, not browser/device guarantees. Both
saved 68 frames and used 1,000 internal FCT substeps. Pure-rotation FCT final
range was [0,0.967], with L² error about 0.114. Refining space alone does not
remove time-discretization error: the baseline still uses Δt=π/100.

Five tests in `tests/sparse.test.mjs` check sparse/dense agreement, nonlocal
entries and structured cloning, n=128 Poisson convergence, retained time steps
with thinned playback, and explicit rejection of an inconsistent singular solve.
The full test suite also covers the existing FEM, reaction, SUPG and AFC paths.


### Consistent browser asset versions

Run `node CDR_visualizer/scripts/version-assets.mjs` after editing runtime files
and before serving them directly. `npm start` and the GitHub Pages workflow do
this automatically. The script stamps a shared content hash onto local module,
worker and stylesheet URLs, so cached older dense-only code cannot be mixed with
the sparse solver. `--check` verifies the stamp without editing files. Reload
an already-open tab after updating; a tab's loaded JavaScript does not change
until navigation or reload. If the HTML itself is cached, use a hard refresh.
