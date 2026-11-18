# Modification rationale, boundary implementation, and assignment notes

[中文版](README_修改原因与Assignment.md) | [Build, run, and plot](README_EN.md)

This document explains the physical setup, boundary implementation, numerical configurations, and reasons for the choices in the final package.

## 1. Scope and sources

The Setup & boundaries task is to adapt `mod_usr.t` for the initial conditions and jet inflow, and to explain boundary conditions and ghost cells. Question 2 of `CMfAA_HPCassignment.pdf` covers geometry, initialization, and normalization; Question 3 asks how ghost cells are filled, how many layers are used on each side, and why. Question 4 covers the numerical discretization and its stencil.

The emphasis here is on initialization, boundaries, ghost cells, and their connection to spatial reconstruction. The physical model follows Massaglia et al. (2016), A&A 596 A12, Section 2.2, equations (5)–(8), and Case B in Table 1.

The original example is AMRVAC's `tests/hd/jet_cloud` (located at `/home/peder/codes/amrvac/tests/hd/jet_cloud` on the development machine). It describes a jet entering along $+x$ and interacting with an off-axis cloud. Although the source contains `IFTHREED` sections, its coordinate-system setting is `Cartesian_2D`.

Adapting it to 3D therefore requires consistent changes to the geometry, third velocity component, nozzle radius, boundary arrays, and parameter files; changing only `setup.pl -d` is insufficient. The adapted files are in the package root, with the originals in `original/` for comparison. This keeps the baseline self-contained and traceable.

## 2. Equations, variables, and physical parameters (assignment Question 1)

The model uses the ideal-gas Euler equations without magnetic fields, gravity, or radiative cooling, plus one passive tracer:

$$
\begin{aligned}
\frac{\partial \rho}{\partial t}+\nabla\cdot(\rho\mathbf{v})&=0,\\
\frac{\partial(\rho\mathbf{v})}{\partial t}
+\nabla\cdot\left(\rho\mathbf{v}\otimes\mathbf{v}+p\mathbf{I}\right)&=0,\\
\frac{\partial E}{\partial t}+\nabla\cdot\left[(E+p)\mathbf{v}\right]&=0,\\
\frac{\partial C}{\partial t}+\nabla\cdot(C\mathbf{v})&=0.
\end{aligned}
$$

Here $\mathbf{I}$ is the identity tensor and $\otimes$ denotes the outer product. The tracer density and total energy density are

$$
C=\rho f,\qquad
E=\frac{p}{\gamma-1}+\frac{\rho}{2}\left(v_x^2+v_y^2+v_z^2\right),
\qquad \gamma=\frac{5}{3}.
$$

The conserved state is

$$
\mathbf{U}=(\rho,m_x,m_y,m_z,E,C)^{\mathsf T},\qquad
\mathbf{m}=\rho\mathbf{v}.
$$

The physical primitive variables are $\rho,v_x,v_y,v_z,p,f$. However, the EOS conversion routines in the target installation do not convert the tracer slot, so that slot must explicitly contain $C$. A tracer is useful because density changes can result from either compression or material mixing; $f$ identifies material of jet origin.

The baseline parameters are

$$
\rho_j=1,\qquad \eta=\frac{\rho_j}{\rho_a(0)}=0.01,\qquad
M=4,\qquad r_c=40,\qquad \alpha=2.
$$

$$
p_j=p_a=\frac{1}{\gamma}=0.6,\qquad
c_{s,j}=\sqrt{\frac{\gamma p_j}{\rho_j}}=1,\qquad
v_j=M c_{s,j}=4.
$$

The Mach number is defined using the jet sound speed. Carrying over the original `jet_cloud` ambient-sound-speed convention would change the meaning of the same numerical value of $M$.

The original cloud parameter $\beta$, cloud-radius setting, and ambient-sound-speed parameter are removed because this model contains no cloud. The variable $r_c$ now denotes the atmospheric core radius, not the original cloud radius; this distinction is also noted in the code.

## 3. Geometry, normalization, and initialization (assignment Question 2)

The 3D Cartesian domain is $x,z\in[-4,4]$, $y\in[0,8]$, with the jet directed along $+y$. The units of length, density, and velocity are $r_j$, $\rho_j$, and $c_{s,j}$. The derived units are

$$
t_0=\frac{r_j}{c_{s,j}},\qquad
p_0=E_0=\rho_j c_{s,j}^{\,2},\qquad
M_0=\rho_j r_j^3.
$$

Here $E_0$ is an energy-density unit and $M_0$ is a mass unit, not the Mach number. The calculation is dimensionless; physical values have not been assigned to AMRVAC's `unit_*` variables. Default unit values printed in the log should not be interpreted as parsecs, years, or kelvin.

If the paper's $r_j=100\, \mathrm{pc}$ is adopted, then $r_c=40r_j=4\, \mathrm{kpc}$. Converting time to years also requires a physical value of $c_{s,j}$.

The spherical radius and initial interior state are

$$
\begin{aligned}
r&=\sqrt{x^2+y^2+z^2},\\
\rho_a(r)&=\frac{100}{1+(r/40)^2},\\
p&=0.6,\qquad v_x=v_y=v_z=0,\qquad C=0.
\end{aligned}
$$

The spherical radius $r$ defines the atmosphere; the transverse radius $R=\sqrt{x^2+z^2}$ defines the nozzle. They must not be interchanged.

The pre-filled jet is removed because the paper initializes the domain with ambient gas and injects the jet through a boundary. The off-axis cloud is removed because the intended problem is propagation through a stratified atmosphere rather than a jet–cloud collision.

Gravity is unnecessary for the specified initial equilibrium: pressure is uniform, so the stationary Euler momentum equation has no pressure-gradient force. Temperature varies as $1/\rho$.

`initonegrid_usr` first assigns density, velocities, and pressure, then calls `eos%to_conserved` once. Initially $C=0$. Do not assume the EOS automatically converts a nonzero fraction $f$ into $\rho f$.

The supplied atmosphere implements Case B with $\alpha=2$. Simply setting $\alpha=0$ would introduce a denominator of 2 in this profile rather than reproduce the intended uniform Case A normalization; that case needs a separate branch.

## 4. The six physical boundaries (assignment Question 3)

For direction $d$, the lower face has boundary index $i_B=2d-1$ and the upper face has index $i_B=2d$. The code variable is named `iB`.

| Domain face | Boundary index | Type |
|---|---|---|
| $x_{\min}$ | $i_B=1$ | `cont` |
| $x_{\max}$ | $i_B=2$ | `cont` |
| $y_{\min}$ | $i_B=3$ | `special`, filled by `specialbound_usr` |
| $y_{\max}$ | $i_B=4$ | `cont` |
| $z_{\min}$ | $i_B=5$ | `cont` |
| $z_{\max}$ | $i_B=6$ | `cont` |

The special boundary changes from `case(1)` to `case(3)` because injection moves from `x_min` to `y_min`; the parameter file must use the corresponding face.

The 3D hydrodynamic variables require five boundary entries: $\rho,m_x,m_y,m_z,E$. In this AMRVAC version the tracer is registered with `need_bc=.false.`, and the built-in boundary type is inherited from density. The special boundary still assigns $C$ explicitly. The parameter arrays therefore use `5*'cont'` or `5*'special'`, rather than the four entries used in 2D or an extra independent entry for the tracer.

### 4.1 Zero-gradient extrapolation: cont

`cont` copies the nearest interior cell into every ghost layer. At the upper $x$ face, with $N$ the last interior cell,

$$
\mathbf{U}_{N+1}=\mathbf{U}_{N},\qquad
\mathbf{U}_{N+2}=\mathbf{U}_{N}.
$$

At the lower $x$ face, with the interior starting at cell 1,

$$
\mathbf{U}_{0}=\mathbf{U}_{1},\qquad
\mathbf{U}_{-1}=\mathbf{U}_{1}.
$$

This imposes a discrete zero normal gradient. It does not mirror successive interior layers. The choice approximates an open boundary and follows the paper's zero-gradient treatment of the other domain faces.

It is not a strictly non-reflecting characteristic boundary and does not prohibit inflow: material can enter if the extrapolated normal velocity points inward. The original example's transverse `noinflow` condition clips inward normal momentum, so it is not equivalent to `cont`. For long-time calculations, domain-size or boundary-sensitivity studies are needed to assess contamination by the outer boundaries.

### 4.2 Reflection outside the nozzle on the inlet plane

The inlet plane is $y=0$. First, a reflected state $Q_r$ is constructed over the entire `y_min` ghost region. For example, if ghost indices in $y$ are 1 and 2 and the interior starts at 3, the outer ghost cell 1 mirrors interior cell 4, and the inner ghost cell 2 mirrors interior cell 3:

```text
... [ghost 1] [ghost 2] | y=0 | [inside 3] [inside 4] ...
```

The general mirror index is

$$
j_{\mathrm{mirror}}=2j_{g,\max}+1-j,
\qquad j_{g,\max}=\texttt{ixOmax2}.
$$

The corresponding source expression is `jmirror=2*ixOmax2+1-j`. Copying interior cell 3 into both ghost layers would instead be constant extrapolation and would change the gradients used in reconstruction.

Density, tangential momenta, energy, and tracer density ($\rho,m_x,m_z,E,C$) are copied with even parity. Only the normal momentum $m_y$ has odd parity:

$$
(m_y)_g=-(m_y)_{\mathrm{mirror}}.
$$

Tangential momenta keep their signs because this is an inviscid slip-wall reflection, not a no-slip wall. Energy is unchanged because reversing the normal velocity does not change its squared contribution to kinetic energy. Mirroring $C$ together with density preserves the material fraction.

The mirrored state already contains conserved variables. Applying `to_conserved` again would incorrectly multiply momenta by density a second time.

### 4.3 Circular nozzle and smooth transition

The prescribed on-axis jet state is

$$
\rho=1,\qquad v_x=v_z=0,\qquad m_y=\rho_jv_j=4,
\qquad p=0.6,\qquad f=1.
$$

The code uses `R=hypot(x,z)` to evaluate $R=\sqrt{x^2+z^2}$ without unnecessarily overflowing the squared coordinates. Following equation (8) of the paper,

$$
Q=Q_r+(Q_j-Q_r)S(R),\qquad
S(R)=\frac{1}{\cosh\!\left[\left(\frac{R}{R_s}\right)^n\right]}.
$$

The smoothing parameters are

$$
(R_s,n)=
\begin{cases}
(1.4,8), & Q=\rho,\\
(1,6), & Q\in\{v_x,m_y,v_z,p,f\}.
\end{cases}
$$

The different density and momentum/velocity profiles follow the paper's inlet prescription. They should not be replaced by a single sharp top-hat profile if that prescription is to be retained.

At $R=0$, $S=1$ and the prescribed jet state is recovered. Far from the nozzle, $S\to0$ and the reflected state is recovered. Thus, “injection inside, reflection outside” describes the limiting behavior; the transition region is a smooth blend.

$Q_r$ is taken from the current interior solution at every boundary update, rather than from the initial atmosphere, because reflection must respond to the evolving flow.

The reconstruction proceeds as follows:

1. Recover $\rho_r,v_{x,r},v_{z,r},m_{y,r},p_r$ and $f_r=C_r/\rho_r$ from the mirrored conserved state $\mathbf{U}$.
2. Apply the smoothing rule to $\rho,v_x,v_z,m_y,p,f$.
3. Rebuild momentum, axial velocity, total energy, and tracer density:

$$
\begin{aligned}
m_x&=\rho v_x,\qquad m_z=\rho v_z,\qquad v_y=\frac{m_y}{\rho},\\
E&=\frac{p}{\gamma-1}
+\frac{1}{2}\left[\rho\left(v_x^2+v_z^2\right)+\frac{m_y^2}{\rho}\right],\\
C&=\rho f.
\end{aligned}
$$

The axial quantity being blended is $m_y$, not $v_y$, as specified by the paper. Directly blending total energy is generally not equivalent to blending pressure and the stated velocity/momentum quantities. The $v_z$ contribution is required to recover the correct pressure in a flow with transverse motion.

The argument is capped before exponentiation to avoid overflow in both the power and `cosh`. The maximum hyperbolic-cosine argument is 80, giving a residual weight of approximately $1/\cosh(80)\approx3.61\times10^{-35}$. This is a negligible numerical tail truncation, not an additional physical edge of the nozzle.

No explicit inlet perturbation is imposed, following the paper. Non-axisymmetric modes can grow from numerical noise, but the supplied short-time configurations do not establish their fully developed evolution.

## 5. Ghost-cell counts and their rationale

The number of ghost layers is determined by the spatial reconstruction stencil, not directly by the number of dimensions or Runge–Kutta stages.

| Parameter files | Reconstruction | Ghost layers per side | Reason |
|---|---|---|---|
| `smoke.par`, `demo.par`, `minmod2.par` | minmod | 2 | Two layers cover the boundary-flux stencil for piecewise-linear reconstruction. |
| `mp5_3.par` | MP5 | 3 | MP5 has a wider stencil; this AMRVAC version requires at least three layers. |

`nghostcells` in the `.par` file's `boundlist` specifies the count. The user module supports both configurations with the same physical inlet and layer-by-layer reflection rule.

For minmod on a uniform grid, the limited increment $\delta q_i$ and slope $s_i$ are

$$
\delta q_i=\operatorname{minmod}\left(q_i-q_{i-1},\ q_{i+1}-q_i\right),
\qquad s_i=\frac{\delta q_i}{\Delta x}.
$$

A cell's slope uses its immediate neighbors, but an interface flux requires reconstructed states from both adjacent cells. Reconstructing the first exterior ghost cell therefore requires the second ghost cell. The complete spatial update stencil for interior cell $i$ can extend from $i-2$ to $i+2$.

MP5 uses a wider reconstruction stencil and consequently needs a third layer of exterior data. It can retain higher-order spatial information in smooth regions; the three-layer choice supplies the required stencil. The number of time-integration stages does not change that spatial requirement.

Physical domain boundaries are filled using the `cont` or `special` rules in Section 4. Internal block halos are filled from adjacent blocks, using MPI communication when necessary. The mirror formula uses the actual ghost-region index bounds, so each of the two or three layers maps to its corresponding interior layer.

Each block contains $8^3$ active cells, allowing block decomposition in a small 3D setup. The rectangular local array, including ghost layers, contains

$$
N_{\mathrm{local}}=
\begin{cases}
(8+2\times2)^3=12^3=1728, & \text{minmod, 2 ghost layers},\\
(8+2\times3)^3=14^3=2744, & \text{MP5, 3 ghost layers}.
\end{cases}
$$

These are local storage counts, not additional physical cells. Setting `refine_max_level=1` keeps a single grid level, so this configuration has no AMR coarse–fine interfaces.

## 6. Numerical method and parameter choices (assignment Question 4 and the CFL part of Question 5)

The baseline combines a fixed-grid finite-volume method, TVDLF fluxes, minmod piecewise-linear reconstruction, and `threestep`. `mp5_3.par` keeps the flux method and time integrator unchanged while selecting MP5 reconstruction and three ghost layers.

In the target installation, `threestep` defaults to `ssprk3`, `ssprk_order` defaults to 3, and the spatial update is unsplit by default. SSPRK(3,3) can be written as

$$
\begin{aligned}
\mathbf{U}^{(1)}&=\mathbf{U}^{n}+\Delta t\,L\!\left(\mathbf{U}^{n}\right),\\
\mathbf{U}^{(2)}&=\frac{3}{4}\mathbf{U}^{n}
+\frac{1}{4}\left[\mathbf{U}^{(1)}+\Delta t\,L\!\left(\mathbf{U}^{(1)}\right)\right],\\
\mathbf{U}^{n+1}&=\frac{1}{3}\mathbf{U}^{n}
+\frac{2}{3}\left[\mathbf{U}^{(2)}+\Delta t\,L\!\left(\mathbf{U}^{(2)}\right)\right].
\end{aligned}
$$

Here $L$ is the spatial discretization operator, $n$ denotes the current time level, and $(1)$ and $(2)$ denote intermediate stages.

The baseline is simple to interpret when checking the inlet and boundaries, making errors easier to isolate than with the original example's MP5/fivestep combination. TVDLF can be understood as an average flux with dissipation controlled by a local maximum signal speed. It is more diffusive at contacts than a contact-resolving HLLC solver, so all smooth mixing visible in the output should not be interpreted as physical turbulent mixing.

Third-order time integration does not make the full method third order. Minmod reconstruction is generally second order in smooth regions and loses order near extrema and discontinuities; spatial errors and solution regularity limit the combined accuracy. There are no additional physical source terms, so no gravitational or cooling source integration or associated extra timestep restriction is present.

The CFL coefficient is $C_{\mathrm{CFL}}=0.3$, with `typecourant='maxsum'`. For a uniform Cartesian grid, the main timestep restriction is approximately

$$
\Delta t\leq
\frac{C_{\mathrm{CFL}}}{
\displaystyle\max_{\mathrm{cells}}\left[
\frac{|v_x|+c_s}{\Delta x}
+\frac{|v_y|+c_s}{\Delta y}
+\frac{|v_z|+c_s}{\Delta z}
\right]}.
$$

The sound speed and normal Euler characteristic speeds are

$$
c_s=\sqrt{\frac{\gamma p}{\rho}},\qquad
\lambda\in\{v_n-c_s,\ v_n,\ v_n+c_s\}.
$$

$v_n$ is the normal velocity; contact, transverse, and tracer disturbances are advected at $v_n$. Output times and the final time can further shorten the actual timestep, so an estimate from the on-axis jet state alone is not the timestep for every update.

The value 0.3 gives a conservative baseline for checking the interfaces; the original value 0.8 is not intrinsically invalid. `small_values_method='error'` avoids silently masking a problem through averaging; the user is responsible for checking finite and positive output values. `dtmin=1e-12` is a stopping threshold for excessively small timesteps, not the timestep prescribed for every iteration.

Using `refine_max_level=1` isolates the setup from AMR effects. A log mentioning a default refinement estimator does not mean multiple refinement levels were used.

| Configuration | Active grid | Cell spacing | Cells per jet radius | Final time | Output interval |
|---|---|---|---|---|---|
| smoke / minmod2 / mp5_3 | $32^3$ | $\Delta x=\Delta y=\Delta z=0.25$ | 4 | $t_{\mathrm{end}}=0.2$ | $0.1$ |
| demo | $64^3$ | $\Delta x=\Delta y=\Delta z=0.125$ | 8 | $t_{\mathrm{end}}=1$ | $0.5$ |

The small domain and short evolution limit local computation and output costs. The demo configuration increases the inlet resolution to show early evolution. Runs produce DAT snapshots for restart/post-processing and VTU files for visualization; neither is preloaded in the package.

Runs with different final times do not constitute a controlled grid-convergence study and should not be used to infer a convergence order.

## 7. Guide to the user routines

| Routine | Implementation and reason |
|---|---|
| `usr_init` | Retains the original five callback registrations and `hd_activate`, and selects `Cartesian_3D` for three spatial and velocity directions. |
| `initglobaldata_usr` | Defines Case B parameters and checks dimensions, gamma, the energy equation, tracer count, inlet position, auxiliary-output count, and ghost-layer count. Inconsistent configurations should stop explicitly. |
| `initonegrid_usr` | Initializes the stationary stratified atmosphere and performs the EOS conversion once, as described in Section 3. |
| `specialbound_usr` | Handles `case(3)`: layer-by-layer reflection, a smooth circular 3D nozzle, and reconstruction of the full energy, as described in Section 4. |
| `extra_var_output` | Keeps the $p/\rho$ temperature proxy. The unused cooling-module import is removed because cooling is not enabled. |
| `extra_var_names_output` | Keeps the output name `T`, consistent with `nwauxio=1`. |

The indices $i,j,k$ correspond to $x,y,z$. `^L` and `^S` are LASY index macros expanded into `mod_usr.f` during compilation.

## Sources and version

The original example is AMRVAC `tests/hd/jet_cloud`; its two input files are retained in `original/`. The physical reference is Massaglia et al. (2016), A&A 596 A12, Section 2.2, equations (5)–(8), and Table 1, Case B.

The course reference is `CMfAA_HPCassignment.pdf`. This contribution focuses on Questions 2 and 3 and their connection to Question 4, with the governing equations and parameters included for context. Compare `original/mod_usr.t` with the adapted root-level `mod_usr.t` to inspect the implementation changes.

The interfaces are based on AMRVAC commit `4ac04b2672f11ba5b8909bed14427aefd1d60425`. In that version, `src/io/mod_input_output.t` sets the two-layer default and raises the requirement to at least three for MP5; `src/mod_boundary_conditions.t` implements the built-in boundary rules; `src/hd/mod_hd_phys.t` and `src/hd/mod_hd_eos.t` define tracer registration and conversion conventions.

The package does not include a copy of the solver source. Compile it against the recipient's AMRVAC installation.
