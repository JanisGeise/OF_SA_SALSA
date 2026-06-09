# Source-Code Consistent Equations

This document summarizes the equations implemented by the OpenFOAM v2506
`SpalartAllmaras`, `SpalartAllmarasDDES`, and the local
`SpalartAllmarasSALSA` / `SpalartAllmarasSALSADDES` models.

The transport equation is written without OpenFOAM's phase-fraction factor.
The density field is denoted by $\rho$; for incompressible usage it is
effectively one. Source terms from `fvOptions` are omitted below.

## Common Definitions

The transported variable is the modified kinematic viscosity $\tilde{\nu}$.
The turbulent viscosity is

$$
\nu_t = \tilde{\nu} f_{v1}.
$$

The common auxiliary functions are

$$
\chi = \frac{\tilde{\nu}}{\nu},
$$

$$
f_{v1} = \frac{\chi^3}{\chi^3 + C_{v1}^3},
$$

$$
f_{v2} = 1 - \frac{\chi}{1 + \chi f_{v1}},
$$

and, if the optional `ft2` switch is enabled,

$$
f_{t2} = C_{t3}\exp(-C_{t4}\chi^2).
$$

If `ft2` is disabled, which is the default, then

$$
f_{t2} = 0.
$$

The effective diffusivity used for the $\tilde{\nu}$ equation is

$$
D_{\tilde{\nu}} = \frac{\tilde{\nu} + \nu}{\sigma_{\nu_t}}.
$$

OpenFOAM uses

$$
\Omega = \sqrt{2}\left|\operatorname{skew}(\nabla \mathbf{U})\right|.
$$

For the standard SA form,

$$
\tilde{S}
= \max\left(
    \Omega
    + f_{v2}\frac{\tilde{\nu}}{(\kappa \tilde{d})^2},
    C_s\Omega
\right).
$$

The wall destruction ratio and damping function are

$$
r = \min\left(
    \frac{\tilde{\nu}}
    {\max(\tilde{S},\epsilon)(\kappa \tilde{d})^2},
    10
\right),
$$

$$
g = r + C_{w2}(r^6 - r),
$$

$$
f_w = g
\left(
    \frac{1 + C_{w3}^6}{g^6 + C_{w3}^6}
\right)^{1/6}.
$$

The standard coefficient relation is

$$
C_{w1} = \frac{C_{b1}}{\kappa^2}
       + \frac{1 + C_{b2}}{\sigma_{\nu_t}}.
$$

The step-like OpenFOAM functions `pos` and `neg` are denoted here by
$H(a)$ and $H(-a)$ respectively. The DDES shielding equations also use
$\nu_{\mathrm{eff}}$, the effective kinematic viscosity returned by
`nuEff()`.

## Standard Spalart-Allmaras RANS

The RANS model uses

$$
\tilde{d} = y,
$$

where $y$ is the wall distance.

The transported equation is

$$
\frac{\partial(\rho\tilde{\nu})}{\partial t}
+ \nabla\cdot(\rho\mathbf{U}\tilde{\nu})
=
\nabla\cdot(\rho D_{\tilde{\nu}}\nabla\tilde{\nu})
+ \frac{C_{b2}}{\sigma_{\nu_t}}\rho
   |\nabla\tilde{\nu}|^2
+ P_{\tilde{\nu}}
- D_{\tilde{\nu},w},
$$

with

$$
P_{\tilde{\nu}}
= C_{b1}\rho\tilde{S}\tilde{\nu}(1 - f_{t2}),
$$

and

$$
D_{\tilde{\nu},w}
=
\left(
    C_{w1}f_w
    - \frac{C_{b1}}{\kappa^2}f_{t2}
\right)
\rho\frac{\tilde{\nu}^2}{\tilde{d}^2}.
$$

Relevant source-code defaults:

$$
\sigma_{\nu_t}=0.66666,\quad
\kappa=0.41,\quad
C_{b1}=0.1355,\quad
C_{b2}=0.622,
$$

$$
C_{w2}=0.3,\quad
C_{w3}=2.0,\quad
C_{v1}=7.1,\quad
C_s=0.3.
$$

## SALSA RANS

The SALSA RANS model also uses

$$
\tilde{d} = y.
$$

The transported equation has the same structure as the SA equation, but the
production and destruction coefficients are made local:

$$
\frac{\partial(\rho\tilde{\nu})}{\partial t}
+ \nabla\cdot(\rho\mathbf{U}\tilde{\nu})
=
\nabla\cdot(\rho D_{\tilde{\nu}}\nabla\tilde{\nu})
+ \frac{C_{b2}}{\sigma_{\nu_t}}\rho
   |\nabla\tilde{\nu}|^2
+ P_{\tilde{\nu}}^{\mathrm{SALSA}}
- D_{\tilde{\nu},w}^{\mathrm{SALSA}}.
$$

The local SALSA strain magnitude is

$$
S^\star
=
\sqrt{2}
\left|
    \operatorname{dev}\left(\operatorname{symm}(\nabla\mathbf{U})\right)
\right|.
$$

With the default `useSmod true`, the source-code variable named `Stilda` is

$$
\tilde{S}_{\mathrm{SALSA}} = S^\star.
$$

With `useSmod false`, the model reverts to the standard SA definition

$$
\tilde{S}_{\mathrm{SALSA}} = \tilde{S}.
$$

The SALSA coefficient multiplier is

$$
r_\gamma
=
\min\left(
    \frac{\tilde{\nu}}
    {\max(\tilde{S}_{\mathrm{SALSA}},\epsilon)(\kappa\tilde{d})^2},
    10
\right),
$$

$$
a_1 = (1.01 r_\gamma)^{0.65},
$$

$$
a_2
=
\left[
    \max\left(0, 1 - \tanh\left(\frac{\chi}{68}\right)\right)
\right]^{0.65},
$$

$$
\gamma = \max(a_1, a_2),
$$

$$
\Gamma_{\mathrm{eff}}
=
\sqrt{
    \min\left(
        1.25,
        \max(\gamma, 0.75)
    \right)
}.
$$

The effective coefficients are

$$
C_{b1,\mathrm{eff}} = C_{b1}\Gamma_{\mathrm{eff}},
$$

$$
C_{w1,\mathrm{eff}}
=
\frac{C_{b1,\mathrm{eff}}}{\kappa^2}
+ \frac{1 + C_{b2}}{\sigma_{\nu_t}}.
$$

The SALSA production term is

$$
P_{\tilde{\nu}}^{\mathrm{SALSA}}
=
C_{b1,\mathrm{eff}}
\rho\tilde{S}_{\mathrm{SALSA}}\tilde{\nu}(1 - f_{t2}),
$$

and the destruction term is

$$
D_{\tilde{\nu},w}^{\mathrm{SALSA}}
=
\left(
    C_{w1,\mathrm{eff}}f_w
    - \frac{C_{b1,\mathrm{eff}}}{\kappa^2}f_{t2}
\right)
\rho\frac{\tilde{\nu}^2}{\tilde{d}^2}.
$$

The local repository default for SALSA is

$$
C_{b2}=0.
$$

All other shared coefficients use the common defaults unless overridden.

### `useRmod`

If `useRmod false`, SALSA uses the standard $r$ and $f_w$ definitions above.

If `useRmod true`, the wall-damping ratio is replaced before computing $f_w$.
The code first forms

$$
\tilde{S}_{r\mathrm{mod}}
=
S_r\left(\frac{1}{\max(\chi,\epsilon)} + f_{v1}\right),
$$

where $S_r$ is the same field currently passed to the wall-damping function.
With `useSmod true`, this is $S^\star$; with `useSmod false`, this is
$\tilde{S}$.

Then

$$
\Psi
=
\sqrt{\frac{\rho_\infty}{\rho}}
\frac{\tilde{\nu}}{(\kappa\tilde{d})^2},
$$

$$
r_{\mathrm{mod}}
=
\min\left(
    1.6\tanh\left(
        \frac{0.7\Psi}{\max(\tilde{S}_{r\mathrm{mod}},\epsilon)}
    \right),
    10
\right).
$$

The model then evaluates

$$
g_{\mathrm{mod}}
=
r_{\mathrm{mod}} + C_{w2}(r_{\mathrm{mod}}^6 - r_{\mathrm{mod}}),
$$

$$
f_{w,\mathrm{mod}}
=
g_{\mathrm{mod}}
\left(
    \frac{1 + C_{w3}^6}{g_{\mathrm{mod}}^6 + C_{w3}^6}
\right)^{1/6}.
$$

In the destruction term, $f_w$ is replaced by $f_{w,\mathrm{mod}}$.

## Standard Spalart-Allmaras DDES

`SpalartAllmarasDDES` uses the same standard SA transport equation as above,
but replaces the RANS length scale by the DDES length scale.

The LES length scale is

$$
\ell_{\mathrm{LES}} = \psi C_{\mathrm{DES}}\Delta,
$$

where $\Delta$ is the selected OpenFOAM LES delta filter, for example
`cubeRootVol`. In a case dictionary this is the `delta` entry in the `LES`
sub-dictionary. The equations below are independent of the particular delta
filter choice except through this value of $\Delta$.

If `lowReCorrection false`, the code leaves

$$
\psi = 0.
$$

If `lowReCorrection true`, which is the source-code default,

$$
\psi
=
\sqrt{
    \min\left[
        100,
        \frac{
            1
            - \frac{C_{b1}}{C_{w1}\kappa^2 f_{w\star}}
              \left(f_{t2} + (1 - f_{t2})f_{v2}\right)
        }{
            \max\left(
                \epsilon,
                f_{v1}\max(10^{-10},1-f_{t2})
            \right)
        }
    \right]
}.
$$

The standard DDES shielding function is

$$
f_d
=
1 - \tanh\left[(C_{d1}r_d)^{C_{d2}}\right],
$$

with

$$
r_d
=
\min\left(
    \frac{\nu_{\mathrm{eff}}}
    {\max(|\nabla\mathbf{U}|,\epsilon)(\kappa y)^2},
    10
\right).
$$

For the standard DDES shielding mode,

$$
\tilde{d}_{\mathrm{DDES}}
=
\max\left[
    y
    - f_d\max(y-\ell_{\mathrm{LES}},0),
    \epsilon_d
\right],
$$

where $\epsilon_d$ is OpenFOAM's small positive length-scale limiter.

The field reported by `LESRegion()` is

$$
\mathrm{LESRegion}
=
H\left(y-\tilde{d}_{\mathrm{DDES}}\right),
$$

implemented as

$$
\operatorname{neg}(\tilde{d}_{\mathrm{DDES}} - y).
$$

### `shielding ZDES2020`

If `shielding ZDES2020`, the model starts from the standard $f_d$ above and
then modifies it.

Let $\mathbf{n}$ be the wall-normal vector from `wallDist`. Define

$$
G_{\tilde{\nu}}
=
\frac{
    C_{d3}\max(\nabla\tilde{\nu}\cdot\mathbf{n},0)
}{
    \max(|\nabla\mathbf{U}|,\epsilon)\kappa y
},
$$

$$
f_{d,G\tilde{\nu}}
=
1 - \tanh\left[(C_{d1}G_{\tilde{\nu}})^{C_{d2}}\right].
$$

Also define

$$
G_\Omega
=
-(\nabla|\nabla\times\mathbf{U}|\cdot\mathbf{n})
\sqrt{
    \frac{\tilde{\nu}}
    {\max(|\nabla\mathbf{U}|^3,\epsilon)}
}.
$$

The code uses

$$
a
=
\frac{\frac{7}{6}C_{d4}-G_\Omega}{C_{d4}/6},
$$

and

$$
f_{RG_\Omega}
=
H(C_{d4}-G_\Omega)
+
\frac{
    H\left(\frac{4}{3}C_{d4}-G_\Omega\right)H(G_\Omega-C_{d4})
}{
    1
    +
    \exp\left[
        \min\left(
            \frac{-6a}{\max(1-a^2,\epsilon)},
            50
        \right)
    \right]
}.
$$

If `usefP2 true`, then

$$
f_{d,G\tilde{\nu}}
\leftarrow
f_{d,G\tilde{\nu}}
\,\frac{
    1-\tanh\left[(C_{d1}\beta_{\mathrm{ZDES}}r_d)^{C_{d2}}\right]
}{
    \max(f_d,\epsilon)
}.
$$

The final ZDES2020 shielding function is

$$
f_d
\leftarrow
f_d
\left(1 - (1-f_{d,G\tilde{\nu}})f_{RG_\Omega}\right).
$$

### `useSigma`

The upstream OpenFOAM `SpalartAllmarasDDES` model inherits the `useSigma`
switch from `SpalartAllmarasDES`.

If `useSigma false`, which is the default, the model uses the standard SA
$\tilde{S}$ with the DDES length scale.

If `useSigma true`, the production strain measure is replaced by

$$
\tilde{S}_\sigma
=
\Omega
- f_d H(y-\ell_{\mathrm{LES}})(\Omega - S_\sigma),
$$

and

$$
\tilde{S}
=
\max\left(
    \tilde{S}_\sigma
    + f_{v2}\frac{\tilde{\nu}}{(\kappa\tilde{d}_{\mathrm{DDES}})^2},
    C_s\tilde{S}_\sigma
\right).
$$

The default shielding coefficient also changes:

$$
C_{d1} =
\begin{cases}
8, & \text{if `useSigma false`},\\
10, & \text{if `useSigma true`, read as `Cd1Sigma`}.
\end{cases}
$$

The sigma invariant is computed from

$$
\mathbf{G} = (\nabla\mathbf{U})^T\nabla\mathbf{U}.
$$

Let

$$
I_1 = \operatorname{tr}(\mathbf{G}),
$$

$$
I_2 = \frac{1}{2}\left[I_1^2-\operatorname{tr}(\mathbf{G}^2)\right],
$$

$$
I_3 = \det(\mathbf{G}),
$$

$$
a_1 = \max\left(\frac{I_1^2}{9}-\frac{I_2}{3},\epsilon\right),
$$

$$
a_2 = \frac{\min(I_1,M)^3}{27}-\frac{I_1I_2}{6}+\frac{I_3}{2},
$$

$$
a_3
=
\frac{1}{3}
\arccos\left[
    \max\left(
        -1+\epsilon,
        \min\left(1-\epsilon,\frac{a_2}{a_1^{3/2}}\right)
    \right)
\right],
$$

where $M$ is OpenFOAM's large limiter for a field with dimensions
$1/t^2$. Then

$$
\sigma_1
=
\sqrt{
    \max\left(
        \frac{I_1}{3}
        +2\sqrt{a_1}\cos(a_3),
        \epsilon
    \right)
},
$$

$$
\sigma_2
=
\sqrt{
    \max\left(
        \frac{I_1}{3}
        -2\sqrt{a_1}\cos\left(\frac{\pi}{3}+a_3\right),
        \epsilon
    \right)
},
$$

$$
\sigma_3
=
\sqrt{
    \max\left(
        \frac{I_1}{3}
        -2\sqrt{a_1}\cos\left(\frac{\pi}{3}-a_3\right),
        \epsilon
    \right)
},
$$

and

$$
S_\sigma
=
C_{\mathrm{trans}}
\,
\frac{
    \sigma_3(\sigma_1-\sigma_2)(\sigma_2-\sigma_3)
}{
    \max(\sigma_1^2,\epsilon)
}.
$$

For SA-DES/DDES, OpenFOAM sets

$$
C_{\mathrm{trans}} = 67.7.
$$

## SALSA DDES

The local `SpalartAllmarasSALSADDES` model uses the SALSA transport equation,
SALSA production/destruction coefficient modifications, and the same DDES
length-scale and shielding definitions as above:

$$
\tilde{d}
=
\tilde{d}_{\mathrm{DDES}}.
$$

Thus the equation is

$$
\frac{\partial(\rho\tilde{\nu})}{\partial t}
+ \nabla\cdot(\rho\mathbf{U}\tilde{\nu})
=
\nabla\cdot(\rho D_{\tilde{\nu}}\nabla\tilde{\nu})
+ \frac{C_{b2}}{\sigma_{\nu_t}}\rho
   |\nabla\tilde{\nu}|^2
+ P_{\tilde{\nu}}^{\mathrm{SALSA}}
- D_{\tilde{\nu},w}^{\mathrm{SALSA}},
$$

where all SALSA definitions are evaluated with
$\tilde{d}=\tilde{d}_{\mathrm{DDES}}$.

The DDES controls implemented by the local model are

$$
C_{\mathrm{DES}},\quad
\mathrm{lowReCorrection},\quad
\mathrm{shielding}\in\{\mathrm{standard},\mathrm{ZDES2020}\},
\quad
C_{d1}, C_{d2}, C_{d3}, C_{d4},
\quad
\beta_{\mathrm{ZDES}},
\quad
\mathrm{usefP2}.
$$

The SALSA controls remain

$$
\mathrm{useSmod},\quad
\mathrm{useRmod},\quad
\rho_\infty,\quad
\mathrm{ft2}.
$$

### Deliberate `useSigma` omission

The local `SpalartAllmarasSALSADDES` model does not expose the upstream
`useSigma` switch. In OpenFOAM's standard `SpalartAllmarasDDES`, `useSigma`
changes the production strain measure to $\tilde{S}_\sigma$. In SALSA,
`useSmod` already changes the source-code `Stilda` field from the vorticity
form to the strain-rate field $S^\star$, and this field is also used in the
SALSA $\Gamma_{\mathrm{eff}}$ and optional `useRmod` path.

For that reason, the first local DDES implementation is conservative:

$$
\text{SALSA source/destruction terms}
\quad + \quad
\text{DDES } \tilde{d}
\quad + \quad
\text{standard or ZDES2020 shielding}.
$$

No sigma-production override is applied.

## Spalart-Allmaras QCR2020

The local QCR2020 variants,
`SpalartAllmarasQCR2020` and `SpalartAllmarasDDESQCR2020`, use the standard
OpenFOAM Spalart-Allmaras RANS and DDES transport equations and replace the
turbulent-stress relation by the QCR2020 nonlinear constitutive relation from
the NASA/TMR Spalart-Allmaras page:

$$
\tau_{ij,QCR2020}
=
\tau_{ij}
- C_{cr1}''(O_{ik}\tau_{jk}+O_{jk}\tau_{ik})
- C_{cr2}''\mu_t\sqrt{2W_{mn}W_{mn}}\delta_{ij}.
$$

The coefficients are

$$
C_{cr1}''=C_{cr1}'(1+C_{fw1}f_w),
\qquad
C_{cr2}''=C_{cr2}'(1+C_{fw2}f_w),
$$

with

$$
C_{cr1}'=0.20,\qquad
C_{cr2}'=\frac{1}{3a_1}=2.150537634408602,\qquad
a_1=0.155,\qquad
C_{fw1}=2.0,\qquad
C_{fw2}=0.3.
$$

The normalized rotation tensor is

$$
O_{ik}
=
\frac{2W_{ik}}
{\sqrt{
    \frac{\partial u_m}{\partial x_n}
    \frac{\partial u_m}{\partial x_n}
}},
\qquad
W_{ik}
=
\frac{1}{2}
\left(
    \frac{\partial u_i}{\partial x_k}
    -
    \frac{\partial u_k}{\partial x_i}
\right),
$$

with a small-gradient guard so QCR has no effect in zero-gradient regions.

For the QCR2020 coefficient functions, the standard SA wall function is used
with the smoother recommended replacement

$$
\Omega_s
=
\left[
    \frac{1}{2}
    \left(
        2W_{ij}W_{ij}+2S_{ij}S_{ij}
    \right)
\right]^{1/2}.
$$

Then

$$
\bar{S}
=
\frac{\tilde{\nu}}{\kappa^2d^2}f_{v2},
$$

and

$$
\hat{S}
=
\Omega_s+\bar{S}
\quad\mathrm{for}\quad
\bar{S}\ge -C_2\Omega_s,
$$

$$
\hat{S}
=
\Omega_s
+
\frac{
    \Omega_s(C_2^2\Omega_s+C_3\bar{S})
}{
    (C_3-2C_2)\Omega_s-\bar{S}
}
\quad\mathrm{for}\quad
\bar{S}< -C_2\Omega_s,
$$

where

$$
C_2=0.7,\qquad C_3=0.9.
$$

The resulting $r$, $g$, and $f_w$ are

$$
r=\min\left[
    \frac{\tilde{\nu}}{\hat{S}\kappa^2d^2},
    10
\right],
\qquad
g=r+c_{w2}(r^6-r),
$$

$$
f_w
=
g
\left[
    \frac{1+c_{w3}^6}{g^6+c_{w3}^6}
\right]^{1/6}.
$$

QCR2020 is independent of the local SALSA equations in this repository. The
implemented DDES variant combines QCR2020 with standard SA-DDES length-scale and
shielding controls because those do not replace the QCR2020 constitutive-stress
relation. QCR2020 is not compatible with another QCR variant in the same model.
If a future base model includes a true modeled $k$ contribution in the
Boussinesq relation, the final isotropic QCR2020 approximation should be omitted
to avoid double counting.
