(theory)=
# Theory 


## Radiative transfer

Radiative transfer is a fundamental problem in astrophysics, particularly in the study of stellar and accretion disk atmospheres. The basic form of the radiative transfer equation in plane-parallel geometry can be expressed as:

$$
\mu\frac{\partial I(z,\mu,\nu)}{\partial z} = \chi(z,\mu,\nu) I(z,\mu,\nu) -\eta(z,\mu,\nu)
$$

where $z$, $\mu$, and $\nu$ denote the specific depth, angle cosine, and frequency (or energy), respectively. Here, $I$ is the specific intensity, $\chi$ is the opacity, and $\eta$ is the emissivity. Photons traveling through the plasma interact with various elements, resulting in absorption and emission.

Leveraging the sophisticated architecture of `XSTAR`, `DAO` accounts for a comprehensive range of atomic physics. It includes nearly all relevant processes, such as free-free, bound-free, three-body recombination, collisional excitation, and radiative recombination. These processes determine the opacity and emissivity terms in the radiative transfer equation, with Compton scattering being treated as a distinct component.

When constructing the model, it is convenient to represent this equation using optical depth, defined as $d\tau(E) = \chi(E)dz$. Consequently, transfer equation can be rewritten as:

$$
\mu \frac{\partial I(\tau_T,\mu,\nu)}{\partial \tau(\nu)} = I(\tau_T,\mu,\nu) - S(\tau_T,\mu,\nu)
$$

(Note: We omit the explicit dependence on optical depth for brevity).

The radiation field depends on the Thomson optical depth $\tau_T$ (where $d\tau_T = n_e\sigma_T dz$), rather than the frequency-dependent optical depth $\tau$. We adopt this convention because the model calculations are performed on a fixed depth grid, and $\tau_T$ provides the most stable reference for determining these grid points. Similar approaches are found in other radiative transfer models, such as `reflionx`.There are various methods to solve the radiative transfer equation, including approximate, direct, and iterative methods. Considering both accuracy and computational efficiency, we employ the Lambda-iteration method combined with the Feautrier method (hereafter FTM).

### Feautrier Method

FTM is based on second order of radiative transfer equation. First if we consider symmetric and antisymmetric averages of the specific intensity at frequency $\nu$ propagating along a ray in the $\pm\mu$ direction are

$$
    u_{\mu\nu} = \frac{1}{2} \left[I_\nu(\mu) + I_\nu(-\mu)\right] \quad \mathrm{and} \quad v_{\mu\nu} = \frac{1}{2} \left[I_\nu(\mu) - I_\nu(-\mu)\right] 
$$

The transfer equation for two directions are 

$$
    \mu \frac{\partial I_\nu(\mu)}{\partial \tau_\nu} = I_\nu(\mu) - S_\nu(\mu) \quad \mathrm{and} \quad  \frac{\partial I_\nu(-\mu)}{\partial \tau_\nu} = I_\nu(-\mu) - S_\nu(-\mu)
$$

Adding above two equations, 

$$
    \mu\frac{\partial v_{\mu\nu}}{\partial\tau_\nu} = u_{\mu\nu}-S_\nu \quad\mathrm{and}\quad \mu\frac{\partial u_{\mu\nu}}{\partial\tau_\nu} = v_{\mu\nu}
$$

Then substituting the expression for $v_{\mu\nu}$, we obtain \textit{Schuster’s second-order transfer equation {cite:p}`1905ApJ....21....1S`.

$$
    \mu^2 \frac{\partial^2 u_{\mu\nu}}{\partial\tau_\nu^2} = u_{\mu\nu}-S_\nu
$$

where $0\leq\mu\leq1$.

To numerically solve that, the discrete form could be get by finite difference. We use m,n,d represent the angle, energy, and depth points.

$$
    \Delta\tau_{d+\frac{1}{2},n}\equiv\tau_{d+1,n}-\tau_{d,n}
$$

Then the transfer funciton:

$$
    \frac{\mu^2 u_{d-1}}{\Delta\tau_{d-\frac{1}{2}}\Delta\tau_{d}} - \frac{\mu^2u_{d}}{\Delta\tau_d}\left(\frac{1}{\Delta\tau_{d-\frac{1}{2}}}+\frac{1}{\Delta\tau_{d+\frac{1}{2}}}\right)+\frac{\mu^2u_{d+1}}{\Delta\tau_{d+\frac{1}{2}}\Delta\tau_{d}} = u_d-S_d \quad (d=2,...,ND-1)
$$

where we omit the subscript of angle and energy for convenient and ND represent the maximum number of depth grids. In this expression, one need to take care for source function:

$$
    S_{n} = \frac{j_{n}}{\alpha_n+\alpha^c_n}\quad,\quad j_n = j^\alpha_n + \alpha^c_nJ^c_n
$$

which angle-independent. In the denominator, $\alpha$ and $\alpha^c$ represent opacity of absorption and scattering, respectively. And $j^\alpha_n $ represent the emissivity include all possible atomic process, $J^c$ represent the scattered mean-intensity, 

$$
    J^c_n = \sum_{n'} R_{n,n'} J_{n'}w_{n'}
$$

where $R_{n,n}$ is redistribution function, $J_{n'}$ is mean-intensity before scattering and $w_{n'}$ is integral weight. $R_{n,n'}$ describe the probablity of energy redistributoin in Compton scattering process. Redistribution we used in `DAO` is calculated by complete quantum mechanism and suitable for relativistic electrons, we will discuss later in [Compton scattering](comptonscattering)

We consider external illumination at both of upper and bottom surface. At these boundaries, the second exact form from {cite:t}`1967ApJ...150L..53A` is used:

$$
    \mu \frac{u_2-u_1}{\Delta \tau_{3/2}} = u_1+\Delta\tau_{3/2}(u_1-S_1)/2\mu +I_{\mathrm{inc}}
$$

and the same expression could be found for lower surface.

The discretized radiative transfer equation forms a block tridiagonal system:

$$
-A_du_{d-1} + B_du_d -C_du_{d+1} = L_d\\
B_1u_1 -C_1u_{2} = L_1\\
-A_{n_\tau} u_{{n_\tau-1}} + B_{n_\tau}u_{n_\tau} = L_{n_\tau}
$$

where A,B,C,L are diagonal in the problem we considering. Equations above could be solved by forward-elimination and back-substitution procedure. Once it's been solved, the new radiation field will been scattered and update the source function, which named $\Lambda$-iteration.

The complete iterative procedure is summarized as follows:

* Obtain the initial estimates of the emissivity and absorption coefficients using a naive radiation field.
* Re-Compute the source function using.
* Apply the Feautrier method to obtain an updated radiation field based on the source function from step~2.
* Repeat steps~2–3 until convergence is achieved.


The convergence criterion adopted in `DAO` is defined as

$$
E_\tau\!\left[\int \frac{S_n(E)-S_{n-1}(E)}{S_{n-1}(E)}\,dE\right] < \varepsilon, \qquad \varepsilon = 10^{-10}
$$

which ensures that the source function ceases to change between successive iterations. 

While more advanced iteration processes such as Accelerated Lambda Iteration (ALI) exist, we did not implement them due to the prohibitive memory costs associated with a large number of energy grids. ALI requires solving linear equations of the form:

$$
\mathbf{T}\cdot \mathbf{x} = \mathbf{b}
$$

In general, $\mathbf{T}$ is a full matrix, while $\mathbf{x}$ and $\mathbf{b}$ are vectors. When the complete redistribution function $R_{\nu,\nu'}$ is considered, the problem becomes considerably more complex, as a $200\times5000\times5000$ matrix is introduced. At each depth point, computing the inverse of a $5000\times5000$ matrix would be necessary, scaling with a computational complexity of $O(n^3)$. In contrast, the method implemented in `DAO` scales linearly with $O(n)$ complexity, although it requires more iterations to achieve convergence.



### Boundary condition

Second-order radiative transfer equation requires the specification of the intensity as a boundary condition.
For the upper surface, the boundary condition can be obtained from the first moment of the radiation field:

$$
F(E) = \int_0^1 u(E,\mu)\mu d\mu
$$

Using $I^+(E) = I_\mathrm{inc}\delta(\mu-\mu_\mathrm{inc})$, one can then derive the incident specific intensity at the incidence angle $\mu_\mathrm{inc}$ as

$$
I_\mathrm{inc}(E) = \frac{2F(E)}{\mu_\mathrm{inc}}.
$$

For the lower surface, the reflection region receives isotropic thermal emission from the disk, which can be simply expressed as $I_\mathrm{bot}(E) = B(E)$, where $B(E)$ is the blackbody spectrum corresponding to the disk temperature. 

In particular, if the diffusion approximation is valid at $\tau_\mathrm{max}$,then

$$
    I_\mathrm{bot}(E,\mu) = B(E) + \mu \frac{\partial B}{\partial \tau}_{\tau=\tau_\mathrm{max}}
$$

## Compton scattering

Essentially, Compton scattering is a stochastic process that alters both the propagation
direction and the energy of photons, and it should be described using the exact quantum
mechanical formalism. In addition, in high temperature gas, motion of relativistic elec-
trons follow the Maxwellian distribution, which need to be taken into account in complete
redistribution function.

Historically, many authors had investigated such unpolarized and quantum mechanism  scattering redistribution function: {citt:t}`1981MNRAS.197..451G` incorporated the relativistic Maxwellian velocity distribution for electrons into the redistribution function, but a computational error was present in the original paper. Subsequently, {cite:t}`1993AstL...19..262N,1996ApJ...470..249P,2010ApJS..189..286P` investigated the  redistribution function between photons and relativistic electrons using quantum electrodynamical methods. Later, {cite:t}`2017MNRAS.469.2032M` corrected the earlier error of Guilbert and demonstrated consistency among these methods. Recently, {cite:t}`2020ApJ...897...67G` showed the large difference between exactly redistribution function and Gaussian-approximated. The former is derived from {cite:t}`1993AstL...19..262N` and the later is been widely use in reflection model {cite:t}`1978ApJ...219..292R,2005MNRAS.358..211R,2000ApJ...537..833N`.

In `DAO` model, the python parameter `kernelpath` (`rfcomp_file` in fortran) allow user to set any type of redistribution function they want. Considering dependence of accuracy for physic probelm, we highly recommend user generate the exact redistribution function by using script in `get_redistribuion`, where we use the code of `2020ApJ...897...67G` to calculate the second-order exact redistribution function as described in `1993AstL...19..262N` and use formula which described in `1978ApJ...219..292R` to caculate the Guassian redistribution function. The calculation routine as follows:

### Exact redistribution

$$
R_E(x,x_1,\mu,\gamma)= \frac{2}{Q} + \frac{u}{v} \left(1-\frac{2}{q}\right) + u \frac{(u^2-Q^2)(u^2+5v)}{2q^2v^3}+u\frac{Q^2}{q^2 v^2},
$$

with $x$ and $x_1$ the dimensionless photon energies before and after scattering, respectively, $\mu$ is the cosine of the scattering angle, and $\gamma$ represents electron Lorentz factor. Here, $q = xx_1 (1-\mu)$ and $Q^2=(x-x_1)^2+2q$, while


$$
a_-^2 = (\gamma-x)^2 + \frac{1+\mu}{1-\mu}, \qquad
a_+^2 = (\gamma+x_1)^2 + \frac{1+\mu}{1-\mu}
$$

$$
v = a_-a_+, \qquad
u = a_+ - a_-
$$

The redistribution function averaged over a relativistic Maxwellian distribution is

$$
R(x,x_1,\mu) = \frac{3}{32\mu\Theta K_2(1/\Theta)}\int _{\left(x-x_1+Q\sqrt{1+2/q}\right)/2}^\infty R(x,x_1,\mu,\gamma)\exp(-\gamma/\Theta)d\gamma
$$

where $K_2$ is the modified Bessel function of the second kind and $\Theta=kT/m_ec^2$. Currently we use angle-averaged redistribution so the resulting function is 

$$
    R(x,x_1) = \int_{\mu_{\mathrm{min}}}^{\mu_\mathrm{max}}R(x,x_1,\mu)d\mu
$$

### Gaussian redistribution function

$$
P(E_i,E_f) = \frac{1}{\sqrt{2\pi}\Sigma}
\exp\left[-\frac{(E_f-E_c)^2}{2\Sigma^2}\right]
$$

with central energy $E_c$ and standard deviation $\Sigma$ defined as

$$
E_c = E_i\left(1 + \frac{4kT}{m_ec^2} - \frac{E_i}{m_ec^2}\right)
$$

$$
\Sigma = E_i\left[\frac{2kT}{m_ec^2} + \frac{2}{5}\left(\frac{E_i}{m_ec^2}\right)^2\right]^{1/2}
$$

```{figure} images/kernel.svg
:name: fig-kernel
:width: 100%
:align: center

Gaussian (red dash line) and Exact redistribution (orange solid line) function with three initial energy: E$_i$ = 1, 10, 100 keV and three different gas temperature: 10$^7$, 10$^8$, 10$^9$ K.
```

In the Figure above, we shows the large deviation between exactly redistribution function and Gaussian-approximated function when gas temperature is high and photons have large energy. Same results have been reported in {cite:t}`2020ApJ...897...67G`. We show here just for convenient for readers to see differences. 

### Cross section

Klein-Nishina cross section is usually used in non-relativistic scattering process, where

$$
    \sigma(x) = \sigma_T \frac{3}{4}\left\{\frac{1+x}{x^3}\left[\frac{2x(1+x)}{1+2x}-\ln{\left(1+2x\right)}\right]+\frac{1}{2x}\ln{\left(1+2x\right)}-\frac{1+3x}{\left(1+2x\right)^2}\right\}
$$
where $x$ is dimensionless energy $x=E/m_ec^2$. However, electrons in hot plasma follow relativistic Maxwellian distribution:

$$
\sigma_{CS}(x) = 
\frac{3\sigma_T}{16x^2\Theta K_2(1/\Theta)}
\int_1^\infty e^{-\gamma/\Theta}
\Biggl\{ &
\left(x\gamma+\frac{9}{2}+\frac{2\gamma}{x}\right)
\ln{\left[\frac{1+2x(\gamma+z)}{1+2x(\gamma-z)}\right]} - 2xz \\
& + z\left(x-\frac{2}{x}\right)\ln(1+4x\gamma+4x^2)
+\frac{4x^2z(\gamma+x)}{1+4x\gamma+4x^2} \\
& -2\int_{x(\gamma-z)}^{x(\gamma+z)} 
\ln(1+2\varepsilon)\frac{d\varepsilon}{\varepsilon}
\Biggr\} d\gamma
$$

In `DAO`, we use Compton scattering cross section averaged over a relativistic Maxwellian electron distribution which following {cite:t}`1996ApJ...470..249P`, instead of the  Klein–Nishina cross section. In certain limiting cases (e.g., $\Theta \ll 1$), the thermal cross section can be approximated using simpler expressions (see {cite:t}`1982ApJ...258..321S` for more details). {cite:t}`2020ApJ...897...67G` showed the difference between the Klein–Nishina cross section and cross section we mentioned here.

We have not yet included the effects of induced scattering, which are planned for future implementation. The electron scattering opacity accounting for the induced scattering is 

$$
    \sigma(x,\mu) = \frac{\sigma_T}{x} \int_0^\infty x_1dx_1\int_{-1}^1 d\mu_1 R(x_1,\mu_1,x,\mu)\left(1+\frac{CI(x_1,\mu_1)}{x_1^3}\right),\quad C=\frac{1}{2m_e}\left(\frac{h}{m_ec^2}\right)^3
$$
which add one more dimension $\mu$ for source function.

### Normalization of Redistribution Functions

Integration of redistribution should give the scattering cross section. Therefore we normalize redistribution function like 

$$
R_n(E_i,E_f) = \sigma_\nu\frac{R(E_i,E_f)}{\int_0^\infty dE_iR(E_i,E_F)}
$$

### Scattered radiation field

In expression of source function, the scattered mean-intensity $J^c$ is 

$$
 J^c(E) = \int_{E_{min}}^{E_{max}} R(E',E)J(E')dE'
$$
Currently, angular redistribution is not included in `DAO`. If it were, scattered intensity would take the form

$$
I(\mu,E) = \int_{E_{min}}^{E_{max}} dE'\int_{-1}^1 d\mu' R(E';\mu',E;\mu)I(\mu',E')
$$

where $E'$ and $\mu'$ denote the photon energy and direction cosine after scattering and $E$, $\mu$ denote the energy and direction cosine before scattering. In this case, the redistribution becomes more complex, as it depends on both the incident and scattered energies and angles. If one considers the scattering angle $\eta$ (defined as $\eta = \boldsymbol{\mu} \cdot \boldsymbol{\mu}'$), the redistribution function becomes a function of $\eta$, $E'$, and $E$, which is more convenient for numerical calculations. In summary, the inclusion of induced scattering and angle-dependent redistribution makes the source function angle-dependent. To date, there has been no research demonstrating whether such angle dependence is important in X-ray reprocessing. However, this aspect is worth investigating, especially when polarization effects are taken into account.

In addition, it is important to note that we use the inverse redistribution function $R(E',E)$ rather than $R(E,E')$. This is because we are interested in the probability that a photon with an initial energy $E'$ is scattered into the current energy $E$. The redistribution function satisfies the detailed balance condition {cite:p}`1973erh..book.....P`:

$$
R(E',E) = R(E,E') \left(\frac{E}{E'}\right)^{2}
\exp\left(\frac{E' - E}{kT_e}\right)
$$

## Atom Physics

`DAO` is fully based on `XSTAR`. All atomic processes, including photoionization, thermal equilibrium, and level population calculations, are derived from `XSTAR`. We employ the latest version (v2.59) together with the updated atomic database released in 2024. For more details, please refer to the [[XSTAR manual](https://heasarc.gsfc.nasa.gov/xstar/xstar.html)

## Algorithm

Convergence must be achieved through an iterative procedure. We begin by assuming that the initial radiation field at each depth is equal to the coronal illumination. Then, we use `XSTAR` to calculate the emissivity and opacity at each depth (excluding scattering in this step). The `XSTAR` calculations enforce thermal equilibrium, satisfying

$$
H=C
$$

where $H$ is the heating rate and $C$ is the cooling rate. Once the plasma reaches thermal equilibrium under the current radiation field, we assume that it remains in this equilibrium state during the radiative-transfer calculation. In other words, the plasma temperature, emissivity, and absorption coefficients are held fixed in this stage.

The solution to the radiative transfer equation is obtained through a standard $\Lambda$-iteration scheme as follows:

* Using the current temperature, emissivity, and absorption, we compute the source function at each Thomson optical depth $\tau_T$.
* We solve the second-order radiative transfer equation using the Feautrier method.
* We update the source function based on the new radiation field obtained at step~2 and use this new source function to update the radiation field.

This iteration loop stops when the source function converges:

$$
\int\frac{S^{n}(E)-S^{n-1}(E)}{S^{n-1}(E)}dE <e_n,\quad e_n = 10^{-10}
$$

After the radiation field converges in the radiative-transfer step, we again use this updated radiation field as the input to `XSTAR` to obtain a new plasma profile. We then solve the radiative transfer equation again. The iterative procedure continues until both the temperature and ionization profiles at each layer no longer change. When `DAO` calculates low-ionization models, photoionization is the dominant process, with $\alpha \gg \alpha^c$. In this regime, convergence are more rapidly than in a scattering-dominated gas. In contrast, for high-ionization models, the gas is scattering-dominated, with $\alpha^c \gg \alpha$, leading to slower convergence.
