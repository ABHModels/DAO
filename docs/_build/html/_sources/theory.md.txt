(theory)=
# Theory 


## Radiative transfer

Radiative transfer is a fundamental problem in astrophysics, particularly in the study of stellar and accretion disk atmospheres. The basic form of the radiative transfer equation in plane-parallel geometry can be expressed as:

$$
\mu\frac{\partial I(z,\mu,\nu)}{\partial z} = \chi(z,\mu,\nu) I(z,\mu,\nu) -\eta(z,\mu,\nu)
$$

where $z$, $\mu$, and $\nu$ denote the specific depth, angle cosine, and frequency (or energy), respectively. Here, $I$ is the specific intensity, $\chi$ is the opacity, and $\eta$ is the emissivity. Photons traveling through the plasma interact with various elements, resulting in absorption and emission.

Leveraging the sophisticated architecture of `XSTAR`, `DAO` accounts for a comprehensive range of atomic physics. It includes nearly all relevant processes, such as free-free, bound-free, three-body recombination, collisional excitation, and radiative recombination. These processes determine the opacity and emissivity terms in the radiative transfer equation, with Compton scattering being treated as a distinct component.

When constructing the model, it is convenient to represent this equation using optical depth, defined as $d\tau(E) = \chi(E)dz$. Consequently, Eq. {eq}eq-rte can be rewritten as:

$$
\mu \frac{\partial I(\tau_T,\mu,\nu)}{\partial \tau(\nu)} = I(\tau_T,\mu,\nu) - S(\tau_T,\mu,\nu)
$$

(Note: We omit the explicit dependence on optical depth for brevity).

It is important to note that in Eq. {eq}`eq-rte-tau`, the radiation field depends on the Thomson optical depth $\tau_T$ (where $d\tau_T = n_e\sigma_T dz$), rather than the frequency-dependent optical depth $\tau$. We adopt this convention because the model calculations are performed on a fixed depth grid, and $\tau_T$ provides the most stable reference for determining these grid points. Similar approaches are found in other radiative transfer models, such as `reflionx`.There are various methods to solve the radiative transfer equation, including approximate, direct, and iterative methods. Considering both accuracy and computational efficiency, we employ the Lambda-iteration method combined with the Feautrier method (hereafter FTM).

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

where A,B,C,L are diagonal in the problem we considering. Equation \ref{eq:tridiagonal}-\ref{eq:tridiagonal_bot} could be solved by forward-elimination and back-substitution procedure. Once it's been solved, the new radiation field will been scattered and update the source function, which named $\Lambda$-iteration.

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

In general, $\mathbf{T}$ is a full matrix, while $\mathbf{x}$ and $\mathbf{b}$ are vectors. When the complete redistribution function $R_{\nu,\nu'}$ is considered, the problem becomes considerably more complex, as a $200\times5000\times5000$ matrix is introduced. At each depth point, computing the inverse of a $5000\times5000$ matrix would be necessary, scaling with a computational complexity of $O(n^3)$. In contrast, the method implemented in DAO scales linearly with $O(n)$ complexity, although it requires more iterations to achieve convergence.



### Boundary condition

It should be emphasized that the second-order radiative transfer equation (Eq.~\ref{eq:second-order rte}) requires the specification of the intensity as a boundary condition.
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