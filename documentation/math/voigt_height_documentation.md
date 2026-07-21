# Height-Normalized Voigt Profile — Mathematical Documentation

## 1. Overview

This document derives and explains the `voigt_height` function used in the
double-Voigt spectral fitting pipeline for exciton/trion PL peaks:

```python
def voigt_height(x, height, cen, sigma, gamma):
    z  = ((x - cen) + 1j*gamma) / (sigma * np.sqrt(2))
    z0 = (1j*gamma) / (sigma * np.sqrt(2))
    return height * np.real(wofz(z)) / np.real(wofz(z0))
```

It computes a **Voigt line shape parameterized directly by peak height**,
rather than by peak area (the more common convention). The derivation below
proceeds from the physical origin of the two component line shapes, through
the convolution integral, to the closed-form expression via the Faddeeva
function, and finally to the height-normalization step implemented in code.

---

## 2. Physical origin of the two components

**Gaussian component** — inhomogeneous broadening. Arises from
sample-to-sample or point-to-point variation in the transition energy
(e.g. local strain, dielectric environment, disorder across a TMD
monolayer). If the local center energy is normally distributed with
standard deviation $\sigma$:

$$
G(x) = \frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{(x-x_0)^2}{2\sigma^2}\right)
$$

**Lorentzian component** — homogeneous broadening. Arises from the finite
coherence lifetime of the transition (radiative + non-radiative decay),
identical for every emitter. With half-width-at-half-maximum $\gamma$:

$$
L(x) = \frac{1}{\pi}\,\frac{\gamma}{(x-x_0)^2+\gamma^2}
$$

---

## 3. The convolution

Each emitter contributes an *identical* Lorentzian, but centered at a
randomly shifted position $x_0 + \delta$, where $\delta$ is drawn from the
Gaussian distribution above. The ensemble-averaged (measured) line shape
is therefore the convolution of the two:

$$
V(x) = (G * L)(x) = \int_{-\infty}^{\infty} G(\delta)\, L(x - x_0 - \delta)\, d\delta
$$

This integral has **no elementary closed form** — it cannot be written in
terms of exponentials, rational functions, or standard error functions
alone. Evaluating it directly would require numerical integration (or an
FFT-based convolution) at every fitting iteration, which is computationally
inefficient and less numerically stable near the peak.

---

## 4. Closed form via the Faddeeva function

Define the complex argument:

$$
z = \frac{(x - x_0) + i\gamma}{\sigma\sqrt{2}}
$$

The **Faddeeva function** (also called the complex probability function /
plasma dispersion function) is defined as:

$$
w(z) = e^{-z^2}\,\operatorname{erfc}(-iz)
$$

It is a classical result (originally from astrophysical and plasma line
theory) that the real part of $w(z)$ evaluated at this specific complex
argument gives the Gaussian–Lorentzian convolution exactly:

$$
V(x) = \frac{1}{\sigma\sqrt{2\pi}}\,\operatorname{Re}\big[w(z)\big]
$$

This is the **area-normalized Voigt profile** (integrates to 1 over all
$x$). The real part of $z$ carries the detuning $(x - x_0)$; the imaginary
part carries the Lorentzian width $\gamma$. The function $w(z)$ jointly
encodes both the Gaussian decay and the Lorentzian's algebraic tails in a
single analytic expression.

In code, `scipy.special.wofz(z)` implements $w(z)$ numerically (via
Humlíček's rational/continued-fraction algorithm), which is accurate and
efficient across the full range of $\sigma/\gamma$ ratios — from
Gaussian-dominated ($\gamma \ll \sigma$) to Lorentzian-dominated
($\gamma \gg \sigma$) regimes.

---

## 5. Height normalization

Area normalization is convenient mathematically, but for fitting purposes
a **peak-height** parameter is often preferable: it decorrelates more
cleanly from $\sigma$ and $\gamma$ during nonlinear least-squares
optimization (`curve_fit`), whereas an area parameter is entangled with
both widths.

To convert from area-normalized to height-normalized, evaluate the
area-normalized profile at the peak center ($x = x_0$, i.e. $\text{Re}(z) = 0$):

$$
z_0 = \frac{i\gamma}{\sigma\sqrt{2}}
$$

$$
V(x_0) = \frac{1}{\sigma\sqrt{2\pi}}\,\operatorname{Re}\big[w(z_0)\big]
$$

Dividing the profile by this peak value and multiplying by the desired
`height` parameter rescales the curve so that $V(x_0) = \text{height}$
exactly, for any $\sigma, \gamma$:

$$
V_{\text{height}}(x) = \text{height} \cdot \frac{\operatorname{Re}[w(z)]}{\operatorname{Re}[w(z_0)]}
$$

This is exactly what the function computes. Note that the
$\dfrac{1}{\sigma\sqrt{2\pi}}$ prefactor cancels between numerator and
denominator — it never needs to be computed explicitly, which is why the
code only calls `wofz` twice and takes a ratio.

---

## 6. Limiting cases (sanity checks)

**Pure Gaussian limit ($\gamma \to 0$):**
$z_0 \to 0$, and $w(0) = 1$, so $\operatorname{Re}[w(z_0)] \to 1$. The
expression reduces to $V_{\text{height}}(x) = \text{height} \cdot
\operatorname{Re}[w(z)]$, which correctly recovers a height-normalized
Gaussian.

**Pure Lorentzian limit ($\sigma \to 0$):** both $z$ and $z_0$ have
magnitude $\to \infty$ along the imaginary axis. $w(z)$ has a known
asymptotic form for large $|z|$ that correctly recovers Lorentzian
behavior, but this limit is numerically delicate — $\sigma$ approaching
zero during optimization can produce very large arguments to `wofz` and
should be guarded against with fitting bounds (e.g. a small positive lower
bound on $\sigma$).

---

## 7. Summary of parameter roles

| Parameter | Physical meaning | Role in $z$ |
|---|---|---|
| `cen` ($x_0$) | Peak center (transition energy) | Real part of numerator |
| `sigma` ($\sigma$) | Gaussian (inhomogeneous) width | Denominator scale |
| `gamma` ($\gamma$) | Lorentzian (homogeneous) HWHM | Imaginary part of numerator |
| `height` | Peak amplitude | Overall scale factor, fixed exactly at $x=x_0$ |

---

## 8. Physical interpretation of the height, cen, sigma, and gamma parameters

Each fit parameter maps onto a distinct physical mechanism governing
exciton photoluminescence in monolayer TMDs (e.g. MoSe₂, WSe₂) at room
temperature.

### `cen` — transition energy

The peak center tracks the free-exciton (or trion) transition energy. At
elevated temperature it red-shifts due to thermal band-gap renormalization
(electron–phonon renormalization plus lattice thermal expansion), commonly
described by a Varshni-type or Bose–Einstein-type temperature dependence
of the gap. This shift is well documented from cryogenic to room
temperature in monolayer MoSe₂ (Tongay et al., 2012; Horng et al., 2018).

### `height` — peak amplitude / relative oscillator strength

The height (in a height-normalized model, as opposed to an integrated-area
model) tracks the peak PL intensity, which is influenced by the exciton
population, radiative quantum yield, and — for the exciton/trion doublet
— the free-carrier density that governs the exciton–trion branching ratio.
At room temperature, thermal population of excitons across the light cone
and strong exciton–phonon scattering reduce the radiative efficiency
relative to cryogenic conditions; first-principles calculations of
intrinsic radiative lifetimes show a few-ps lifetime at low temperature
increasing to a few-ns lifetime at room temperature in monolayer TMDs,
consistent with a thermally averaged population spread outside the light
cone (Palummo, Bernardi & Grossman, *Nano Lett.* 2015, "Exciton Radiative
Lifetimes in Two-Dimensional Transition Metal Dichalcogenides"). The
trion peak height additionally encodes free-carrier (doping) density,
since trion formation requires a resident charge to bind to the exciton
(Mak et al., *Nat. Mater.* 2013, "Tightly bound trions in monolayer
MoS₂").

### `sigma` — inhomogeneous (Gaussian) broadening

σ represents *inhomogeneous* broadening: a static, ensemble/spatial
distribution of local transition energies, caused by substrate-induced
disorder, dielectric inhomogeneity, local strain fields, and defect
density variations across the flake. Encapsulation studies show that this
contribution can be strongly suppressed with high-quality, flat
hBN-encapsulated samples (down to ~2 meV FWHM at 4 K), directly
identifying substrate/surface disorder as the dominant source of σ
(Cadiz et al., *Phys. Rev. X* 2017, "Excitonic Linewidth Approaching the
Homogeneous Limit in MoS₂-Based van der Waals Heterostructures"). This is
the mechanism most directly responsible for the spatial inhomogeneity you
observe in σ across the hyperspectral PL map — it reflects real,
point-to-point variation in local dielectric/strain environment rather
than a measurement artifact.

### `gamma` — homogeneous (Lorentzian) broadening

γ represents *homogeneous* broadening: the intrinsic dephasing rate of
the exciton coherence, set by population relaxation (radiative +
non-radiative decay) and pure dephasing from exciton–phonon and
exciton–exciton scattering. At low temperature the homogeneous linewidth
in WSe₂ has been measured directly (via 2D Fourier-transform spectroscopy)
to be nearly two orders of magnitude narrower than the inhomogeneous
width, with a residual (population-relaxation-limited) linewidth of
~1.6 meV corresponding to a coherence time T₂ ≈ 0.4 ps (Moody et al.,
*Nat. Commun.* 2015, "Intrinsic homogeneous linewidth and broadening
mechanisms of excitons in monolayer transition metal dichalcogenides").
Crucially, γ grows strongly with temperature because acoustic- and
optical-phonon scattering rates increase with phonon occupation: at room
temperature, exciton–phonon coupling dominates the total linewidth, with
measured total linewidths increasing from ~19 meV at 4 K to ~33 meV at
300 K in monolayer MoSe₂ (Shree et al., *Phys. Rev. B* 2018, "Observation
of exciton-phonon coupling in MoSe₂ monolayers"; see also Selig et al.,
*Nat. Commun.* 2016, "Excitonic linewidth and coherence lifetime in
monolayer transition metal dichalcogenides", for microscopic theory of
the exciton–phonon dephasing rate and its temperature/valley dependence).
This is why, at room temperature, the Lorentzian (homogeneous) component
of the Voigt profile is expected to be significant or even dominant
relative to the Gaussian component — the reverse of the low-temperature
regime, where inhomogeneous broadening dominates.

### Summary table

| Parameter | Broadening type | Dominant physical mechanism at RT | Key reference |
|---|---|---|---|
| `cen` | — | Band-gap renormalization / thermal shift | Tongay et al. 2012; Horng et al. 2018 |
| `height` | — | Radiative yield, exciton/trion population, doping | Palummo et al. 2015; Mak et al. 2013 |
| `sigma` | Inhomogeneous (Gaussian) | Static disorder: strain, dielectric environment, substrate roughness | Cadiz et al. 2017 |
| `gamma` | Homogeneous (Lorentzian) | Exciton–phonon scattering, population relaxation | Moody et al. 2015; Shree et al. 2018; Selig et al. 2016 |

### References

1. Moody, G. et al. "Intrinsic homogeneous linewidth and broadening mechanisms of excitons in monolayer transition metal dichalcogenides." *Nat. Commun.* 6, 8315 (2015), DOI: [10.1038/ncomms9315](https://www.nature.com/articles/ncomms9315).
2. Cadiz, F. et al. "Excitonic Linewidth Approaching the Homogeneous Limit in MoS₂-Based van der Waals Heterostructures." *Phys. Rev. X* 7, 021026 (2017), DOI: [10.1103/PhysRevX.7.021026](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.7.021026).
3. Shree, S. et al. "Observation of exciton-phonon coupling in MoSe₂ monolayers." *Phys. Rev. B* 98, 035302 (2018), DOI: [10.1103/PhysRevB.98.035302](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.98.035302).
4. Selig, M. et al. "Excitonic linewidth and coherence lifetime in monolayer transition metal dichalcogenides." *Nat. Commun.* 7, 13279 (2016), DOI: [10.1038/ncomms13279](https://www.nature.com/articles/ncomms13279).
5. Palummo, M., Bernardi, M. & Grossman, J. C. "Exciton Radiative Lifetimes in Two-Dimensional Transition Metal Dichalcogenides." *Nano Lett.* 15, 2794–2800 (2015), DOI: [10.1021/nl503799t](https://pubs.acs.org/doi/10.1021/nl503799t).
6. Mak, K. F. et al. "Tightly bound trions in monolayer MoS₂." *Nat. Mater.* 12, 207–211 (2013), DOI: [10.1038/nmat3505](https://www.nature.com/articles/nmat3505).
7. Tongay, S. et al. "Thermally driven crossover from indirect toward direct bandgap in 2D semiconductors: MoSe₂ versus MoS₂." *Nano Lett.* 12, 5576–5580 (2012), DOI: [10.1021/nl3026357](https://pubs.acs.org/doi/10.1021/nl3026357).
8. Christiansen, D. et al. "Phonon sidebands in monolayer transition metal dichalcogenides." *Phys. Rev. Lett.* 119, 187402 (2017), DOI: [10.1103/PhysRevLett.119.187402](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.187402).

---
