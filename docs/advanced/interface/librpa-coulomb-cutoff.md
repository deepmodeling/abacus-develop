# Smooth Coulomb cutoff for LibRPA exports

ABACUS can export a smooth truncated Coulomb interaction for post-SCF RPA/GW
workflows. The default remains the existing hard spherical (Spencer) cutoff
with the legacy radius convention.

## Usage

In an existing LCAO RPA export calculation, use:

```text
rpa                     1
rpa_out_vel             1
rpa_outdir              ./OUT.librpa/
out_unshrinked_v        1

rpa_cut_coulomb         log_erfc
rpa_cut_radius          legacy
rpa_cut_width           1.091868692359
rpa_ccp_rmesh_times     10
```

Keep the structure, pseudopotentials, orbital/auxiliary bases and k mesh from
the intended calculation. The width and radial-range setting above were used
for the silicon comparison; they are not universal convergence settings.
The width must be finite and satisfy `1 < rpa_cut_width <= 1.25`.

Check standard output for `RPA cut Coulomb: mode=log_erfc radius=legacy` and
the actual `Rcut_bohr` and `width`. Then regenerate the downstream LibRPA
response and screened interaction. Consumers must use matching RI coefficients,
cut Coulomb matrices and screened interaction from the same cutoff workflow.
Changing only a downstream BSE input does not convert previously exported data.

The existing LibRPA-interface options in the example request velocity and
uncompressed Coulomb output for downstream optical and matrix diagnostics.
They do not enable smoothing by themselves.

To recover the hard-cutoff control, set `rpa_cut_coulomb hard` and
`rpa_cut_radius legacy`. Use separate output directories for each variant.

## Kernel and radius

The smooth radial kernel is

$$
v_{\mathrm{cut}}(r)=\frac{1}{2r}
\operatorname{erfc}\!\left(\frac{\ln(r/R_c)}{\ln w}\right).
$$

The implementation evaluates its Fourier factor as a positive Gaussian
average of hard-sphere factors. With $\sigma=\ln w$,

$$
\frac{k^2 v_{\mathrm{cut}}(k)}{4\pi}
=\frac{1}{\sqrt{\pi}}\int_{-\infty}^{\infty}e^{-t^2}
\left[1-\cos\!\left(kR_c e^{\sigma t}\right)\right]dt.
$$

This tends to the hard-sphere factor as $w\to1^+$. Its finite zero-wavevector
kernel limit is $2\pi R_c^2\exp(\sigma^2)$. The stored multiplicative factor
is zero at exactly $k=0$ because it represents $k^2v(k)/(4\pi)$.

`legacy` retains the existing radius rule. The optional `inscribed` convention
uses half the shortest nonzero Born-von Karman superlattice vector:

$$
R_c=\frac12\min_{\mathbf n\in\mathbb Z^3\setminus\{0\}}
\left|\sum_i n_i N_i\mathbf a_i\right|.
$$

Lattice scale is included and lengths are in Bohr. The shortest-vector search
includes integer combinations, not just the three supplied lattice vectors.
Radius selection and kernel smoothing are independent controls. The smooth
kernel has a tail beyond $R_c$; the inscribed radius is not a strict support
boundary. `rpa_ccp_rmesh_times` controls a numerical radial range, not this
physical radius. Check radial-range convergence for the smooth tail.

## Scope

The option is passed explicitly into the post-SCF truncated-Coulomb path,
including compressed and uncompressed cut exports. The full-Coulomb kernel
definition is unchanged. RI fits using the cut metric also change, so this
is not generally a fixed-response replacement of the outer Coulomb factors
in the screened interaction. Downstream Fourier folding and interpolation
rules are unchanged.

The default width was taken from a silicon FHI-aims log; this implementation
does not reproduce the complete FHI-aims short-/long-range split. Forces and
stresses for the smooth kernel have not been validated.

Focused kernel tests are registered as `MODULE_RI_smooth_cut_coulomb` in an
MPI-enabled LCAO/LibRI test build. They cover the small-wavevector moment,
hard-cutoff limit, positivity, invalid arguments and cubic/fcc/skew-cell
radius selection. INPUT validation tests also cover the new switches.
