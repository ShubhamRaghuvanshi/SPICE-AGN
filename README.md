# SPICE-AGN

**SPICE-AGN** extends the SPICE simulation framework by incorporating supermassive black hole (SMBH) seeding, growth, and AGN feedback into the RAMSES-RT infrastructure used for epoch-of-reionization (EoR) studies.

The project is designed to investigate how AGN mechanical feedback reshapes the ISM/CGM of high-redshift galaxies and regulates the escape of ionizing radiation during cosmic dawn.

---

## Core Goals

- Integrate OBELISK-style BH physics into SPICE.
- Isolate AGN-driven structural changes via controlled ablation experiments.
- Quantify porosity-driven escape of stellar ionizing radiation.
- Maintain strict numerical consistency with the SPICE reference branch.

---

## Build

Requirements:
- RAMSES-RT (SPICE-compatible branch)
- MPI toolchain
- Python for analysis

Clone and compile:
git clone https://github.com/ShubhamRaghuvanshi/SPICE-AGN.git
cd SPICE-AGN
make clean
make


Optional Makefile flags:
- `SINKFLAG = TRUE` (enable SMBH module)
- `AGNRT = TRUE` (future AGN radiation support)

---

## Experimental Design

| Run | BH | AGN FB | Purpose |
|-----|----|--------|---------|
| G | No | No | Stellar-only baseline |
| G+BH | Yes | No | BH presence only |
| G+BH+AGNfb | Yes | Yes | Full AGN mechanical pathway |

Escape diagnostics use neutral column ray-casting and open-sky fraction \(f_\Omega\) as a structural proxy.

---

## Status

SPICE-AGN preserves numerical identity with SPICE when SMBH physics is disabled. Any differences arise from physical model activation.

---

## Citation

Raghuvanshi et al., *SPICE-AGN I: AGN-regulated ISM/CGM porosity and stellar ionizing escape at cosmic dawn*, ApJ (submitted).

---

Maintainer: Shubham Raghuvanshi  
