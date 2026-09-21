# Encoding a compartmental AgNP–exosome–Raman theranostic concept as a multi-observation Disease Profile research object

**Thesis #15** (series label R4). Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

Nanobiocomposite proposals often fuse photothermal/ROS, exosome delivery, and Raman sensing into one platform story; without named observation channels and refusal rules, those layers collapse into a single "nanotherapy works" claim.

The deposit encodes that concept as a Disease Profile with three channels (Raman reporter, exosome compartment, photothermal or ROS assay) and an empty therapeutic parameter list. Schema `1.1.0-multiobs` extends the Thesis #3 contract locally. A validator accepts the worked example and refuses a file that writes sensing and killing into one symbol θ.

On a four-preparation linear map the legal Fisher rank is 3 of 3. Each channel alone is rank 1 of 3, and the unused coefficients have flat profiles. The illegal scalar has rank 1 of 1 and a closed profile centred at θ̂ = 1.060, against generating values ρ = 1.70, λ = 0.85, and κ = 0.42. A product of the reporter gain and the damage coefficient leaves a flat profile along the hyperbola. Seed 20260921. The loadings are synthetic.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered. It does not repeat Thesis #8's assay-to-ODE channel, and it does not re-tabulate the 2022 papaya AgNP project.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `schema/multiobservation_profile.schema.json` | JSON Schema for the extended profile |
| `profiles/agnp_exosome_raman.profile.yaml` | Legal worked example |
| `profiles/agnp_exosome_raman.profile.json` | JSON twin of the legal example |
| `profiles/illegal_merged_theta.profile.yaml` | File the validator refuses |
| `sim/channel_profile.py` | Seeded ranks, profiles, and refusal checks (seed 20260921) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Channel means, profiles, and Fisher spectra |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/channel_profile.py
python3 build_pdf.py
```

NumPy, Matplotlib, and PyYAML are required for the sketches. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json`, `profiles/agnp_exosome_raman.profile.json`, and `sim/figures/`.

## Cite

Ogbonna KE. Encoding a compartmental AgNP–exosome–Raman theranostic concept as a multi-observation Disease Profile research object [Internet]. Thesis #15 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-15-nanobiocomposite-multiobservation-profile

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
