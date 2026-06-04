# Data sources, provenance and licenses

The repository's top-level `LICENSE` (**Apache-2.0**) applies to the **source code only**
(the `code/` directory and `run_host_parasite_interactions.m`).

The files under `data/` and `models/` are **data**, not code. Some of them are
third-party datasets redistributed here for reproducibility and carry their **own
licenses**, which are listed below. Apache-2.0 does **not** apply to these files.
Several upstream sources are **non-commercial** — reuse of the corresponding files
(and of any model derived from them) is restricted accordingly.

This file documents each dataset's origin, license and citation, in line with
FAIR principle **R1.1** (data released with a clear, accessible usage license).

> ⚠️ **Items marked “(confirm)” should be verified against the upstream source
> before relying on them.** In particular, read the current Recon3D/VMH terms on
> <https://www.vmh.life> and the BiGG terms on <http://bigg.ucsd.edu/license>.

## `data/`

| File | Origin | License | Cite |
|---|---|---|---|
| `Recon3D.xls` | Human metabolic reconstruction Recon3D, via the Virtual Metabolic Human (VMH, <https://www.vmh.life>) | **Non-commercial Creative Commons** (CC BY-NC 4.0, *confirm exact variant on vmh.life*) — **not** Apache-2.0 | Brunk et al., *Nat Biotechnol* **36**, 272 (2018), doi:10.1038/nbt.4072 |
| `bigg_models_metabolites.txt` | Verbatim BiGG Models universal-metabolites table (<http://bigg.ucsd.edu/data_access>) | **BiGG / UC San Diego license — non-commercial**; redistribution allowed for research/education with copyright & license retained; commercial use requires UCSD permission. **Not** Apache-2.0 / not CC BY | King et al., *Nucleic Acids Res* **44**, D515 (2016), doi:10.1093/nar/gkv1049 |
| `transcriptomics_hepatocyte.mat` | Liver transcriptomics derived from the Human Protein Atlas (<https://www.proteinatlas.org>) | **CC BY 4.0** (attribution required) | Uhlén et al., *Science* **347**, 1260419 (2015), doi:10.1126/science.1260419 |
| `humanGenome_hgnc_entrez.mat` | HGNC ↔ Entrez gene-ID mapping (HGNC, <https://www.genenames.org>) | Open / no usage restriction (attribution requested) — *confirm* | Tweedie et al., *Nucleic Acids Res* **49**, D939 (2021) |
| `thermoData_Recon3.mat` | Thermodynamic data (group-contribution) used by matTFA, mapped to Recon3D | Inherits matTFA terms (Apache-2.0) + group-contribution data terms — *confirm* | Salvy et al. (matTFA), *Bioinformatics* **35**, 167 (2019), doi:10.1093/bioinformatics/bty499 |
| `ipfaMedium.mat` | *P. falciparum* growth-medium definition (derived from the iPfa model, see `models/iPfa.mat`) | Derived from iPfa — inherits upstream terms (*confirm*) | see `models/iPfa.mat` below |
| `MetabolomicsData.csv` | This study — measured/curated metabolite concentration bounds | **CC BY 4.0** (original to this work) | this work (cite the paper) |
| `Uptakes.csv` | This study — curated exchange/uptake bounds (uses BiGG/SEED identifiers) | **CC BY 4.0** (original to this work; identifiers only) | this work (cite the paper) |
| `PmetsRecon3_redHUMAN.mat` | This study — *P. falciparum*-relevant metabolite set mapped onto Recon3D | Derived from Recon3D → **non-commercial** constraints apply | this work + Recon3D |
| `matches_table_para_host.mat` | This study — host↔parasite metabolite matching table | **CC BY 4.0** (original to this work) | this work |

## `models/`

The generated models are original outputs of this study but are **derived from**
Recon3D (and, for the integrated models, from the *P. falciparum* model). They
therefore **inherit the non-commercial restrictions** of Recon3D/VMH (and any
upstream iPfa terms). Recommended license for the genuinely original model
artefacts: **CC BY 4.0**, *subject to the upstream non-commercial constraints noted above*.

| File | Origin | License |
|---|---|---|
| `iPfa.mat` | *Plasmodium falciparum* genome-scale model (iPfa) | Upstream iPfa terms — **(confirm original source & citation)** |
| `GEM_Recon3_thermocurated_redHUMAN.mat` | Thermodynamically curated, reduced human GEM (from Recon3D) | Derived from Recon3D → non-commercial |
| `hepatocyte_thermocurated.mat` | Context-specific hepatocyte model (from Recon3D) | Derived from Recon3D → non-commercial |
| `hepatocyte_with_drains_for_para_as_auxotrophic.mat` | Hepatocyte model with parasite drains (auxotrophic) | Derived from Recon3D → non-commercial |
| `hepatocyte_with_drains_for_para_as_prototrophic.mat` | Hepatocyte model with parasite drains (prototrophic) | Derived from Recon3D → non-commercial |
| `integrated_hepatocyte_model_for_para_as_auxotrophic.mat` | Integrated host–parasite model (auxotrophic) | Derived from Recon3D + iPfa → non-commercial |
| `integrated_hepatocyte_model_for_para_as_prototrophic.mat` | Integrated host–parasite model (prototrophic) | Derived from Recon3D + iPfa → non-commercial |

## A note on KEGG

This repository does **not** redistribute bulk KEGG content. It stores only KEGG
**identifiers** and queries the public KEGG REST API at run time
(`code/database_query/queryKEGGIDS.m`). KEGG content itself remains © Kanehisa
Laboratories; users must obtain it from <https://www.kegg.jp> under KEGG's terms.

## Workspaces

`workspaces/*.mat` are intermediate MATLAB workspaces saved at each pipeline step.
They are convenience snapshots of the outputs above and inherit the same provenance
and license constraints as the data/models they contain.
