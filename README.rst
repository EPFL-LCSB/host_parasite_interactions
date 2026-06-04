Modeling host-parasite interactions
====================================

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.20543332.svg
   :target: https://doi.org/10.5281/zenodo.20543332
   :alt: Zenodo DOI (concept — all versions)

Generation of models for human hepatocyte and plasmodium falciparum contextualized to the nutrient environment on the host hepatocyte. Analysis of nutritional interactions between host and parasite and host gene essentiality for parasite survival. Companion code for: Maurizio M, Masid M, *et al.* "Host cell CRISPR genomics and modelling reveal shared metabolic vulnerabilities in the intracellular development of *Plasmodium falciparum* and related hemoparasites." *Nature Communications* **15**, 6145 (2024). doi:10.1038/s41467-024-50405-x 


Requirements
------------

You will need to have `Git LFS <https://git-lfs.github.com/>`_ in order to properly download the binary files:

.. code:: bash

    git clone https://github.com/EPFL-LCSB/host_parasite_interactions.git /path/to/host_parasite_interactions
    cd /path/to/host_parasite_interactions
    git lfs install
    git lfs pull

The scripts have been developed with Matlab 2021b, and CPLEX 12.10 (freely downloadable with the `IBM Academic initiative <https://developer.ibm.com/academic/>`_), and successfully ran on several other versions of both softwares. However, it is important to respect the IBM compatibility specs sheets between Matlab, CPLEX, and the computer OS - available `on IBM's website <https://www.ibm.com/software/reports/compatibility/clarity/index.html>`_.

This module requires `matTFA <https://github.com/EPFL-LCSB/mattfa/>`_, `redGEM <https://github.com/EPFL-LCSB/redgem/>`_


Generating host, parasite and host-parasite metabolic models and performing gene essentiality analysis
-------------------------------------------------------------------------------------------------------
Run the script `run_host_parasite_interactions.m <https://github.com/EPFL-LCSB/host_parasite_interactions/blob/main/run_host_parasite_interactions.m>`_ from the main folder.
There are four steps:
1. Generation of a hepatocyte model from the human Recon 3D model
2. Generation of parasitosomes for Plasmodium falciparum
3. Reconstruction of an integrated host-parasite model
4. Performing gene essentiality analysis

Associated data
===============
The CRISPR screen and RNA-seq data underlying the paper are deposited in the
European Nucleotide Archive (ENA) at EMBL-EBI under BioProject
`PRJEB73635 <https://www.ebi.ac.uk/ena/browser/view/PRJEB73635>`_
(RNA-seq samples: ERS16273013, ERS16273014, ERS16273015, ERS18408949, ERS18408950, ERS18408951).

Citation
========
If you use this software, please cite **both** the archived software and the paper.

* Software archive (Zenodo): concept DOI ``10.5281/zenodo.20543332`` (all versions; resolves to the latest). This release v1.0.0: ``10.5281/zenodo.20543333``.
* Paper: Maurizio M, Masid M, Woods K, *et al.* *Host cell CRISPR genomics and modelling reveal shared metabolic vulnerabilities in the intracellular development of Plasmodium falciparum and related hemoparasites.* **Nature Communications** 15, 6145 (2024). `doi:10.1038/s41467-024-50405-x <https://doi.org/10.1038/s41467-024-50405-x>`_

A machine-readable ``CITATION.cff`` is provided at the repository root.

License
=======
The **source code** in this repository (the ``code/`` directory and
``run_host_parasite_interactions.m``) is released under the Apache-2.0 license —
see the `LICENSE <https://github.com/EPFL-LCSB/host_parasite_interactions/blob/main/LICENSE>`_ file.

The **data and models** under ``data/`` and ``models/`` are *not* covered by
Apache-2.0. Several are third-party datasets (Recon3D, BiGG, the Human Protein Atlas, …)
redistributed under their own — in some cases **non-commercial** — licenses. See
`DATA_LICENSES.md <https://github.com/EPFL-LCSB/host_parasite_interactions/blob/main/DATA_LICENSES.md>`_
for per-file provenance, license and citation.

