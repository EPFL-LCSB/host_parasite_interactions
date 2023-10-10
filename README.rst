Modeling host-parasite interactions
====================================

Generation of models for human hepatocyte and plasmodium falciparum contextualized to the nutrient environment on the host hepatocyte. Analysis of nutritional interactions between host and parasite and host gene essentiality for parasite survival. Paper: Marina Maurizio, Maria Masid et al. "Host metabolic pathways essential for malaria and related hemoparasites in the infection of nucleated cells" 


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
Run the scripts `run_host_parasite_interactions <https://github.com/EPFL-LCSB/redhuman/run_host_parasite_interactions>`_ from the main folder.
There are four steps:
1. Generation of a hepatocyte model from the human Recon 3D model
2. Generation of parasitosomes for Plasmodium falciparum
3. Reconstruction of an integrated host-parasite model
4. Performing gene essentiality analysis

License
=======
The software in this repository is put under an APACHE licensing scheme - please see the `LICENSE <https://github.com/EPFL-LCSB/host_parasite_interactions/blob/master/LICENSE>`_ file for more details.

