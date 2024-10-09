.. ABACUS documentation master file, created by
   sphinx-quickstart on Fri Mar 11 10:42:27 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=================================================
ABACUS Documentation
=================================================

ABACUS (Atomic-orbital Based Ab-initio Computation at UStc) is
an open-source computer code package based on density functional
theory (DFT). The package utilizes both plane wave and numerical
atomic basis sets with the usage of norm-conserving pseudopotentials
to describe the interactions between nuclear ions and valence electrons.
ABACUS supports LDA, GGA, meta-GGA, and hybrid functionals. Apart from
single-point calculations, the package allows geometry optimizations
and ab-initio molecular dynamics with various ensembles. The package
also provides a variety of advanced functionalities for simulating materials,
including the DFT+U, VdW corrections, and implicit solvation model, etc.
In addition, ABACUS strives to provide a general infrastructure to facilitate
the developments and applications of novel machine-learning-assisted DFT methods
(DeePKS, DP-GEN, DeepH, etc.) in molecular and material simulations.

.. toctree::
   :maxdepth: 2
   :caption: Quick Start

   quick_start/easy_install
   quick_start/hands_on
   quick_start/input
   quick_start/output

.. toctree::
   :maxdepth: 2
   :caption: Advanced

   advanced/install
   advanced/scf/index
   advanced/pp_orb
   advanced/opt
   advanced/md
   advanced/acceleration/index
   advanced/elec_properties/index
   advanced/interface/index
   advanced/input_files/index

.. toctree::
   :maxdepth: 2
   :caption: Citing ABACUS

   CITATIONS

.. toctree::
   :maxdepth: 2
   :caption: Developing Team

   DevelopingTeam

.. toctree::
   :maxdepth: 2
   :caption: Community

   community/contribution_guide
   CONTRIBUTING
   community/cicd.md

.. toctree::
   :glob:
   :titlesonly:

   community/faq
