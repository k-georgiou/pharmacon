Pharmacon Documentation
========================

**Pharmacon — A Molecular Dynamics Simulation Analysis Toolkit**

A unified command-line toolkit for analyzing molecular-dynamics trajectories
and static structures, persisting results as signed, self-describing HDF5
artifacts, and producing publication-ready plots.

.. list-table::
   :widths: 30 70
   :stub-columns: 0

   * - **Author**
     - Kyriakos Georgiou
   * - **Affiliation**
     - Department of Pharmacy, University of Athens
   * - **Version**
     - 1.0.1
   * - **License**
     - GPLv3-only (GNU General Public License v3.0 only)
   * - **Python**
     - ≥ 3.12
   * - **Manuscript**
     - *Journal of Chemical Information and Modeling* (2026) — DOI: 10.1021/acs.jcim.6c00837

.. note::

   If you use Pharmacon in your research, please cite:

   Georgiou, K.; Kolocouris, A. *Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit.*
   Journal of Chemical Information and Modeling (2026).
   DOI: ``10.1021/acs.jcim.6c00837``

----

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Command Reference

   commands/index
   commands/trajectory
   commands/structure
   commands/dump
   commands/export
   commands/plot
   commands/merge

.. toctree::
   :maxdepth: 2
   :caption: Reference

   architecture
   analyzer
   pbc
   file_formats
   plot_configuration
   logging
