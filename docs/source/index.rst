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
     - Kyriakos Georgiou ()
   * - **Affiliation**
     - Department of Pharmacy, University of Athens
   * - **Version**
     - 1.0.0
   * - **License**
     - Apache 2.0
   * - **Python**
     - ≥ 3.12
   * - **Manuscript**
     - `Nature Scientific Data (2023) <https://www.nature.com/articles/s41597-023-00972-3>`_

.. note::

   If you use Pharmacon in your research, please cite:

   Georgiou, K. *Pharmacon: a unified command-line suite for molecular-dynamics
   trajectory and structure analysis.* Nature Scientific Data (2023).
   DOI: ``[TO BE UPDATED]``

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
