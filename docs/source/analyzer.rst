Analyzer Module
===============

The analyzer module (``pharmacon.analyzer``) is the computational core of
Pharmacon.  It contains all analysis algorithms as standalone, stateless
routines that operate on MDAnalysis selections and return raw numerical
results.

.. contents:: Analysis types
   :local:
   :depth: 1

----

RMSD
----

Computes the Root Mean Square Deviation of atomic positions relative to a
reference frame.  Supports multiple atom-group selections per run with a
dedicated fitting group for structural alignment.  Results are stored
per-frame in the ``.pta`` file and plotted as time series.

CLI: ``pharmacon trajectory rmsd``

----

Distances
---------

Computes pairwise or group-to-group distance time series between user-defined
atom selections.  Four computation modes are supported:

- **MIN** — minimum distance between any atom pair across the two groups
- **MAX** — maximum distance between any atom pair across the two groups
- **COM** — distance between the centres of mass
- **COG** — distance between the centres of geometry

All calculations are PBC-aware.  Multiple distance pairs can be calculated in
a single run.

CLI: ``pharmacon trajectory distances``

----

Angles
------

Measures geometric angles across a trajectory.  Three angle types are
supported, selected by the format of the MDAnalysis selection string:

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Type
     - Selection syntax
     - Example
   * - Three-atom angle
     - Select exactly 3 atoms; angle is measured at the central atom
     - ``"index 10 or index 11 or index 12"``
   * - Vector angle
     - Two atom pairs separated by ``->``
     - ``"index 10 or index 100 -> index 300 or index 301"``
   * - Dihedral / torsion
     - Select exactly 4 atoms; torsion measured around the central bond
     - ``"index 10 or index 11 or index 12 or index 13"``

All calculations are PBC-aware when box dimensions are available.

CLI: ``pharmacon trajectory angles``

----

Non-Covalent Interactions
--------------------------

Detects and quantifies **nine types** of non-covalent interactions on a
per-frame basis.  Atom-type classification uses RDKit SMARTS pattern matching
for accuracy.  Available for both protein–ligand and protein–protein pairs.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Interaction
     - Detection criterion
   * - **Hydrophobic contacts**
     - Distance-based between non-polar atoms
   * - **Hydrogen bonds**
     - Donor–hydrogen–acceptor geometry: distance + D–H···A angle
   * - **Ionic / salt-bridge**
     - Oppositely charged groups within a distance cutoff
   * - **Halogen bonds**
     - C–X···A geometry with angle validation
   * - **Metal contacts**
     - Metal ion coordination to nearby electronegative atoms
   * - **Water bridges (1st degree)**
     - Single bridging water mediating an H-bond chain between two groups
   * - **Water bridges (2nd degree)**
     - Two bridging waters forming an extended H-bond network *(not yet implemented)*
   * - **Pi–cation**
     - Aromatic ring to charged group: distance + angle threshold
   * - **Pi–stacking**
     - Parallel or T-shaped aromatic ring pairs classified by geometry

Each interaction record includes:

- Full atom metadata for both partners (name, residue, chain, segment)
- Interaction-specific geometric details (distances, angles, orientation)

See :ref:`interaction-schema` in :doc:`file_formats` for the full column schema.

CLI: ``pharmacon trajectory pl-interactions`` / ``pharmacon trajectory pp-interactions``

----

Hydrogen Bonds (standalone)
----------------------------

Dedicated hydrogen-bond analysis that can be run independently from the full
interaction profiler.  Uses the same donor–acceptor detection and geometric
criteria but produces a focused output containing only hydrogen-bond records.
Useful when the hydrogen-bonding network is the sole interest.

CLI: ``pharmacon trajectory h-bonds``

----

Principal Component Analysis (PCA)
------------------------------------

Performs PCA on atomic coordinates across a trajectory to identify dominant
modes of conformational motion.

**Outputs stored in the** ``.pta`` **file:**

- Eigenvalues (explained variance ratios)
- Eigenvectors (principal components)
- Per-frame projections onto each selected component

**Plot types produced by** ``pharmacon plot pta``:

- PC projections vs time (time series)
- PC1 vs PC2 scatter plot
- Variance-ratio scree plot
- Free-energy surface (FES) heatmap
- Probability density heatmap

CLI: ``pharmacon trajectory pca``

----

Average Structure
-----------------

Computes per-atom mean coordinates over a specified frame range of a
trajectory, with optional alignment to a reference frame.  The resulting
average structure is written as a coordinate file for visualisation or further
analysis.

CLI: ``pharmacon trajectory average-st``

----

Sequence Extraction
-------------------

Extracts amino-acid sequences from a topology file, organised by chain.

**Outputs:**

- Single-letter sequence (``aa1_seq``)
- Three-letter residue list (``aa3_list``)
- Residue ID mapping (``resid_seq``)

Supports FASTA export via ``pharmacon export psa -f fasta``.

CLI: ``pharmacon structure sequence``

----

Molecular Properties
--------------------

Computes structural and chemical descriptors for small molecules using RDKit.

**Scalar properties:**

- Molecular weight, LogP, TPSA
- Rotatable bond count, ring count, aromaticity flags
- Stereocenter count, net formal charge
- Element composition, molecular volume, fragment counts

**Fingerprints:**

- Morgan / ECFP (circular fingerprint)
- Topological torsion
- Atom pair
- MACCS keys

Results support CSV export and ML-ready parquet fingerprint dumps.

CLI: ``pharmacon structure properties``
