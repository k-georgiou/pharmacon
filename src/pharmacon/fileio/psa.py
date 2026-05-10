"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Handles PharmaconPSA specific HDF5 file operations.

This module contains the `PharmaconPSAFile` class, a specialized file handler
for PSA (Pharmacon Structure Analysis) data files. It extends
`PharmaconPTAFile` and provides dedicated methods for reading and writing
sequence and molecular-property data, as well as ML-ready exporters.

Classes:
    PharmaconPSAFile: Manages PharmaconPSA-specific HDF5 file operations.
"""


import csv
import json
import numpy as np
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Union

from pharmacon.fileio.pta import PharmaconPTAFile




__all__ = [
    "PharmaconPSAFile",
]


class PharmaconPSAFile(PharmaconPTAFile):
    """
    Represents a file handler for PSA (Pharmacon Structure Analysis) data files.

    This class extends functionality from the `PharmaconPTAFile` base class and
    provides mechanisms for managing PSA-related data files. Specific options
    like overwrite behavior and file mode can be assigned during initialization.

    :ivar path: The file path where PSA data is stored.
    :type path: Union[str, Path]
    :ivar overwrite: Determines whether to overwrite the existing file on save.
    :type overwrite: bool
    :ivar mode: Defines the file access mode such as 'read', 'write', or 'append'.
    :type mode: str
    """

    # ------------------------------------------------------------------ #
    # Schema constants                                                   #
    # ------------------------------------------------------------------ #
    SCALAR_PROPERTY_COLUMNS: List[str] = [
        "total_atoms",
        "molecular_weight",
        "logP",
        "net_charge",
        "volume",
        "rotatable_bonds",
        "tpsa",
        "stereo_centers",
        "rings",
        "aromatic_rings",
    ]

    INT_SCALAR_COLUMNS: frozenset[str] = frozenset({
        "total_atoms",
        "rotatable_bonds",
        "stereo_centers",
        "rings",
        "aromatic_rings",
    })

    DICT_PROPERTY_COLUMNS: List[str] = [
        "elements_dictionary",
        "fragments_dictionary",
    ]

    FINGERPRINT_FIELDS: Dict[str, str] = {
        "morgan": "morgan_fingerprint",
        "topological_torsion": "topological_torsion_fingerprint",
        "maccs": "maccs_keys",
        "atom_pair": "atom_pair_fingerprint",
    }

    DEFAULT_FINGERPRINT_PARAMS: Dict[str, Dict[str, Any]] = {
        "morgan": {"radius": 2, "n_bits": 2048, "include_chirality": True},
        "topological_torsion": {"n_bits": 2048},
        "atom_pair": {"n_bits": 2048},
        "maccs": {"n_bits": 167},
    }

    def __init__(self, path: Union[str, Path], *, overwrite: bool = False, mode: str = "a",
                 command: str = "", subcommand: str = "") -> None:
        super().__init__(path=path, overwrite=overwrite, mode=mode,
                         command=command, subcommand=subcommand)

    # ------------------------------------------------------------------ #
    # Sequence I/O                                                       #
    # ------------------------------------------------------------------ #
    def write_sequence(
        self,
        *,
        sequences: Dict[str, Dict],
        group_name: str = "sequence",
        overwrite: bool = False,
    ) -> None:
        """
        Writes sequence data to an HDF5-like group.

        Expected structure of *sequences*::

            {
              chain_id: {
                "aa1_seq": str,
                "aa3_list": list[str],
                "resids_topology": np.ndarray,
                "resid_seq": np.ndarray,
              }
            }

        :param sequences: A dictionary keyed by chain ID.
        :param group_name: Group name. Defaults to ``"sequence"``.
        :param overwrite: Whether to overwrite an existing group.
        :raises ValueError: If the group already exists and `overwrite` is False.
        """
        if self.group_exists(group_name):
            if not overwrite:
                raise ValueError(f"Group already exists: {group_name}")
            self.delete_group(group_name)

        root = self.create_group(group_name)

        for chain_id, data in sequences.items():
            grp = root.create_group(str(chain_id))

            aa1_bytes = data["aa1_seq"].encode("utf-8")
            aa3_bytes = np.asarray(
                [x.encode("utf-8") for x in data["aa3_list"]],
                dtype="S",
            )

            grp.create_dataset("aa1_seq", data=aa1_bytes)
            grp.create_dataset("aa3_list", data=aa3_bytes)
            grp.create_dataset(
                "resids_topology",
                data=np.asarray(data["resids_topology"], dtype=int),
            )
            grp.create_dataset(
                "resid_seq",
                data=np.asarray(data["resid_seq"], dtype=int),
            )

            grp.attrs["chain_id"] = str(chain_id)
            grp.attrs["n_residues"] = len(data["aa1_seq"])

    def read_sequence(self, *, group_name: str = "sequence") -> Dict[str, Dict]:
        """
        Reads a sequence group from the file.

        :param group_name: Group name. Defaults to ``"sequence"``.
        :return: A dictionary keyed by chain ID.
        :raises KeyError: If the sequence group is not found.
        """
        if not self.group_exists(group_name):
            raise KeyError(f"Sequence group not found: {group_name}")

        grp = self.file[group_name]
        results: Dict[str, Dict] = {}

        for chain_id, chain_grp in grp.items():
            aa1_seq = chain_grp["aa1_seq"][()].decode()
            aa3_list = [x.decode() for x in chain_grp["aa3_list"][()]]
            resids_topology = chain_grp["resids_topology"][()].astype(int)
            resid_seq = chain_grp["resid_seq"][()].astype(int)

            results[str(chain_id)] = {
                "aa1_seq": aa1_seq,
                "aa3_list": aa3_list,
                "resids_topology": resids_topology,
                "resid_seq": resid_seq,
            }

        return results

    def write_sequence_fasta(
        self,
        path: Union[str, Path],
        *,
        group_name: str = "sequence",
        overwrite: bool = False,
    ) -> Path:
        """Writes sequences in FASTA format to the specified file path."""
        path = Path(path).expanduser().resolve()

        if path.exists() and not overwrite:
            raise FileExistsError(path)

        seqs = self.read_sequence(group_name=group_name)

        with path.open("w") as f:
            for chain_id, data in seqs.items():
                f.write(f">{chain_id}\n")
                f.write(f"{data['aa1_seq']}\n")

        return path

    def write_sequence_to_csv(
        self,
        path: Union[str, Path],
        *,
        group_name: str = "sequence",
        overwrite: bool = False,
    ) -> Path:
        """Writes sequence data to a CSV file at the specified path."""
        return self._write_sequence_delimited(
            path=path, group_name=group_name, overwrite=overwrite, delimiter=",",
        )

    def write_sequence_to_tsv(
        self,
        path: Union[str, Path],
        *,
        group_name: str = "sequence",
        overwrite: bool = False,
    ) -> Path:
        """Writes sequence data to a TSV file at the specified path."""
        return self._write_sequence_delimited(
            path=path, group_name=group_name, overwrite=overwrite, delimiter="\t",
        )

    def _write_sequence_delimited(
        self,
        *,
        path: Union[str, Path],
        group_name: str,
        overwrite: bool,
        delimiter: str,
    ) -> Path:
        path = Path(path).expanduser().resolve()

        if path.exists() and not overwrite:
            raise FileExistsError(path)

        seqs = self.read_sequence(group_name=group_name)

        with path.open("w", newline="") as f:
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(["chain_id", "resid_seq", "resid_topology", "aa1", "aa3"])

            for chain_id, data in seqs.items():
                for i in range(len(data["aa1_seq"])):
                    writer.writerow([
                        chain_id,
                        int(data["resid_seq"][i]),
                        int(data["resids_topology"][i]),
                        data["aa1_seq"][i],
                        data["aa3_list"][i],
                    ])

        return path

    # ------------------------------------------------------------------ #
    # Properties I/O — raw reader                                        #
    # ------------------------------------------------------------------ #
    def read_properties_data(
        self,
        *,
        group_name: str = "properties",
    ) -> Dict[str, Dict[str, Dict[str, str]]]:
        """
        Reads property data from a specified group and returns a nested dict.

        Returned structure::

            {
              source_file: {
                molecule_name: {
                  property_name: value (str),
                  ...
                }
              }
            }

        :param group_name: The name of the group to read data from. Default is
            ``"properties"``.
        :return: A nested dictionary of source files → molecules → properties.
        :raises KeyError: If the specified group is not found.
        """
        if not self.group_exists(group_name):
            raise KeyError(f"Properties group not found: {group_name}")

        root = self.file[group_name]
        out: Dict[str, Dict[str, Dict[str, str]]] = {}

        for source_file, src_grp in root.items():
            # Skip the group-level `completed` attribute container
            if not hasattr(src_grp, "items"):
                continue
            src_data: Dict[str, Dict[str, str]] = {}

            for mol_name, mol_grp in src_grp.items():
                props: Dict[str, str] = {}
                for key, value in mol_grp.attrs.items():
                    props[str(key)] = str(value)
                src_data[str(mol_name)] = props

            out[str(source_file)] = src_data

        return out

    # ------------------------------------------------------------------ #
    # Properties I/O — typed reader                                      #
    # ------------------------------------------------------------------ #
    def read_properties_typed(
        self,
        *,
        group_name: str = "properties",
    ) -> Dict[str, Dict[str, Dict[str, Any]]]:
        """
        Reads property data and coerces scalars to numeric types, parsing
        ``elements_dictionary`` and ``fragments_dictionary`` from their JSON
        string representations.

        Scalar columns in :attr:`INT_SCALAR_COLUMNS` are cast to ``int``; all
        other :attr:`SCALAR_PROPERTY_COLUMNS` are cast to ``float``. Failed
        coercions yield ``None`` (which CSV writers emit as an empty cell,
        read back as ``NaN`` by pandas).

        :param group_name: The name of the group to read data from.
        :return: Nested dict with numeric values and parsed dicts.
        """
        raw = self.read_properties_data(group_name=group_name)
        out: Dict[str, Dict[str, Dict[str, Any]]] = {}

        for src, mols in raw.items():
            src_out: Dict[str, Dict[str, Any]] = {}
            for mol_name, props in mols.items():
                typed: Dict[str, Any] = {}

                for col in self.SCALAR_PROPERTY_COLUMNS:
                    typed[col] = self._coerce_numeric(
                        props.get(col, ""),
                        is_int=col in self.INT_SCALAR_COLUMNS,
                    )

                for col in self.DICT_PROPERTY_COLUMNS:
                    raw_val = props.get(col, "")
                    try:
                        parsed = json.loads(raw_val) if raw_val else {}
                    except (ValueError, TypeError):
                        parsed = {}
                    typed[col] = parsed

                for kind, field in self.FINGERPRINT_FIELDS.items():
                    typed[field] = props.get(field, "")

                src_out[mol_name] = typed
            out[src] = src_out

        return out

    @staticmethod
    def _coerce_numeric(value: Any, *, is_int: bool) -> Optional[Union[int, float]]:
        if value is None or value == "" or value == "None":
            return None
        try:
            f = float(value)
        except (ValueError, TypeError):
            return None
        if is_int:
            try:
                return int(f)
            except (ValueError, OverflowError):
                return None
        return f

    # ------------------------------------------------------------------ #
    # Properties I/O — ML-ready scalar table                             #
    # ------------------------------------------------------------------ #
    def write_properties_to_csv(
        self,
        path: Union[str, Path],
        *,
        group_name: str = "properties",
        overwrite: bool = False,
    ) -> Path:
        """
        Writes molecular properties to a ML-ready CSV file.

        Emits one row per molecule with the following columns:

        - ``source_file``, ``molecule_name`` (join keys)
        - All scalar descriptors in :attr:`SCALAR_PROPERTY_COLUMNS`, numerically
          typed (int / float / empty for NaN).
        - Dynamic ``elem_<symbol>`` columns, one per element observed across
          the entire dataset. Missing elements default to ``0``.
        - ``frag_<name>`` columns for every RDKit ``fr_*`` fragment present in
          the data. Missing fragments default to ``0``.

        Fingerprints are **not** written here — use
        :meth:`write_fingerprints_to_parquet`.

        Rows are sorted by ``(source_file, molecule_name)`` for deterministic
        output.
        """
        return self._write_properties_delimited(
            path=path, group_name=group_name, overwrite=overwrite, delimiter=",",
        )

    def write_properties_to_tsv(
        self,
        path: Union[str, Path],
        *,
        group_name: str = "properties",
        overwrite: bool = False,
    ) -> Path:
        """TSV counterpart of :meth:`write_properties_to_csv`."""
        return self._write_properties_delimited(
            path=path, group_name=group_name, overwrite=overwrite, delimiter="\t",
        )

    def _write_properties_delimited(
        self,
        *,
        path: Union[str, Path],
        group_name: str,
        overwrite: bool,
        delimiter: str,
    ) -> Path:
        path = Path(path).expanduser().resolve()
        if path.exists() and not overwrite:
            raise FileExistsError(path)

        data = self.read_properties_typed(group_name=group_name)

        elements_sorted, fragments_sorted = self._collect_dynamic_schema(data)
        elem_cols = [f"elem_{e}" for e in elements_sorted]
        frag_cols = [f"frag_{f}" for f in fragments_sorted]

        header = [
            "source_file",
            "molecule_name",
            *self.SCALAR_PROPERTY_COLUMNS,
            *elem_cols,
            *frag_cols,
        ]

        with path.open("w", newline="") as f:
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(header)

            for src in sorted(data.keys()):
                for mol_name in sorted(data[src].keys()):
                    props = data[src][mol_name]
                    row: List[Any] = [src, mol_name]

                    for col in self.SCALAR_PROPERTY_COLUMNS:
                        v = props.get(col)
                        row.append("" if v is None else v)

                    edict = props.get("elements_dictionary", {}) or {}
                    for e in elements_sorted:
                        row.append(int(edict.get(e, 0) or 0))

                    fdict = props.get("fragments_dictionary", {}) or {}
                    for fr in fragments_sorted:
                        row.append(int(fdict.get(fr, 0) or 0))

                    writer.writerow(row)

        return path

    @staticmethod
    def _collect_dynamic_schema(
        data: Dict[str, Dict[str, Dict[str, Any]]],
    ) -> tuple[List[str], List[str]]:
        """Scan all molecules and return (sorted elements, sorted fragments)."""
        elements: set[str] = set()
        fragments: set[str] = set()
        for mols in data.values():
            for props in mols.values():
                elements.update((props.get("elements_dictionary") or {}).keys())
                fragments.update((props.get("fragments_dictionary") or {}).keys())
        return sorted(elements), sorted(fragments)

    # ------------------------------------------------------------------ #
    # Properties I/O — fingerprints → Parquet                            #
    # ------------------------------------------------------------------ #
    def write_fingerprints_to_parquet(
        self,
        path: Union[str, Path],
        *,
        kind: str,
        group_name: str = "properties",
        overwrite: bool = False,
        compression: str = "snappy",
    ) -> Path:
        """
        Writes a single fingerprint family to a Parquet file.

        Produces a wide matrix with columns ``source_file``, ``molecule_name``,
        ``bit_0``..``bit_{N-1}`` where N is the native bit width of the
        fingerprint (auto-detected from the data). Bits are stored as
        ``uint8``.

        :param path: Output Parquet path.
        :param kind: One of :attr:`FINGERPRINT_FIELDS` keys
            (``"morgan"``, ``"maccs"``, ``"atom_pair"``,
            ``"topological_torsion"``).
        :param group_name: Source PSA group. Defaults to ``"properties"``.
        :param overwrite: Overwrite behavior.
        :param compression: Parquet compression codec. Defaults to ``"snappy"``.
        :raises ImportError: If ``pyarrow`` is not installed.
        :raises ValueError: If the fingerprint kind is unknown, the group is
            empty, or fingerprint widths are inconsistent across molecules.
        """
        if kind not in self.FINGERPRINT_FIELDS:
            raise ValueError(
                f"Unknown fingerprint kind: {kind!r}. "
                f"Expected one of: {sorted(self.FINGERPRINT_FIELDS)}"
            )

        path = Path(path).expanduser().resolve()
        if path.exists() and not overwrite:
            raise FileExistsError(path)

        try:
            import pyarrow as pa
            import pyarrow.parquet as pq
        except ImportError as exc:
            raise ImportError(
                "write_fingerprints_to_parquet requires pyarrow. "
                "Install with: pip install pyarrow"
            ) from exc

        field = self.FINGERPRINT_FIELDS[kind]
        data = self.read_properties_data(group_name=group_name)

        rows: List[tuple[str, str, str]] = []
        for src in sorted(data.keys()):
            for mol_name in sorted(data[src].keys()):
                bitstr = data[src][mol_name].get(field, "") or ""
                rows.append((src, mol_name, bitstr))

        if not rows:
            raise ValueError(f"No molecules found in group {group_name!r}")

        widths = {len(r[2]) for r in rows if r[2]}
        if not widths:
            raise ValueError(
                f"No {kind!r} fingerprints found in group {group_name!r}"
            )
        if len(widths) > 1:
            raise ValueError(
                f"Inconsistent {kind!r} fingerprint widths across molecules: "
                f"{sorted(widths)}"
            )
        n_bits = widths.pop()

        n_rows = len(rows)
        bit_matrix = np.zeros((n_rows, n_bits), dtype=np.uint8)
        for i, (_, _, bs) in enumerate(rows):
            if len(bs) != n_bits:
                continue  # leave zeros for missing-width rows
            bit_matrix[i] = np.frombuffer(bs.encode("ascii"), dtype=np.uint8) - ord("0")

        source_files = [r[0] for r in rows]
        mol_names = [r[1] for r in rows]

        arrays: List[Any] = [pa.array(source_files), pa.array(mol_names)]
        names: List[str] = ["source_file", "molecule_name"]
        for j in range(n_bits):
            arrays.append(pa.array(bit_matrix[:, j]))
            names.append(f"bit_{j}")

        table = pa.table(arrays, names=names)
        pq.write_table(table, str(path), compression=compression)
        return path

    # ------------------------------------------------------------------ #
    # Properties I/O — metadata sidecar                                  #
    # ------------------------------------------------------------------ #
    def write_properties_metadata_sidecar(
        self,
        path: Union[str, Path],
        *,
        group_name: str = "properties",
        overwrite: bool = False,
        fingerprint_params: Optional[Dict[str, Dict[str, Any]]] = None,
        extra: Optional[Dict[str, Any]] = None,
    ) -> Path:
        """
        Writes a JSON sidecar describing the ML-ready export schema.

        The sidecar enables downstream pipelines to reconstruct feature groups,
        dtypes, and join keys without inspecting the HDF5 file.

        :param path: Output JSON path.
        :param group_name: PSA group to inspect.
        :param overwrite: Overwrite behavior.
        :param fingerprint_params: Optional override of the fingerprint
            parameter block (defaults to :attr:`DEFAULT_FINGERPRINT_PARAMS`).
        :param extra: Optional additional top-level fields to merge into the
            sidecar.
        """
        path = Path(path).expanduser().resolve()
        if path.exists() and not overwrite:
            raise FileExistsError(path)

        data = self.read_properties_typed(group_name=group_name)
        elements_sorted, fragments_sorted = self._collect_dynamic_schema(data)

        n_molecules = sum(len(mols) for mols in data.values())
        source_files = sorted(data.keys())

        file_attrs = {str(k): str(v) for k, v in self.file.attrs.items()}

        dtypes: Dict[str, str] = {}
        for col in self.SCALAR_PROPERTY_COLUMNS:
            dtypes[col] = "int" if col in self.INT_SCALAR_COLUMNS else "float"
        for e in elements_sorted:
            dtypes[f"elem_{e}"] = "int"
        for fr in fragments_sorted:
            dtypes[f"frag_{fr}"] = "int"

        meta: Dict[str, Any] = {
            "pharmacon_file": str(self.path),
            "pharmacon_version": file_attrs.get("pharmacon_version", ""),
            "command": file_attrs.get("command", ""),
            "subcommand": file_attrs.get("subcommand", ""),
            "input_file": file_attrs.get("input_file", ""),
            "n_molecules": n_molecules,
            "source_files": source_files,
            "row_key": ["source_file", "molecule_name"],
            "feature_groups": {
                "scalars": list(self.SCALAR_PROPERTY_COLUMNS),
                "elements": [f"elem_{e}" for e in elements_sorted],
                "fragments": [f"frag_{fr}" for fr in fragments_sorted],
            },
            "dtypes": dtypes,
            "fingerprint_params": fingerprint_params or self.DEFAULT_FINGERPRINT_PARAMS,
        }

        if extra:
            meta.update(extra)

        with path.open("w") as f:
            json.dump(meta, f, indent=2, sort_keys=False)

        return path

    # ------------------------------------------------------------------ #
    # Properties I/O — orchestrator                                      #
    # ------------------------------------------------------------------ #
    def write_ml_ready_export(
        self,
        base_path: Union[str, Path],
        *,
        group_name: str = "properties",
        overwrite: bool = False,
        include_fingerprints: bool = True,
        fingerprint_kinds: Iterable[str] = ("morgan", "maccs", "atom_pair", "topological_torsion"),
        delimiter: str = ",",
    ) -> Dict[str, Path]:
        """
        One-shot ML-ready export.

        Produces, next to *base_path*:

        - ``<stem>_scalars.csv``  (or ``.tsv`` if ``delimiter='\\t'``)
        - ``<stem>_fp_<kind>.parquet`` for each requested fingerprint kind
        - ``<stem>.meta.json`` sidecar

        Returns a dict mapping logical output name → resolved ``Path``.

        :param base_path: Base output path. The stem is reused for every
            generated file; the parent directory is created if missing.
        :param group_name: PSA group to read.
        :param overwrite: Overwrite behavior for every emitted file.
        :param include_fingerprints: Toggle fingerprint Parquet emission.
        :param fingerprint_kinds: Which fingerprint families to export.
        :param delimiter: ``","`` for CSV, ``"\\t"`` for TSV.
        """
        base_path = Path(base_path).expanduser().resolve()
        base_path.parent.mkdir(parents=True, exist_ok=True)

        stem = base_path.name
        parent = base_path.parent
        ext = "tsv" if delimiter == "\t" else "csv"

        outputs: Dict[str, Path] = {}

        scalars_path = parent / f"{stem}_scalars.{ext}"
        outputs["scalars"] = self._write_properties_delimited(
            path=scalars_path,
            group_name=group_name,
            overwrite=overwrite,
            delimiter=delimiter,
        )

        if include_fingerprints:
            for kind in fingerprint_kinds:
                fp_path = parent / f"{stem}_fp_{kind}.parquet"
                outputs[f"fp_{kind}"] = self.write_fingerprints_to_parquet(
                    fp_path,
                    kind=kind,
                    group_name=group_name,
                    overwrite=overwrite,
                )

        meta_path = parent / f"{stem}.meta.json"
        outputs["meta"] = self.write_properties_metadata_sidecar(
            meta_path,
            group_name=group_name,
            overwrite=overwrite,
        )

        return outputs
