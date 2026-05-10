"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Module :mod:`pharmacon.fileio.base`.
"""
import re
from collections import defaultdict
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Union, Final, Tuple, Dict, List

import h5py
from rich.align import Align
from rich.box import HEAVY, ROUNDED, SIMPLE_HEAVY
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.tree import Tree


from pharmacon.constants import __version__, BASE_PHARMACON_META
from pharmacon.logger import get_logger, PharmaconLogger
from pharmacon.utils.fingerprint import create_pharmacon_signature




__all__ = [
    "logger",
    "PharmaconHDF5Types",
    "PharmaconHDF5File",
]


logger: PharmaconLogger = get_logger(__name__)


class PharmaconHDF5Types(Enum):
    """
    Defines an enumeration for different types of HDF5 datasets related to Pharmacon
    analysis.

    This class provides a way to categorize datasets based on their purpose within
    Pharmacon's workflow. Each member of the Enum corresponds to a specific type
    of analysis, such as trajectory or structure analysis.

    :ivar TRAJECTORY_ANALYSIS: Represents HDF5 datasets related to trajectory
        analysis in the Pharmacon workflow.
    :type TRAJECTORY_ANALYSIS: str
    :ivar STRUCTURE_ANALYSIS: Represents HDF5 datasets related to structure
        analysis in the Pharmacon workflow.
    :type STRUCTURE_ANALYSIS: str
    """
    TRAJECTORY_ANALYSIS = "pta"
    STRUCTURE_ANALYSIS  = "psa"


class PharmaconHDF5File:
    """
    Represents an interface for handling Pharmacon-compatible HDF5 files.

    This class provides functionalities for working with HDF5 files that conform
    to Pharmacon's specifications. It enables file creation, opening, modification,
    and metadata management, while also ensuring compliance with required formats
    and extensions.

    :ivar SUPPORTED_EXTENSIONS: Tuple containing supported file extensions for
        Pharmacon HDF5 files.
    :type SUPPORTED_EXTENSIONS: Final[Tuple[str, ...]]
    :ivar file_meta: Metadata dictionary containing predefined and current metadata
        for the Pharmacon file.
    :type file_meta: Dict[str, str]
    """

    SUPPORTED_EXTENSIONS: Final[Tuple[str, ...]] = (".pta", ".psa")
    FILE_METADATA_CASES: Dict[str, List[str]] = {
        "command": [
            "command",
            "subcommand",
            "arguments",
        ],

        "core": [
            "file_type",
            "pharmacon_version",
            "created_at",
        ],

        "analysis": [
            "labels",
            "reference_frame",
            "begin",
            "end",
            "step",
        ],

        "input": [
            "topology_file",
            "trajectory_file",
            "fitting_group",
            "calculation_groups",
        ],
    }
    _FRAME_RE = re.compile(r"(.+?)_(\d+)$")

    def __init__(self, path: Union[str, Path], *, overwrite: bool = False, mode: str = "a",
                 command: str = "", subcommand: str = "") -> None:

        if mode == "r" and overwrite:
            raise ValueError("overwrite=True is incompatible with mode='r'")

        self.overwrite = overwrite
        self.mode = mode
        self.path = path

        self.file_type = PharmaconHDF5Types.STRUCTURE_ANALYSIS if ".psa" in str(self.path) \
            else PharmaconHDF5Types.TRAJECTORY_ANALYSIS

        # If overwrite is enabled and the file already exists, remove it first
        # so we start with a clean HDF5 file rather than appending to stale data.
        if overwrite and self.mode in ("a", "w") and self.path.exists():
            self.path.unlink()

        is_new_file = (
                self.mode == "w"
                or (self.mode == "a" and not self.path.exists())
        )

        self.file = h5py.File(self.path, self.mode)


        self.file_meta: Dict[str, str] = dict(BASE_PHARMACON_META)
        self.file_meta["file_type"] = str(self.file_type)
        self.file_meta["command"] = command
        self.file_meta["subcommand"] = subcommand

        if is_new_file:
            self.file_meta["created_at"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

        if is_new_file and command and subcommand:
            sig = create_pharmacon_signature(
                format_name=self.file_type.value,
                command=command,
                subcommand=subcommand,
            )
            self.file_meta["signature"] = sig.signature
            self.file_meta["fingerprint"] = sig.fingerprint

        if is_new_file:
            self.add_file_metadata(self.file_meta)

    @property
    def overwrite(self) -> bool:
        """
        Indicates whether an overwrite operation is enabled.

        This property retrieves the state of the `_overwrite` attribute, which specifies
        if overwrite behavior is enabled or not.

        :return: A boolean value indicating if overwrite is enabled.
        :rtype: bool
        """
        return self._overwrite

    @overwrite.setter
    def overwrite(self, value: bool) -> None:
        """
        Sets the `overwrite` property for the instance. This property determines whether
        a specific operation should override existing data or behavior. Ensures that the
        value assigned is of type `bool`.

        :param value: Specifies whether to enable or disable the overwrite behavior.
        :type value: bool
        :return: None
        """
        assert isinstance(value, bool), "overwrite must be boolean"
        self._overwrite = value

    @property
    def mode(self) -> str:
        """
        Represents the mode attribute as a property.

        This property method allows access to the private `_mode` attribute.

        :rtype: str
        :return: The current value of the `_mode` attribute.
        """
        return self._mode

    @mode.setter
    def mode(self, value: str) -> None:
        """
        Sets the mode attribute for HDF5 functionality.

        This property setter enforces the input mode to be a valid value among
        "a", "w", or "r". It processes the input value by converting it to
        lowercase and stripping leading/trailing whitespace. If the processed
        value is not one of the accepted modes, it raises a ValueError.

        :param value: The mode to be set for HDF5 operations. Needs to be one of
            {"a", "w", "r"}.
        :type value: str
        :raises ValueError: If the provided mode is not a valid HDF5 mode.
        """
        value = value.lower().strip()
        if value not in {"a", "w", "r"}:
            raise ValueError(f"Invalid HDF5 mode: {value}")
        self._mode = value

    @property
    def path(self) -> Path:
        """
        Provides access to the file path associated with the instance. This property
        returns the path object, retrieved from an internal attribute, which represents
        the location of a file or directory.

        :rtype: Path
        :return: The file system path as a Path object.
        """
        return self._path

    @path.setter
    def path(self, value: Union[str, Path]) -> None:
        """
        Sets the file path while ensuring it meets specific conditions based on the
        provided mode and configurations. The function validates the path-type, resolves
        absolute paths, checks for supported file extensions, and enforces rules for
        existing files.

        :param value: The path to be set. It can be a string or a Path object.
        :raises TypeError: If the provided value is not of type str or Path.
        :raises ValueError: If the file extension is not supported or if overwrite=True
            is used with mode='r'.
        :raises FileExistsError: If the file already exists and overwrite is set to False.
        """
        if not isinstance(value, (str, Path)):
            raise TypeError(f"Invalid path type: {type(value)}")

        value = Path(value).expanduser().resolve()

        if value.suffix not in self.SUPPORTED_EXTENSIONS:
            raise ValueError(f"Unsupported file extension: {value.suffix}")

        if value.exists():
            if self.mode == "r":
                # Read-only mode: existing file is REQUIRED
                pass
            else:
                # Write / append modes
                if not self.overwrite:
                    raise FileExistsError(f"File already exists: {value}")
                value.unlink()

        self._path = value

    def validate_version(self, *, strict: bool = False) -> None:
        """
        Validates the Pharmacon version compatibility between the metadata in a file and
        the runtime version. Ensures that the major version is consistent and provides
        appropriate warnings or raises errors based on the strictness parameter.

        :param strict: If True, any validation failure will raise a ValueError. If False,
            warnings will be logged instead of raising exceptions.
        :return: None
        """
        file_v = self.file_meta.get("pharmacon_version")
        runtime_v = str(__version__)

        if not file_v:
            msg = "File is missing 'pharmacon_version' metadata"
            if strict:
                raise ValueError(msg)
            logger.warning(msg)
            return

        try:
            f_major, _, _ = map(int, file_v.split("."))
            r_major, _, _ = map(int, runtime_v.split("."))
        except Exception:
            msg = (
                "Invalid Pharmacon version format:\n"
                f"  File version    : {file_v}\n"
                f"  Runtime version : {runtime_v}"
            )
            if strict:
                raise ValueError(msg)
            logger.warning(msg)
            return

        if f_major != r_major:
            msg = (
                "Pharmacon major version mismatch detected:\n"
                f"  File version    : {file_v}\n"
                f"  Runtime version : {runtime_v}"
            )
            if strict:
                raise ValueError(msg)
            logger.warning(msg)

    def close(self) -> None:
        """
        Closes the file resource if it exists.

        This method checks whether the instance has an attribute named `file`.
        If the attribute exists, it invokes the `close` method on the `file`
        object to release the associated resource.

        :return: None
        """
        if hasattr(self, "file"):
            self.file.close()

    def group_exists(self, group_name: str) -> bool:
        """
        Checks if a given group exists in the HDF5 file and is of type h5py.Group.

        :param group_name: The name of the group to check.
        :type group_name: str
        :return: True if the group exists and is of the correct type; otherwise, False.
        :rtype: bool
        """
        return group_name in self.file and isinstance(self.file[group_name], h5py.Group)

    def delete_group(self, group_name: str) -> None:
        """
        Deletes a group with the specified name from the current file.

        This function checks whether the group exists in the file and deletes it
        if found.

        :param group_name: The name of the group to be deleted.
        :type group_name: str
        :return: None
        """
        if group_name in self.file:
            del self.file[group_name]

    def create_group(self, group_name: str) -> h5py.Group:
        """
        Creates a new group within the HDF5 file with the specified name. If the group
        already exists, it will only be overwritten if the `overwrite` attribute is
        set to True.

        :param group_name: Name of the group to be created.
        :type group_name: str
        :return: The created HDF5 group.
        :rtype: h5py.Group
        :raises ValueError: If the group already exists and `overwrite` is False.
        """
        if self.group_exists(group_name):
            if not self.overwrite:
                raise ValueError(f"Group already exists: {group_name}")
            self.delete_group(group_name)

        return self.file.create_group(group_name)

    def create_dataset(self, *, group_name: str, dataset_name: str, data, metadata: Dict[str, str] | None = None):
        """
        Creates a new dataset within the specified group in an HDF5 file. If a dataset
        with the same name already exists and `overwrite` is not enabled, an error is
        raised. Optionally, metadata can be added to the created dataset.

        :param group_name: Name of the group in the HDF5 file where the dataset should
                           be created. Must already exist.
        :type group_name: str
        :param dataset_name: Name of the dataset to create within the group.
        :type dataset_name: str
        :param data: Data to store in the dataset.
        :type data: Any
        :param metadata: Optional metadata to associate with the dataset as key-value
                         pairs. Metadata will be added only if specified.
                         Defaults to None.
        :type metadata: dict[str, str] | None
        :return: The created HDF5 dataset.
        :rtype: h5py.Dataset
        :raises KeyError: If the specified group does not exist.
        :raises ValueError: If the dataset already exists and overwriting is not
                            enabled.
        """
        if not self.group_exists(group_name):
            raise KeyError(f"Group does not exist: {group_name}")

        group = self.file[group_name]

        if dataset_name in group:
            if not self.overwrite:
                raise ValueError(f"Dataset already exists: {group_name}/{dataset_name}")
            del group[dataset_name]

        dset = group.create_dataset(dataset_name, data=data)

        if metadata:
            self.add_dataset_metadata(
                group_name=group_name,
                dataset_name=dataset_name,
                metadata=metadata,
                overwrite=True,
            )

        return dset

    def add_group_metadata(self, *, group_name: str, metadata: Dict[str, str], overwrite: bool = False) -> None:
        """
        Adds metadata attributes to a specified group if it exists in the file. Each metadata key-value pair will
        be added to the group's attributes. If the key already exists and `overwrite` is False, an error will be
        raised. If `overwrite` is True, existing attribute values will be replaced.

        :param group_name: The name of the group to which metadata will be added.
        :param metadata: A dictionary containing key-value pairs representing metadata attributes.
        :param overwrite: A boolean flag to indicate whether to overwrite existing attributes with the same key.
                          Defaults to False. If False and a key already exists, an error will be raised.
        :raises KeyError: If the specified group does not exist in the file.
        :raises ValueError: If a metadata key already exists in the group and `overwrite` is False.
        :return: None
        """
        if not self.group_exists(group_name):
            raise KeyError(f"Group does not exist: {group_name}")

        group = self.file[group_name]

        for key, value in metadata.items():
            if key in group.attrs and not overwrite:
                raise ValueError(
                    f"Group attribute '{key}' already exists in {group_name}"
                )
            group.attrs[key] = str(value)

    def add_dataset_metadata(self, *, group_name: str, dataset_name: str, metadata: Dict[str, str],overwrite: bool = False) -> None:
        """
        Adds metadata to a specified dataset within a group in the file.

        This method allows adding key-value pair metadata to a specific dataset. The
        metadata is stored as attributes of the dataset. If a key already exists and
        the `overwrite` flag is not set, an error will be raised.

        :param group_name: The name of the group containing the dataset.
        :param dataset_name: The name of the dataset to which metadata is added.
        :param metadata: A dictionary containing metadata as key-value pairs.
        :param overwrite: Whether to overwrite existing metadata attributes. Defaults
            to False.
        :return: None
        :raises KeyError: If the specified group or dataset does not exist.
        :raises ValueError: If a metadata key already exists and overwrite is set to
            False.
        """
        if not self.group_exists(group_name):
            raise KeyError(f"Group does not exist: {group_name}")

        group = self.file[group_name]

        if dataset_name not in group:
            raise KeyError(f"Dataset does not exist: {group_name}/{dataset_name}")

        dset = group[dataset_name]

        for key, value in metadata.items():
            if key in dset.attrs and not overwrite:
                raise ValueError(
                    f"Dataset attribute '{key}' already exists in {group_name}/{dataset_name}"
                )
            dset.attrs[key] = str(value)

    def add_file_metadata(self, metadata: Dict[str, str], *, overwrite: bool = False) -> None:
        """
        Adds metadata to a file. This function assigns key-value pairs from the provided
        metadata dictionary to the file's attributes. If the key already exists in the
        file attributes and the `overwrite` flag is set to `False`, a KeyError will be raised.

        :param metadata: Dictionary containing key-value metadata entries to be added.
                         Keys and values must be strings.
        :type metadata: Dict[str, str]
        :param overwrite: Flag indicating whether to overwrite existing metadata in the file.
        :type overwrite: bool
        :return: None
        """
        for key, value in metadata.items():
            if key in self.file.attrs and not overwrite:
                raise KeyError(f"File metadata '{key}' already exists")
            self.file.attrs[key] = str(value)

    def get_groups(self, group_name: str | None = None) -> List[str]:
        """
        Retrieve the list of group names from the dataset or the specific group, if provided.

        This method checks the dataset structure to identify all groups it contains.
        If a specific group name is supplied as an argument, it attempts to retrieve
        the sub-groups from the specified group. If the provided group name does
        not exist within the dataset, an empty list is returned.

        :param group_name: Optional name of a group to retrieve its sub-groups.
                           If None, retrieve top-level groups from the dataset.
                           Defaults to None.
        :return: A list of names of groups in the dataset or in the specified group.
        :rtype: List[str]
        """
        if group_name is None:
            grp = self.file
        else:
            if group_name not in self.file:
                return []
            grp = self.file[group_name]

        return [
            key for key, obj in grp.items()
            if isinstance(obj, h5py.Group)
        ]

    def print_file_metadata(self) -> None:
        """
        Prints metadata related to a Pharmacon file in a tabular format.

        This function retrieves file metadata from the `self.file.attrs` dictionary, processes
        the key-value pairs to differentiate between core and additional metadata, and displays
        both types of metadata in separate tables using the `Console` class. Core metadata is
        mapped to a user-friendly format, whereas additional metadata is presented in a more
        compact view. Specific keys are excluded from the display if they match the
        predefined exclusion list.

        :raises KeyError: if required keys are missing from the `attrs` data structure.
        """

        console = Console()
        attrs = dict(self.file.attrs)

        CORE_KEYS = {
            "command": "Command",
            "subcommand": "Subcommand",
            "date_created": "Date Created",
            "file_type": "File Type",
            "hostname": "Hostname",
            "username": "Username",
            "pharmacon_version": "Pharmacon Version",
            "is_merged": "Is Merged",
            "uuid": "UUID",
            "blueprint": "Blueprint",
            "description": "Description",
            "signature": "Signature",
            "fingerprint": "Fingerprint",
        }

        EXCLUDED_KEYS = {
            "artifact_status",
            "artifact_status_code",
            "artifact_token",
            "artifact_token_version",
        }

        core_rows: List[Tuple[str, str]] = []
        extra_rows: List[Tuple[str, str]] = []

        for key, value in attrs.items():
            k = str(key).strip().lower()
            v = str(value)

            if k in CORE_KEYS:
                core_rows.append((CORE_KEYS[k], v))
            elif k not in EXCLUDED_KEYS:
                extra_rows.append((str(key), v))

        # Core metadata
        core_table = Table(
            title="[bold]Pharmacon File Metadata[/bold]",
            show_header=True,
            header_style="bold cyan",
            box=ROUNDED,
            title_style="bold white",
            border_style="bright_blue",
            pad_edge=True,
            padding=(0, 1),
        )

        core_table.add_column("Property", style="bold magenta", no_wrap=True, min_width=20)
        core_table.add_column("Value", style="green", overflow="fold")

        for k, v in core_rows:
            core_table.add_row(k, v)

        console.print()
        console.print(core_table)

        # Additional metadata (compact)
        if extra_rows:
            extra_table = Table(
                title="[bold]Arguments Passed[/bold]",
                show_header=True,
                header_style="bold cyan",
                box=SIMPLE_HEAVY,
                title_style="bold white",
                border_style="dim",
                pad_edge=True,
                padding=(0, 1),
            )

            extra_table.add_column(
                "Key",
                style="cyan",
                no_wrap=True,
                min_width=24,
                max_width=32,
            )
            extra_table.add_column(
                "Value",
                style="white",
                overflow="fold",
                max_width=max(40, console.width - 38),
            )

            for k, v in sorted(extra_rows):
                extra_table.add_row(k, v)

            console.print()
            console.print(extra_table)

    def _group_indexed_children(self, h5group):
        groups = defaultdict(list)
        others = []

        for name, obj in h5group.items():
            if isinstance(obj, h5py.Group):
                m = self._FRAME_RE.match(name)
                if m:
                    groups[m.group(1)].append(int(m.group(2)))
                else:
                    others.append(name)
            else:
                others.append(name)

        return groups, others

    def print_tree(self, *, group_meta: bool = True, dataset_meta: bool = True, compact: bool = False,
                   max_items: int = 10, file_attrs: bool = True) -> None:
        """
        Prints the hierarchical structure of an HDF5 file or group as a tree.

        This method creates a visual representation of the HDF5 structure using a tree
        format. It displays HDF5 groups as folders and datasets as files, with
        additional metadata like attributes, shape, and datatype optionally included.

        :param group_meta: Whether to display the attributes of HDF5 groups. Defaults to True.
        :type group_meta: bool
        :param dataset_meta: Whether to display the attributes of HDF5 datasets. Defaults to True.
        :type dataset_meta: bool
        :param compact: Whether to display indexed groups (e.g., frames, steps) in a compact manner
            when the count exceeds `max_items`. Defaults to False.
        :type compact: bool
        :param max_items: Maximum number of indexed group entries to display before switching
            to a compact representation. Defaults to 10.
        :type max_items: int
        :param file_attrs: Whether to display root-level (file) attributes in the tree.
            Set to False when file metadata is already shown elsewhere. Defaults to True.
        :type file_attrs: bool
        :return: None
        """


        console = Console()

        root = Tree(
            f"[bold magenta]{self.path.name}[/bold magenta]  "
            f"[dim italic]{self.path.parent}[/dim italic]",
            guide_style="bold bright_blue",
        )

        def add_attrs(tree, attrs):
            for k, v in sorted(attrs.items()):
                val_str = str(v)
                if len(val_str) > 80:
                    val_str = val_str[:77] + "..."
                tree.add(f"[dim cyan]@{k}[/dim cyan] [dim]=[/dim] [dim white]{val_str}[/dim white]")

        def walk(h5group, tree, depth=0, is_root=False):
            show_attrs = (file_attrs or not is_root) and group_meta and h5group.attrs
            if show_attrs:
                add_attrs(tree, h5group.attrs)

            indexed, others = self._group_indexed_children(h5group)

            # --- indexed groups (frames, steps, etc.) ---
            for prefix, indices in indexed.items():
                indices.sort()
                count = len(indices)

                if compact and count >= max_items:
                    tree.add(
                        f"[bold green]{prefix}_*[/bold green]  "
                        f"[dim bright_black]{count} groups[/dim bright_black] "
                        f"[dim]([/dim][bright_cyan]{indices[0]}[/bright_cyan]"
                        f"[dim] .. [/dim]"
                        f"[bright_cyan]{indices[-1]}[/bright_cyan][dim])[/dim]"
                    )
                else:
                    for i in indices[:max_items]:
                        gname = f"{prefix}_{i}"
                        branch = tree.add(f"[green]{gname}[/green]")
                        walk(h5group[gname], branch, depth + 1, is_root=False)

                    if count > max_items:
                        tree.add(f"[dim italic]+ {count - max_items} more ...[/dim italic]")

            # --- non-indexed entries ---
            for name in sorted(others):
                obj = h5group[name]

                if isinstance(obj, h5py.Group):
                    n_children = len(obj)
                    branch = tree.add(
                        f"[bold green]{name}[/bold green]  "
                        f"[dim bright_black]{n_children} item{'s' if n_children != 1 else ''}[/dim bright_black]"
                    )
                    walk(obj, branch, depth + 1, is_root=False)

                elif isinstance(obj, h5py.Dataset):
                    shape_str = "x".join(str(d) for d in obj.shape) if obj.shape else "scalar"
                    dnode = tree.add(
                        f"[yellow]{name}[/yellow]  "
                        f"[dim bright_black]{shape_str}[/dim bright_black] "
                        f"[dim]{obj.dtype}[/dim]"
                    )
                    if dataset_meta and obj.attrs:
                        add_attrs(dnode, obj.attrs)

        walk(self.file, root, is_root=True)
        console.print()
        console.print(root)
        console.print()

    def print_file_validity_status(self) -> None:
        """
        Prints the validity status of a file based on its artifact status attribute.

        This method checks the "artifact_status" attribute of the file. If the status
        is "SUCCESS", it displays a success panel indicating that the file is valid.
        Otherwise, it displays an error panel indicating that the file is corrupted.

        When the file is valid but was created by a newer version of Pharmacon than
        the current runtime, a warning panel is shown instead.

        :raises AttributeError: If the "artifact_status" attribute is not found.
        """

        console = Console()

        status = str(self.file.attrs.get("artifact_status", "")).strip().upper()
        file_version = str(self.file.attrs.get("pharmacon_version", "")).strip()
        file_type = str(self.file.attrs.get("file_type", ""))
        runtime_version = str(__version__)

        # ── Determine if the runtime is older than the file ──────────────────
        version_warning = False
        if file_version and runtime_version:
            try:
                f_parts = tuple(map(int, file_version.split(".")))
                r_parts = tuple(map(int, runtime_version.split(".")))
                if f_parts > r_parts:
                    version_warning = True
            except (ValueError, AttributeError):
                pass

        if status != "SUCCESS":
            body = (
                "[bold red]PHARMACON FILE IS CORRUPTED[/bold red]\n"
                f"[dim]{self.path.name}[/dim]"
            )
            if status:
                body += f"\n[dim red]artifact_status: {status}[/dim red]"

            panel = Panel(
                Align.center(body, vertical="middle"),
                border_style="red",
                title="[bold red]Status[/bold red]",
                title_align="left",
                box=HEAVY,
                padding=(1, 2),
                width=min(80, console.width),
            )
        elif version_warning:
            body = (
                "[bold yellow]PHARMACON FILE IS VALID[/bold yellow]\n"
                f"[dim]{self.path.name}[/dim]"
            )
            if file_version:
                body += f"  [dim]v{file_version}[/dim]"
            if file_type:
                body += f"  [dim]({file_type})[/dim]"
            body += (
                f"\n\n[bold yellow]Warning:[/bold yellow] "
                f"[yellow]File was created with Pharmacon v{file_version} "
                f"but the current runtime is v{runtime_version}. "
                f"Some features may not be supported.[/yellow]"
            )

            panel = Panel(
                Align.center(body, vertical="middle"),
                border_style="yellow",
                title="[bold yellow]Status[/bold yellow]",
                title_align="left",
                box=HEAVY,
                padding=(1, 2),
                width=min(80, console.width),
            )
        else:
            body = (
                "[bold green]PHARMACON FILE IS VALID[/bold green]\n"
                f"[dim]{self.path.name}[/dim]"
            )
            if file_version:
                body += f"  [dim]v{file_version}[/dim]"
            if file_type:
                body += f"  [dim]({file_type})[/dim]"

            panel = Panel(
                Align.center(body, vertical="middle"),
                border_style="green",
                title="[bold green]Status[/bold green]",
                title_align="left",
                box=HEAVY,
                padding=(1, 2),
                width=min(80, console.width),
            )

        console.print()
        console.print(panel)

    def __enter__(self):
        """
        Provides the ability to use the object as a context manager.

        This method allows the object to be used with the `with` statement, ensuring
        that proper setup and teardown can occur within the context.

        :return: The current instance of the object
        :rtype: self
        """
        return self

    def __exit__(self, exc_type, exc, tb):
        """
        Handles cleanup operations required when exiting a context managed by the object. This method
        is typically invoked at the end of a `with` statement to ensure that necessary resource
        management is performed correctly, especially for releasing resources.

        :param exc_type: The exception type, if an exception occurred, otherwise None.
        :type exc_type: Optional[Type[BaseException]]
        :param exc: The exception instance, if an exception occurred, otherwise None.
        :type exc: Optional[BaseException]
        :param tb: The traceback object corresponding to the exception, if applicable, otherwise None.
        :type tb: Optional[TracebackType]
        :return: Indicates whether the exception, if any, should be suppressed.
        :rtype: Optional[bool]
        """
        self.close()

    def __del__(self):
        """
        Safely handles the deletion of the instance by attempting to close resources that
        might have been opened during its lifecycle. Ensures robustness by catching any
        exceptions that occur during the cleanup process, avoiding unhandled errors
        during object destruction.

        :return: None
        """
        try:
            self.close()
        except Exception:
            pass

    def __repr__(self):
        """
        Provides a string representation of the object for debugging and logging
        purposes. This method returns a concise and comprehensive string that
        represents the state of the current instance, including its important
        attributes.

        :return: A string representation of the object.
        :rtype: str
        """
        return (
            f"{self.__class__.__name__}("
            f"path={self.path!r}, overwrite={self.overwrite}, mode={self.mode!r})"
        )

    def __str__(self):
        """
        Generate a string representation of the PharmaconFile object.

        This method constructs and returns a human-readable string
        representation of the PharmaconFile instance by displaying its
        associated file path.

        :return: A string in the format "PharmaconFile -> {self.path}" where
            `self.path` is replaced with the file path of the object.
        :rtype: str
        """
        return f"PharmaconFile -> {self.path}"

    def __getitem__(self, key):
        """
        Retrieve an item from the file using the specified key.

        This method allows indexing into the object to retrieve specific items stored
        in the internal file attribute. It behaves like a dictionary or similar
        mapping, returning the value associated with the specified key.

        :param key: The key associated with the value to retrieve.
        :return: The value corresponding to the specified key.
        """
        return self.file[key]

    def __contains__(self, key: str) -> bool:
        """
        Checks whether a given key exists in the file.

        This method determines if the specified key is present in the file.

        :param key: The string key to check for existence in the file.
        :type key: str
        :return: Returns True if the key exists in the file, otherwise False.
        :rtype: bool
        """
        return key in self.file

    def __iter__(self):
        """
        Provides an iterator for the `file` attribute of the class.

        :return: An iterator over the `file` attribute.
        :rtype: iterator
        """
        return iter(self.file)

    def __len__(self):
        """
        Provides the length of the object it represents. This method allows the
        object to be used with the built-in `len()` function, returning the length
        of an associated resource, typically a file or collection.

        :return: The length of the associated resource
        :rtype: int
        """
        return len(self.file)

    def __bool__(self):
        """
        Evaluates whether the current instance represents a valid file.

        This method checks if both the `file` attribute is not `None` and whether
        its `id` is valid. If both conditions are satisfied, it returns `True`,
        otherwise `False`.

        :return: A boolean indicating the validity of the current instance.
        :rtype: bool
        """
        return self.file is not None and self.file.id.valid
