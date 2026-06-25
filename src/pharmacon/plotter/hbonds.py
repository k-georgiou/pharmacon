"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Hydrogen-bond-specific plots.

The residue-pair plots (heatmap, pair timeline) reuse the protein–protein
plotters; this module adds the H-bond-specific summaries:

- :func:`plot_hbonds_count_per_frame_from_file` — number of H-bonds vs frame.
- :func:`plot_hbonds_occupancy_from_file`        — ranked per-pair occupancy.
- :func:`plot_hbonds_network_from_file`          — residue H-bond network graph.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pharmacon.logger import get_logger
from pharmacon.constants.smarts import AA3_to_AA1
from pharmacon.plotter.universal import _get_x_value


__all__ = [
    "plot_hbonds_count_per_frame_from_file",
    "plot_hbonds_occupancy_from_file",
    "plot_hbonds_network_from_file",
]

logger = get_logger(__name__)


def plot_hbonds_count_per_frame_from_file(pta_file, *, group_name: str,
                                          settings, out_dir: Path,
                                          is_merged: bool = False) -> None:
    """
    Plot the number of hydrogen bonds detected in each frame.

    Reads the per-frame ``frame_<N>/interactions`` datasets (using the stored
    ``n_interactions`` attribute, falling back to the row count) and draws a
    single time series of H-bond counts.
    """
    if group_name not in pta_file.file:
        raise RuntimeError(f"Group not found: {group_name}")

    x_axis = str(getattr(settings, "x_axis", "frame_index"))

    # Auto-label the x-axis from the chosen unit when the user kept the default
    # "Frame" label (mirrors choosing time_ns/time_ps on the RMSD plot).
    _AXIS_LABELS = {"frame_index": "Frame", "time_ps": "Time (ps)",
                    "time_ns": "Time (ns)", "time_us": "Time (µs)"}
    x_label = settings.x_label
    if x_label == "Frame" and x_axis in _AXIS_LABELS:
        x_label = _AXIS_LABELS[x_axis]

    xs: list[float] = []
    ys: list[int] = []
    for frame in pta_file._iter_frames(group_name):
        dpath = f"{group_name}/frame_{frame}/interactions"
        if dpath not in pta_file.file:
            continue
        dset = pta_file.file[dpath]
        attrs = dset.attrs
        try:
            count = int(attrs["n_interactions"])
        except (KeyError, ValueError):
            count = int(getattr(dset, "size", 0))
        xs.append(_get_x_value(attrs, frame, x_axis))
        ys.append(count)

    if not xs:
        raise RuntimeError("No frames with interaction data found.")

    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    order = np.argsort(x)
    x, y = x[order], y[order]

    fig, ax = plt.subplots(
        figsize=(settings.fig_size_width, settings.fig_size_height),
        dpi=settings.fig_dpi,
        facecolor=settings.bg_color,
    )
    plt.rcParams["font.family"] = settings.font_family

    color = settings.line_colors[0] if settings.line_colors else "#1f77b4"
    ax.plot(x, y, color=color, linewidth=settings.line_width,
            alpha=settings.line_alpha, linestyle=settings.line_style)

    if not settings.disable_title:
        ax.set_title(settings.fig_title, fontsize=settings.font_size_title,
                     fontweight=settings.font_weight_title)
    ax.set_xlabel(x_label, fontsize=settings.font_size_label,
                  fontweight=settings.font_weight_label)
    ax.set_ylabel(settings.y_label, fontsize=settings.font_size_label,
                  fontweight=settings.font_weight_label)
    ax.tick_params(labelsize=settings.font_size_ticks)

    if settings.enable_grid:
        ax.grid(True, linestyle=settings.grid_style, alpha=settings.grid_alpha)

    if settings.tight_layout:
        fig.tight_layout()

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = (out_dir / f"{group_name}_count_per_frame").with_suffix(f".{settings.fig_format}")
    fig.savefig(out_path, dpi=settings.fig_dpi, transparent=settings.fig_transparent)
    plt.close(fig)
    logger.info("Saved H-bond count-per-frame plot → %s", out_path)


def _pair_label(r1, r2, settings) -> str:
    """Build a 'resA–resB' label for a residue pair using the representation."""
    def one(res) -> str:
        resname, resid, chainid, segid = res[0], res[1], res[2], res[3]
        if settings.aa3_to_aa1:
            resname = AA3_to_AA1.get(str(resname).upper(), resname)
        label = settings.representation
        label = label.replace("resname", str(resname))
        label = label.replace("resid", str(resid))
        label = label.replace("chainid", str(chainid))
        label = label.replace("segid", str(segid))
        return label
    return f"{one(r1)}–{one(r2)}"


def plot_hbonds_occupancy_from_file(pta_file, *, group_name: str, mode_name: str,
                                    settings, out_dir: Path,
                                    attach_to_name: str = "",
                                    is_merged: bool = False) -> None:
    """
    Ranked horizontal bar of per-pair hydrogen-bond occupancy.

    Reads ``modes/<mode_name>/table`` (use mode2, the once-per-frame mode, so
    ``frequency`` is a true 0–1 occupancy), keeps the top-N pairs above the
    occupancy threshold, and draws them strongest-at-top.
    """
    sub = "modes_merged" if is_merged else "modes"
    dset_path = f"{group_name}/{sub}/{mode_name}/table"
    if dset_path not in pta_file.file:
        raise RuntimeError(f"Dataset not found: {dset_path}")

    dset = pta_file.file[dset_path]
    if getattr(dset, "size", 0) == 0:
        raise RuntimeError(f"Dataset empty: {dset_path}")

    freq_field = "mean_frequency" if is_merged else "frequency"

    pairs: dict[str, float] = {}
    for row in dset:
        rec = json.loads(row)
        r1, r2, _label = rec["key"]
        occ = float(rec[freq_field])
        if occ < settings.threshold:
            continue
        label = _pair_label(r1, r2, settings)
        # Same residue pair may recur per interaction sub-type — keep the max.
        pairs[label] = max(pairs.get(label, 0.0), occ)

    if not pairs:
        raise RuntimeError("No pairs above occupancy threshold.")

    items = sorted(pairs.items(), key=lambda kv: kv[1], reverse=True)
    if settings.top_n > 0:
        items = items[:settings.top_n]

    labels = [k for k, _ in items]
    values = [v for _, v in items]

    fig, ax = plt.subplots(
        figsize=(settings.fig_size_width, settings.fig_size_height),
        dpi=settings.fig_dpi,
        facecolor=settings.bg_color,
    )
    plt.rcParams["font.family"] = settings.font_family

    ypos = np.arange(len(labels))
    # Draw strongest at the top.
    ax.barh(ypos, values, color=settings.bar_color, alpha=settings.bar_alpha,
            edgecolor=settings.bar_edge_color, linewidth=settings.bar_edge_width)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=settings.font_size_ticks)
    ax.invert_yaxis()

    if not settings.disable_title:
        ax.set_title(settings.fig_title, fontsize=settings.font_size_title,
                     fontweight=settings.font_weight_title)
    ax.set_xlabel(settings.x_label, fontsize=settings.font_size_label,
                  fontweight=settings.font_weight_label)
    ax.set_ylabel(settings.y_label, fontsize=settings.font_size_label,
                  fontweight=settings.font_weight_label)
    ax.tick_params(labelsize=settings.font_size_ticks)

    if settings.enable_grid:
        ax.grid(True, axis="x", linestyle=settings.grid_style, alpha=settings.grid_alpha)

    if settings.tight_layout:
        fig.tight_layout()

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = (out_dir / f"{attach_to_name}{settings.fig_basename}").with_suffix(f".{settings.fig_format}")
    fig.savefig(out_path, dpi=settings.fig_dpi, transparent=settings.fig_transparent)
    plt.close(fig)
    logger.info("Saved H-bond occupancy plot → %s", out_path)


def _residue_label(res, settings) -> str:
    """Build a node label for a single residue using the representation."""
    resname, resid, chainid, segid = res[0], res[1], res[2], res[3]
    if settings.aa3_to_aa1:
        resname = AA3_to_AA1.get(str(resname).upper(), resname)
    label = settings.representation
    label = label.replace("resname", str(resname))
    label = label.replace("resid", str(resid))
    label = label.replace("chainid", str(chainid))
    label = label.replace("segid", str(segid))
    return label


def plot_hbonds_network_from_file(pta_file, *, group_name: str,
                                  settings, out_dir: Path,
                                  is_merged: bool = False) -> None:
    """
    Draw the hydrogen-bond network: nodes are residues, edges are H-bonds
    weighted by occupancy (read from ``settings.mode``, default mode2 = 0–1).

    Tuned to stay readable on large proteins: an occupancy ``threshold`` keeps
    only persistent bonds, ``top_n`` caps the edges, ``min_seq_sep`` can drop
    near-sequence backbone bonds, isolated residues are never drawn (only nodes
    that have an edge appear), and only the ``label_top_n`` highest-degree hubs
    are labelled.
    """
    import matplotlib as mpl
    import networkx as nx

    sub = "modes_merged" if is_merged else "modes"
    dset_path = f"{group_name}/{sub}/{settings.mode}/table"
    if dset_path not in pta_file.file:
        raise RuntimeError(f"Dataset not found: {dset_path}")

    dset = pta_file.file[dset_path]
    if getattr(dset, "size", 0) == 0:
        raise RuntimeError(f"Dataset empty: {dset_path}")

    freq_field = "mean_frequency" if is_merged else "frequency"

    edges: dict = {}            # (node_i, node_j) -> occupancy (max)
    node_labels: dict = {}      # node key -> display label
    for row in dset:
        rec = json.loads(row)
        r1, r2, _label = rec["key"]
        occ = float(rec[freq_field])
        if occ < settings.threshold:
            continue
        k1 = (str(r1[2]), int(r1[1]))   # (chain, resid)
        k2 = (str(r2[2]), int(r2[1]))
        if k1 == k2:
            continue
        if (settings.min_seq_sep > 0 and k1[0] == k2[0]
                and abs(k1[1] - k2[1]) <= settings.min_seq_sep):
            continue
        ekey = tuple(sorted((k1, k2)))
        edges[ekey] = max(edges.get(ekey, 0.0), occ)
        node_labels[k1] = _residue_label(r1, settings)
        node_labels[k2] = _residue_label(r2, settings)

    if not edges:
        raise RuntimeError("No H-bonds above the occupancy threshold.")

    items = sorted(edges.items(), key=lambda kv: kv[1], reverse=True)
    if settings.top_n > 0:
        items = items[:settings.top_n]

    graph = nx.Graph()
    for (a, b), w in items:
        graph.add_edge(a, b, weight=w)

    if settings.largest_component and graph.number_of_nodes() > 0:
        largest = max(nx.connected_components(graph), key=len)
        graph = graph.subgraph(largest).copy()

    if graph.number_of_edges() == 0:
        raise RuntimeError("Network empty after filtering.")

    # Layout
    if settings.layout == "circular":
        pos = nx.circular_layout(graph)
    elif settings.layout == "kamada_kawai":
        pos = nx.kamada_kawai_layout(graph, weight="weight")
    else:
        pos = nx.spring_layout(graph, weight="weight", seed=settings.seed)

    # Node sizes scale with degree.
    degrees = dict(graph.degree())
    d_lo = min(degrees.values())
    d_hi = max(degrees.values())

    def _nsize(d: int) -> float:
        if d_hi == d_lo:
            return 0.5 * (settings.node_size_min + settings.node_size_max)
        t = (d - d_lo) / (d_hi - d_lo)
        return settings.node_size_min + t * (settings.node_size_max - settings.node_size_min)

    nodes = list(graph.nodes())
    node_sizes = [_nsize(degrees[n]) for n in nodes]

    # Edge occupancy → color + width.
    e_occ = [graph[u][v]["weight"] for u, v in graph.edges()]
    e_widths = [settings.edge_width_min
                + w * (settings.edge_width_max - settings.edge_width_min)
                for w in e_occ]
    cmap = mpl.colormaps[settings.edge_cmap]

    fig, ax = plt.subplots(
        figsize=(settings.fig_size_width, settings.fig_size_height),
        dpi=settings.fig_dpi,
        facecolor=settings.bg_color,
    )
    plt.rcParams["font.family"] = settings.font_family

    nx.draw_networkx_edges(
        graph, pos, ax=ax, edge_color=e_occ, edge_cmap=cmap,
        edge_vmin=0.0, edge_vmax=1.0, width=e_widths, alpha=settings.edge_alpha,
    )
    nx.draw_networkx_nodes(
        graph, pos, ax=ax, nodelist=nodes, node_size=node_sizes,
        node_color=settings.node_color, edgecolors=settings.node_edge_color,
        linewidths=0.5,
    )

    # Label only the highest-degree hubs (or all / none).
    if settings.label_top_n != 0:
        if settings.label_top_n < 0:
            label_nodes = nodes
        else:
            label_nodes = [n for n, _ in sorted(
                degrees.items(), key=lambda kv: kv[1], reverse=True
            )[:settings.label_top_n]]
        nx.draw_networkx_labels(
            graph, pos, ax=ax,
            labels={n: node_labels.get(n, "") for n in label_nodes},
            font_size=settings.font_size_labels,
        )

    if not settings.disable_colorbar:
        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.Normalize(0.0, 1.0))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label("Occupancy", fontsize=settings.font_size_cbar)
        cbar.ax.tick_params(labelsize=settings.font_size_cbar)

    ax.set_axis_off()
    if not settings.disable_title:
        ax.set_title(settings.fig_title, fontsize=settings.font_size_title,
                     fontweight=settings.font_weight_title)

    if settings.tight_layout:
        fig.tight_layout()

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = (out_dir / settings.fig_basename).with_suffix(f".{settings.fig_format}")
    fig.savefig(out_path, dpi=settings.fig_dpi, transparent=settings.fig_transparent)
    plt.close(fig)
    logger.info("Saved H-bond network plot → %s (%d nodes, %d edges)",
                out_path, graph.number_of_nodes(), graph.number_of_edges())
