"""Generic Scanpy clustering pipeline: HVG selection -> PCA -> neighbors -> UMAP -> Leiden.

Extracted from the pseudobulk overclustering step in cas.py and parameterized
for reuse on any AnnData with raw integer counts in `.X` or a named layer.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import scanpy as sc
from anndata import AnnData


def cluster_adata(
    adata: AnnData,
    *,
    counts_layer: str | None = None,
    n_top_genes: int = 2000,
    hvg_flavor: str = "seurat_v3",
    n_comps: int = 10,
    neighbors_metric: str = "cosine",
    resolution: float = 1.0,
    cluster_key: str = "leiden",
    random_state: int = 0,
    zero_center: bool = True,
    n_iterations: int = 2,
    compute_umap: bool = True,
    copy: bool = False,
) -> AnnData:
    """Cluster cells with a standard HVG -> PCA -> neighbors -> UMAP -> Leiden pipeline.

    Expects raw integer counts in `adata.X`, or pass `counts_layer` to copy
    counts from `adata.layers[counts_layer]` first. Writes a neighbor graph,
    a `X_umap` embedding (unless `compute_umap=False`), and a `cluster_key`
    column in `.obs`. The intermediate `X_pca` embedding is dropped from
    `obsm` after neighbors are computed, since SCP's ingest validator treats
    every obsm entry as a 2D/3D visualization coordinate array.

    UMAP only affects the `X_umap` embedding used for visualization; it has
    no effect on the `cluster_key` assignments, which are computed from the
    neighbor graph directly. Set `compute_umap=False` to skip it when only
    the cluster labels are needed.

    Set `zero_center=False` for large datasets: mean-centering forces the
    scaling step to densify a sparse matrix, which can exhaust memory well
    before `n_obs` reaches even the low hundreds of thousands.

    `n_iterations` caps how many full Leiden passes are run (each pass keeps
    the best-quality partition found so far); `-1` reruns until an iteration
    fails to improve quality, which is unbounded and can take a very long
    time on large graphs. `2` is the standard value used in current Scanpy
    tutorials and the sc-best-practices docs.
    """

    adata = adata.copy() if copy else adata

    if counts_layer is not None:
        adata.X = adata.layers[counts_layer].copy()

    sc.pp.highly_variable_genes(adata, flavor=hvg_flavor, n_top_genes=n_top_genes)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10, zero_center=zero_center)
    sc.pp.pca(adata, use_highly_variable=True, n_comps=n_comps, random_state=random_state)
    sc.pp.neighbors(adata, method="umap", metric=neighbors_metric, random_state=random_state)
    if compute_umap:
        # sc.tl.umap() re-derives the neighbors' representation from
        # obsm["X_pca"] internally, so it must still be present here.
        sc.tl.umap(adata, random_state=random_state)
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=cluster_key,
        random_state=random_state,
        n_iterations=n_iterations,
    )

    # SCP's ingest validator treats every obsm entry as a visualization
    # coordinate array and rejects anything but 2 or 3 columns, so the
    # n_comps-dimensional PCA embedding can't be left in obsm. Everything
    # that needs it (neighbors, UMAP) has already run by this point.
    del adata.obsm["X_pca"]

    return adata


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Cluster an H5AD with a standard HVG -> PCA -> neighbors -> Leiden pipeline."
    )
    parser.add_argument("h5ad_path", help="Path to an H5AD with raw integer counts")
    parser.add_argument(
        "--counts-layer",
        default=None,
        help="Layer to copy raw counts from before clustering (default: use .X as-is)",
    )
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--hvg-flavor", default="seurat_v3")
    parser.add_argument("--n-comps", type=int, default=10)
    parser.add_argument("--neighbors-metric", default="cosine")
    parser.add_argument("--resolution", type=float, default=1.0)
    parser.add_argument("--cluster-key", default="leiden")
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument(
        "--n-iterations",
        type=int,
        default=2,
        help="Leiden passes to run; -1 reruns until quality stops improving "
        "(unbounded, can be very slow on large graphs)",
    )
    parser.add_argument(
        "--no-zero-center",
        dest="zero_center",
        action="store_false",
        help="Skip mean-centering during scaling, keeping the matrix sparse "
        "(recommended for large datasets, e.g. hundreds of thousands+ of cells)",
    )
    parser.set_defaults(zero_center=True)
    parser.add_argument(
        "--no-umap",
        dest="compute_umap",
        action="store_false",
        help="Skip computing the UMAP embedding, since it only affects "
        "visualization and not the cluster_key assignments",
    )
    parser.set_defaults(compute_umap=True)
    parser.add_argument(
        "--output",
        default=None,
        help="Output H5AD path (default: <input>.clustered.h5ad)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    h5ad_path = Path(args.h5ad_path)
    if not h5ad_path.exists():
        print(f"error: H5AD not found: {h5ad_path}", file=sys.stderr)
        return 1

    output_path = (
        Path(args.output)
        if args.output
        else h5ad_path.with_suffix("").with_suffix(".clustered.h5ad")
    )

    adata = sc.read(h5ad_path)

    # cluster_adata() overwrites .X with normalized/log1p/scaled values, so
    # stash the raw counts first and restore them before writing; otherwise
    # the output H5AD loses raw counts and can densify under zero-centering.
    counts_layer = args.counts_layer
    stashed_counts = counts_layer is None
    if stashed_counts:
        counts_layer = "counts"
        adata.layers[counts_layer] = adata.X.copy()

    cluster_adata(
        adata,
        counts_layer=counts_layer,
        n_top_genes=args.n_top_genes,
        hvg_flavor=args.hvg_flavor,
        n_comps=args.n_comps,
        neighbors_metric=args.neighbors_metric,
        resolution=args.resolution,
        cluster_key=args.cluster_key,
        random_state=args.random_state,
        zero_center=args.zero_center,
        n_iterations=args.n_iterations,
        compute_umap=args.compute_umap,
    )

    adata.X = adata.layers[counts_layer]
    if stashed_counts:
        del adata.layers[counts_layer]
    adata.write(output_path)

    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
