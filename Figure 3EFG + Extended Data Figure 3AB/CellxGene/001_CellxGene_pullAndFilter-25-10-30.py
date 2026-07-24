#!/usr/bin/env python3
import numpy as np
import pandas as pd
import cellxgene_census
import os
import math
import json

############################################
# CONFIG
############################################

OUT_TSV = "census_v4_summary_incremental.tsv"

# ----------------------------
# Gene list
# ----------------------------
GENE_ENSEMBL_IDS = [
"ENSG00000198089", #SFI1
    "ENSG00000147862", #NFIB
    "ENSG00000105669", #COPE
    "ENSG00000118689", #FOXO3
    "ENSG00000169554", #ZEB2
    "ENSG00000123983", #ACSL3
    "ENSG00000068366", #ACSL4
    "ENSG00000137936", #BCAR3
    "ENSG00000091831", #ESR1
    "ENSG00000114353", #GNAI2
    "ENSG00000151623", #NR3C2
    "ENSG00000160113", #NR2F6
    "ENSG00000175745", #NR2F1
    "ENSG00000157933", #SKI
    "ENSG00000172819", #RARG
    "ENSG00000112561", #TFEB
    "ENSG00000186350", #RXRA
    "ENSG00000143171", #RXRG
    "ENSG00000120068", #HOXB8
    "ENSG00000037965", #HOXC8
    "ENSG00000010818", #HIVEP2
    "ENSG00000132170", #PPARG
    "ENSG00000148516", #ZEB1
    "ENSG00000185630", #PBX1
    "ENSG00000167081", #PBX3
    "ENSG00000244405", #ETV5
    "ENSG00000163219", #ARHGAP25
    "ENSG00000213923", #CSNK1E
    "ENSG00000117650", #NEK2
    "ENSG00000173744", #AGFG1
    "ENSG00000182809", #CRIP2
    "ENSG00000160199", #PKNOX1
    "ENSG00000165495", #PKNOX2
    "ENSG00000104518", #GSDMD
    "ENSG00000113430", #IRX4
    "ENSG00000182568", #SATB1
    "ENSG00000119042", #SATB2
    "ENSG00000085276", #MECOM
    "ENSG00000142611", #PRDM16
    "ENSG00000117399", #CDC20
    "ENSG00000054598", #FOXC1
    "ENSG00000167670", #CHAF1A
    "ENSG00000145604", #SKP2
    "ENSG00000013619", #MAMLD1
    "ENSG00000102554", #KLF5
    "ENSG00000007968", #E2F2
    "ENSG00000112242", #E2F3
    "ENSG00000187079", #TEAD1
    "ENSG00000074219", #TEAD2
    "ENSG00000007866", #TEAD3
    "ENSG00000197905", #TEAD4
    "ENSG00000160051", #IQCC
    "ENSG00000184661", #CDCA2
    "ENSG00000083307", #GRHL2
    "ENSG00000160551", #TAOK1
    "ENSG00000123607", #TTC21B
    "ENSG00000184384", #MAML2
    "ENSG00000196323", #ZBTB44
    "ENSG00000137812", #KNL1
    "ENSG00000158411", #MITD1
    "ENSG00000101493", #ZNF516
    "ENSG00000087510", #TFAP2C
    "ENSG00000137203", #TFAP2A
    "ENSG00000156787", #TBC1D31
    "ENSG00000196233", #LCOR
    "ENSG00000120647", #CCDC77
    "ENSG00000100505", #TRIM9
    "ENSG00000018408", #WWTR1
    "ENSG00000119866", #BCL11A
    "ENSG00000151461", #UPF2
    "ENSG00000151849", #CPAP
    "ENSG00000163257", #DCAF16
    "ENSG00000075218", #GTSE1
    "ENSG00000134317", #GRHL1
    "ENSG00000122008", #POLK
    "ENSG00000104447", #TRPS1
    "ENSG00000161405", #IKZF3
    "ENSG00000108175", #ZMIZ1
    "ENSG00000143842", #SOX13
    "ENSG00000124226", #RNF114
    "ENSG00000169083", #AR
    "ENSG00000082175", #PGR
    "ENSG00000106004", #HOXA5
    "ENSG00000049768", #FOXP3
    "ENSG00000082014", #SMARCD3
    "ENSG00000137818", #RPLP1
    "ENSG00000131469", #RPLP27
    "ENSG00000075624", #ACTB
    "ENSG00000111640", #GAPDH
]

# ----------------------------
# Tissues and filters
# ----------------------------
TISSUES_VALUEFILTER = {
    "breast": "tissue_general == 'breast'",
    "blood":  "tissue_general == 'blood'",
    "ureter": "tissue_general == 'ureter'",
}

# ----------------------------
# Assay whitelist
# ----------------------------
ALLOWED_ASSAY_EFO = {
    "EFO:0010550", "EFO:0009901", "EFO:0011025", "EFO:0009899", "EFO:0009900",
    "EFO:0009922", "EFO:0030003", "EFO:0030004", "EFO:0008995", "EFO:0008919",
    "EFO:0008722", "EFO:0010010",
}

# ----------------------------
# Gene-length-normalised assays (Smart-seq, Smart-seq2)
# ----------------------------
ASSAYS_FOR_GENE_LENGTH_NORMALIZATION = {
    "EFO:0008930",  # Smart-seq
    "EFO:0008853",  # Smart-seq2
}

CHUNK_SIZE = 25000
UMI_MIN_FOR_POSITIVE = 2.0  # raw_count >= 2
MIN_GENES_PER_CELL = 500
TARGET_LIBRARY_SIZE = 10_000.0

# ontology rollup 
CELL_TYPE_CL_GROUPS = {
    "luminal_mammary_epithelial_cell": [
        "CL:0002326",
        "CL:4033057",
        "CL:4033058",
    ],
    "urothelial_cell": [
        "CL:0000731",
        "CL:1001428",
        "CL:1000486",
        "CL:1000703",
        "CL:4030055",
        "CL:4030056",
    ],
    "cd4_tcell": [
        "CL:0000624",
        "CL:0000492",
        "CL:0000896",
        "CL:0000792",
        "CL:0000895",
        "CL:0000897",
        "CL:0000934",
        "CL:0001043",
        "CL:0001044",
        "CL:0000545",
        "CL:0000546",
        "CL:0000899",
        "CL:0001042",
        "CL:0002038",
        "CL:0000903",
        "CL:0000904",
        "CL:0000905",
        "CL:0002677",
        "CL:0002678",
        "CL:4033038",
    ],
}


############################################
# HELPERS
############################################

def open_census():
    return cellxgene_census.open_soma()

def load_var_hs(census):
    """Load gene metadata for Homo sapiens including feature_length."""
    return cellxgene_census.get_var(
        census,
        "homo_sapiens",
        column_names=["soma_joinid", "feature_id", "feature_name", "feature_length"],
    )

def build_gene_index(var, ensembl_ids):
    """Map ENSG prefix → soma_joinid"""
    mapper = {}
    for ens in ensembl_ids:
        hit = var[var["feature_id"].str.startswith(ens)]
        if not hit.empty:
            mapper[ens] = int(hit["soma_joinid"].iloc[0])
    return mapper

def load_obs_for_tissue(census, tissue_filter_expr):
    """Load QC-filtered obs for a tissue."""
    obs_df = cellxgene_census.get_obs(
        census,
        "homo_sapiens",
        value_filter=f"({tissue_filter_expr}) and is_primary_data==True and disease == 'normal'",
        column_names=[
            "soma_joinid", "dataset_id", "tissue_general",
            "cell_type_ontology_term_id", "assay_ontology_term_id",
            "raw_sum", "nnz", "disease", "is_primary_data"
        ],
    )
    obs_df = obs_df[obs_df["assay_ontology_term_id"].isin(ALLOWED_ASSAY_EFO)].copy()
    obs_df = obs_df[obs_df["nnz"] >= MIN_GENES_PER_CELL].copy()
    obs_df = obs_df[obs_df["raw_sum"].notna() & (obs_df["raw_sum"] > 0)].copy()
    return obs_df

############################################
# NORMALISATION — matches CellxGene pipeline
############################################

def compute_logCPTT(raw_count, libsize, assay_id, gene_length):
    """
    Implements the exact per-cell normalisation:
        if assay in Smart-seq/Smart-seq2:
            expr = (raw / gene_length) * TARGET_LIBRARY_SIZE / libsize
        else:
            expr = raw * TARGET_LIBRARY_SIZE / libsize
        log_expr = ln(expr + 1)
    """
    if libsize <= 0:
        return 0.0

    if assay_id in ASSAYS_FOR_GENE_LENGTH_NORMALIZATION and gene_length > 0:
        expr = (raw_count / gene_length) * TARGET_LIBRARY_SIZE / libsize
    else:
        expr = raw_count * TARGET_LIBRARY_SIZE / libsize

    return math.log(expr + 1.0)

############################################
# Summarisation
############################################

def summarize_gene_for_cells_per_dataset_datasetavg(
    raw_X, cells_df, gene_joinid, gene_length,
    chunk_size=25000, umi_min_for_positive=2.0,
):
    """Same logic as your script, now with CellxGene normalisation."""
    if cells_df.empty:
        return dict(
            total_cells=0, pos_cells=0, pct_pos=0.0,
            mean_logCPTT_pos_dataset_avg=0.0,
            mean_logCPTT_allcells_dataset_avg=0.0,
        )

    total_cells_all = 0
    pos_cells_all = 0
    dataset_pos_means, dataset_allcell_means = [], []

    for ds_id, df_ds in cells_df.groupby("dataset_id", observed=False):
        soma_ids_ds = df_ds["soma_joinid"].to_numpy()
        raw_sum_ds = df_ds["raw_sum"].to_numpy()
        assay_ds = df_ds["assay_ontology_term_id"].to_numpy()
        soma_to_rawsum = dict(zip(soma_ids_ds, raw_sum_ds))
        soma_to_assay = dict(zip(soma_ids_ds, assay_ds))
        n_cells_ds = len(soma_ids_ds)
        total_cells_all += n_cells_ds

        pos_cells_ds = 0
        sum_log_expr_pos_ds = 0.0
        sum_log_expr_all_ds = 0.0

        for start in range(0, n_cells_ds, chunk_size):
            end = min(start + chunk_size, n_cells_ds)
            soma_chunk = soma_ids_ds[start:end]
            rr = raw_X.read(coords=(soma_chunk, np.array([gene_joinid])))

            chunk_counts = {cid: 0.0 for cid in soma_chunk}
            for tbl in rr.tables():
                if tbl.num_rows == 0:
                    continue
                pdf = tbl.to_pandas()
                for cell_id, val in zip(pdf["soma_dim_0"], pdf["soma_data"]):
                    chunk_counts[cell_id] = float(val)

            for cell_id, raw_count in chunk_counts.items():
                lib = soma_to_rawsum[cell_id]
                assay = soma_to_assay[cell_id]
                log_expr = compute_logCPTT(raw_count, lib, assay, gene_length)

                sum_log_expr_all_ds += log_expr
                if raw_count >= umi_min_for_positive:
                    pos_cells_ds += 1
                    sum_log_expr_pos_ds += log_expr

        # dataset-level means
        if pos_cells_ds > 0:
            dataset_pos_means.append(sum_log_expr_pos_ds / pos_cells_ds)
        if n_cells_ds > 0:
            dataset_allcell_means.append(sum_log_expr_all_ds / n_cells_ds)
        pos_cells_all += pos_cells_ds

    pct_pos = (pos_cells_all / total_cells_all * 100.0) if total_cells_all else 0.0
    mean_logCPTT_pos_dataset_avg = np.mean(dataset_pos_means) if dataset_pos_means else 0.0
    mean_logCPTT_allcells_dataset_avg = np.mean(dataset_allcell_means) if dataset_allcell_means else 0.0

    return dict(
        total_cells=int(total_cells_all),
        pos_cells=int(pos_cells_all),
        pct_pos=float(pct_pos),
        mean_logCPTT_pos_dataset_avg=float(mean_logCPTT_pos_dataset_avg),
        mean_logCPTT_allcells_dataset_avg=float(mean_logCPTT_allcells_dataset_avg),
    )

############################################
# MAIN
############################################

def main():
    census = open_census()
    var = load_var_hs(census)
    raw_X = census["census_data"]["homo_sapiens"]["ms"]["RNA"]["X"]["raw"]

    gene_map = build_gene_index(var, GENE_ENSEMBL_IDS)

    print("Resolved gene indices:")
    for ens, gj in gene_map.items():
        print(f"  {ens} -> soma_joinid {gj}")

    write_header = not os.path.exists(OUT_TSV)
    with open(OUT_TSV, "a") as fh_out:
        if write_header:
            fh_out.write("\t".join([
                "tissue", "cell_type_group", "cl_ids_grouped",
                "gene_ensembl", "gene_soma_joinid",
                "total_cells", "pos_cells", "pct_pos",
                "mean_logCPTT_pos_dataset_avg",
                "mean_logCPTT_allcells_dataset_avg",
            ]) + "\n")
            fh_out.flush()

        for ens, gj in gene_map.items():
            gene_len = float(var.loc[var["soma_joinid"] == gj, "feature_length"].iloc[0])
            print(f"\n=== GENE {ens} (soma_joinid={gj}, len={gene_len}) ===")
            for tissue_label, tissue_filter_expr in TISSUES_VALUEFILTER.items():
                print(f"   Tissue: {tissue_label}")
                obs_tissue = load_obs_for_tissue(census, tissue_filter_expr)
                print(f"      rows after filter: {obs_tissue.shape[0]}")

                for group_label, cl_list in CELL_TYPE_CL_GROUPS.items():
                    cells_df = obs_tissue[
                        obs_tissue["cell_type_ontology_term_id"].isin(cl_list)
                    ].copy()
                    n_cells = cells_df.shape[0]
                    print(f"      {group_label} -> {n_cells} cells (CLs: {','.join(cl_list)})")

                    if n_cells == 0:
                        stats = dict(
                            total_cells=0, pos_cells=0, pct_pos=0.0,
                            mean_logCPTT_pos_dataset_avg=0.0,
                            mean_logCPTT_allcells_dataset_avg=0.0,
                        )
                    else:
                        stats = summarize_gene_for_cells_per_dataset_datasetavg(
                            raw_X, cells_df, gj, gene_len,
                            chunk_size=CHUNK_SIZE,
                            umi_min_for_positive=UMI_MIN_FOR_POSITIVE,
                        )

                    fh_out.write("\t".join([
                        tissue_label, group_label, ",".join(cl_list),
                        ens, str(gj),
                        str(stats["total_cells"]), str(stats["pos_cells"]),
                        f"{stats['pct_pos']:.6f}",
                        f"{stats['mean_logCPTT_pos_dataset_avg']:.6f}",
                        f"{stats['mean_logCPTT_allcells_dataset_avg']:.6f}",
                    ]) + "\n")
                    fh_out.flush()

    print(f"\nDone. Wrote incremental results to {OUT_TSV}")

if __name__ == "__main__":
    main()

