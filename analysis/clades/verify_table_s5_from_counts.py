#!/usr/bin/env python3
"""Count-only reproduction of current SEN genomics main:15_amr_byclade.py.
Source blob: a937dee96017627756145183f0e26abf203fedfb.
The exact omnibus/posthoc/sorting and helper blocks are copied from main.
Only sample-dependent outputs are omitted; 42-gene-only BH scope is explicit.
"""
import argparse
import itertools
import os
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact, MonteCarloMethod
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
MIN_TOTAL_PRESENT = 10
ALPHA = 0.05
CLADE_PALETTE_SPLIT = {"Clade 1A":"#FF7F0E","Clade 1B":"#FF0000","Clade 2":"#0000FF","Clade 3":"#008000","Unassigned":"#808080"}
CLADE_PALETTE_COMBINED = {"Clade 1":"#FF0000","Clade 2":"#0000FF","Clade 3":"#008000","Unassigned":"#808080"}
N = {"Clade 1A":169,"Clade 1B":759,"Clade 2":894,"Clade 3":1357,"Unassigned":127}

SUPERSCRIPT = {
    "a": "ᵃ", "b": "ᵇ", "c": "ᶜ", "d": "ᵈ", "e": "ᵉ", "f": "ᶠ", "g": "ᵍ", "h": "ʰ", "i": "ᶦ",
    "j": "ʲ", "k": "ᵏ", "l": "ˡ", "m": "ᵐ", "n": "ⁿ", "o": "ᵒ", "p": "ᵖ", "r": "ʳ", "s": "ˢ",
    "t": "ᵗ", "u": "ᵘ", "v": "ᵛ", "w": "ʷ", "x": "ˣ", "y": "ʸ", "z": "ᶻ"
}

def letters_to_sup(letters: str) -> str:
    return "".join(SUPERSCRIPT.get(ch, ch) for ch in letters)

def pct_format(count: int, denom: int) -> str:
    if denom <= 0:
        return ""
    if count == 0:
        return "0%"
    pct = 100.0 * count / denom
    if 0 < pct < 0.1:
        return "<0.1%"
    if abs(pct - round(pct)) < 1e-12:
        return f"{int(round(pct))}%"
    return f"{pct:.1f}%"

def format_cell(count: int, denom: int, sup_letters: str = "") -> str:
    p = pct_format(count, denom)
    if p == "":
        return ""
    if sup_letters:
        return f"{count} ({p}){letters_to_sup(sup_letters)}"
    return f"{count} ({p})"

def chi2_with_optional_montecarlo(table: np.ndarray):
    result = chi2_contingency(table, correction=False)
    expected = result.expected_freq

    if np.any(expected < 5):
        # Fixed seed + fixed number of resamples makes sparse-count handling
        # reproducible. SciPy >=1.15 requires a MonteCarloMethod instance.
        mc = MonteCarloMethod(
            n_resamples=5000,
            rng=np.random.default_rng(12345),
        )
        mc_result = chi2_contingency(
            table,
            correction=False,
            method=mc,
        )
        return (
            float(mc_result.pvalue),
            "chi-square (montecarlo, 5000 resamples, seed=12345)",
            float(mc_result.statistic),
            np.nan,
        )

    return (
        float(result.pvalue),
        "chi-square",
        float(result.statistic),
        int(result.dof),
    )

def odds_ratio_with_haldene(a, b, c, d, add=0.5):
    return ((a + add) / (b + add)) / ((c + add) / (d + add))

def make_cld_letters(groups, p_adj_map, alpha=0.05):
    if not p_adj_map or len(groups) == 0:
        return {g: "" for g in groups}

    def different(g1, g2):
        if g1 == g2:
            return False
        key = tuple(sorted((g1, g2)))
        p = p_adj_map.get(key, np.nan)
        return (pd.notna(p) and p <= alpha)

    letter_pool = list("abcdefghijklmnopqrstuvwxyz")
    assign = {g: set() for g in groups}
    remaining = set(groups)

    while remaining and letter_pool:
        L = letter_pool.pop(0)
        block = []
        for g in groups:
            if g not in remaining:
                continue
            ok = True
            for h in block:
                if different(g, h):
                    ok = False
                    break
            if ok:
                block.append(g)

        if not block:
            g = next(iter(remaining))
            block = [g]

        for g in block:
            assign[g].add(L)
            if g in remaining:
                remaining.remove(g)

    return {g: "".join(sorted(assign[g])) for g in groups}

def get_scheme_config(scheme: str, mode: str):
    include_unassigned = (mode == "With_Unassigned")

    if scheme == "Combined_Clade1":
        clade_order = ["Clade 1", "Clade 2", "Clade 3"]
        palette = CLADE_PALETTE_COMBINED.copy()
    elif scheme == "Split_Clade1A1B":
        clade_order = ["Clade 1A", "Clade 1B", "Clade 2", "Clade 3"]
        palette = CLADE_PALETTE_SPLIT.copy()
    else:
        raise ValueError(f"Unknown scheme: {scheme}")

    if include_unassigned:
        clade_order = clade_order + ["Unassigned"]

    return clade_order, palette, include_unassigned

def build_output_paths(base_out: str):
    return {
        "excel": os.path.join(base_out, "AMRFinder_ByClade_Summary.xlsx"),
        "table": os.path.join(base_out, "AMRFinder_ByClade_ElementCounts.tsv"),
        "matrix": os.path.join(base_out, "AMRFinder_PresenceAbsence_Matrix.tsv"),
        "badlines": os.path.join(base_out, "master_bad_lines.tsv"),
        "heat_png": os.path.join(base_out, "AMRFinder_Heatmap_ByType.png"),
        "heat_pdf": os.path.join(base_out, "AMRFinder_Heatmap_ByType.pdf"),
        "omnibus": os.path.join(base_out, "AMRFinder_Omnibus_Stats.tsv"),
        "pairwise": os.path.join(base_out, "AMRFinder_Pairwise_Fisher.tsv"),
        "summary": os.path.join(base_out, "AMRFinder_ByClade_Summary.tsv"),
    }



def manuscript_style(ws, table):
    """Readable publication layout. Display only; no statistical calculations."""
    from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
    from openpyxl.utils import get_column_letter

    ws.sheet_view.showGridLines = False
    ws.freeze_panes = "C2"
    ws.auto_filter.ref = f"A1:{get_column_letter(table.shape[1])}{len(table)+1}"
    ws.row_dimensions[1].height = 36
    header_fill = PatternFill("solid", fgColor="203451")
    for cell in ws[1]:
        cell.fill = header_fill
        cell.font = Font(name="Aptos", size=11, bold=True, color="FFFFFF")
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    widths = {"Type":17, "Gene":21, "Total (n)":14,
              "Stats p-value":18, "Stats FDR":18}
    for j, column in enumerate(table.columns, 1):
        width = widths.get(column, 23 if "(n=" in column else 19)
        ws.column_dimensions[get_column_letter(j)].width = width

    prior_type = None
    for i, category in enumerate(table["Type"], 2):
        base_fill = ("F1F5FA" if i%2==0 else "FFFFFF")
        if category != prior_type:
            base_fill = {"AMR":"E6EFFB", "STRESS":"E8F3EA",
                         "VIRULENCE":"F6EDF7"}.get(category, "E6EFFB")
        for cell in ws[i]:
            cell.fill = PatternFill("solid", fgColor=base_fill)
            cell.font = Font(name="Aptos", size=10,
                             bold=(cell.column<=2),
                             color="203451" if cell.column<=2 else "1F2937")
            cell.alignment = Alignment(
                vertical="center",
                horizontal="left" if cell.column<=2 else "center"
            )
            if category != prior_type:
                cell.border = Border(top=Side(style="medium",color="53647B"))
        ws.row_dimensions[i].height = 23
        prior_type = category
    # Data range remains exactly the manuscript rows (no trailing note
    # rows that could be misinterpreted as extra genes by Excel/pandas).
    from openpyxl.comments import Comment
    ws.cell(1, table.shape[1]).comment = Comment(
        "Superscript letters reproduce the original SEN Holm-adjusted "
        "Fisher groupings only for tested variable clades. Blank letters "
        "can mean FDR nonsignificance or original 0%/100% clade exclusion. "
        "BH FDR is limited to available published genes, not the "
        "unavailable complete AMRFinder master.", "SEN pipeline"
    )


def run():
    p = argparse.ArgumentParser()
    p.add_argument("--counts",type=Path,default=ROOT/"data/table_s5/Table_S5_published_pure_counts.tsv")
    p.add_argument("--out",type=Path,default=ROOT/"results/table_s5_count_only_main_replica")
    args = p.parse_args()
    df = pd.read_csv(args.counts,sep="\t")
    if list(df.columns)!=["Type","Gene",*N] or len(df)!=42 or not df.Gene.is_unique:
        raise ValueError("Expected 42 unique Table S5 rows with original five-clade counts")
    for c,n in N.items():
        if not pd.api.types.is_integer_dtype(df[c]) or not df[c].between(0,n).all():
            raise ValueError("Invalid count: "+c)
    if df[["Type","Gene"]].isna().any().any():
        raise ValueError("Missing element identity")
    lookup = df.set_index("Gene")
    args.out.mkdir(parents=True, exist_ok=True)
    df.to_excel(args.out / "Table_S5_published_pure_counts.xlsx", sheet_name="Pure_Counts", index=False)
    excel_results = []
    print("[REFERENCE] SEN_genomics_pipeline main 15_amr_byclade.py blob a937dee96017627756145183f0e26abf203fedfb")
    print("[LIMIT] BH FDR is computed over the 42 displayed genes only")
    print("[LIMIT] Only group-level data and output tables are reproducible")
    for scheme in ["Combined_Clade1","Split_Clade1A1B"]:
        for mode in ["With_Unassigned","Without_Unassigned"]:
            outdir = args.out/scheme/mode
            outdir.mkdir(parents=True,exist_ok=True)
            out = build_output_paths(str(outdir))
            clade_order, clade_palette, include_unassigned = get_scheme_config(scheme,mode)
            clade_sizes = {c: N["Clade 1A"]+N["Clade 1B"] if c=="Clade 1" else N[c] for c in clade_order}
            def get_n(row,clade):
                return int(row["Clade 1A"])+int(row["Clade 1B"]) if clade=="Clade 1" else int(row[clade])
            symbols = sorted(g for g,r in lookup.iterrows() if sum(get_n(r,c) for c in clade_order)>0)
            # Original core uses pa.columns only; we do not fabricate sample rows.
            pa = pd.DataFrame(columns=symbols)
            counts_by_clade = {c:pd.Series({g:get_n(lookup.loc[g],c) for g in symbols},dtype=int) for c in clade_order}
            total_present = pd.Series({g:sum(int(counts_by_clade[c].loc[g]) for c in clade_order) for g in symbols},dtype=int)
            sym_to_type = lookup["Type"].to_dict()
            sym_to_name = {g:"NA" for g in symbols}
            def clade_header(c):
                return f"{c} (n={clade_sizes.get(c,0)})"
            print(f"[MODE] {scheme}/{mode} : {sum(clade_sizes.values())} genomes, {len(symbols)} present elements")
            omnibus_rows = []
            used_clades_by_element = {}

            for symbol in pa.columns:
                pres = {c: int(counts_by_clade[c].loc[symbol]) for c in clade_order}
                ns   = {c: int(clade_sizes.get(c, 0)) for c in clade_order}
                absn = {c: ns[c] - pres[c] for c in clade_order}
                tot_pres = int(total_present.loc[symbol])

                variable_clades = [c for c in clade_order if ns[c] > 0 and (0 < pres[c] < ns[c])]
                used_clades_by_element[symbol] = variable_clades

                pval = np.nan
                method = ""
                chi2 = np.nan
                dof = np.nan
                reason = ""

                if tot_pres < MIN_TOTAL_PRESENT:
                    reason = f"rare (<{MIN_TOTAL_PRESENT} total present)"
                elif len(variable_clades) < 2:
                    reason = "insufficient variable clades after zero/full exclusion"
                else:
                    table = np.array([[pres[c], absn[c]] for c in variable_clades], dtype=int)
                    try:
                        pval, method, chi2, dof = chi2_with_optional_montecarlo(table)
                    except Exception as e:
                        reason = f"chi2 failed: {e}"

                omnibus_rows.append({
                    "Element symbol": symbol,
                    "Total_Present": tot_pres,
                    "Tested_Clades": ",".join(variable_clades),
                    "Test_Method": method,
                    "Chi2": chi2,
                    "DoF": dof,
                    "Stats p-value": pval,
                    "Stats FDR": np.nan,
                    "Reason_Not_Tested": reason,
                })

            omnibus_df = pd.DataFrame(omnibus_rows)

            p_num = pd.to_numeric(omnibus_df["Stats p-value"], errors="coerce")
            tested_mask = p_num.notna()
            if tested_mask.sum() > 0:
                _, fdr_vals, _, _ = multipletests(p_num.loc[tested_mask].values.astype(float), method="fdr_bh")
                omnibus_df.loc[tested_mask, "Stats FDR"] = fdr_vals

            omnibus_df.to_csv(out["omnibus"], sep="\t", index=False)
            print(f"[OK] Omnibus stats: {out['omnibus']}")

            print("📊 Running post-hoc pairwise Fisher (Holm within element) for significant elements...")

            fdr_num = pd.to_numeric(omnibus_df["Stats FDR"], errors="coerce")
            sig_elements = omnibus_df.loc[tested_mask & (fdr_num <= ALPHA), "Element symbol"].tolist()

            pairwise_rows = []
            cld_letters_by_element = {e: {c: "" for c in clade_order} for e in pa.columns}

            for symbol in sig_elements:
                pres = {c: int(counts_by_clade[c].loc[symbol]) for c in clade_order}
                ns   = {c: int(clade_sizes.get(c, 0)) for c in clade_order}
                absn = {c: ns[c] - pres[c] for c in clade_order}

                variable_clades = used_clades_by_element.get(symbol, [])
                if len(variable_clades) < 2:
                    continue

                pairs = list(itertools.combinations(variable_clades, 2))
                pvals = []
                tmp_rows = []

                for a, b in pairs:
                    a1, a0 = pres[a], absn[a]
                    b1, b0 = pres[b], absn[b]

                    _, p_raw = fisher_exact([[a1, a0], [b1, b0]], alternative="two-sided")
                    or05 = odds_ratio_with_haldene(a1, a0, b1, b0, add=0.5)

                    tmp_rows.append({
                        "Element symbol": symbol,
                        "Group1": a,
                        "Group2": b,
                        "G1_present": a1,
                        "G1_absent": a0,
                        "G2_present": b1,
                        "G2_absent": b0,
                        "OR_0.5": or05,
                        "p_raw": float(p_raw),
                        "p_holm": np.nan,
                    })
                    pvals.append(float(p_raw))

                if pvals:
                    _, p_holm, _, _ = multipletests(pvals, method="holm")
                    for i in range(len(tmp_rows)):
                        tmp_rows[i]["p_holm"] = float(p_holm[i])

                pairwise_rows.extend(tmp_rows)

                p_adj_map = {}
                for r in tmp_rows:
                    key = tuple(sorted((r["Group1"], r["Group2"])))
                    p_adj_map[key] = r["p_holm"]

                cld_var = make_cld_letters(variable_clades, p_adj_map, alpha=ALPHA)
                for c in variable_clades:
                    cld_letters_by_element[symbol][c] = cld_var.get(c, "")

            pairwise_df = pd.DataFrame(pairwise_rows)
            pairwise_df.to_csv(out["pairwise"], sep="\t", index=False)
            print(f"[OK] Pairwise stats: {out['pairwise']}")

            summary_rows = []
            p_map = dict(zip(omnibus_df["Element symbol"], omnibus_df["Stats p-value"]))
            f_map = dict(zip(omnibus_df["Element symbol"], omnibus_df["Stats FDR"]))

            for symbol in pa.columns:
                t = sym_to_type.get(symbol, "NA")
                nm = sym_to_name.get(symbol, "NA")

                letters = cld_letters_by_element.get(symbol, {c: "" for c in clade_order})

                row = {
                    "Type": t,
                    "Element symbol": symbol,
                    "Element name": nm,
                    "Total (n)": int(total_present.loc[symbol]),
                }
                for c in clade_order:
                    n = int(clade_sizes.get(c, 0))
                    cnt = int(counts_by_clade[c].loc[symbol]) if n > 0 else 0
                    row[clade_header(c)] = format_cell(cnt, n, letters.get(c, ""))

                pval = p_map.get(symbol, np.nan)
                fdr  = f_map.get(symbol, np.nan)
                row["Stats p-value"] = "" if pd.isna(pval) else float(pval)
                row["Stats FDR"]     = "" if pd.isna(fdr) else float(fdr)

                summary_rows.append(row)

            summary_df = pd.DataFrame(summary_rows)

            _tested = pd.to_numeric(summary_df["Stats p-value"], errors="coerce").notna()
            _total  = pd.to_numeric(summary_df["Total (n)"], errors="coerce").fillna(0).astype(int)

            summary_df["_tested_sort"] = _tested
            summary_df["_total_sort"] = _total

            summary_df = (
                summary_df.sort_values(["_tested_sort", "_total_sort", "Type", "Element symbol"],
                                       ascending=[False, False, True, True])
                         .drop(columns=["_tested_sort", "_total_sort"])
                         .reset_index(drop=True)
            )

            summary_df.to_csv(out["summary"], sep="\t", index=False)
            summary_df.to_csv(out["table"], sep="\t", index=False)
            print(f"[OK] Summary (paper-style): {out['summary']}")
            print(f"[OK] ElementCounts TSV (same content): {out['table']}")


            # The original workbook is preserved. Excel mirrors of every
            # individual TSV are also written for direct inspection on AjuloPC.
            # Reorder the original pipeline's formatted table to the exact
            # published manuscript gene order and hide the unavailable
            # Element name field. Superscript letters remain unchanged.
            manuscript_columns = (
                ["Type", "Gene", "Total (n)"]
                + [clade_header(c) for c in clade_order]
                + ["Stats p-value", "Stats FDR"]
            )
            paper_by_gene = summary_df.set_index("Element symbol",drop=False)
            published_order = [
                symbol for symbol in df["Gene"].tolist() if symbol in paper_by_gene.index
            ]
            manuscript_rows = paper_by_gene.loc[published_order]
            manuscript_df = manuscript_rows.rename(
                columns={"Element symbol":"Gene"}
            )[manuscript_columns].reset_index(drop=True)
            ms_tsv = outdir / "Table_S5_Manuscript_Ready.tsv"
            ms_xlsx = outdir / "Table_S5_Manuscript_Ready.xlsx"
            manuscript_df.to_csv(ms_tsv,sep="\t",index=False)
            with pd.ExcelWriter(ms_xlsx,engine="openpyxl") as paper_writer:
                manuscript_df.to_excel(
                    paper_writer,sheet_name="Manuscript_Table_S5",index=False
                )
                manuscript_style(
                    paper_writer.sheets["Manuscript_Table_S5"],manuscript_df
                )
            with pd.ExcelWriter(out["excel"],engine="openpyxl") as xw:
                manuscript_df.to_excel(xw,sheet_name="Manuscript_Table_S5",index=False)
                manuscript_style(xw.sheets["Manuscript_Table_S5"],manuscript_df)
                summary_df.to_excel(xw,sheet_name="Summary_ForPaper",index=False)
                omnibus_df.to_excel(xw,sheet_name="Omnibus_Stats",index=False)
                pairwise_df.to_excel(xw,sheet_name="Pairwise_Fisher_Holm",index=False)
            # Summary.xlsx above already mirrors Summary.tsv (and includes
            # omnibus/pairwise sheets). Do not overwrite that workbook.
            omnibus_df.to_excel(Path(out["omnibus"]).with_suffix(".xlsx"),
                                sheet_name="Omnibus_Stats",index=False)
            pairwise_df.to_excel(Path(out["pairwise"]).with_suffix(".xlsx"),
                                 sheet_name="Pairwise_Fisher_Holm",index=False)
            summary_df.to_excel(Path(out["table"]).with_suffix(".xlsx"),
                                sheet_name="ElementCounts",index=False)
            excel_results.append((scheme, mode, summary_df.copy(),
                                  omnibus_df.copy(), pairwise_df.copy(),
                                  manuscript_df.copy()))
            (outdir/"COUNT_ONLY_LIMITATIONS.txt").write_text(
                "Source: SEN_genomics_pipeline main 15_amr_byclade.py (a937dee96017627756145183f0e26abf203fedfb)\n"
                "Original tests, parameters, exclusions, Holm/BH and sort copied unchanged.\n"
                "BH scope: 42 published genes only (not the complete master catalog).\n"
                "Element name: NA, unavailable in pure counts.\n"
                "Cannot generate original Sample_Clade_Metadata, PresenceAbsence_Matrix, "
                "AMRFinder_PresenceAbsence_Matrix.tsv, sample QC or heatmap PNG/PDF.\n"
            )
            print(f"[OK] Excel and count-derived results: {outdir}")
    # All four analysis configurations in one workbook, produced locally
    # whenever this script is run. This does not recompute or modify any test.
    combined_path = args.out / "Table_S5_All_Configurations.xlsx"
    with pd.ExcelWriter(combined_path, engine="openpyxl") as writer:
        df.to_excel(writer,sheet_name="Original_Pure_Counts",index=False)
        pd.DataFrame({
            "Important notes": [
                "Statistics copied from SEN genomics main 15_amr_byclade.py",
                "BH-FDR family is limited to published Table S5 determinants",
                "Source input is aggregate clade counts; no sample-level matrix is reconstructed",
                "Combined and split schemes are analysed separately, with and without Unassigned",
            ]
        }).to_excel(writer,sheet_name="ReadMe",index=False)
        for scheme, mode, summary, omnibus, pairwise, manuscript in excel_results:
            prefix = ("C1" if scheme=="Combined_Clade1" else "Split") + (
                "_With" if mode=="With_Unassigned" else "_Without"
            )
            manuscript_sheet = prefix+"_Manuscript"
            manuscript.to_excel(writer,sheet_name=manuscript_sheet,index=False)
            manuscript_style(writer.sheets[manuscript_sheet],manuscript)
            summary.to_excel(writer,sheet_name=prefix+"_Summary",index=False)
            omnibus.to_excel(writer,sheet_name=prefix+"_Omnibus",index=False)
            pairwise.to_excel(writer,sheet_name=prefix+"_Pairwise",index=False)
    print(f"[OK] Combined Excel workbook: {combined_path}")
    print("[OK] Manuscript-ready Table S5 as separate TSV/XLSX plus first worksheet in each summary workbook")
    print("[PASS] Four original main-branch configurations complete")
if __name__=="__main__":
    run()
