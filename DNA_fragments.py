"""
DNA Fragment Grouper – GUI wrapper
"""

from __future__ import annotations
import os
import threading
import warnings
from pathlib import Path
from typing import Callable, List, Set

import numpy as np
import pandas as pd
import random
import distance
from sklearn.cluster import AffinityPropagation
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# ----------------------------------------------------------------------
BsmBI      = "CGTCTC";  BsmBI_rev  = "GAGACG"
BsaI       = "GGTCTC";  BsaI_rev   = "GAGACC"
BbsI       = "GAAGACT"; BbsI_rev   = "AGTCTTC"
BspMI      = "ACCTGCTTA"; BspMI_rev = "TAAGCAGGT"

NUCL       = ["A", "T", "C", "G"]
MIN_LENGTH = 301
MAX_LENGTH = 499

# ----------------------------------------------------------------------

def nth_repl(s: str, old: str, new: str, n: int) -> str:
    """Replace the n‑th occurrence of *old* in *s* by *new*."""
    find = s.find(old)
    i = int(find != -1)
    while find != -1 and i != n:
        find = s.find(old, find + 1)
        i += 1
    if i == n and i <= len(s.split(old)) - 1:
        return s[:find] + new + s[find + len(old) :]
    return s

# ----------------------------------------------------------------------
# Pipeline functions

def compute_distance_matrix(
    fragments: np.ndarray,
    progress: Callable[[float], None] | None = None,
) -> np.ndarray:
    """Compute Levenshtein similarity matrix, emitting ≤100 progress updates."""
    n = len(fragments)
    matrix = np.empty((n, n), dtype=float)
    step = max(1, n // 100)
    for i, f2 in enumerate(fragments):
        if i % step == 0 and progress:
            progress((i / n) * 80 + 2)
        matrix[i] = [-distance.levenshtein(f1, f2) for f1 in fragments]
    if progress:
        progress(90)
    return matrix


def group_fragments(
    df: pd.DataFrame,
    *,
    max_length: int = MAX_LENGTH,
    progress: Callable[[float], None] | None = None,
) -> pd.DataFrame:
    fragments = df["Sequence"].tolist()
    used_fragments: Set[str] = set()
    group_labels = np.empty(len(fragments))
    group_counter = 0
    iter_counter = 0

    while len(used_fragments) < len(fragments):
        group: List[str] = []
        total_length = 0
        clusters: Set[int] = set()

        for idx, row in df.iterrows():
            frag = row["Sequence"]
            cluster = row["Cluster"]
            if frag in used_fragments:
                continue
            frag_length = len(frag)
            if (
                len(group) < 3 and total_length + frag_length <= max_length and cluster not in clusters
            ):
                group.append(frag)
                clusters.add(cluster)
                total_length += frag_length
                used_fragments.add(frag)
                group_labels[idx] = group_counter
            elif len(group) == 0 and frag_length >= 500:
                group_labels[idx] = group_counter
                group_counter += 1
                used_fragments.add(frag)
        group_counter += 1
        iter_counter += 1
        if progress:
            progress(90 + (len(used_fragments) / len(fragments)) * 5)  # 90–95 %
        if iter_counter > len(fragments) * 3:
            break
    out = df.copy()
    out["Group"] = group_labels.astype(int)
    return out


def apply_aggressive_grouping(df1: pd.DataFrame) -> pd.DataFrame:
    df = df1.copy()
    singleton_groups = [g for g, sub in df.groupby("Group") if len(sub) == 1]
    for i, gid in enumerate(singleton_groups):
        if i % 3 == 0:
            new = gid
        else:
            df.loc[df["Group"] == gid, "Group"] = new
    return df


def replace_cut_sites_and_pad(grouped_df: pd.DataFrame) -> pd.DataFrame:
    df = grouped_df.copy()
    for idx, row in df.iterrows():
        group, sequence, name, length = row
        new_seq = nth_repl(sequence, BsmBI, BbsI, 2)
        new_seq = nth_repl(new_seq, BsmBI_rev, BbsI_rev, 2)
        new_seq = nth_repl(new_seq, BsmBI, BspMI, 2)
        new_seq = nth_repl(new_seq, BsmBI_rev, BspMI_rev, 2)
        if len(new_seq) < 300:
            padding = "".join(random.choices(NUCL, k=MIN_LENGTH - len(new_seq)))
            for cut in (
                BsmBI, BsaI, BbsI, BspMI,
                BsmBI_rev, BsaI_rev, BbsI_rev, BspMI_rev,
            ):
                padding = padding.replace(cut, "ATCCGATGGTC")
            new_seq += padding
        df.at[idx, "Sequence"] = new_seq
        df.at[idx, "Length"] = len(new_seq)
    return df


def pipeline(
    csv_path: os.PathLike,
    *,
    aggressive: bool,
    progress: Callable[[float], None] | None = None,
) -> Path:
    if progress:
        progress(0)
    df = pd.read_csv(csv_path, sep=";")
    if progress:
        progress(2)

    df = df.apply(lambda x: x.astype(str).str.upper())
    df_large = df[df["Sequence"].apply(len) >= 500]
    # filter short fragments and reset index so numpy arrays align
    df = df[df["Sequence"].apply(len) < MAX_LENGTH].reset_index(drop=True)

    # distance matrix + clustering
    fragments = df["Sequence"].to_numpy()
    sim = compute_distance_matrix(fragments, progress)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        affprop = AffinityPropagation(affinity="precomputed", random_state=0)
        affprop.fit(sim)

    clusters = np.empty(len(df), dtype=int)
    for cid in np.unique(affprop.labels_):
        clusters[affprop.labels_ == cid] = cid
    df["Cluster"] = clusters
    if progress:
        progress(90)

    df1 = group_fragments(df, progress=progress)
    if aggressive:
        df1 = apply_aggressive_grouping(df1)

    grouped = (
        df1.groupby("Group").agg({"Sequence": "".join, "Name": list}).reset_index()
    )
    grouped["Length"] = grouped["Sequence"].str.len()

    # add back large fragments
    max_group = grouped["Group"].max() if not grouped.empty else 0
    for _, row in df_large.iterrows():
        name, author, frag = row
        max_group += 1
        grouped = pd.concat(
            [
                grouped,
                pd.DataFrame(
                    {
                        "Group": [max_group],
                        "Sequence": [frag],
                        "Name": [[name]],
                        "Length": [len(frag)],
                    }
                ),
            ],
            ignore_index=True,
        )

    grouped = replace_cut_sites_and_pad(grouped)
    if progress:
        progress(98)

    out_path = Path(csv_path).with_name(Path(csv_path).stem + "_grouped.csv")
    grouped.to_csv(out_path, index=False, sep=";")
    if progress:
        progress(100)
    return out_path

# ----------------------------------------------------------------------
# GUI

class GrouperApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("DNA Fragment Grouper")
        self.geometry("500x210")
        # bring to foreground immediately (mac‑friendly)
        self.lift()
        self.attributes("-topmost", True)
        self.after(100, self.attributes, "-topmost", False)
        self.after_idle(self.focus_force)

        self.csv_path = tk.StringVar()
        self.aggressive = tk.BooleanVar()
        self.build_ui()

    # ---------------- UI -----------------
    def build_ui(self) -> None:
        pad = {"padx": 10, "pady": 5}
        ttk.Label(self, text="CSV file:").grid(row=0, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.csv_path, width=42).grid(row=0, column=1, sticky="we", **pad)
        ttk.Button(self, text="Browse…", command=self.browse).grid(row=0, column=2, **pad)

        ttk.Checkbutton(
            self,
            text="Aggressive grouping (combine singletons)",
            variable=self.aggressive,
        ).grid(row=1, column=1, columnspan=2, sticky="w", **pad)

        self.run_button = ttk.Button(self, text="Run", command=self.start)
        self.run_button.grid(row=2, column=2, sticky="e", **pad)

        self.prog_var = tk.DoubleVar()
        self.prog = ttk.Progressbar(self, maximum=100, variable=self.prog_var)
        self.prog.grid(row=3, column=0, columnspan=3, sticky="we", **pad)

        self.prog_label = ttk.Label(self, text="Idle (0 %)", anchor="w")
        self.prog_label.grid(row=4, column=0, columnspan=3, sticky="we", **pad)

        self.grid_columnconfigure(1, weight=1)

    # ------------- Callbacks -------------
    def browse(self) -> None:
        file = filedialog.askopenfilename(
            title="Select CSV file",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if file:
            self.csv_path.set(file)

    def start(self) -> None:
        if not self.csv_path.get():
            messagebox.showerror("Error", "Please choose a CSV file first.")
            return
        self.prog_var.set(0)
        self.prog_label.config(text="Queued… (0 %)")
        self.run_button.configure(state="disabled")
        threading.Thread(target=self._worker, daemon=True).start()

    def update_progress(self, percent: float) -> None:
        """Update bar and label, ignoring out‑of‑order callbacks."""
        def _set() -> None:
            current = self.prog_var.get()
            if percent < current:
                return  # skip stale update
            self.prog_var.set(percent)
            if percent < 2:
                phase = "Loading data…"
            elif percent < 90:
                phase = "Computing distance matrix & clustering…"
            elif percent < 95:
                phase = "Grouping fragments…"
            elif percent < 98:
                phase = "Final processing…"
            elif percent < 100:
                phase = "Saving output…"
            else:
                phase = "Done."
            self.prog_label.config(text=f"{phase} ({percent:.0f} %)")
        self.after(0, _set)

    # ---------------- Worker thread ----------------
    def _worker(self) -> None:
        try:
            out = pipeline(
                Path(self.csv_path.get()),
                aggressive=self.aggressive.get(),
                progress=self.update_progress,
            )
            def _show(path=out):
                self.prog_var.set(100)
                self.prog_label.config(text="Done. (100 %)")
                messagebox.showinfo("Finished", f"Grouped CSV written to:\n{path}")
            self.after_idle(_show)
        except Exception as exc:
            msg = str(exc)
            self.after(0, lambda m=msg: messagebox.showerror("Error", f"An error occurred:\n{m}"))
        finally:
            self.after(0, lambda: self.run_button.configure(state="normal"))


def main() -> None:
    GrouperApp().mainloop()


if __name__ == "__main__":
    main()
