#!/usr/bin/env python3
"""birth_genealogy.py — render the cluster birth tree from a trace file.

Each non-Initial cluster has a `birth_from` field pointing at the
cluster it was budded from; chasing those edges across all clusters
gives a forest of trees rooted at the Initial cluster(s).

This script renders that forest as a Graphviz DOT file. Pipe through
`dot -Tpdf -o tree.pdf` to render. Cluster nodes are labeled with
their ID and abundance; edges are colored by birth p-value (darker =
more significant).

Usage:
  python3 birth_genealogy.py < trace.json > tree.dot
  python3 birth_genealogy.py trace.json | dot -Tpdf -o tree.pdf

Pure stdlib; no networkx dependency.
"""

from __future__ import annotations
import json
import math
import sys
from pathlib import Path


def load(src: str | Path | None) -> dict:
    text = Path(src).read_text() if src else sys.stdin.read()
    obj = json.loads(text)
    if obj.get("dada2_rs_command") != "cluster-trace":
        sys.exit(f"not a cluster-trace file (tag={obj.get('dada2_rs_command')!r})")
    return obj


def edge_color(pval: float | None) -> str:
    """Map birth pval to a Graphviz greyscale shade. Smaller pval → darker."""
    if pval is None or pval <= 0:
        return "black"
    # log10(pval) ranges roughly [-300, 0]; clip to [-40, 0] for color.
    lp = max(math.log10(pval), -40.0)
    # 0 (white-ish) at lp=0, 0.9 (darkest) at lp=-40
    shade = min(0.9, max(0.0, -lp / 40.0))
    grey = int((1.0 - shade) * 255)
    return f'"#{grey:02x}{grey:02x}{grey:02x}"'


def render(trace: dict, out=sys.stdout) -> None:
    clusters = trace["clusters"]
    sample = trace.get("sample", "?")
    iter_str = (
        f"iter {trace['iteration']}" if trace.get("iteration") is not None else "final"
    )

    print('digraph cluster_genealogy {', file=out)
    print('  rankdir=TB;', file=out)
    print(f'  label="{sample} ({iter_str}, {len(clusters)} clusters)";', file=out)
    print('  labelloc=t;', file=out)
    print('  node [shape=ellipse, fontsize=10];', file=out)
    print('  edge [fontsize=8];', file=out)

    for c in clusters:
        cid = c["id"]
        n = c["abundance"]
        bt = c["birth_type"]
        shape = "ellipse" if bt == "Initial" else "box"
        fill = {
            "Initial": "lightblue",
            "Abundance": "lightyellow",
            "Prior": "lightgreen",
            "Singleton": "lightpink",
        }.get(bt, "white")
        label = f"#{cid}\\n{n} reads\\n{bt}"
        print(f'  c{cid} [label="{label}", style=filled, fillcolor="{fill}", shape={shape}];',
              file=out)

    for c in clusters:
        if c["birth_type"] == "Initial":
            continue
        parent = c["birth_from"]
        col = edge_color(c.get("birth_pval"))
        h = c.get("birth_hamming", 0)
        elabel = f'h={h}'
        if c.get("birth_pval") is not None:
            elabel += f"\\np={c['birth_pval']:.1e}"
        print(f'  c{parent} -> c{c["id"]} [color={col}, label="{elabel}"];', file=out)

    print('}', file=out)


def main(argv: list[str]) -> int:
    src = argv[1] if len(argv) > 1 else None
    trace = load(src)
    render(trace)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
