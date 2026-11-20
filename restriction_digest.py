import math
import os
import sys
import urllib.request

import matplotlib.pyplot as plt

ENZYMES = {
    "EcoRI":   ("GAATTC", 1),
    "BamHI":   ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
    "TaqI":    ("TCGA",   1),
    "HaeIII":  ("GGCC",   2),
}

NCBI_ACCESSION = "NC_000913.3"
SEQ_START, SEQ_STOP = 1, 10000
CACHE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequence.fasta")


def read_fasta(path):
    with open(path, encoding="utf-8") as f:
        lines = f.read().splitlines()
    header = lines[0][1:].strip() if lines and lines[0].startswith(">") else "sequence"
    seq = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
    return header, seq


def fetch_from_ncbi(accession, start, stop):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={accession}&rettype=fasta&retmode=text"
        f"&seq_start={start}&seq_stop={stop}"
    )
    req = urllib.request.Request(url, headers={"User-Agent": "bioinfo-lab7"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        text = resp.read().decode("utf-8")
    if not text.startswith(">"):
        raise ValueError("Unexpected NCBI response")
    with open(CACHE_FILE, "w", encoding="utf-8") as f:
        f.write(text)
    return read_fasta(CACHE_FILE)


def load_sequence():
    if len(sys.argv) > 1:
        return read_fasta(sys.argv[1])
    if os.path.exists(CACHE_FILE):
        return read_fasta(CACHE_FILE)
    print(f"Downloading {NCBI_ACCESSION} [{SEQ_START}-{SEQ_STOP}] from NCBI...")
    return fetch_from_ncbi(NCBI_ACCESSION, SEQ_START, SEQ_STOP)


def find_cut_positions(seq, site, offset):
    cuts = []
    i = seq.find(site)
    while i != -1:
        cuts.append(i + offset)
        i = seq.find(site, i + 1)
    return cuts


def fragments_from_cuts(seq_len, cuts):
    bounds = [0] + sorted(set(cuts)) + [seq_len]
    return [bounds[k + 1] - bounds[k] for k in range(len(bounds) - 1)]


def digest(seq, enzymes):
    results = {}
    for name, (site, offset) in enzymes.items():
        cuts = find_cut_positions(seq, site, offset)
        results[name] = {
            "site": site,
            "offset": offset,
            "cuts": cuts,
            "fragments": fragments_from_cuts(len(seq), cuts),
        }
    all_cuts = sorted(set(c for r in results.values() for c in r["cuts"]))
    results["All enzymes"] = {
        "site": "-",
        "offset": None,
        "cuts": all_cuts,
        "fragments": fragments_from_cuts(len(seq), all_cuts),
    }
    return results


def format_site(site, offset):
    return f"5'-{site[:offset]}^{site[offset:]}-3'"


def print_report(header, seq, results):
    print("=" * 70)
    print(f"Sequence: {header}")
    print(f"Length:   {len(seq)} bp")
    print("=" * 70)
    for name, r in results.items():
        print(f"\n{name}", end="")
        if r["offset"] is not None:
            print(f"  recognition site {r['site']}, cleavage {format_site(r['site'], r['offset'])}")
        else:
            print("  (combined digest)")
        print(f"  Number of cleavages: {len(r['cuts'])}")
        if r["cuts"]:
            print(f"  Cleavage positions (cut after nt): {r['cuts']}")
        frags = r["fragments"]
        print(f"  Number of fragments: {len(frags)}")
        print(f"  Fragment lengths (5'->3' order): {frags}")
        print(f"  Fragment lengths (sorted, largest first): {sorted(frags, reverse=True)}")


LADDER = [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500,
          1000, 750, 500, 250, 100]
MIN_VISIBLE_BP = 50


def migration(length, min_bp=50, max_bp=12000):
    length = min(max(length, min_bp), max_bp)
    return (math.log10(max_bp) - math.log10(length)) / (math.log10(max_bp) - math.log10(min_bp))


def simulate_gel(results, seq_len, out_file="gel.png"):
    lanes = [("Ladder", LADDER)] + [(n, r["fragments"]) for n, r in results.items()]
    n = len(lanes)

    fig, ax = plt.subplots(figsize=(1.4 * n + 1.5, 8))
    ax.set_facecolor("#111111")
    fig.patch.set_facecolor("white")
    lane_w = 0.6

    for i, (name, frags) in enumerate(lanes):
        x = i + 1
        ax.add_patch(plt.Rectangle((x - lane_w / 2, -0.04), lane_w, 0.02, color="#555555"))
        counts = {}
        for f in frags:
            counts[f] = counts.get(f, 0) + 1
        for size, cnt in counts.items():
            if size < MIN_VISIBLE_BP:
                continue
            y = migration(size)
            mass = size * cnt
            alpha = 0.35 + 0.65 * min(1.0, mass / (0.25 * max(seq_len, max(LADDER))))
            color = "#ff9966" if name == "Ladder" else "#ffffff"
            ax.add_patch(plt.Rectangle((x - lane_w / 2, y - 0.004), lane_w, 0.008,
                                       color=color, alpha=alpha, linewidth=0))
            if name == "Ladder":
                ax.text(x - lane_w / 2 - 0.08, y, f"{size}", color="#dddddd",
                        ha="right", va="center", fontsize=7, clip_on=False)

    ax.set_xlim(0.0, n + 0.8)
    ax.set_ylim(1.03, -0.06)
    ax.set_xticks(range(1, n + 1))
    ax.set_xticklabels([name for name, _ in lanes], rotation=30, ha="right")
    ax.set_yticks([])
    ax.set_title(f"Simulated agarose gel electrophoresis ({seq_len} bp, {len(results) - 1} enzymes)")
    ax.text(n + 0.7, -0.03, "(-)", color="white", ha="right", va="center")
    ax.text(n + 0.7, 1.01, "(+)", color="white", ha="right", va="center")
    plt.tight_layout()
    plt.savefig(out_file, dpi=150)
    print(f"\nGel image saved to {out_file}")
    plt.show()


def main():
    header, seq = load_sequence()
    seq = "".join(c for c in seq if c in "ACGT")
    results = digest(seq, ENZYMES)
    print_report(header, seq, results)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gel.png")
    simulate_gel(results, len(seq), out)


if __name__ == "__main__":
    main()
