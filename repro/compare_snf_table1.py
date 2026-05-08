"""Compare repro SNF metrics (repro_out/<CANCER>.snf.metrics.json) to paper Table 1 SNF row."""

from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
_REF = REPO_ROOT / "data" / "paper_table1_reference.json"


def _load_paper_snf() -> dict[str, tuple[float, int]]:
    if not _REF.is_file():
        raise FileNotFoundError(f"Missing {_REF} (expected data/paper_table1_reference.json)")
    data = json.loads(_REF.read_text(encoding="utf-8"))
    snf = data.get("snf", {})
    out: dict[str, tuple[float, int]] = {}
    for k, v in snf.items():
        if isinstance(v, (list, tuple)) and len(v) >= 2:
            out[str(k).upper()] = (float(v[0]), int(v[1]))
    return out


def main() -> None:
    paper_snf = _load_paper_snf()
    out: list[dict] = []
    for cancer, (paper_lp, paper_n) in paper_snf.items():
        p = REPO_ROOT / "repro_out" / f"{cancer}.snf.metrics.json"
        if not p.is_file():
            out.append({"cancer": cancer, "status": "missing", "path": str(p)})
            continue
        m = json.loads(p.read_text(encoding="utf-8"))
        lp = m.get("neg_log10_logrank_p")
        n = m.get("n_significant_clinical")
        out.append(
            {
                "cancer": cancer,
                "paper_neg_log10_logrank_p": paper_lp,
                "paper_n_significant_clinical": paper_n,
                "repro_neg_log10_logrank_p": lp,
                "repro_n_significant_clinical": n,
                "n_samples_eval": m.get("n_samples_eval"),
            }
        )
    dest = REPO_ROOT / "repro_out" / "snf_table1_compare.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(dest)


if __name__ == "__main__":
    main()
