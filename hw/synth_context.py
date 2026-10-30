#!/usr/bin/env python3
"""Sensitivity of the synthesised areas to the synthesis context.

The committed reports rpt_<unit>.txt are synthesised from route_units.v, which holds every unit. Yosys and ABC are
deterministic, but their result can depend on the rest of the source file (module order, generated identifiers).
For every unit this script synthesises, with the commands of synth.sh,
  (a) "alone":    a file that contains only the unit and the modules it instantiates, in their original order;
  (b) "reversed": the whole of route_units.v with its modules in reverse order;
and compares the areas with the committed reports. It writes synth_context.json, read by the article's generator
(paper3/sec_results.py, D8), and prints a table.
Usage (from hw/, with yowasp-yosys and NangateOpenCellLibrary_typical.lib):
    python3 synth_context.py > synth_context.log
"""
import json
import re
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
LIB = HERE / "NangateOpenCellLibrary_typical.lib"
UNITS = ["rc_xyz", "rc_zxy", "rc_nf", "rc_dyad", "rc_ham", "rc_o1t", "rc_twophase", "rc_fa", "rca_select", "rca_agg",
         "qr_select", "qr_update", "rc_ta", "rc_ef", "rc_fl", "adele_sel", "vc_fifo", "xbar7"]
AREA = re.compile(r"Chip area for module '\\?(\w+)': ([0-9.]+)")


def modules(src: str) -> list:
    """(name, text) of every module of the file, in order."""
    return [(m.group(1), m.group(0)) for m in re.finditer(r"(?ms)^module\s+(\w+)\b.*?^endmodule[^\n]*\n", src)]


def deps(name: str, mods: dict) -> set:
    """The modules instantiated by `name`, recursively."""
    out, todo = set(), [name]
    while todo:
        body = mods[todo.pop()]
        for other in mods:
            if other not in out and other != name and re.search(rf"(?m)^\s*{other}\b", body):
                out.add(other)
                todo.append(other)
    return out


def area(verilog: str, top: str) -> float:
    # YoWASP sees only the working directory, so the temporary files live under hw/ and paths are relative.
    with tempfile.TemporaryDirectory(dir=HERE, prefix=".ctx_") as tmp:
        d = Path(tmp).name
        (Path(tmp) / "units.v").write_text(verilog)
        lib = LIB.name
        script = (f"read_verilog {d}/units.v; synth -top {top} -flatten; dfflibmap -liberty {lib}; "
                  f"abc -liberty {lib}; opt_clean; tee -o {d}/rpt.txt stat -liberty {lib}")
        subprocess.run(["yowasp-yosys", "-q", "-p", script], check=True, capture_output=True, cwd=HERE)
        return float(AREA.search((Path(tmp) / "rpt.txt").read_text()).group(2))


def main() -> None:
    src = (HERE / "route_units.v").read_text()
    order = modules(src)
    mods = dict(order)
    reversed_src = "".join(text for _, text in reversed(order))
    res = {}
    print(f"{'unit':12s} {'committed':>10s} {'alone':>10s} {'reversed':>10s} {'max dev.':>9s}")
    for u in UNITS:
        committed = float(AREA.search((HERE / f"rpt_{u}.txt").read_text()).group(2))
        keep = deps(u, mods) | {u}
        alone = area("".join(text for n, text in order if n in keep), u)
        rev = area(reversed_src, u)
        res[u] = {"committed": committed, "alone": alone, "reversed": rev}
        dev = max(abs(alone / committed - 1), abs(rev / committed - 1))
        print(f"{u:12s} {committed:10.3f} {alone:10.3f} {rev:10.3f} {100 * dev:8.1f}%", flush=True)
    (HERE / "synth_context.json").write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
