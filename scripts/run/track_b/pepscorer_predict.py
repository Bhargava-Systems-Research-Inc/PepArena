#!/usr/bin/env python
import sys
import runpy
from pathlib import Path

from mordred import Calculator

PS = Path("/12TBDrive1/mega_pep_bench/methods/external/PepScorer/PepScorerRMSD/PepScorerRMSD")

_orig_pandas = Calculator.pandas

def _pandas_fill_missing(self, mols, **kw):
    return _orig_pandas(self, mols, **kw).fill_missing()

Calculator.pandas = _pandas_fill_missing

if __name__ == "__main__":
    sys.argv = ["predict.py"] + sys.argv[1:]
    runpy.run_path(str(PS / "predict.py"), run_name="__main__")
