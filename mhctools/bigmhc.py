# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Wrapper for BigMHC (https://github.com/KarchinLab/bigmhc).

BigMHC offers two models:

* **BigMHC_EL** — pMHC presentation (eluted-ligand trained)
* **BigMHC_IM** — immunogenicity (transfer-learned on IFN-γ / MANAFEST)

Both are pan-allelic deep neural networks that accept peptides of any
length and produce a probability score (0–1).
"""

import csv
import os
import subprocess
import tempfile

from .pred import Kind, Pred, PeptidePreds


def _find_bigmhc_dir(bigmhc_path=None):
    """Resolve the BigMHC installation directory.

    Checks, in order:
    1. The *bigmhc_path* argument
    2. The ``BIGMHC_DIR`` environment variable
    3. ``~/bigmhc``
    """
    if bigmhc_path:
        return bigmhc_path
    env = os.environ.get("BIGMHC_DIR")
    if env:
        return env
    home = os.path.join(os.path.expanduser("~"), "bigmhc")
    if os.path.isdir(home):
        return home
    raise FileNotFoundError(
        "BigMHC not found. Set BIGMHC_DIR or pass bigmhc_path= to the constructor. "
        "Clone from https://github.com/KarchinLab/bigmhc")


class BigMHC:
    """Wrapper for BigMHC presentation and immunogenicity predictions.

    Parameters
    ----------
    alleles : list of str
        HLA alleles (e.g. ``["HLA-A*02:01", "HLA-B*07:02"]``).
    mode : str
        ``"el"`` for presentation or ``"im"`` for immunogenicity.
    bigmhc_path : str, optional
        Path to the cloned BigMHC repository root.
    device : str
        PyTorch device string. Default ``"cpu"``.
    """

    VALID_MODES = ("el", "im")

    def __init__(self, alleles, mode="el", bigmhc_path=None, device="cpu"):
        if isinstance(alleles, str):
            alleles = [alleles]
        if mode not in self.VALID_MODES:
            raise ValueError(
                "mode must be one of %s, got %r" % (self.VALID_MODES, mode))
        self.alleles = list(alleles)
        self.mode = mode
        self.bigmhc_dir = _find_bigmhc_dir(bigmhc_path)
        self.device = device

        predict_script = os.path.join(self.bigmhc_dir, "src", "predict.py")
        if not os.path.isfile(predict_script):
            raise FileNotFoundError(
                "predict.py not found at %s" % predict_script)

    def __str__(self):
        return "BigMHC(mode=%s, alleles=%s)" % (self.mode, self.alleles)

    def __repr__(self):
        return str(self)

    def _pred_kind(self):
        if self.mode == "im":
            return Kind.immunogenicity
        return Kind.pMHC_presentation

    def _predictor_name(self):
        return "bigmhc_%s" % self.mode

    # ------------------------------------------------------------------
    # Core prediction
    # ------------------------------------------------------------------

    def _run_bigmhc(self, peptides, alleles):
        """Write CSV, run predict.py, return parsed rows."""
        predict_script = os.path.join(self.bigmhc_dir, "src", "predict.py")

        with tempfile.TemporaryDirectory() as tmpdir:
            input_csv = os.path.join(tmpdir, "input.csv")
            output_csv = input_csv + ".prd"

            with open(input_csv, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["mhc", "pep"])
                for pep, allele in zip(peptides, alleles):
                    writer.writerow([allele, pep])

            result = subprocess.run(
                ["python", predict_script,
                 "-i", input_csv,
                 "-m", self.mode,
                 "-d", self.device],
                capture_output=True,
                cwd=os.path.join(self.bigmhc_dir, "src"),
                timeout=600,
            )
            if result.returncode != 0:
                raise RuntimeError(
                    "BigMHC failed (exit %d):\n%s" % (
                        result.returncode,
                        result.stderr.decode("utf-8", "replace")))

            rows = []
            with open(output_csv) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    rows.append(row)
            return rows

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def predict(self, peptides):
        """Predict scores for peptides against all alleles.

        Parameters
        ----------
        peptides : list of str

        Returns
        -------
        list of PeptidePreds
            One entry per peptide; each contains one Pred per allele.
        """
        if isinstance(peptides, str):
            peptides = [peptides]

        # Build all (peptide, allele) combinations
        all_peptides = []
        all_alleles = []
        for pep in peptides:
            for allele in self.alleles:
                all_peptides.append(pep)
                all_alleles.append(allele)

        rows = self._run_bigmhc(all_peptides, all_alleles)

        # Parse into Pred objects, grouped by peptide
        score_col = [c for c in rows[0] if c.startswith("BigMHC")][0]
        kind = self._pred_kind()
        name = self._predictor_name()

        # Group by peptide in input order
        idx = 0
        results = []
        for pep in peptides:
            preds = []
            for allele in self.alleles:
                score = float(rows[idx][score_col])
                preds.append(Pred(
                    kind=kind,
                    score=score,
                    peptide=pep,
                    allele=allele,
                    predictor_name=name,
                ))
                idx += 1
            results.append(PeptidePreds(preds=tuple(preds)))
        return results

    def predict_dataframe(self, peptides, sample_name=""):
        """``predict()`` flattened to a DataFrame."""
        import pandas as pd
        from .pred import COLUMNS
        dfs = [pp.to_dataframe(sample_name) for pp in self.predict(peptides)]
        if not dfs:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)
