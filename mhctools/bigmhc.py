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

Models are loaded lazily on first prediction and kept in memory for
subsequent calls.
"""

import os
import sys

import pandas as pd
import torch

from .pred import Kind, Pred, PeptideResult, COLUMNS


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


def _import_bigmhc_modules(bigmhc_dir):
    """Import BigMHC's internal modules by temporarily adding src/ to sys.path."""
    src_dir = os.path.join(bigmhc_dir, "src")
    added = src_dir not in sys.path
    if added:
        sys.path.insert(0, src_dir)
    try:
        # BigMHC modules use bare imports (e.g. "from bigmhc import BigMHC")
        import bigmhc as bigmhc_mod
        import mhcenc as mhcenc_mod
        import dataset as dataset_mod
        import mhcuid as mhcuid_mod
        import encseq as encseq_mod
        return bigmhc_mod, mhcenc_mod, dataset_mod, mhcuid_mod, encseq_mod
    finally:
        if added:
            sys.path.remove(src_dir)


class BigMHC:
    """Wrapper for BigMHC presentation and immunogenicity predictions.

    Models are loaded lazily on the first call to :meth:`predict` and
    kept in memory for all subsequent calls.

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
    max_batch_size : int
        Maximum peptides per inference batch. Default 4096.
    """

    VALID_MODES = ("el", "im")

    def __init__(self, alleles, mode="el", bigmhc_path=None, device="cpu",
                 max_batch_size=4096):
        if isinstance(alleles, str):
            alleles = [alleles]
        if mode not in self.VALID_MODES:
            raise ValueError(
                "mode must be one of %s, got %r" % (self.VALID_MODES, mode))
        self.alleles = list(alleles)
        self.mode = mode
        self.bigmhc_dir = _find_bigmhc_dir(bigmhc_path)
        self.device = device
        self.max_batch_size = max_batch_size

        # Validate installation
        src_dir = os.path.join(self.bigmhc_dir, "src")
        if not os.path.isfile(os.path.join(src_dir, "bigmhc.py")):
            raise FileNotFoundError(
                "bigmhc.py not found in %s" % src_dir)

        # Lazy-loaded state
        self._models = None
        self._mhcenc = None
        self._modules = None

    def __str__(self):
        loaded = "loaded" if self._models is not None else "not loaded"
        return "BigMHC(mode=%s, alleles=%s, %s)" % (
            self.mode, self.alleles, loaded)

    def __repr__(self):
        return str(self)

    def _pred_kind(self):
        if self.mode == "im":
            return Kind.immunogenicity
        return Kind.pMHC_presentation

    def _predictor_name(self):
        return "bigmhc_%s" % self.mode

    @property
    def _model_name(self):
        return "BigMHC_EL" if self.mode == "el" else "BigMHC_IM"

    # ------------------------------------------------------------------
    # Model loading (lazy, cached)
    # ------------------------------------------------------------------

    def _ensure_loaded(self):
        """Load models and MHC encoder on first use."""
        if self._models is not None:
            return

        self._modules = _import_bigmhc_modules(self.bigmhc_dir)
        bigmhc_mod, mhcenc_mod, _, _, _ = self._modules

        # Load MHC pseudo-sequence encoder
        pseudoseqs = os.path.join(self.bigmhc_dir, "data", "pseudoseqs.csv")
        self._mhcenc = mhcenc_mod.MHCEncoder.read(pseudoseqs)

        # Resolve model paths (same logic as cli._parseModel)
        models_dir = os.path.join(self.bigmhc_dir, "models")
        el_dirs = [
            os.path.join(models_dir, "bat%d" % (2 ** x))
            for x in range(9, 16)
        ]
        if self.mode == "im":
            model_dirs = [os.path.join(d, "im") for d in el_dirs]
        else:
            model_dirs = el_dirs

        # Load ensemble
        self._models = []
        for d in model_dirs:
            model = bigmhc_mod.BigMHC.load(d)
            model.eval()
            model = bigmhc_mod.BigMHC.accelerate(model, devices=[])
            self._models.append(model)

    # ------------------------------------------------------------------
    # Core prediction (in-process)
    # ------------------------------------------------------------------

    def _predict_raw(self, peptides, alleles):
        """Run BigMHC in-process. Returns list of float scores."""
        self._ensure_loaded()
        _, _, dataset_mod, mhcuid_mod, _ = self._modules

        # Build DataFrame in the format BigMHC expects
        df = pd.DataFrame({"mhc": alleles, "pep": peptides})
        df["pep"] = df["pep"].str.upper()
        df["tgt"] = None
        df["len"] = df["pep"].apply(len)
        df.insert(
            0, "uid",
            df["mhc"].apply(mhcuid_mod.mhcuid) + "_" + df["pep"])
        df.set_index("uid", drop=True, inplace=True)
        df = df.astype({
            "mhc": "string", "pep": "string",
            "tgt": "float32", "len": "int8"})

        # Create Dataset and batch
        data = dataset_mod.Dataset(pmhcs=df, mhcenc=self._mhcenc)
        data.makebats(self.max_batch_size)
        loader = torch.utils.data.DataLoader(
            data, batch_size=None, shuffle=False,
            num_workers=0)

        # Run inference
        preds_parts = []
        with torch.no_grad():
            for idx, bat in enumerate(loader):
                out_stack = []
                for model in self._models:
                    dev = next(model.parameters()).device
                    out, _ = model(
                        mhc=bat.mhc.to(dev),
                        pep=bat.pep.to(dev))
                    out_stack.append(torch.sigmoid(out))
                ensemble_out = torch.mean(torch.stack(out_stack), dim=0)
                rawbat = data.getbat(idx=idx, enc=False)
                rawbat[self._model_name] = ensemble_out.cpu().numpy()
                preds_parts.append(rawbat)

        result_df = pd.concat(preds_parts).sort_index()
        return result_df[self._model_name].tolist()

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
        list of PeptideResult
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

        scores = self._predict_raw(all_peptides, all_alleles)

        kind = self._pred_kind()
        name = self._predictor_name()

        idx = 0
        results = []
        for pep in peptides:
            preds = []
            for allele in self.alleles:
                preds.append(Pred(
                    kind=kind,
                    score=float(scores[idx]),
                    peptide=pep,
                    allele=allele,
                    predictor_name=name,
                ))
                idx += 1
            results.append(PeptideResult(preds=tuple(preds)))
        return results

    def predict_dataframe(self, peptides, sample_name=""):
        """``predict()`` flattened to a DataFrame."""
        dfs = [pp.to_dataframe(sample_name) for pp in self.predict(peptides)]
        if not dfs:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)


class BigMHC_EL(BigMHC):
    """BigMHC in presentation (eluted-ligand) mode."""
    def __init__(self, alleles, bigmhc_path=None, device="cpu",
                 max_batch_size=4096):
        BigMHC.__init__(
            self, alleles, mode="el", bigmhc_path=bigmhc_path,
            device=device, max_batch_size=max_batch_size)


class BigMHC_IM(BigMHC):
    """BigMHC in immunogenicity mode."""
    def __init__(self, alleles, bigmhc_path=None, device="cpu",
                 max_batch_size=4096):
        BigMHC.__init__(
            self, alleles, mode="im", bigmhc_path=bigmhc_path,
            device=device, max_batch_size=max_batch_size)
