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

"""T-cell receptor input type for pMHC:TCR binding predictors.

mhctools' predictors have historically taken ``(peptide, allele)`` inputs
and scored MHC-ligand binding/presentation/processing. TCR specificity
predictors (e.g. NetTCR) add a new input modality — the receptor itself,
represented by its complementarity-determining regions (CDRs) — and a new
output kind (:attr:`~mhctools.pred.Kind.pMHC_TCR_binding`).

A :class:`TCR` is a paired αβ receptor described by its six CDR loops.
The field names use the biological convention (``cdr1a`` = CDR1 of the α
chain, ``cdr3b`` = CDR3 of the β chain, ...); the ``a1``/``b3`` short
forms used by NetTCR are exposed as properties.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields


_SHORT_CDR_ALIASES = {
    "a1": "cdr1a",
    "a2": "cdr2a",
    "a3": "cdr3a",
    "b1": "cdr1b",
    "b2": "cdr2b",
    "b3": "cdr3b",
}


@dataclass(frozen=True)
class TCR:
    """A paired αβ T-cell receptor, described by its six CDR loops.

    Parameters
    ----------
    cdr1a, cdr2a, cdr3a : str
        CDR1/2/3 amino-acid sequences of the α chain (NetTCR ``A1``/``A2``/``A3``).
    cdr1b, cdr2b, cdr3b : str
        CDR1/2/3 amino-acid sequences of the β chain (NetTCR ``B1``/``B2``/``B3``).
    name : str, optional
        Human-readable identifier (e.g. a clonotype name). If omitted,
        :attr:`identifier` falls back to the CDR3α/CDR3β pair.

    Notes
    -----
    CDR3β (``cdr3b``) carries most of the antigen-contact specificity, but
    NetTCR-2.2's pan model expects all six loops; leaving loops empty will
    encode as fully-padded (uninformative) features.
    """

    cdr1a: str = ""
    cdr2a: str = ""
    cdr3a: str = ""
    cdr1b: str = ""
    cdr2b: str = ""
    cdr3b: str = ""
    name: str = ""

    # NetTCR-style short aliases -------------------------------------------
    @property
    def a1(self) -> str:
        return self.cdr1a

    @property
    def a2(self) -> str:
        return self.cdr2a

    @property
    def a3(self) -> str:
        return self.cdr3a

    @property
    def b1(self) -> str:
        return self.cdr1b

    @property
    def b2(self) -> str:
        return self.cdr2b

    @property
    def b3(self) -> str:
        return self.cdr3b

    @property
    def identifier(self) -> str:
        """Stable string identifier for this receptor.

        Uses :attr:`name` if set, otherwise the ``CDR3α/CDR3β`` pair, which
        together define the clonotype for most practical purposes.
        """
        if self.name:
            return self.name
        return "%s/%s" % (self.cdr3a, self.cdr3b)

    def cdr_dict(self) -> dict:
        """Return the six CDRs keyed by NetTCR feature name (``a1``..``b3``)."""
        return {
            "a1": self.cdr1a,
            "a2": self.cdr2a,
            "a3": self.cdr3a,
            "b1": self.cdr1b,
            "b2": self.cdr2b,
            "b3": self.cdr3b,
        }

    def __str__(self):
        return "TCR(%s)" % self.identifier

    def __repr__(self):
        return str(self)

    def to_dict(self):
        """Serialize to a JSON-friendly dict."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        """Deserialize from a dict.

        Accepts canonical ``cdr1a``..``cdr3b`` keys plus NetTCR/IMMREP-style
        short aliases ``a1``..``b3`` / ``A1``..``B3``. Canonical field names
        win if both a canonical key and an alias are present.
        """
        valid = {f.name for f in fields(cls)}
        values = {}
        for key, value in d.items():
            lower = str(key).lower()
            canonical = lower if lower in valid else _SHORT_CDR_ALIASES.get(lower)
            if canonical is None:
                continue
            if canonical not in values or lower in valid:
                values[canonical] = value
        return cls(**values)

    @classmethod
    def from_series(cls, row):
        """Deserialize from a pandas Series or other mapping-like row."""
        if hasattr(row, "to_dict"):
            row = row.to_dict()
        return cls.from_dict(row)
