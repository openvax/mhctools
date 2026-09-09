# Peptidase cleavage evidence

```python
from mhctools import CleavageInput, DPP4qPISA

peptide = CleavageInput("HAEGTFTSD", source_id="GLP-1 fragment")
result = DPP4qPISA().predict(peptide)
print(result.sites[0].bond)   # 2: HA | EGTFTSD
print(result.sites[0].score)  # 2.1694, native qPISA score
print(result.to_dict())      # input, model, assay, limitations and source bonds
```

This API describes individual peptide bonds. Bond `b` splits
`sequence[:b] | sequence[b:]`; `source_bond` adds the zero-based
`source_start` offset. Existing proteasomal scoring APIs remain available.

Inputs describe canonical, linear L-peptides. Strings assume free N and C
termini. Use `CleavageInput` to explicitly record `n_term="acetylated"` or
`c_term="amidated"`, or `unknown`. These forms are recorded but the initial
qPISA implementation abstains because they are outside its supported domain.
Other modifications, D-residues, cyclization and conjugates are unsupported;
do not strip their chemistry to obtain a score.

Exopeptidases assess exposed termini. An internal `HAE` is not an immediate
DPP4 site. To ask about a hypothesized product, use
`parent.fragment(start, end, n_term="free", c_term="free")` and predict that
fragment. The source offset is retained. This conditional analysis does not
predict formation of the fragment, cleavage order, rates or competition.

## Human DPP4 qPISA

The [Gudipati et al. 2024 paper](https://doi.org/10.1038/s44320-024-00071-4)
reports a model fitted to substrate depletion by purified human DPP4. The
source assay used tryptic HeLa peptides in HEPES pH 7.4, at 21 C for 4 hours.
The implementation independently evaluates the three rearranged terms in
Dataset EV2: `P1 + P2:P1 + P1:P1-prime`, where the first three peptide
residues are P2, P1 and P1-prime. Only bond 2 is assessed.

Higher scores indicate greater predicted log2 depletion relative to buffer
control in that assay. Negative values are retained. Scores are **not
cleavage probabilities, serum half-lives or stability ranks across enzymes**.
Compartment metadata records where the enzyme can act, not where the model
has been calibrated. Physiological exposure, structure, concentration and
competing enzymes are not modeled.

All triplets with complete coefficients can be evaluated, including those
without Pro or Ala at P1. Of 8,000 canonical triplets, 6,420 have complete
coefficients; the other 1,580 return an explicit missing-coefficient reason.
Evaluability does not establish that an individual triplet occurred in the
training data. The related C. elegans DPF-3 model is not used for human DPP8/9.

### Parameter provenance

`mhctools/data/dpp4_qpisa.json` contains the numeric cells from
`44320_2024_71_MOESM3_ESM.xlsx`, sheet `dpp4_modelParams`. Published `NA`
cells become JSON `null`; no numerical imputation or refitting is performed.
The source workbook SHA-256 is
`ee449da13b5ec66fd6fb08c16203da8c44f2e9c10572a355abdcb985694c6ddc`.
It was retrieved via the [Europe PMC supplementary archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11612144/supplementaryFiles).
The article assigns associated data to [CC0](https://creativecommons.org/publicdomain/zero/1.0/),
unless otherwise credited; Dataset EV2 has no separate credit restriction.
Attribution: Rajani Kanth Gudipati and colleagues, 2024, DOI above. No figures
or upstream R source code are redistributed.

## Interpreting rule-based evidence

The shared result contract also supports `motif_rule` evidence, reported as
`matched` or `not_matched` with no numerical score. A rule match indicates a
known recognition pattern; a non-match does not establish resistance.
`unsupported_reason` means the input could not be assessed and carries no
score. These states must remain separate in downstream displays and ranking.

## Model panel and JSON command

```sh
mhctools cleavage --list-models
mhctools cleavage --sequence RPPGFSPFR --model app2-xp --model cpn-basic
mhctools cleavage --sequence VPYGSFKHV --compartment cytosol --out cleavage.json
mhctools cleavage --sequence HAEGTFTSD --model dpp4-qpisa --n-term acetylated
```

```python
from mhctools import predict_cleavage
results = predict_cleavage("TSGPNQ", models=["fap-endo-gp", "prep-pro"])
```

The default panel evaluates all 11 built-in models and returns separate
results. `--model` and `--sequence` can be repeated. The JSON output retains
unmatched and unsupported results, source coordinates, chemistry and model
provenance. `--source-start` is a zero-based offset shared by the supplied
inputs; use separate calls when fragments have different offsets.

Compartment filtering uses exact, conservative enzyme-location annotations.
`serum`, `plasma`, `extracellular`, `cytosol`, and `er` are distinct. The
`extracellular` filter is useful for broader candidate screening; `serum`
is not an exhaustive inventory of everything potentially present in a serum
sample. Presence, concentration, activation, inhibitors and exposure are not
inferred. No combination is converted into overall stability.

| Model | Assessed recognition pattern | Main scope and primary evidence |
| --- | --- | --- |
| `dpp4-qpisa` | N-terminal P2-P1\|P1′ score | Human DPP4; quantitative model described above |
| `cpn-basic` | C-terminal \|Lys/Arg | Human plasma CPN; [Oshima et al. 1975](https://doi.org/10.1016/0003-9861(75)90104-6) |
| `app2-xp` | N-terminal X\|Pro | Human XPNPEP2; [Molinaro et al.](https://pubmed.ncbi.nlm.nih.gov/15361070/) |
| `fap-dipeptidyl` | N-terminal X-Pro\|non-Pro | Human FAP; [Edosada et al.](https://pubmed.ncbi.nlm.nih.gov/16410248/) |
| `fap-endo-gp` | Gly-Pro\|non-Pro | FAP endopeptidase; [substrate profiling](https://pubmed.ncbi.nlm.nih.gov/16480718/), [prime-side constraint](https://pubmed.ncbi.nlm.nih.gov/22750443/) |
| `enpep-acidic` | N-terminal Asp/Glu\|X | Aminopeptidase A; [human specificity study](https://pubmed.ncbi.nlm.nih.gov/23888046/); calcium and sequence affect activity |
| `anpep-ala` | N-terminal Ala\|X preference | Aminopeptidase N; [human structure/biochemistry](https://pubmed.ncbi.nlm.nih.gov/22932899/); many other substrates omitted |
| `dpp8-xp-xa` | N-terminal X-(Pro/Ala)\|non-Pro | Cytosolic DPP8; [characterization](https://pubmed.ncbi.nlm.nih.gov/11012666/), [degradomics](https://pmc.ncbi.nlm.nih.gov/articles/PMC3656252/) |
| `dpp9-xp-xa` | N-terminal X-(Pro/Ala)\|non-Pro | Cytosolic DPP9; [antigen processing](https://pubmed.ncbi.nlm.nih.gov/19667070/), same degradomics study |
| `prep-pro` | Internal X-Pro\|X | PREP/POP, conservative 4–30-residue domain; [human profiling](https://pubmed.ncbi.nlm.nih.gov/22750443/); flanking preferences omitted |
| `erap2-basic` | N-terminal Arg/Lys\|X preference | ERAP2 in ER; [biochemistry](https://pubmed.ncbi.nlm.nih.gov/12799365/), [peptide structures](https://pubmed.ncbi.nlm.nih.gov/26381406/); not a full-context predictor |
| `eramer-step` (optional) | Initial N-terminal bond; length-specific PWM score | ERAP1 in ER; [ERAMER](https://pubmed.ncbi.nlm.nih.gov/38925438/), 9–16 residues |

The motif models return decisions without numerical scores. For example,
APP removes the first residue of `RPPGFSPFR` at `R|PPGFSPFR`; DPP-like
activity would remove two residues and is a separate assessment. FAP's
endopeptidase rule can also assess an N-acetylated input, supported by its
blocked-substrate activity. The other initial rules conservatively accept
free termini only; an unsupported modified input does not establish that
its bonds resist enzymatic cleavage.

### ERAP1 and existing processing models

```sh
mhctools fetch eramer
# Install openpyxl in the environment if it is not already available.
mhctools cleavage --sequence LAAAFGAAA --model eramer-step
```

Alternatively, construct `ERAMERCleavage(pwm_path="/path/to/PWM.xlsx")` or set
`ERAMER_HOME`. This optional model is listed without loading assets and is
excluded from the default panel. Explicit selection fails clearly if the
external asset or its runtime is missing.

`eramer-step` computes the existing ERAMER intermediate PWM specificity for
one exposed 9–16-residue precursor, assigning it to bond 1. It does not report
the average of later trimming intermediates. Its version contains the SHA-256
of the actual workbook snapshot used for inference. The GPL-licensed workbook
is loaded at runtime and is not included in the mhctools distribution.

The existing `ERAMER` cascade API, `NetChop`, `Pepsickle` and other proteasome
predictors remain available through their existing interfaces. Their output
scales must not be mixed with qPISA scores or motif decisions.

## Wider candidate inventory and next PRs

This is a starting panel, not a complete inventory of human proteolysis.
The next work is prioritized by useful substrate coverage and evidence:

1. **Benchmark the shipped models** on held-out, assay-specific substrates:
   [#291](https://github.com/openvax/mhctools/issues/291). Include observed
   non-cleavages, terminal modifications, homologous-sequence leakage checks,
   and coverage/abstention. Purified-enzyme turnover and disappearance of intact
   peptide in serum are separate endpoints.
2. **ACE, MME/neprilysin and activated CPB2/TAFI**:
   [#326](https://github.com/openvax/mhctools/issues/326). ACE includes both
   terminal dipeptide processing and exceptional substrate behavior. MME needs
   length and local sequence context. CPB2 requires activation state. Human
   plasma bradykinin data demonstrate [concentration-dependent enzyme contributions](https://pubmed.ncbi.nlm.nih.gov/10749699/);
   serum carboxypeptidase activity also [differs from plasma](https://pubmed.ncbi.nlm.nih.gov/2760564/).
3. **TPP2, NPEPPS, THOP1, NLN and XPNPEP1**, plus endosomal **LNPEP/IRAP**:
   [#327](https://github.com/openvax/mhctools/issues/327). N-terminal removal
   of three residues describes TPP2 topology, not a selective substrate model.
   [TPP2 precursor processing](https://pubmed.ncbi.nlm.nih.gov/16849449/)
   and [THOP1/NLN substrate studies](https://pubmed.ncbi.nlm.nih.gov/11284698/)
   are starting evidence. IRAP belongs in an endosomal context, not a default
   ER panel.
4. **Activated blood and inflammation-associated proteases**: F2/thrombin,
   PLG/plasmin, KLKB1, ELANE, CTSG, PRTN3 and relevant extracellular cathepsins
   and matrix metalloproteases. Curate activation, inhibitors, tissue exposure,
   and assay conditions before choosing a predictor. A generic Arg/Lys or
   hydrophobic-residue scan would not distinguish these enzymes.
5. **Other exopeptidases/compartments**: DPP7, DPP2-family annotations,
   cathepsins C/H, ACE2, carboxypeptidases M/E and lysosomal/endosomal processing
   merit separate domain reviews. Do not infer ER or serum activity merely
   from membership in a peptidase family.

Downstream sequence overlays and fragment/construct projections are tracked
in [topiary #288](https://github.com/openvax/topiary/issues/288) and
[vaxrank #422](https://github.com/openvax/vaxrank/issues/422). They should retain
exact fragment chemistry and parent coordinates rather than treating internal
exo motifs as immediate cleavage sites. The broader validated serum-cleavage
work remains tracked in [#278](https://github.com/openvax/mhctools/issues/278).
