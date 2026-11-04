# NetMHCpan 4.1b coordinate fixtures

Captured from the installed native NetMHCpan 4.1b on 2026-11-03 with:

```sh
netMHCpan-4.1 -f input.fasta -a HLA-A02:01 -l 9 -BA
netMHCpan-4.1 -f input.fasta -a HLA-A02:01 -l 9
```

Input is `>offset_probe` followed by `ACDEFGHIKLMN`. Only the temporary
directory diagnostic was removed. The first peptide `ACDEFGHIK` starts at
Python offset 0 but native `Pos` 1; the remaining positions are 2, 3 and 4.

The [DTU output example](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1a/3-Output_format.php)
also places `HQAAMQMLK` at `Pos` 15 in `TPQDLNTMLNTVGGHQAAMQMLKETINEEA`
(Python offset 14), despite the prose describing zero-based positions.
The executable output and sequence slices determine this regression.
