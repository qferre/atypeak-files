# atypeak-files

ChIP-seq data treated by the atyPeak method (https://github.com/qferre/atypeak). Data drawn from selected datasets and transcriptional regulators from ReMap 2018.

Scores are given to each peak. A higher score means the peak has, broadly speaking, more of its usual correlators present and is less likely to be an anomalous peak.

## Contents

- *bed*: BED files giving scores for the selected peaks. Data files are given after all normalizations and corrections. Raw atyPeak scores before any correction are provided in the *raw* subdirectory.
- *model*: architecture of the atyPeak models.
- *diagnostic*: partial diagnostic data as produced by the atyPeak method per cell lines.
- *script*: other scripts used in results.