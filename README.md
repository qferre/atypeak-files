# atypeak-files

ChIP-seq data treated by the atyPeak method. The data was drawn from selected datasets and transcriptional regulators from ReMap 2018.

Scores are given to each peak. A higher score means the peak has, broadly speaking, more of its usual correlators present and is less likely to be an anomalous peak.

The source code for the atyPeak method can be found at: https://github.com/qferre/atypeak


## Principle

ChIP-seq can be obscured by data anomalies and biases. *atyPeak* is a deep-learning based method to identify ChIP-seq peaks that have “atypical profiles”, meaning that those peaks are found without their usual collaborators. "Collaborators" are defined as the other Transcriptional Regulators (or corroborating datasets) which are usually found in the same neighborhood, in the Cis Regulatory Elements of this cell line. Peaks get a higher score when more of their usual collaborators are present.

Each file in the *bed* directory is a BED file. Analysis was restricted to the densest 65K CRM in the human genome, and to selected TFs and selected datasets. Each file gives, for a given cell line, the atyPeak score in the “score” column.

Anomaly score thresholds and interpretation are at the user’s discretion. Anomalies usually represent noisy peaks. However, a focused study of a single experimental series may rely on low-scoring peaks as they might be caused by certain events of interest.
Regions with a high density of high-scoring peaks are the strongest candidate CREs. Furthermore, when comparing your own ChIP-seq data to those tracks, similarity (meaning peaks for the same TFs in the same cell line at the same position) to high scoring peaks suggests this is a typical event. But if the data matches a low scoring peak or no peak at all, it is atypical. What this means precisely depends on the data being studied: it can be noise, or if you are confident in your data it could be an interesting anomaly worth studying.

For more information, please see *Ferré et al.* <DOI: TBD>


## Contents

This data is also available at: http://remap2020.univ-amu.fr/download_page

- *bed*: BED files giving scores for the selected peaks. Data files are given after all normalizations and corrections. Raw atyPeak scores before any correction are provided in the *raw* subdirectory.
- *model*: architecture of the atyPeak models.
- *diagnostic*: partial diagnostic data as produced by the atyPeak method per cell lines.
- *script*: other scripts used in results.
