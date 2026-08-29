Metabolic model
This folder contains `Lactobacillus_plantarum_WCFS1.xml`, the AGORA2
genome-scale metabolic reconstruction of Lactiplantibacillus plantarum
WCFS1 (Heinken et al., 2023; AGORA2 collection), obtained from the
Virtual Metabolic Human (VMH) database
(https://www.vmh.life/#microbe/Lactobacillus_plantarum_WCFS1).
Model specifications (verified directly from this file)
Parameter	Value
Reactions	1,245
Metabolites	1,083
Genes	902
Objective function	`biomass524`
Format	SBML
AGORA2 version	2.0
Citation
Heinken, A., Hertel, J., Acharya, G., Ravcheev, D.A., Nyga, M., Okpala,
O.E., Hogan, M., Magnusdottir, S., Martinelli, F., Nap, B., Preciat, G.,
Edirisinghe, J.N., Henry, C.S., Fleming, R.M.T., Thiele, I. "AGORA2:
Knowledge-driven genome-scale reconstruction of 7,302 human microbes for
personalised medicine." Nature Biotechnology (2023).
https://doi.org/10.1038/s41587-022-01628-0
License
AGORA/AGORA2 reconstructions are distributed by VMH under a
Creative Commons Attribution-NonCommercial 2.0 Generic license
(CC BY-NC 2.0). This model file is redistributed here under those
terms — attribution given above, and this repository and its associated
academic publication are non-commercial. If you reuse this model file
outside this repository, retain the citation above and do not use it for
commercial purposes without separate permission from VMH/the Thiele Lab.
Genome source used for annotation (not the metabolic model itself)
Genome annotation (GFF3) used for the manual iron-gene curation step was
downloaded from NCBI RefSeq, not re-hosted here — regenerate via:
```bash
datasets download genome accession GCF_000203855.3 --include gff3,protein,genome
```
Reference genomes used in this project:
Organism	Assembly accession
L. plantarum WCFS1	GCF_000203855.3
E. coli K-12 MG1655	GCF_000005845.2
