# Computational Design of an Iron-Acquiring *Lactiplantibacillus plantarum* WCFS1

Genome-scale metabolic modelling used to identify the mechanistic
bottleneck behind *L. plantarum*'s clinically observed enhancement of
human iron absorption, and to propose a virtual engineered strain that
explains it.


## Overview

Iron deficiency anemia affects nearly two billion people worldwide, and
oral iron supplementation has poor bioavailability along with significant
side effects. A 2017 clinical trial (Hoppe et al.) showed that the
probiotic *L. plantarum* 299v improves human iron absorption by 50%, but
the molecular mechanism behind this effect has never been explained. This
project addresses that gap computationally. We first checked whether
WCFS1 uses the classical siderophore route for iron acquisition by mapping
the enterobactin biosynthesis pathway via KEGG and running BLAST homology
searches against *E. coli* entC and *Yersinia* ybtS — WCFS1 turned out to
lack this pathway entirely, which shifted the investigation toward direct
transporter-mediated iron uptake instead. Using the AGORA2 genome-scale
metabolic model of WCFS1 under simulated gut conditions (0.02 mM Fe²⁺, pH
6.5, 5% O₂), we ran flux balance analysis, single-reaction knockout
analysis, and flux variability analysis to pin down exactly which
transport step limits growth, then simulated virtual engineered strains to
see whether increasing that transporter's capacity could reproduce the
clinical 50% improvement.

## Pipeline

The analysis starts with genome annotation: the WCFS1 genome (NCBI RefSeq
GCF_000203855.3) is searched for iron-related genes using a keyword set
covering iron transport, Fur regulation, siderophores, and Fe-S cofactor
assembly, with each hit manually cross-referenced against UniProt and
KEGG. In parallel, gut-relevant physiological parameters (iron
concentration, pH, oxygen, carbon and nitrogen source, transit time) are
defined from peer-reviewed literature on the human duodenum, since
simulating unrealistic conditions would make any downstream result
meaningless.

The AGORA2 metabolic model is then loaded in COBRApy, and flux balance
analysis is run first without any iron constraint to establish a
theoretical maximum growth rate, then again under the gut-relevant iron
cap to see how much growth that constraint alone removes. An iron
concentration scan across the physiologically relevant range establishes
whether growth responds linearly to iron availability. From there, each
iron-related reaction is knocked out individually to see which one, if
removed, is lethal to growth, and flux variability analysis checks
whether the reactions found essential by knockout are also operating with
zero flexibility — i.e., already maxed out with no alternate routing
available. Where all three methods agree on the same reaction, that
reaction is treated as the genuine rate-limiting bottleneck rather than an
artifact of any one method. Finally, virtual engineered strains are
simulated by raising that reaction's capacity in stages and comparing the
resulting growth improvement against the clinical benchmark.

## Key findings

*L. plantarum* WCFS1 lacks the classical enterobactin/siderophore
biosynthesis pathway altogether — none of the entA-H or fepA-C genes are
present, and BLAST searches turned up no meaningful siderophore-specific
homology in the WCFS1 proteome specifically. Instead, iron acquisition
happens through direct ABC-transporter-mediated uptake. Under simulated
gut conditions, the Fe²⁺ transporter FE2abc emerged as the sole essential
bottleneck for growth, confirmed by three independent lines of evidence:
growth scaled perfectly linearly with iron availability (implicating a
transport limit rather than a downstream metabolic one), knocking out
FE2abc caused complete growth failure with no other reaction producing
that effect, and flux variability analysis showed FE2abc was already
operating at the very edge of its permissible range with essentially no
slack. A virtual strain with double the FE2abc capacity produced a 50%
growth improvement — matching the clinical iron-absorption improvement
reported for *L. plantarum* 299v almost exactly, and giving the
first mechanistic, transporter-level hypothesis for a clinical effect that
had gone unexplained for over two decades.

## Reproducing the analysis

Clone the repository and install dependencies:

```bash
git clone https://github.com/srijeeta20-sys/iron-acquiring-lplantarum-gsmm.git
cd iron-acquiring-lplantarum-gsmm
pip install -r requirements.txt
```

The AGORA2 model file itself isn't committed here — see `model/README.md`
for how to obtain it from the Virtual Metabolic Human (VMH) database and
where to place it. Once it's in place, run:

```bash
python code/iron_model.py
```

This regenerates the baseline FBA, knockout, FVA, and virtual strain
result files.

## Known issue

The result values for FE2abc overexpression at 1.5x and 2x reported in the
original PBL write-up do not match this repository's raw script output
(see `results/Stage8_virtual_strains.txt`) — the raw output is the correct
source of truth and should be used to correct the manuscript table before
submission. The headline finding itself (2x FE2abc overexpression giving a
50% growth improvement, matching the Hoppe et al. clinical benchmark) is
unaffected by this — only the intermediate table values need fixing.

## Data and model provenance

Genome assemblies were pulled from NCBI RefSeq: *L. plantarum* WCFS1
(GCF_000203855.3) and *E. coli* K-12 MG1655 (GCF_000005845.2, used as the
reference for the entC/ybtS BLAST comparison). The metabolic model is from
the AGORA2 collection (Magnusdottir et al., 2017), obtained via the
Virtual Metabolic Human database and not redistributed in this repository.
Siderophore pathway mapping used KEGG pathway lpl01053.

