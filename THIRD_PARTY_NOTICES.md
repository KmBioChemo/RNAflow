# Third-party data and resources

RNAflow bundles a small number of third-party datasets and biological
networks so that its demo loaders and offline activity inference work
without a network connection. **The RNAflow MIT license (see `LICENSE`)
covers only RNAflow’s own source code, not these resources.** Each
remains under its own upstream license and terms — please cite the
original sources when you use them, and check the upstream terms before
any commercial use.

## Bundled demo datasets (`inst/extdata/`)

### `demo_airway_*.csv` — airway

- **Source:** airway smooth-muscle RNA-seq (Himes *et al.* 2014), GEO
  accession GSE52778, distributed as the Bioconductor `airway`
  experiment-data package.
- **Reference:** Himes BE, *et al.* “RNA-Seq transcriptome profiling
  identifies CRISPLD2 as a glucocorticoid responsive gene…” *PLoS
  ONE* 2014. <doi:10.1371/journal.pone.0099625>
- A subset re-exported as a plain counts + metadata CSV for
  demonstration. Regenerate with `dev/make_demo_airway.R`.

### `demo_tcga_*.csv` — TCGA pan-cancer via GSE62944

- **Source:** TCGA RNA-seq reprocessed by Rahman *et al.* 2015, GEO
  accession GSE62944, distributed via the Bioconductor `GSE62944`
  package. Underlying TCGA data are open-access.
- **Reference:** Rahman M, *et al.* “Alternative preprocessing of
  RNA-Sequencing data in The Cancer Genome Atlas leads to improved
  analysis results.” *Bioinformatics* 2015.
  <doi:10.1093/bioinformatics/btv377>
- A small many-cancer subset re-exported for demonstration. Regenerate
  with `dev/make_demo_tcga.R`.

## Bundled activity networks (`inst/extdata/`)

Both are convenience snapshots retrieved from
[OmniPath](https://omnipathdb.org/); see `data-raw/make_networks.R` for
the exact queries. Human only.

### `collectri_human.rds` — CollecTRI transcription-factor regulons

- **Source:** CollecTRI, retrieved from OmniPath with the
  `license=academic` filter.
- **Reference:** Müller-Dott S, *et al.* “Expanding the coverage of
  regulons from high-confidence prior knowledge for accurate estimation
  of transcription factor activities.” *Nucleic Acids Research* 2023.
  <doi:10.1093/nar/gkad841>
- CollecTRI is a **composite** resource; its constituent databases
  retain their own licenses. The `license=academic` filter restricts the
  download to resources that permit academic use. Commercial users must
  review the upstream terms of each constituent resource.

### `progeny_human.rds` — PROGENy pathway footprints

- **Source:** PROGENy model, retrieved from OmniPath.
- **Reference:** Schubert M, *et al.* “Perturbation-response genes
  reveal signaling footprints in cancer gene expression.” *Nature
  Communications*
  2018. <doi:10.1038/s41467-017-02391-6>
- PROGENy is distributed under the **Apache License 2.0**.

## Notes

- These files are snapshots for convenience and reproducibility. To
  rebuild them from the upstream sources, run `data-raw/make_networks.R`
  (networks) and `dev/make_demo_airway.R` / `dev/make_demo_tcga.R`
  (demos).
- Mouse and rat activity inference still fetch networks live from
  OmniPath; only the human snapshots are bundled.
- If you redistribute RNAflow, keep this notice and consult each
  resource’s current upstream license before commercial use.
