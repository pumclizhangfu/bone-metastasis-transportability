# Pan-cancer bone-metastasis single-cell and spatial transcriptomics

This directory is the audited code release for the manuscript-associated
reproducibility materials. The public repository is
`https://github.com/pumclizhangfu/bone-metastasis-transportability`. Original
project code is released under the MIT License. The directory contains sanitized
analysis/rendering code, frozen configuration, provenance ledgers, and an offline
verifier. The separate local `companion/` directory contains accepted static
figures and source-data payload but is outside the present public-code
authorization and is not part of this GitHub release.

## Verified in Gate12BM

- frozen Gate12BL figures and source-data files are present and hash-verified;
- allow-listed public scripts parse successfully;
- manifests can be verified offline with Python's standard library;
- private paths, personal identifiers, credentials, disallowed raw objects, and
  repository files larger than 5 MiB are rejected.

Run from this directory:

```bash
workflow/run.sh verify
```

When using the two candidate archives, extract both into the same parent
directory; they create sibling `repository/` and `companion/` directories.

## Not yet validated

`render`, `analysis-frozen`, and `full-raw` intentionally return
`DOCUMENTED_NOT_YET_VALIDATED`. The project currently lacks a validated public
lockfile/container, distributable discovery SCE objects, distributable NicheNet
priors, and an end-to-end clean-environment reconstruction receipt.

## Author actions before public deposit

1. verify the public GitHub landing page and its release checksum after push;
2. insert the real GitHub URL in the manuscript Data Availability Statement;
3. separately authorize the companion/source-data payload before any Zenodo
   deposit; this has not been inferred from the code-only authorization;
4. retain all third-party exclusions and the disclosed raw-to-final replay limit.

The GitHub URL and MIT licence were author-authorized. No DOI is claimed. Author
metadata and Zhangfu Li's ORCID were supplied by the author; the ORCID checksum
passed validation.
