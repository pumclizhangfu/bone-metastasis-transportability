# Gate12BM: sanitized public reproducibility candidate

## Purpose

Gate12BM builds a local, journal-neutral candidate for public deposition. It does
not upload files, create a GitHub repository, mint a DOI, choose a software
license, or insert author metadata. Those actions remain under author control.

The candidate is split into two deliberately different payloads:

1. `repository/` is the GitHub-oriented code, configuration, documentation, and
   offline-verification payload. It is restricted to a frozen allow-list and is
   scanned for credentials, personal paths, personal identifiers, symlinks,
   redistributability risks, and unexpectedly large files.
2. `companion/` is the Zenodo-oriented static evidence payload copied from the
   accepted Gate12BL release: final figures, source-data tables, supplementary
   tables, and the active manuscript. Its files are independently rehashed and
   remain subject to an explicit author redistribution review before deposit.

## Scientific and reproducibility boundary

Gate12BM may establish that the accepted static publication payload is intact,
that its source-data files are present, and that the selected analysis/rendering
scripts are syntactically valid after path sanitization. It must not claim that:

- the raw data can currently be rebuilt end to end in a clean environment;
- every figure can currently be rerendered from the deposited source data;
- cross-platform graphics are byte-identical;
- restricted NODE/OEP matrices can be redistributed;
- third-party reference objects are licensed for redistribution;
- a DOI, repository URL, or software/data license already exists.

The public interface therefore has four modes:

| mode | Gate12BM target status | meaning |
|---|---|---|
| `verify` | validated offline | verify manifests, hashes, file counts, security rules, and static evidence contracts |
| `render` | documented, not validated | future figure rerendering entry point; no success claim in Gate12BM |
| `analysis-frozen` | documented, not validated | future rerun from frozen processed inputs; required inputs are not all distributable |
| `full-raw` | documented, not validated | future raw-to-final reconstruction; clean-environment validation is outside Gate12BM |

## Authoritative upstream release

- Gate12BL output: `results/gate12bl_reproducibility_release/run_v1`
- Gate12BL receipt SHA-256:
  `43f0aebc80882a03e262f4ef51d58768cefaec92569eeefe8679347dee903f2b`
- Active manuscript SHA-256:
  `9e8d70f511815be1338a0d2fc3d7c6fc20942dc5a6fd6ee462f8741b98406ac5`
- Supplementary workbook SHA-256:
  `b6a8072977325ba1101278c7cc6d86a952c7ecc369978878a09d5d34f9734398`
- Expected final figures: 42 files (21 PNG/PDF pairs)
- Expected source-data files: 85

Gate12BM refuses to run if the frozen receipt or active manuscript hash differs.

## Frozen repository allow-list

Only files listed in the following manifests can enter `repository/`:

- `config/gate12bm_script_whitelist_v1.tsv`
- `config/gate12bm_config_whitelist_v1.tsv`

The script manifest assigns each file to acquisition, discovery, UMAP,
communication, external validation, spatial analysis, rendering, or audit/release.
Files retain their original flat `scripts/` and `config/` relative paths so that
existing project-root discovery and configuration lookups are not broken merely
for visual organization. The config manifest records the parameters needed by those stages. Source and
candidate SHA-256 values are recorded separately because private absolute-path
defaults are sanitized in the candidate copy without modifying the canonical
project scripts.

## Excluded material

The following material is excluded by design:

- passwords, tokens, private keys, SSH configuration, server host names, private
  IP addresses, personal email addresses, and personal filesystem paths;
- patient-level or quasi-identifying clinical spreadsheets and audit tables;
- restricted NODE/OEP expression matrices and SFTP credentials/paths;
- third-party RDS/RData objects, raw SOFT archives, raw sequencing archives, and
  cached bytecode;
- local logs, temporary files, environment receipts containing personal paths,
  and historical superseded release directories;
- any invented `LICENSE`, GitHub URL, Zenodo DOI, ORCID, affiliation, or author
  list.

Public accession identifiers and coded public sample identifiers may be retained
when necessary for scientific traceability. Their role is documented in the data
ledger and must not be presented as newly de-identified clinical data.

## Sanitization rules

Candidate script/config copies are transformed only by explicit literal rules:

- known macOS and Ubuntu project roots become `${PROJECT_ROOT}` placeholders in
  documentation or runtime root discovery in code;
- the private bundled Python path becomes `python3`;
- the private bundled LibreOffice path becomes `soffice`;
- personal email addresses and user names are rejected, not silently retained;
- unresolved absolute paths are a hard failure.

Each transformed file receives both a source hash and a candidate hash plus a
machine-readable transformation record.

## Data-role ledger

The candidate must list at least the following public-accession roles:

| accession | role in the manuscript/repository |
|---|---|
| GSE143791 | discovery single-cell cohort |
| GSE202813 | external single-cell cohort |
| GSE266330 | external paired/clinical single-cell support |
| GSE323357 | spatial transcriptomic cohort |
| OEP005136 | NODE/BMDC external bone-metastasis cohort; access-controlled matrix |
| OEZ00021715 / OED01122886 | NODE analysis/data identifiers for the paired mBone archive |
| GSE225209 | supportive communication co-detection dataset used in the active analysis |
| GSE190772 | prespecified supportive dataset retained as non-evaluable under the primary eligibility rule |

## New disclosure blocker

The Gate12BL manuscript uses GSE225209 and records GSE190772 in the communication
support analysis, but its current Data availability paragraph lists only five
core accessions. Gate12BM must add the following blocker without silently editing
the manuscript:

- `BI-18 DATASET_DISCLOSURE_GAP` (`BLOCKING`, `UNRESOLVED`): add GSE225209 and
  GSE190772 to the study/data-role and Data availability disclosures, explicitly
  stating that GSE190772 was non-evaluable under the prespecified eligibility
  rule.

This blocker is expected to be repaired in a later manuscript-alignment gate.

## Candidate contents

### `repository/`

- `README.md`, `CHANGELOG.md`, and `CITATION.cff.in`;
- license-decision placeholder and third-party reuse matrix, but no asserted
  license;
- frozen script/config allow-lists and sanitized copies;
- accession, authoritative-run, method-to-code, claim-to-code, stage-map, and
  artifact-DAG ledgers;
- environment history and explicit missing-lockfile/container limitations;
- portable offline verifier and data-contract tests;
- manifest of every repository file.

Repository and companion manifests deliberately exclude their own manifest file
to avoid a self-hash cycle. The top-level output manifest similarly excludes
itself and `GATE12BM_RECEIPT.json`; the receipt stores the output-manifest hash,
and its final SHA-256 is reported externally after execution.

### `companion/`

- exact copies of the 42 accepted Gate12BL figure files;
- exact copies of the 85 accepted Gate12BL source-data files;
- accepted supplementary workbook, PDF, and table legends;
- active manuscript and figure legends;
- a sanitized Gate12BL source receipt and full companion manifest.

### `admin/`

- Gate12BM receipt;
- acceptance matrix;
- blocker ledger including inherited Gate12BL blockers and BI-18;
- security/privacy scan;
- source-to-candidate transformation ledger;
- package and archive hashes;
- smoke-test log and independent verification results.

## Acceptance matrix

Formal execution must report every dimension independently:

| dimension | pass criterion |
|---|---|
| upstream identity | frozen Gate12BL receipt/manuscript hashes match |
| package structure | all required repository, companion, and admin paths exist |
| manifest integrity | all files independently rehash to their manifests |
| static publication payload | 42 figure files and 85 source-data files match Gate12BL |
| code syntax | all copied Python/R/shell scripts parse without writing cache files |
| privacy/security | no credential, private-infrastructure, personal identifier, or private absolute-path hit |
| repository size/type policy | no disallowed binary/raw object; no repository file over 5 MiB; no symlink |
| companion boundary | all larger/static files are isolated and tagged for author redistribution review |
| accession ledger | all eight accession/data-identifier roles above are present |
| offline verification | `repository/workflow/run.sh verify` succeeds on the built candidate |
| figure replay | `NOT_VALIDATED` unless a clean rerender is actually run |
| clean-environment smoke | `NOT_VALIDATED` unless a lock/container is created and tested |
| full raw reconstruction | `NOT_VALIDATED` unless all access-controlled inputs are supplied and the DAG succeeds |
| redistribution clearance | `PENDING_AUTHOR_APPROVAL` |
| repository license | `NOT_SELECTED` |
| persistent identifier | `NOT_ASSIGNED` |
| manuscript disclosure | `FAIL` while BI-18 remains unresolved |
| submission readiness | `NOT_READY` while any blocking item remains |

## Formal verdict vocabulary

Gate12BM may return one of:

- `CANDIDATE_BUILT_WITH_BLOCKERS`
- `DEPOSITION_READY_PENDING_AUTHOR_ACTION`
- `FAILED`

For the frozen input state described above, the expected honest verdict is
`CANDIDATE_BUILT_WITH_BLOCKERS`, because BI-18, license selection,
redistribution clearance, persistent identifier assignment, and clean
raw-to-final validation remain unresolved.

## Execution policy

Preparation and smoke testing use a disposable temporary directory only. Formal
execution writes once to:

`results/gate12bm_public_repository_candidate/run_v1`

The formal builder refuses to overwrite an existing run. External upload,
repository creation, DOI minting, license selection, and manuscript modification
are all outside Gate12BM and require separate author confirmation.
