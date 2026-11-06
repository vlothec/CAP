# Dependency source-obligation review

This review applies to packages in the **final Linux runtime image**, not the
larger development environment. Package metadata is evidence for review, not a
legal conclusion.

## Workflow

1. Generate a compliance bundle from the final image.
2. Classify GPL/LGPL-family candidates:

   ```bash
   scripts/compliance/classify-source-obligations.py \
     BUNDLE/inventory/conda-source-urls.csv \
     BUNDLE/source-obligations.csv
   ```

3. Apply CAP's reviewed conservative policy. It fails if the candidate set has
   changed, forcing explicit review of new or removed dependencies:

   ```bash
   scripts/compliance/apply-source-review-policy.py \
     BUNDLE/source-obligations.csv \
     BUNDLE/source-obligations-reviewed.csv
   ```

   The current policy marks all 53 candidates `yes`. This preserves evidence for
   GPL/LGPL packages, permissive alternatives, and exception-bearing runtimes;
   it does not claim that every package legally requires source distribution.
4. Collect exact binary evidence and verified upstream archives:

   ```bash
   scripts/compliance/collect-conda-sources.py \
     BUNDLE/source-obligations-reviewed.csv \
     BUNDLE/source/conda
   ```

   During review, `--conservative` collects all candidates, including unresolved
   `REVIEW` rows. This preserves source proactively but does not resolve their
   legal classification.
5. For every `yes`, preserve the exact upstream source archive and checksum,
   plus the exact package recipe, patches, and build scripts used for the
   distributed binary.
6. Collect exact embedded Conda notice evidence for all runtime packages:

   ```bash
   scripts/compliance/collect-conda-notices.py \
     BUNDLE/inventory/conda-source-urls.csv \
     BUNDLE/notices/conda
   ```

   The current 186-package runtime has embedded license or notice files in 154
   exact artifacts when both Conda metadata and payload archives are inspected.
   The other 32 require source-level review; do not treat a declared license
   expression as a complete notice. Collect the reviewed supplemental source,
   extract source notice evidence, and apply the exact review policy:

   ```bash
   scripts/compliance/collect-conda-notice-sources.py \
     licenses/conda-notice-supplemental-sources.csv \
     BUNDLE/source/conda-notice-supplemental
   scripts/compliance/collect-conda-source-notices.py \
     BUNDLE/notices/conda/conda-notice-manifest.csv \
     BUNDLE/source/conda-all/source-to-binary.csv \
     BUNDLE/source/conda-all BUNDLE/notices/conda-source \
     --shared-license-dir BUNDLE/notices/R-shared-licenses
   scripts/compliance/finalize-conda-source-notices.py \
     BUNDLE/notices/conda-source/source-notice-manifest.csv \
     licenses/conda-source-notice-review.csv \
     BUNDLE/notices/R-shared-licenses \
     BUNDLE/notices/conda-source/final-review.csv
   ```

   Merge the primary and supplemental source mappings/corpora into
   `BUNDLE/source/conda-all` before extraction. Export the canonical license
   directory from the exact R runtime. Validation resolved the 32 gaps as 11
   source license sets, 20 R `DESCRIPTION` declarations backed by those exact
   canonical texts, and one packaging exception whose exact recipe declares MIT
   but whose two-file helper payload omits the license text. Preserve that exact
   recipe and payload as evidence; do not claim a bundled notice for it.
   Nine Bioconductor recipes provide only MD5 for their exact source inputs;
   the supplemental collector verifies those recipe digests and records
   SHA-256 for preservation.
7. Export the installed Debian binary-to-source map and copyright files, then
   collect exact Debian source:

   ```bash
   scripts/compliance/export-debian-compliance-inputs.sh IMAGE BUNDLE/debian-inputs
   scripts/compliance/collect-debian-sources.py \
     BUNDLE/debian-inputs/binary-to-source.csv \
     BUNDLE/debian-inputs/installed-copyright.tar.gz \
     BUNDLE/source/debian
   ```

   Validation mapped 91 binary packages to 65 exact source packages and
   checksum-verified 209 source files. Current Bookworm indices supplied 64
   versions; Debian Snapshot supplied the unavailable patched `pcre2` version.
8. Preserve required copyright and license notices for every package regardless
   of whether corresponding source is required.

## Review classes

- `gpl-source-review`: GPL-family metadata without a listed alternative or
  exception. Normally collect corresponding source.
- `lgpl-source-review`: LGPL-family metadata. Review dynamic/static linkage and
  relinking requirements; collecting source is the conservative default.
- `alternative-license-review`: metadata contains `OR`; determine whether the
  distributed binary can be used under the non-copyleft alternative.
- `exception-review`: metadata names an exception, such as the GCC Runtime
  Library Exception or Classpath Exception. Preserve notices and verify that the
  distributed use satisfies the exception.

A missing feedstock commit in older package metadata does not prevent
collection: the extracted package evidence remains the authoritative exact
recipe snapshot.

## Exact package evidence

Conda packages commonly embed their rendered recipe, patches, build scripts,
and license files under `info/`. Use those files from the exact binary artifact
listed in `conda-source-urls.csv`; do not use the current recipe repository
branch as proof of how an older artifact was built.

A Conda binary URL is not corresponding source. The final external compliance
bundle must map each distributed copyleft binary to exact upstream source and
its recipe/patch material. `collect-conda-sources.py --conservative` was
validated against all 53 current runtime candidates: it produced 64 mappings to
48 unique checksum-verified source archives. Fifty artifacts record an exact
feedstock commit; the older `bioconductor-msa`, `glpk`, and `gsl` artifacts omit
that commit metadata but still embed their exact rendered recipes, build
scripts, patches, source URLs, and checksums. Debian packages require the
analogous mapping to the matching Debian source package and version.
