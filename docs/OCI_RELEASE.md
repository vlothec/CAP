# OCI image release and GHCR publication

GitHub Container Registry (GHCR) is an OCI registry. CAP images do not need to
be published to Docker Hub. The public image repository is:

```text
ghcr.io/mellifluousjdg/cap
```

CAP 1.0.3 is published under immutable release tags and digests:

```text
OCI image:
ghcr.io/mellifluousjdg/cap:1.0.3
ghcr.io/mellifluousjdg/cap@sha256:2563036a264a3d7b72bba57215229df599af15a36ae0c4d7155b70a3127b6f95

Corresponding source:
ghcr.io/mellifluousjdg/cap-source:1.0.3
ghcr.io/mellifluousjdg/cap-source@sha256:33daf3e3b1b90011201c528d93b3bd86052258cf9a410bb1d1376c4a2a0699e0
```

The source archive is 3,206,034,233 bytes with SHA-256
`eeb0e951ca253dff3fc13569ac1a5e4945f8586b7681ee495f341c79a16cc57d`.

## Tag and digest policy

Published version tags are immutable. For each release, publish one full version
tag such as `1.1.0`; do not replace an existing version tag or rely on mutable
`latest`, major, or minor aliases. Documentation, release notes, and archived
validation records must also record the immutable manifest digest.

## Manual publication

CAP release images are currently published manually. The commands below
document the manual release procedure.
Build the release image from the exact clean release checkout using the image
build command documented in `docs/OCI_IMAGE.md`. Keep the canonical local name
`cap:local` for the test and publication commands below. OCI labels are useful
metadata but are not themselves a release gate; establish provenance from the
exact checkout, tested artifact, immutable tag, and recorded registry digest.

Before assigning or pushing the immutable release tag, run the integration test
inside the exact final local image:

```bash
podman run --rm --entrypoint /bin/bash \
  cap:local /opt/CAP/scripts/test-container-image.sh
```

Do not publish if this test fails. Podman and ORAS can use different default
credential stores. Authenticate with the client that will publish the artifact:

```bash
# For OCI image publication with Podman:
podman login ghcr.io --username MellifluousJDG

# For corresponding-source publication with ORAS:
oras login ghcr.io --username MellifluousJDG
```

Use a classic GitHub personal access token with `read:packages` and
`write:packages` when prompted for the password. Interactive prompting keeps the
token out of command-line arguments and shell history. After `oras login`,
invoke ORAS normally and let it use its default credential store. Pass
`--registry-config` only when intentionally supplying an alternate credential
file; doing so overrides the default ORAS credentials.

Tag and push the tested release image, replacing `<version>` with the new version:

```bash
version=<version>
podman tag cap:local "ghcr.io/mellifluousjdg/cap:$version"
podman push "ghcr.io/mellifluousjdg/cap:$version"
```

After pushing, display the registry descriptor and record its `digest` field.
Do not confuse a local image ID with the registry digest:

```bash
oras manifest fetch --descriptor \
  "ghcr.io/mellifluousjdg/cap:$version"
```

Set the GHCR package visibility deliberately in GitHub's package settings and
link the package to the CAP repository.

## Corresponding-source bundle

Generate the initial bundle from the final tested image as described in
`docs/OCI_IMAGE.md`, then complete the reviewed source and notice collection.
By default, generation verifies image revision labels against the checkout. For
an image already proven to come from an exact clean release checkout, explicit
expected revisions may be supplied instead:

```bash
CAP_EXPECTED_REVISION=$(git rev-parse HEAD) \
CAP_EXPECTED_TRASH2_REVISION=$(git rev-parse HEAD:modules/TRASH_2) \
make oci-compliance IMAGE="$image" OUTDIR="$bundle"
```

Both explicit revisions are required together and must match the checkout. The
following commands are intentionally explicit: generation requires network
access and can be lengthy, and each collector verifies downloaded checksums.

```bash
bundle=../cap-<version>-compliance-bundle
image=cap:local

scripts/compliance/apply-source-review-policy.py \
  "$bundle/source-obligations.csv" \
  "$bundle/source-obligations-reviewed.csv"
scripts/compliance/collect-conda-sources.py \
  "$bundle/source-obligations-reviewed.csv" \
  "$bundle/source/conda"
scripts/compliance/collect-conda-notices.py \
  "$bundle/inventory/conda-source-urls.csv" \
  "$bundle/notices/conda"
scripts/compliance/collect-conda-notice-sources.py \
  licenses/conda-notice-supplemental-sources.csv \
  "$bundle/source/conda-notice-supplemental"
```

Merge the two verified Conda corpora and their mappings into one canonical
corpus. The command consumes the two input directories after collection so the
same large binary and source archives are not duplicated in the release bundle:

```bash
scripts/compliance/merge-conda-source-corpora.py \
  "$bundle/source/conda" \
  "$bundle/source/conda-notice-supplemental" \
  "$bundle/source/conda-all"
```

Export the canonical R license texts from the exact image, resolve all reviewed
Conda notice gaps, and collect exact Debian sources:

```bash
mkdir "$bundle/notices/R-shared-licenses"
podman run --rm --entrypoint /bin/bash "$image" -c \
  'tar -C /opt/cap-env/lib/R/share/licenses -cf - .' | \
  tar -C "$bundle/notices/R-shared-licenses" -xf -
scripts/compliance/collect-conda-source-notices.py \
  "$bundle/notices/conda/conda-notice-manifest.csv" \
  "$bundle/source/conda-all/source-to-binary.csv" \
  "$bundle/source/conda-all" "$bundle/notices/conda-source" \
  --shared-license-dir "$bundle/notices/R-shared-licenses"
scripts/compliance/finalize-conda-source-notices.py \
  "$bundle/notices/conda-source/source-notice-manifest.csv" \
  licenses/conda-source-notice-review.csv \
  "$bundle/notices/R-shared-licenses" \
  "$bundle/notices/conda-source/final-review.csv"
scripts/compliance/export-debian-compliance-inputs.sh \
  "$image" "$bundle/debian-inputs"
scripts/compliance/collect-debian-sources.py \
  "$bundle/debian-inputs/binary-to-source.csv" \
  "$bundle/debian-inputs/installed-copyright.tar.gz" \
  "$bundle/source/debian"
scripts/compliance/finalize-compliance-checksums.sh "$bundle"
(
  cd "$bundle"
  sha256sum --check --strict SHA256SUMS
)
```

The existing reviewed-policy commands reject changed dependency sets and
unresolved notice decisions. The checksum finalizer removes obsolete nested
manifests, creates the single whole-bundle manifest, and the final command
verifies every recorded payload before archiving.

Create one deterministic Zstandard archive outside the checkout:

```bash
tar --sort=name --mtime='UTC 1970-01-01' --owner=0 --group=0 \
  --numeric-owner -C "$(dirname "$bundle")" \
  -cf - "$(basename "$bundle")" | \
  zstd -T0 -19 -o ../cap-<version>-compliance-bundle.tar.zst
sha256sum ../cap-<version>-compliance-bundle.tar.zst
```

## Publish the corresponding source artifact

Authenticate with `oras login` as described above, then check whether the
release source tag already exists:

```bash
version=<version>
oras manifest fetch "ghcr.io/mellifluousjdg/cap-source:$version"
```

Proceed only if GHCR reports that the tag was not found. If the command returns
a manifest or any other error, stop and investigate. Never replace an existing
release tag. Publish the single tested archive:

```bash
oras push \
  --artifact-type application/vnd.cap.compliance-bundle.v1 \
  "ghcr.io/mellifluousjdg/cap-source:$version" \
  "../cap-$version-compliance-bundle.tar.zst:application/zstd"
```

Display the published descriptor and record its `digest` field:

```bash
oras manifest fetch --descriptor \
  "ghcr.io/mellifluousjdg/cap-source:$version"
```

Record the returned manifest digest and archive SHA-256. Pull by digest into a
new directory, verify that the downloaded archive is byte-identical, and test
its complete checksum manifest before considering publication complete.

## Release artifacts

In addition to the registry image, preserve:

- the CAP Git tag and recursive submodule revision;
- the immutable GHCR digest;
- an immutable external compliance bundle from
  `scripts/compliance/generate-oci-compliance-bundle.sh`, including inventory files,
  revision records, BCT corresponding source, notices, and checksums;
- a CycloneDX 1.6 SBOM generated from the exact Conda, Debian, and
  embedded-component image inventories;
- exact corresponding-source archives and source-to-binary mappings required
  by bundled dependencies, following `docs/DEPENDENCY_SOURCE_REVIEW.md`;
- the validation log;
- optionally, a compressed OCI archive and SHA-256 checksum.

An OCI archive protects against registry loss but may be too large for GitHub's
normal release-asset limits. If so, store it in institutional archival storage
and publish its checksum and location in the release notes.

## Redistribution requirements

- Include CAP's root MIT license.
- Include the `bin/src/BCT` corresponding source and README, plus the GPL
  texts, with the compiled `ctw-calc`. The standalone BCT C++ repository has no
  separate explicit license notice; the related BCT R package published by the
  same author declares `GPL (>= 2)` in its CRAN metadata. On that basis, CAP
  applies `GPL-2.0-or-later` to the BCT-derived source and executable.
- Include TRASH_2's MIT-style `license.txt`. Its tracked Windows-only
  MAFFT/HMMER bundles and historical `temp/` outputs are excluded from the Linux
  image, which uses the locked Conda MAFFT/HMMER packages instead.
- Generate dependency inventories from the final image. If they differ from the
  reviewed candidate inventories, review the changed Conda and Debian packages
  and collect any additional notices or corresponding source before publication.
- Preserve the reviewed Conda notices and source mappings, plus the Debian
  copyright files, binary-to-source mappings, and corresponding source, in the
  external compliance bundle.
