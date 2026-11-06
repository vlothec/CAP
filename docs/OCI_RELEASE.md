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

