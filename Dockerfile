# syntax=docker/dockerfile:1

ARG DEBIAN_IMAGE=docker.io/library/debian@sha256:3783cc01769c7b2b1b83a5c5ad96c815348e28ed7da68e2e3687004faa906251
FROM ${DEBIAN_IMAGE} AS env-builder

ARG MINICONDA_VERSION=py313_25.7.0-2
ARG MINICONDA_SHA256=dda3629462ba1cfa72eb74535214c2e315c77f1cfb0f02046537e99f1bf64abc
ARG TARGETARCH

ENV DEBIAN_FRONTEND=noninteractive \
    PATH=/opt/conda/bin:/opt/cap-env/bin:${PATH}

RUN apt-get update \
    && apt-get install --yes --no-install-recommends ca-certificates curl bzip2 \
    && rm -rf /var/lib/apt/lists/*

RUN test "${TARGETARCH}" = "amd64" \
    && test -n "${MINICONDA_SHA256}" \
    && curl --fail --location --retry 3 \
        "https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh" \
        --output /tmp/miniconda.sh \
    && echo "${MINICONDA_SHA256}  /tmp/miniconda.sh" | sha256sum --check --strict \
    && bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh

COPY conda-runtime-linux-64.lock /tmp/conda-runtime-linux-64.lock
RUN conda create --yes --prefix /opt/cap-env \
        --override-channels \
        --channel conda-forge \
        --channel bioconda \
        --file /tmp/conda-runtime-linux-64.lock \
    && conda clean --all --yes \
    && rm /tmp/conda-runtime-linux-64.lock

FROM ${DEBIAN_IMAGE} AS ctw-builder
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get install --yes --no-install-recommends g++ make \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/CAP
COPY bin/src/BCT bin/src/BCT
RUN make -C bin/src/BCT \
    && test -x bin/ctw-calc

FROM ${DEBIAN_IMAGE} AS runtime

ARG CAP_REVISION=UNSPECIFIED
ARG TRASH_2_REVISION=UNSPECIFIED
LABEL org.opencontainers.image.title="CAP" \
      org.opencontainers.image.description="Centromere Analysis Pipeline with its locked Conda environment" \
      org.opencontainers.image.licenses="MIT AND GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/MellifluousJDG/CAP" \
      org.opencontainers.image.revision="${CAP_REVISION}" \
      org.cap-project.trash-2.revision="${TRASH_2_REVISION}"

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH=/opt/cap-env/bin:${PATH} \
    NXF_HOME=/home/cap/.nextflow \
    NXF_ANSI_LOG=false

RUN apt-get update \
    && apt-get install --yes --no-install-recommends bash ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd --gid 1000 cap \
    && useradd --uid 1000 --gid cap --create-home --shell /bin/bash cap \
    && mkdir -p /data /results /work \
    && chown cap:cap /data /results /work

COPY --from=env-builder /opt/cap-env /opt/cap-env
COPY --chown=cap:cap . /opt/CAP
COPY --from=ctw-builder --chown=cap:cap /opt/CAP/bin/ctw-calc /opt/CAP/bin/ctw-calc
COPY --chown=cap:cap scripts/cap-container-entrypoint.sh /usr/local/bin/cap

RUN test -x /opt/CAP/bin/ctw-calc \
    && test -f /opt/CAP/modules/TRASH_2/src/TRASH.R \
    && test -f /opt/CAP/modules/TRASH_2/dep/HOR.V3.3 \
    && test ! -e /opt/CAP/modules/TRASH_2/dep/hmmer \
    && test ! -e /opt/CAP/modules/TRASH_2/dep/mafft-7.520-win64-signed \
    && test ! -e /opt/CAP/modules/TRASH_2/dep/new.hor \
    && test ! -e /opt/CAP/modules/TRASH_2/temp \
    && chmod 0555 /usr/local/bin/cap

USER cap
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/cap"]
CMD ["--help"]
