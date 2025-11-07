# Dockerfile
FROM mambaorg/micromamba:1.5-alpine

# --- Copy and install Conda env ---
COPY environment.yml /tmp/environment.yml
RUN micromamba create -y -f /tmp/environment.yml && \
    micromamba clean --all --yes

# --- Activate env ---
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# --- Switch to root ---
USER root

# --- Install system packages ---
RUN apk add --no-cache git bash curl openjdk17

# --- Clone TRASH_2 ---
RUN git clone https://github.com/vlothec/TRASH_2.git /opt/TRASH_2 && \
    chmod +x /opt/TRASH_2/src/TRASH.R

# --- Install Nextflow ---
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# --- Set PATH ---
ENV PATH="/opt/TRASH_2/src:$PATH"

# --- Copy pipeline files INTO /work ---
COPY bin/ /usr/local/bin/
COPY model/ /model/
COPY main.nf nextflow.config /work/

WORKDIR /work

CMD ["nextflow", "run", "main.nf", "--help"]