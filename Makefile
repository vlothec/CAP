ENV_NAME ?= cap-pipeline
CONDA ?= conda
ASSEMBLY ?=
OUTDIR ?=
CORES ?= 1
MAX_REP_SIZE ?= 1000
MIN_REP_SIZE ?= 7
TEMPLATES ?=

.PHONY: install compile run run-trash2 lint test-ctw-build test-conda test-trash2-standalone test-container test-slurm-wrapper test-appimage-slurm-wrapper appimage test-appimage clean-appimage oci-inventory oci-compliance remove-env clean-build clean

install:
	@echo "Setting up Conda environment '$(ENV_NAME)'..."
	CAP_ENV_NAME="$(ENV_NAME)" bash setup_conda.sh

compile:
	CAP_ENV_NAME="$(ENV_NAME)" CAP_CONDA_EXE="$(CONDA)" scripts/build-ctw.sh

run:
	@test -n "$(ASSEMBLY)" || { \
		echo "Usage: make run ASSEMBLY=path/to/genome.fasta" >&2; \
		exit 1; \
	}
	conda run -n "$(ENV_NAME)" nextflow run . --assembly "$(ASSEMBLY)"

run-trash2:
	@test -n "$(ASSEMBLY)" && test -n "$(OUTDIR)" || { \
		echo "Usage: make run-trash2 ASSEMBLY=path/to/genome.fasta OUTDIR=path/to/output [ENV_NAME=cap-pipeline]" >&2; \
		exit 1; \
	}
	CAP_ENV_NAME="$(ENV_NAME)" CAP_CONDA_EXE="$(CONDA)" scripts/run-trash2.sh \
		--assembly "$(ASSEMBLY)" \
		--outdir "$(OUTDIR)" \
		--cores "$(CORES)" \
		--max-rep-size "$(MAX_REP_SIZE)" \
		--min-rep-size "$(MIN_REP_SIZE)" \
		$(if $(TEMPLATES),--templates "$(TEMPLATES)")

lint:
	bash scripts/lint.sh

test-ctw-build:
	bash scripts/test-ctw-build.sh

test-conda:
	CAP_ENV_NAME="$(ENV_NAME)" CAP_CONDA_EXE="$(CONDA)" bash scripts/test-conda.sh

test-trash2-standalone:
	CAP_TEST_CONDA="$(CONDA)" CAP_TEST_ENV_NAME="$(ENV_NAME)" bash scripts/test-trash2-standalone.sh

test-container:
	bash scripts/test-container.sh

test-slurm-wrapper:
	bash scripts/test-slurm-wrapper.sh

test-appimage-slurm-wrapper:
	bash scripts/test-appimage-slurm-wrapper.sh

appimage:
	bash scripts/appimage/build-appimage.sh

test-appimage:
	@test -n "$(APPIMAGE)" || { echo "Usage: make test-appimage APPIMAGE=dist/CAP-...-x86_64.AppImage" >&2; exit 1; }
	bash scripts/appimage/test-built-appimage.sh "$(APPIMAGE)"

clean-appimage:
	rm -rf build/appimage dist/*.AppImage dist/*.AppImage.sha256

oci-inventory:
	@test -n "$(IMAGE)" || { echo "Usage: make oci-inventory IMAGE=cap:local [OUTDIR=oci-inventory]" >&2; exit 1; }
	scripts/compliance/generate-oci-inventory.sh "$(IMAGE)" "$(or $(OUTDIR),oci-inventory)"

oci-compliance:
	@test -n "$(IMAGE)" || { echo "Usage: make oci-compliance IMAGE=cap:local [OUTDIR=../cap-compliance-bundle]" >&2; exit 1; }
	scripts/compliance/generate-oci-compliance-bundle.sh "$(IMAGE)" "$(or $(OUTDIR),../cap-compliance-bundle)"

remove-env:
	@echo "Removing Conda environment '$(ENV_NAME)'..."
	conda env remove -n "$(ENV_NAME)"

clean-build:
	@echo "Removing compiled CTW/BCT binary..."
	rm -f bin/ctw-calc

clean:
	@echo "ERROR: 'make clean' is intentionally disabled because cleanup is ambiguous." >&2
	@echo "Use 'make remove-env [ENV_NAME=name]' to remove a Conda environment." >&2
	@echo "Use 'make clean-build' to remove the compiled CTW/BCT binary." >&2
	@echo "Use 'make clean-appimage' to remove AppImage build output and artifacts." >&2
	@exit 2
