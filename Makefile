PYTHON ?= python3
JUPYTER ?= jupyter
CFFCONVERT ?= cffconvert
IPYTHON_DIR ?= outputs/ipython
JUPYTER_RUNTIME ?= outputs/jupyter-runtime
PYTHON_JUPYTER_PATH ?= $(shell $(PYTHON) -c 'import sys; print(sys.prefix + "/share/jupyter")')
INPUT_ROOT ?= data/private
OUTPUT_ROOT ?= outputs/stata
ANALYSIS_ROOT ?= .
STATA_BIN ?=
STATA_MODE ?= auto
RUN_ID ?=
EXPECTED_INPUT_SHA256 ?=
BASELINE_RUN ?=
CANDIDATE_RUN_1 ?=
CANDIDATE_RUN_2 ?=
COMPARISON_REPORT ?= outputs/validation/comparison_report.json

STATA_RUN_ARGS = --input-root "$(INPUT_ROOT)" --output-root "$(OUTPUT_ROOT)" --analysis-root "$(ANALYSIS_ROOT)" --stata-mode "$(STATA_MODE)"
ifneq ($(strip $(STATA_BIN)),)
STATA_RUN_ARGS += --stata-bin "$(STATA_BIN)"
endif
ifneq ($(strip $(RUN_ID)),)
STATA_RUN_ARGS += --run-id "$(RUN_ID)"
endif
ifneq ($(strip $(EXPECTED_INPUT_SHA256)),)
STATA_RUN_ARGS += --expected-input-sha256 "$(EXPECTED_INPUT_SHA256)"
endif

.PHONY: help check diagram-smoke stata-run stata-compare

help:
	@echo "Public, data-free targets:"
	@echo "  make check          Validate metadata, safety rules, and tests"
	@echo "  make diagram-smoke  Execute the diagram notebook into ignored outputs"
	@echo "Restricted-data targets:"
	@echo "  make stata-run      Run guarded Stata analysis into a unique folder"
	@echo "  make stata-compare  Compare one baseline and two candidate runs"

check:
	$(PYTHON) -m unittest discover -s tests -p 'test_*.py'
	$(PYTHON) scripts/check_public_surface.py
	$(CFFCONVERT) --validate
	git diff --check

diagram-smoke:
	@command -v dot >/dev/null || { echo "Graphviz dot is required." >&2; exit 1; }
	dot -V
	mkdir -p outputs/notebooks outputs/figures $(IPYTHON_DIR) $(JUPYTER_RUNTIME)
	IPYTHONDIR="$(IPYTHON_DIR)" JUPYTER_PATH="$(PYTHON_JUPYTER_PATH)" JUPYTER_RUNTIME_DIR="$(JUPYTER_RUNTIME)" $(JUPYTER) nbconvert \
		--to notebook \
		--execute "Case Definitions Consort.ipynb" \
		--output "Case Definitions Consort.executed.ipynb" \
		--output-dir outputs/notebooks \
		--ExecutePreprocessor.timeout=120
	test -s outputs/figures/consort_diagram.tiff
	git diff --exit-code -- "Case Definitions Consort.ipynb"

stata-run:
	PYTHON="$(PYTHON)" scripts/run_stata.sh $(STATA_RUN_ARGS)

stata-compare:
	@test -n "$(BASELINE_RUN)" || { echo "BASELINE_RUN is required." >&2; exit 2; }
	@test -n "$(CANDIDATE_RUN_1)" || { echo "CANDIDATE_RUN_1 is required." >&2; exit 2; }
	@test -n "$(CANDIDATE_RUN_2)" || { echo "CANDIDATE_RUN_2 is required." >&2; exit 2; }
	$(PYTHON) scripts/compare_stata_runs.py \
		--baseline-run "$(BASELINE_RUN)" \
		--candidate-run-1 "$(CANDIDATE_RUN_1)" \
		--candidate-run-2 "$(CANDIDATE_RUN_2)" \
		--input-file "$(INPUT_ROOT)/full_db.dta" \
		$(if $(strip $(EXPECTED_INPUT_SHA256)),--expected-input-sha256 "$(EXPECTED_INPUT_SHA256)",) \
		--report "$(COMPARISON_REPORT)"
