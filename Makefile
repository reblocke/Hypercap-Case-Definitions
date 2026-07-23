PYTHON ?= python3
JUPYTER ?= jupyter
CFFCONVERT ?= cffconvert
IPYTHON_DIR ?= outputs/ipython
JUPYTER_RUNTIME ?= outputs/jupyter-runtime
PYTHON_JUPYTER_PATH ?= $(shell $(PYTHON) -c 'import sys; print(sys.prefix + "/share/jupyter")')

.PHONY: help check diagram-smoke

help:
	@echo "Public, data-free targets:"
	@echo "  make check          Validate metadata, safety rules, and tests"
	@echo "  make diagram-smoke  Execute the diagram notebook into ignored outputs"

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
