# One command from pinned data to REPORT.md and CONTROLS.md:
#   make reproduce
# Requires: rustup (toolchain pinned in rust-toolchain.toml) and uv.

QG      := target/release/qg
UV      := uv run
ORDERS  := 1 2 3 4 5 6
DER     := data/derived

.PHONY: reproduce build rust-test runs analysis residual beyond report py-test clean-derived paper oeis-fingerprints independent-checks

reproduce: build rust-test runs analysis residual beyond report py-test

build:
	cargo build --release --locked

rust-test:
	cargo test --release --locked

runs: build
	mkdir -p $(DER) data/floor
	$(QG) laws > $(DER)/laws_survey.tsv
	$(QG) parastrophe $(DER)/parastrophe.tsv > $(DER)/parastrophe_summary.tsv
	: > $(DER)/enumerate.tsv
	for n in $(ORDERS); do $(QG) enumerate $$n data/floor/order$$n.tsv.gz >> $(DER)/enumerate.tsv; done
	: > $(DER)/count_latin.tsv
	for n in $(ORDERS); do $(QG) count-latin $$n >> $(DER)/count_latin.tsv; done
	$(QG) g4 5 > $(DER)/g4.tsv
	$(QG) controls $(foreach n,$(ORDERS),data/floor/order$(n).tsv.gz) > $(DER)/controls.tsv
	{ printf 'cpu\t%s\n' "$$(sysctl -n machdep.cpu.brand_string 2>/dev/null || uname -m)"; \
	  printf 'logical_cores\t%s\n' "$$(sysctl -n hw.ncpu 2>/dev/null || nproc)"; \
	  printf 'memory_bytes\t%s\n' "$$(sysctl -n hw.memsize 2>/dev/null || echo unknown)"; \
	  printf 'rustc\t%s\n' "$$(rustc --version)"; \
	  printf 'python\t%s\n' "$$($(UV) python --version)"; } > $(DER)/machine.tsv

analysis:
	$(UV) python analysis/compute.py

# Needs mace4 on PATH (not exhaustive; see REPORT.md). Skips cleanly without it.
residual: build
	$(UV) python analysis/residual.py

# Exhaustive models, at order k = 7..9, of every variety still unrealised at order k-1.
# No timeout: a result is exact only if its run finishes. Order 9 takes hours.
beyond: build analysis
	mkdir -p data/beyond
	rm -f data/beyond/models_*.tsv.gz data/beyond/open_*.txt
	: > data/beyond/runs.tsv
	set -e; for n in 7 8 9; do \
	  $(UV) python analysis/beyond.py --open-before $$n > data/beyond/open_$$n.txt; \
	  for g in $$(cat data/beyond/open_$$n.txt); do \
	  f=data/beyond/models_o$${n}_$$(echo $$g | tr , -).tsv.gz; \
	  if ! $(QG) models $$n $$g $$f >> data/beyond/runs.tsv; then \
	    printf 'order\t%s\tlaws\t%s\tstatus\terror\n' $$n $$g >> data/beyond/runs.tsv; rm -f $$f; fi; \
	done; done
	$(UV) python analysis/beyond.py

# Needs network; pins OEIS matches for class model counts. Not part of reproduce.
oeis-fingerprints:
	$(UV) python analysis/oeis_fingerprints.py

report:
	$(UV) python analysis/report.py

# Independent re-checks: own evaluator, raw ground truth, SAT via kissat. Not part of reproduce.
independent-checks:
	$(UV) python analysis/independent_checks.py --cover --models --sat

# Needs TeX (latexmk). Not part of reproduce.
paper:
	$(MAKE) -C paper clean all

py-test:
	$(UV) pytest -q tests/test_gates.py

clean-derived:
	rm -rf $(DER) data/floor data/signatures data/cover
