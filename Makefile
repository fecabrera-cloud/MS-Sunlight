.PHONY: all validate

all:
	bash scripts/run_all.sh

validate:
	python scripts/validate_release.py
