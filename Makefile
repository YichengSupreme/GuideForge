.PHONY: help install run test clean dev-check

help:  ## Show this help message
	@echo "GuideForge - Available Commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install:  ## Install dependencies
	pip install -r requirements.txt

run:  ## Run full CRISPR pipeline (UCSC → PAM → QC)
	python run_pipeline.py targets.txt --scan-pam --qc

test:  ## Run a small integration test (no IDT scoring)
	@echo "chr17:7668402-7668521:+" > test_targets.txt
	python utils/get_ucsc_sequences.py test_targets.txt --scan-pam --qc
	@echo "✅ Test completed successfully."
	@rm -f test_targets.txt CRISPR_candidates*.csv Upstream_sequences.txt Downstream_sequences.txt

clean:  ## Remove generated files (keeps manifests for reproducibility)
	rm -f CRISPR_candidates*.csv CRISPR_candidates.txt Upstream_sequences.txt Downstream_sequences.txt top_CRISPR_candidates*.csv *.log

dev-check:  ## Optional developer checks (lint, YAML validation)
	@echo "🔍 Running developer checks..."
	@echo "📝 Checking YAML syntax..."
	@python -c "import yaml; yaml.safe_load(open('config.yaml')); yaml.safe_load(open('policy.yaml')); print('✅ YAML syntax valid')"
	@echo "🔍 Checking Python syntax..."
	@python -m py_compile utils/*.py run_pipeline.py test_pipeline.py 2>/dev/null && echo "✅ Python syntax valid" || echo "⚠️  Python syntax issues found"
	@echo "✅ Developer checks completed"
