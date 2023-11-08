
.PHONY: Makefile setup-docs preview docs clear-docs watch

preview:
	quarto preview

watch:
	python -m quartodoc build --watch

docs:
	python -m quartodoc build 
	python -m quartodoc interlinks
	quarto render

setup-docs: clear-docs
	quarto add --no-prompt machow/quartodoc
	python -m quartodoc build --verbose
	python -m quartodoc interlinks
	quarto render

clear-docs:
	rm -rf docs
	rm -rf .quarto
	rm -rf _inv
	rm -rf _reference
	rm -f ./objects.json
	rm -f _sidebar.yml

format:
	black src tests
