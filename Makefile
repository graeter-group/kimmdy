
.PHONY: Makefile quartodoc preview

preview:
	python -m quartodoc build 
	quarto preview

quartodoc:
	quarto add --no-prompt machow/quartodoc
	python -m quartodoc build --verbose
	python -m quartodoc interlinks
	quarto render

