
.PHONY: Makefile quartodoc preview docs
	

preview:
	python -m quartodoc build 
	python -m quartodoc interlinks
	quarto preview

docs:
	python -m quartodoc build 
	python -m quartodoc interlinks
	quarto render


quartodoc:
	quarto add --no-prompt machow/quartodoc
	python -m quartodoc build --verbose
	python -m quartodoc interlinks
	quarto render

