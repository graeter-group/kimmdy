
.PHONY: Makefile setup-docs preview docs
	

preview:
	python -m quartodoc build --watch &
	python -m quartodoc interlinks
	quarto preview

docs:
	python -m quartodoc build 
	python -m quartodoc interlinks
	quarto render


setup-docs:
	quarto add --no-prompt machow/quartodoc
	python -m quartodoc build --verbose
	python -m quartodoc interlinks
	quarto render

