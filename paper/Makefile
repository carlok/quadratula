# Build the draft. Run `make reproduce` at the repository root first.
DOC := quadratula

all: $(DOC).pdf

numbers.tex:
	python3 ../analysis/tex_numbers.py numbers.tex

$(DOC).pdf: $(DOC).tex numbers.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error $(DOC).tex

clean:
	latexmk -C $(DOC).tex
	rm -f numbers.tex

.PHONY: all clean numbers.tex
