.PHONY: all clean

all:
	latexmk -pdf -interaction=nonstopmode -halt-on-error talk.tex
	$(MAKE) clean

clean:
	latexmk -c -e '$$clean_ext .= " nav snm";' talk.tex
