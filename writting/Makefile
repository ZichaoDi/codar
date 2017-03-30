paper = abstract_0327

# latex & bibtex each file
abstract_0327.pdf: $(paper).tex
	pdflatex $(paper)
# 	bibtex $(paper)
# 	latex $(paper)
# 	latex $(paper)
# 	dvipdf $(paper).dvi


# create pictures, if needed 
.SUFFIXES: .fig .pdf .eps

clean:
	rm -f *~ *.lof *.fff *.aux *.bbl *.blg *.toc *.log *.bak *.out *.gz *.dvi

reallyclean: clean
	rm -f fig*.pdf fig*.eps
