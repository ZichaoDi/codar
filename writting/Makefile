paper = report_060517
presentation = slide0907

# latex & bibtex each file

report_060517.pdf: $(paper).tex
	latex $(paper)
#	bibtex $(paper)
	latex $(paper)
	latex $(paper)
	dvipdf $(paper).dvi

slide0907.pdf: $(presentation).tex
	pdflatex $(presentation)

# create pictures, if needed 
.SUFFIXES: .fig .pdf .eps

clean:
	rm -f *~ *.lof *.fff *.aux *.bbl *.blg *.toc *.log *.bak *.out *.gz *.dvi *.snm *.nav slide0907.pdf

reallyclean: clean
	rm -f fig*.pdf fig*.eps
