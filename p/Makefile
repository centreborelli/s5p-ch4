PDF = p.pdf

all: $(PDF)

%.pdf : %.tex
	grep ^%SCRIPT $< | cut -c9- | sh
	cat $< | ./pdflatexfilter > $@

clean :
	$(RM) $(PDF) *.png
