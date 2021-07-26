PDF = p.pdf pa.pdf unbiasing.pdf

all: $(PDF)

%.pdf : %.tex
	grep ^%SCRIPT $< | cut -c9- | sh
	cat $< | ./pdflatexfilter > $@

clean :
	$(RM) $(PDF) *.png
