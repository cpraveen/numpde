TEX = $(wildcard *.tex)
PDF = $(TEX:.tex=.pdf)
TARGET = $(PDF)

all: $(TARGET)

%.pdf: %.tex
	pdflatex $<
	pdflatex $<

clean:
	rm -rf $(PDF)
