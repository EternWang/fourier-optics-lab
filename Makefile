.PHONY: all analysis test report clean

all: analysis

analysis:
	python analysis/analyze.py

test:
	python -m unittest discover -s tests

report:
	cd report && pdflatex -interaction=nonstopmode main.tex && pdflatex -interaction=nonstopmode main.tex

clean:
	find report -maxdepth 1 -type f \( -name "*.aux" -o -name "*.log" -o -name "*.out" -o -name "*.toc" -o -name "*.synctex.gz" \) -delete
