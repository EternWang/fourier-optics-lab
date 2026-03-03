.PHONY: analysis report clean

analysis:
	python analysis/analyze.py

report:
	cd report && pdflatex -interaction=nonstopmode main.tex && pdflatex -interaction=nonstopmode main.tex

clean:
	find report -maxdepth 1 -type f \( -name "*.aux" -o -name "*.log" -o -name "*.out" -o -name "*.toc" -o -name "*.synctex.gz" \) -delete
