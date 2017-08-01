
all: test package

test:
	Rscript --vanilla -e "devtools::test()"

package: test document
	R CMD build .
	ls bdist_*.tar.gz | xargs -n1 R CMD check

document:
	Rscript --vanilla -e "devtools::document()"
	Rscript --vanilla -e "rmarkdown::render('README.Rmd')"
