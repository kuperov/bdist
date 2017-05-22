
all: test package

test:
	Rscript --vanilla -e "devtools::test()"

package: test document
	Rscript --vanilla -e "devtools::build()"

document:
	Rscript --vanilla -e "devtools::document()"
