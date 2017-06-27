.PHONY: MatrixFactorizationR_install MatrixFactorizationR_test MatrixFactorizationR_document MatrixFactorizationR_check

## install on krakenator
krakenator_install_MatrixFactorizationR:
	ssh -t cayek@krakenator.imag.fr "cd ~/Projects/Thesis/MatrixFactorizationR/; git pull; make MatrixFactorizationR_install"

## Rpackage
MatrixFactorizationR_install:
	R -e 'devtools::install(pkg = ".")'

MatrixFactorizationR_test:
	R -e 'devtools::test(pkg = ".")'

MatrixFactorizationR_document:
	R -e 'devtools::document(pkg = ".")'

MatrixFactorizationR_check:
	R -e 'devtools::check(pkg = ".")'

MatrixFactorizationR_clean:
	rm -f R/RcppExports.R src/RcppExports.cpp src/*.o src/*.so

