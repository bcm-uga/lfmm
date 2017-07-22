.PHONY: MatrixFactorizationR_install MatrixFactorizationR_test MatrixFactorizationR_document MatrixFactorizationR_check

## krak
krakenator_deploy:
	git status
## git commit --allow-empty -am "deploy on krakenator"
	git push krakenator master

krakenator_push_hook:
	scp ./hooks/post-receive.sh cayek@krakenator:/home/cayek/Gits/2017/MatrixFactorizationR.git/hooks/post-receive

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

