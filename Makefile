.PHONY: lfmm_testonR_install lfmm_test lfmm_document lfmm_check

## krak
krakenator_deploy:
	git status
## git commit --allow-empty -am "deploy on krakenator"
	git push krakenator master

krakenator_push_hook:
	scp ./hooks/post-receive.sh cayek@krakenator:/home/cayek/Gits/2017/lfmm.git/hooks/post-receive

## Rpackage
lfmm_install:
	R -e 'devtools::install(pkg = ".")'

lfmm_test:
	R -e 'devtools::test(pkg = ".")'

lfmm_document:
	R -e 'devtools::document(pkg = ".")'

lfmm_check:
	R -e 'devtools::check(pkg = ".")'

lfmm_clean:
	rm -f R/RcppExports.R src/RcppExports.cpp src/*.o src/*.so

