PKG := lava lava.tobit lava.mixture twinlm lava.nlin mets
RCMD := R
getver := ./getversion.pl
package :=  $(foreach var,$(PKG),$(shell $(getver) $(var)))
version :=  $(foreach var,$(PKG),$(shell $(getver) $(var) version))
winpkg :=  $(foreach var,$(PKG),$(var).zip)
chkpkg :=  $(foreach var,$(PKG),$(var).Rcheck)
##$(shell perl getversion.pl $(PKG))
##version := $(shell perl getversion.pl $(PKG) version)
libpath := $(shell echo "cat(.libPaths()[1],\"\n\",sep=\"\")" | R --no-save --slave | tail -n 1)
curpath := $(shell pwd)

default: install

win: $(winpkg)

index:
	@echo "Package: $(PKG)" > current/PACKAGES
	@echo "Version: $(version)" >> current/PACKAGES
	@echo "Depends: graph" >> current/PACKAGES
	@echo "" >> current/PACKAGES

all: build install buildwin
##	cp $(package) $(PKG).zip versions/
##	@rm -f current/*.zip current/*.tar.gz
##	mv $(package) $(PKG).zip current/
##	@cp current/$(package) current/$(PKG).tar.gz 
####	@cp -f current/$(PKG).zip current/$(PKG)_$(version).zip
	@ln -f -s $(curpath)/current/* ~/public_html/sem/

srcbuild: $(RCMD) CMD build --library=$(libpath) $(PKG)

build: 
	$(RCMD) CMD INSTALL --build --library=$(libpath) $(PKG)

##buildwin: 
%.zip:
	cd $(libpath); zip -r $@ $*/* && mv $@ $(curpath)

install: 
	$(RCMD) CMD INSTALL --no-multiarch $(PKG)

R:
	$(RCMD) CMD SHLIB m.cpp

qcheck: 
	$(RCMD) CMD check $(PKG) --no-tests --no-install --library=$(libpath)

check: 
	$(RCMD) CMD check $(PKG) --no-multiarch --library=$(libpath)

fcheck: 
	$(RCMD) CMD check $(PKG) --as-cran --no-multiarch --library=$(libpath)

clean:
	@rm -Rf $(chkpkg)
	@rm -f *~ *.vrb *.ilg *idx *.lof *.brf *.ind  weaver_debug* _region_* *.pdfsync *.aux *.bbl *.out *.nav *.snm *.toc *aux *.bak *.blg *dvi *.log *tmp

.PHONY: build check R index clean
