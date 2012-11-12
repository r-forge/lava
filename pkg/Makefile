PKG := lava lava.tobit lava.mixture lava.nlin mets
R := R
RCMD := $(R) CMD
RSCRIPT := Rscript
SVNCMD := svn
getver := ./getversion.pl
package :=  $(foreach var,$(PKG),$(shell $(getver) $(var)))
version :=  $(foreach var,$(PKG),$(shell $(getver) $(var) version))
cranpkg := $(foreach var,$(PKG), CRAN/$(var))

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

export:
	rm -Rf $(cranpkg)
	$(SVNCMD) export $(PKG) $(cranpkg)
	cd CRAN; $(RCMD) build $(PKG)

exportcheck: export 
	cd CRAN; $(MAKE) fcheck $(PKG)

srcbuild: 
	$(RCMD) build --library=$(libpath) $(PKG)

wincheck: 
	$(RCMD) build $(PKG)
	ftp -a -u ftp://win-builder.r-project.org/R-release/ $(package)

build: 
	$(RCMD) INSTALL --build --library=$(libpath) $(PKG)

roxygen: 
	echo 'library(roxygen2); roxygenize("$(PKG)")' | $(R) --library=$(libpath) --no-save 

##buildwin: 
%.zip:
	cd $(libpath); zip -r $@ $*/* && mv $@ $(curpath)

install: 
	$(RCMD) INSTALL --no-multiarch $(PKG)

R:
	$(RCMD) SHLIB m.cpp

qcheck: 
	$(RCMD) check $(PKG) --no-tests --no-install --library=$(libpath)

check: 
	$(RCMD) check $(PKG) --timings --no-multiarch --library=$(libpath)

fcheck: 
	$(RCMD) check $(PKG) --timings --as-cran --no-multiarch --library=$(libpath)

cranprep: 
	rm -Rf ~/temp/cran/$(PKG)
	mkdir -p ~/temp/cran
	svn export $(PKG) ~/temp/cran/$(PKG)
	cd ~/temp/cran; R CMD build $(PKG)
	cd ~/temp/cran; R CMD check $(package) --as-cran --library=$(libpath)
	@echo "Upload to:      anonymous@cran.r-project.org"
	@echo "Prepare email:"
	@echo "   To: CRAN@R-project.org"
	@echo "   Subject: CRAN submission $(package))"
	@echo "   Dear CRAN administrators,"
	@echo ""
	@echo "   An updated version of the package $(PKG) has been uploaded to the incoming directory."
	@echo "   Thanks for putting my package on CRAN."
	@echo "   Sincerely,"

srcclean: 
	@rm -Rf $(PKG)/src-i386
	@rm -Rf $(PKG)/src-x86_64
	@rm -f $(PKG)/src/*.o
	@rm -f $(PKG)/src/*.so

clean:
	@rm -Rf $(chkpkg)
	@rm -f *~ *.vrb *.ilg *idx *.lof *.brf *.ind  weaver_debug* _region_* *.pdfsync *.aux *.bbl *.out *.nav *.snm *.toc *aux *.bak *.blg *dvi *.log *tmp

.PHONY: build check R index clean
