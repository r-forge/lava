##pkg := lava lava.tobit lava.mixture lava.nlin mets
pkg := mets
pkgs := lava lava.tobit lava.mixture lava.nlin mets
R := R
RCMD := $(R) CMD
RSCRIPT := Rscript
GITCMD := git
getver := ./getversion.pl
package :=  $(foreach var,$(pkg),$(shell $(getver) $(var)))
version :=  $(foreach var,$(pkg),$(shell $(getver) $(var) version))
cranpkg := $(foreach var,$(pkg), CRAN/$(var))
ifeq ($(wildcard ~/tmp),) 
    tempdir = /tmp
else 
    tempdir = ~/tmp
endif


winpkg :=  $(foreach var,$(pkg),$(var).zip)
chkpkg :=  $(foreach var,$(pkg),$(var).Rcheck)
##$(shell perl getversion.pl $(pkg))
##version := $(shell perl getversion.pl $(pkg) version)
libpath := $(shell echo "cat(.libPaths()[1],\"\n\",sep=\"\")" | R --no-save --slave | tail -n 1)
curpath := $(shell pwd)

default: install

win: $(winpkg)

tmp:
	echo $(tempdir)

index:
	@echo "Package: $(pkg)" > current/PACKAGES
	@echo "Version: $(version)" >> current/PACKAGES
	@echo "Depends: graph" >> current/PACKAGES
	@echo "" >> current/PACKAGES

all: build install buildwin
##	cp $(package) $(pkg).zip versions/
##	@rm -f current/*.zip current/*.tar.gz
##	mv $(package) $(pkg).zip current/
##	@cp current/$(package) current/$(pkg).tar.gz 
####	@cp -f current/$(pkg).zip current/$(pkg)_$(version).zip
	@ln -f -s $(curpath)/current/* ~/public_html/sem/

export:
	@rm -Rf $(cranpkg)
	@mkdir -p $(cranpkg)
	$(GITCMD) archive master | (cd $(cranpkg); tar x)
	cd CRAN; $(RCMD) build $(pkg)

exportcheck: export 
	cd CRAN; $(MAKE) fcheck $(pkg)

srcbuild: 
	$(RCMD) build --library=$(libpath) $(pkg)

wincheck: 
	$(RCMD) build $(pkg)
	ftp -a -u ftp://win-builder.r-project.org/R-release/ $(package)

build: 
	$(RCMD) INSTALL --build --library=$(libpath) $(pkg)

roxygen: roxy
roxy: 
	echo 'library(roxygen2); roxygenize("$(pkg)")' | $(R) --library=$(libpath) --no-save 

##buildwin: 
%.zip:
	cd $(libpath); zip -r $@ $*/* && mv $@ $(curpath)


q:
	$(RCMD) INSTALL --no-multiarch --no-data --no-demo --no-help  --configure-args="CXXFLAGS=-O0" $(pkg)

install: in
in: 
	$(RCMD) INSTALL --no-multiarch $(pkg)

rin: roxy in 

R:
	$(RCMD) SHLIB m.cpp

qcheck: 
	$(RCMD) check $(pkg) --no-tests --no-install --library=$(libpath)

check: 
	$(RCMD) check $(pkg) --timings --no-multiarch --library=$(libpath)

fcheck: 
	$(RCMD) check $(pkg) --timings --as-cran --no-multiarch --library=$(libpath)

cranprep: 
	rm -Rf $(tempdir)/cran/$(pkg)
	mkdir -p $(tempdir)/cran
	svn export $(pkg) $(tempdir)/cran/$(pkg)
	cd $(tempdir)/cran; R CMD build $(pkg)
	cd $(tempdir)/cran; R CMD check $(package) --as-cran --library=$(libpath)
	@echo "Upload to:      anonymous@cran.r-project.org"
	@echo "Prepare email:"
	@echo "   To: CRAN@R-project.org"
	@echo "   Subject: CRAN submission $(package))"
	@echo "   Dear CRAN administrators,"
	@echo ""
	@echo "   An updated version of the package $(pkg) has been uploaded to the incoming directory."
	@echo "   Thanks for putting my package on CRAN."
	@echo "   Sincerely,"

srcclean: cc
cc: 
	@rm -Rf $(pkg)/src-i386
	@rm -Rf $(pkg)/src-x86_64
	@rm -f $(pkg)/src/*.o
	@rm -f $(pkg)/src/*.so

clean:
	@rm -Rf $(chkpkg)
	@rm -f *~ *.vrb *.ilg *idx *.lof *.brf *.ind  weaver_debug* _region_* *.pdfsync *.aux *.bbl *.out *.nav *.snm *.toc *aux *.bak *.blg *dvi *.log *tmp

.PHONY: build check R index clean
