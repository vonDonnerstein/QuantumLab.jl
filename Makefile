srcfiles := $(shell find src/ -type f -name '*.jl')
cov_files_ifexist := $(wildcard src/QuantumLab.jl.*.cov)
ifdef cov_files_ifexist
  cov_files=$(cov_files_ifexist)
else
  cov_files=src/QuantumLab.jl.cov
endif






default:
	@echo "There is no sensible default target (yet). Look into Makefile to see your options."

# MAIN TARGETS
libint2jl:
	cd deps/usr/src && $(MAKE) 

libint2jl-clean:
	cd deps/usr/src && $(MAKE) clean

tags:
	ctags -R

coverage: 
	$(MAKE) clean
	$(MAKE) index.html
	xdg-open index.html

show_coverage: index.html
	xdg-open index.html

full_coverage:
	rm -f test/STO-3G.tx93
	rm -rf test/hackbusch_pretables/
	$(MAKE) coverage

clean: 
	find . -name *.cov -exec rm {} \;
	find . -name index.html -exec rm {} \;
	find . -name index-sort-f.html -exec rm {} \;
	find . -name index-sort-l.html -exec rm {} \;
	find . -name *.jl.*.html -exec rm {} \;
	rm -f lcov.info





# partial generators (helpers)
$(cov_files): 
	find . -name *.cov -exec rm {} \;
	@echo "INFO: The next step may take a while:"
	cd test; julia --code-coverage=user --inline=no runtests.jl &>coverage_run.out

lcov.info: $(cov_files)
	julia -e 'using Coverage; LCOV.writefile("lcov.info",process_folder())'

index.html: lcov.info
	@genhtml -v || (echo "genhtml seems not to be available. Try installing lcov or adding it to PATH and rerun make." && exit 1)
	genhtml lcov.info

.PHONY: default clean coverage show_coverage
