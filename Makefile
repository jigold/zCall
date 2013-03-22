# makefile for zCall by Iain Bancarz, ib5@sanger.ac.uk

PREFIX = 	PREFIX_DIRECTORY    # dummy value as default
DEST =		$(PREFIX)/zCall
SCRIPTS =	src/zcall
ETC =       src/etc

usage:
	@echo -e "Usage: make install PREFIX=<destination directory>\nWill install to the zCall subdirectory of PREFIX.\nPREFIX must exist, zCall subdirectory will be created if necessary."

install: $(PREFIX)
	@echo -e "Installing scripts..."
	install -d $(DEST) $(DEST)/zcall $(DEST)/etc $(DEST)/doc
	@rm -f $(DEST)/zcall/*.pyc    # force recompiling of any .pyc files
	install $(SCRIPTS)/*.py $(SCRIPTS)/*.r $(SCRIPTS)/findMeanSD $(SCRIPTS)/findThresholds $(DEST)/zcall
	install $(ETC)/*.ini $(DEST)/etc
	@echo -e "Writing documentation..."
	$(SCRIPTS)/createDocs.py --out  $(DEST)/doc --recursive
	@echo -e "Installation complete. See $(DEST)/doc/zcall.html for class documentation."
