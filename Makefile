# makefile for zCall by Iain Bancarz, ib5@sanger.ac.uk

PREFIX = 	PREFIX_DIRECTORY    # dummy value as default
DEST =		$(PREFIX)/zCall
SCRIPTS =	src/zcall
CREATE_DOCS = src/create_docs
ETC =       src/etc

usage:
	@echo -e "Usage: make install PREFIX=<destination directory>\nWill install to the zCall subdirectory of PREFIX."

install: $(PREFIX)
	@echo -e "Installing scripts..."
	install -d $(DEST) $(DEST)/zcall $(DEST)/etc $(DEST)/doc
	@rm -f $(DEST)/zcall/*.pyc    # force recompiling of any .pyc files
	install $(SCRIPTS)/*.py $(SCRIPTS)/*.r $(DEST)/zcall
	install $(ETC)/*.ini $(DEST)/etc
	@echo -e "Writing documentation..."
	$(CREATE_DOCS)/createDocs.py --recursive
	install $(CREATE_DOCS)/*.html $(DEST)/doc
	@echo -e "Installation complete. See $(DEST)/doc/zcall.html for class documentation."
