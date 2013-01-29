# makefile for zCall by Iain Bancarz, ib5@sanger.ac.uk

PREFIX = 	PREFIX_DIRECTORY    # dummy value as default
DEST =		$(PREFIX)/zCall
SCRIPTS =	src/scripts
ETC =       src/etc

usage:
	@echo -e "Usage: make install PREFIX=<destination directory>\nWill install to the zCall subdirectory of PREFIX."

install: $(PREFIX)
	install -d $(DEST) $(DEST)/scripts $(DEST)/etc
	rm -f $(DEST)/scripts/*.pyc    # force recompiling of any .pyc files
	install $(SCRIPTS)/*.py $(SCRIPTS)/*.r $(DEST)/scripts
	install $(ETC)/*.ini $(DEST)/etc
