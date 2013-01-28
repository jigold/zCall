# makefile for zCall by Iain Bancarz, ib5@sanger.ac.uk

PREFIX = 	.
DEST =		$(PREFIX)/zCall
SCRIPTS =	src/scripts

usage:
	@echo -e "Usage: make install PREFIX=<destination directory>\nWill install to the zCall subdirectory of PREFIX."

install:
	install -d $(DEST) $(DEST)/scripts 
	install $(SCRIPTS)/*.py $(SCRIPTS)/*.r $(DEST)/scripts
