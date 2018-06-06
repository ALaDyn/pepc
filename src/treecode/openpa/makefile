include $(ROOTDIR)/makefile.defs
include $(ROOTDIR)/tools/build/makefile.paths

OPA_DIR       = openpa
OPA_MAKEFILE  = $(OPA_DIR)/Makefile
COPY_FILES    = $(addprefix $(OPA_DIR)/src/, opa_config.h opa_primitives.h \
		opa_queue.h opa_util.h)
LIBTOOL_FILES = $(OPA_DIR)/src/libopa.la

default: makeopa
	@cp $(COPY_FILES) $(LIBDIR)
	@mkdir -p $(LIBDIR)/primitives
	@cp $(OPA_DIR)/src/primitives/* $(LIBDIR)/primitives
	@$(OPA_DIR)/libtool --mode=install cp $(LIBTOOL_FILES) $(LIBDIR)

makeopa: $(OPA_MAKEFILE)
	@echo "     calling make for openpa"
	@make -C $(OPA_DIR)

clean: $(OPA_MAKEFILE)
	@echo "     calling make clean for openpa"
	@make -C $(OPA_DIR) clean
	@$(RM) $(OPA_MAKEFILE)

$(OPA_MAKEFILE):
	@echo "     calling ./configure for openpa"
	@cd $(OPA_DIR) && ./configure -q --disable-shared --enable-static $(OPENPA_CONFIGFLAGS) AR=$(AR) RANLIB=$(RANLIB) CC=$(CC) CFLAGS="$(CFLAGS)"
