# tests/Makefile -- standalone BPSD test harness.
# Independent of TASK's make.header so it works in a bare BPSD checkout.
# Uses gfortran with -fcheck=all -Wall -Wextra to catch bounds errors.

FC      = gfortran
FFLAGS  = -ffree-form -g -O0 -std=legacy \
          -fcheck=all -Wall -Wextra \
          -ffpe-trap=invalid,zero,overflow -fbacktrace
MODFLAG = -Jmod -Imod

OBJDIR  = obj
MODDIR  = mod
BPSD    = ..

# Order matters: each module must be compiled after its USEd modules.
BPSD_SRCS = \
    bpsd_kinds.f90 \
    bpsd_constants.f90 \
    bpsd_flags.f90 \
    bpsd_libchar.f90 \
    bpsd_libfio.f90 \
    bpsd_libspl.f90 \
    bpsd_types.f90 \
    bpsd_types_internal.f90 \
    bpsd_subs.f90 \
    bpsd_shot.f90 \
    bpsd_device.f90 \
    bpsd_species.f90 \
    bpsd_equ1D.f90 \
    bpsd_metric1D.f90 \
    bpsd_plasmaf.f90 \
    bpsd_trmatrix.f90 \
    bpsd_trsource.f90 \
    bpsd_base.f90

BPSD_OBJS = $(addprefix $(OBJDIR)/, $(BPSD_SRCS:.f90=.o))

TESTS = \
    test_put_get_equ1D \
    test_put_get_metric1D \
    test_put_get_plasmaf \
    test_put_get_species \
    test_put_get_trmatrix \
    test_put_get_trsource \
    test_save_load_roundtrip \
    test_get_device_no_data

TEST_BINS = $(addprefix $(OBJDIR)/, $(TESTS))

# Force serial build to avoid module-file race between BPSD sources.
.NOTPARALLEL:

.PHONY: all test clean

all: $(TEST_BINS)

test: $(TEST_BINS)
	@fail=0; \
	for t in $(TESTS); do \
	  echo "==== $$t ===="; \
	  ./$(OBJDIR)/$$t || fail=1; \
	done; \
	if [ $$fail -ne 0 ]; then echo "TEST FAILURE"; exit 1; fi; \
	echo "ALL TESTS PASSED"

$(OBJDIR) $(MODDIR):
	mkdir -p $@

# BPSD sources live in $(BPSD) (parent dir).
$(OBJDIR)/%.o: $(BPSD)/%.f90 | $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

# Tests live alongside this Makefile.
$(OBJDIR)/test_%.o: test_%.f90 $(BPSD_OBJS) | $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) $(MODFLAG) -c $< -o $@

# Each test links against the full BPSD object set.
$(OBJDIR)/test_%: $(OBJDIR)/test_%.o $(BPSD_OBJS)
	$(FC) $(FFLAGS) -o $@ $< $(BPSD_OBJS)

clean:
	rm -rf $(OBJDIR) $(MODDIR) bpsd.data
