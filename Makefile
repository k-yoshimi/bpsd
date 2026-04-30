# Makfile

include ../task/make.header

.SUFFIXES:
.SUFFIXES: .f90 .f .mod .o .a

FFLAGS=$(OFLAGS)
#FFLAGS=$(DFLAGS)

OBJDIR=./obj

MODINCLUDE= -I./$(MOD)

SRCS = bpsd_kinds.f90 bpsd_constants.f90 bpsd_flags.f90 \
       bpsd_libchar.f90 bpsd_libfio.f90 bpsd_libspl.f90 \
       bpsd_types.f90 bpsd_types_internal.f90 bpsd_subs.f90 \
       bpsd_shot.f90 bpsd_device.f90 bpsd_species.f90 \
       bpsd_equ1D.f90 bpsd_metric1D.f90 bpsd_plasmaf.f90 \
       bpsd_trmatrix.f90 bpsd_trsource.f90 \
       bpsd_base.f90

OBJS=$(addprefix $(OBJDIR)/, $(SRCS:.f90=.o))

LIBS=libbpsd.a

all : libbpsd.a

# Pattern + directory rules below `all` so they don't become the default goal.
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(MOD):
	mkdir -p $(MOD)

${OBJDIR}/%.o:%.f90 | $(OBJDIR) $(MOD)
	$(FCFREE) $(FFLAGS) -c $< -o $@ $(MODDIR) $(MODINCLUDE)

libbpsd.a: $(OBJS)
	$(LD) $(LDFLAGS) $@ $(OBJS)

# Standalone BPSD tests with gfortran -fcheck=all -Wall -Wextra.
# Independent of TASK's make.header so it works in a bare BPSD checkout.
test:
	$(MAKE) -C tests test

testclean:
	$(MAKE) -C tests clean

clean:
	-rm -f core a.out *.o *.mod ./*~ ./#* *.a $(OBJDIR)/*.o $(MOD)/*.mod

veryclean: clean testclean

BPSD_COMMON = bpsd_subs.f90 bpsd_types_internal.f90 bpsd_types.f90 \
	      bpsd_flags.f90 bpsd_constants.f90 bpsd_kinds.f90 
BPSD_TYPES  = bpsd_shot.f90 bpsd_device.f90 bpsd_species.f90 bpsd_equ1D.f90 \
	      bpsd_metric1D.f90 bpsd_plasmaf.f90

$(OBJDIR)/bpsd_kinds.o:		bpsd_kinds.f90
$(OBJDIR)/bpsd_constants.o: 	bpsd_constants.f90 bpsd_kinds.f90
$(OBJDIR)/bpsd_flags.o:		bpsd_flags.f90
$(OBJDIR)/bpsd_libchar.o:	bpsd_libchar.f90
$(OBJDIR)/bpsd_libfio.o:	bpsd_libfio.f90
$(OBJDIR)/bpsd_libspl.o:	bpsd_libspl.f90
$(OBJDIR)/bpsd_types.o:		bpsd_types.f90 bpsd_kinds.f90
$(OBJDIR)/bpsd_types_internal.o:bpsd_types_internal.f90 bpsd_kinds.f90
$(OBJDIR)/bpsd_subs.o:		bpsd_subs.f90 bpsd_types_internal.f90 bpsd_kinds.f90
$(OBJDIR)/bpsd_shot.o:		bpsd_shot.f90 $(BPSD_COMMON)
$(OBJDIR)/bpsd_device.o:	bpsd_device.f90 $(BPSD_COMMON)
$(OBJDIR)/bpsd_species.o:	bpsd_species.f90 $(BPSD_COMMON)
$(OBJDIR)/bpsd_equ1D.o:		bpsd_equ1D.f90 $(BPSD_COMMON)
$(OBJDIR)/bpsd_metric1D.o:	bpsd_metric1D.f90 $(BPSD_COMMON)
$(OBJDIR)/bpsd_plasmaf.o:	bpsd_plasmaf.f90 $(BPSD_COMMON)
$(OBJDIR)/bpsd_trmatrix.o:	bpsd_trmatrix.f90 $(BPSD_COMMON)
$(OBJDIR)/bpsd_trsource.o:	bpsd_trsource.f90 $(BPSD_COMMON)
$(OBJDIR)/bpsd_base.o:		bpsd_base.f90 bpsd_subs.f90 $(BPSD_TYPES) $(BPSD_COMMON)
