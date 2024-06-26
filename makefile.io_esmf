# these settings for compiling standalone on Compaq. Type "make -r"
CPP = /lib/cpp -P
#FC  = mpif90 -r4 -i4  -w -Mfree -byteswapio 
SFC         =       gfortran
FC          =       mpif90 -f90=$(SFC)
ESMF_DIR=/usr2/luca/Devs/esmf-8.6.0
ESMF_INCLUDE  = ${ESMF_DIR}/mod/modO/Linux.gfortran.64.openmpi.default
WRF_DIR=/usr2/luca/Devs/shyfem-wrf-esmf/wrf/WRF
NETCDF_DIR=/home/model/lib/lib
CURRENT_DIR=/usr2/luca/Devs/shyfem-wrf-esmf/wrf/WRF/external/io_esmf
NUOPCMOD      = atmosphere_wrf
NUOPC_OBJ     = nuopc_wrf_cap.o
INCLUDE_MODULES = -I$(CURRENT_DIR)/../../ \
			      -I$(CURRENT_DIR)/../../main \
                  -I$(CURRENT_DIR)/../../external/io_netcdf \
                  -I$(CURRENT_DIR)/../../external/io_int \
                  -I$(CURRENT_DIR)/../../frame \
                  -I$(CURRENT_DIR)/../../share \
                  -I$(CURRENT_DIR)/../../phys \
                  -I$(CURRENT_DIR)/../../wrftladj \
                  -I$(CURRENT_DIR)/../../chem \
                  -I$(CURRENT_DIR)/../../inc

.SUFFIXES: .F90 .o

AR = ar
RANLIB	= echo

OBJS = module_symbols_util.o \
       module_esmf_extensions.o \
       module_utility.o \
       io_esmf.o \
       ext_esmf_open_for_read.o \
       ext_esmf_open_for_write.o \
       ext_esmf_read_field.o \
       ext_esmf_write_field.o

TARGET = libwrfio_esmf.a

nuopc: nuopc_clean nuopc_wrf.mk $(TARGET)

nuopc_clean:
	-rm -f *.mk

nuopc_wrf.mk:
	@echo "creating nuopc makefile fragment"
	@echo "# ESMF self-describing build dependency makefile fragment for WRF model" > $@
	@echo >> $@
	@echo "WRF_DEP_FRONT     = "$(NUOPCMOD)                   >> $@        
	@echo "WRF_DEP_INCPATH   = "$(WRF_DIR)/dyn_em                \
                                    $(WRF_DIR)/main                  \
                                    $(WRF_DIR)/external/io_esmf      \
                                    $(WRF_DIR)/external/io_netcdf    \
                                    $(WRF_DIR)/external/io_int       \
                                    $(WRF_DIR)/frame                 \
                                    $(WRF_DIR)/share                 \
                                    $(WRF_DIR)/phys                  \
                                    $(WRF_DIR)/chem                  \
                                    $(WRF_DIR)/inc                >> $@
	@echo "WRF_DEP_CMPL_OBJS = "$(WRF_DIR)/main/$(NUOPC_OBJ)  >> $@
	@echo "WRF_DEP_LINK_OBJS = "$(WRF_DIR)/main/nuopc_wrf_cap.o                     \
	                            $(WRF_DIR)/main/module_wrf_top.o                     \
                                    $(WRF_DIR)/main/libwrflib.a                          \
                                    $(WRF_DIR)/external/fftpack/fftpack5/libfftpack.a    \
                                    $(WRF_DIR)/external/io_grib1/libio_grib1.a           \
                                    $(WRF_DIR)/external/io_grib_share/libio_grib_share.a \
                                    $(WRF_DIR)/external/io_int/libwrfio_int.a            \
                                    $(WRF_DIR)/external/io_int/libwrfio_int.a            \
                                    $(WRF_DIR)/external/RSL_LITE/librsl_lite.a           \
                                    $(WRF_DIR)/frame/module_internal_header_util.o       \
                                    $(WRF_DIR)/frame/pack_utils.o >> $@
	@echo "WRF_DEP_SHRD_PATH = "$(NETCDF_DIR)                                \
                                    $(WRF_DIR)/external/io_netcdf                \
                                    $(WRF_DIR)/external/io_esmf   >> $@
	@echo "WRF_DEP_SHRD_LIBS = "netcdff netcdf wrfio_nf wrfio_esmf >> $@

$(TARGET) : $(OBJS)
	$(AR) ru $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)

.F90.o :
	$(CPP) -I../ioapi_share -P -traditional -DESMF_COUPLING $*.F90 > $*.f
	$(FC) -I$(ESMF_INCLUDE) -I$(CURRENT_DIR) $(INCLUDE_MODULES) -c -g -I../ioapi_share $*.f

superclean: 
	@/bin/rm -f *.f *.o $(TARGET) *.mod

# DEPENDENCIES : only dependencies after this line 

module_symbols_util.o :
module_esmf_extensions.o :
module_utility.o : module_symbols_util.o module_esmf_extensions.o
io_esmf.o : module_esmf_extensions.o
ext_esmf_open_for_read.o : io_esmf.o
ext_esmf_open_for_write.o : io_esmf.o
ext_esmf_read_field.o : io_esmf.o
ext_esmf_write_field.o : io_esmf.o

