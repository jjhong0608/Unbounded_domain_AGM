CC = icc
AGMS=integrated_AGM.cpp
ALGDIR=ALG_Header
AGMDIR=AGM_Header
ALGH=$(ALGDIR)/AXLGEN.h $(ALGDIR)/Basic++.h $(ALGDIR)/Class++.h $(ALGDIR)/Read.h $(ALGDIR)/Stuff++.h
AGMH=$(AGMDIR)/LAGM.h $(AGMDIR)/AxialData.h $(AGMDIR)/CalcRepresenSol.h $(AGMDIR)/Class.h $(AGMDIR)/ControlData.h $(AGMDIR)/Ftns.h $(AGMDIR)/MatrixProcess.h $(AGMDIR)/Point.h $(AGMDIR)/util.h
SOURCE=AGM

$(SOURCE) : $(AGMS) $(ALGH) $(AGMH)
	$(CC)  -m64  -w -I"/home/jhjo/intel/compilers_and_libraries_2017.4.196/linux/mkl/include" \
	$(AGMS) -Wl,--start-group \
	"/home/jhjo/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_intel_lp64.a" \
	"/home/jhjo/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_intel_thread.a" \
	"/home/jhjo/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_core.a" \
	-Wl,--end-group -L"/home/jhjo/intel/compilers_and_libraries_2017.4.196/linux/mkl/../compiler/lib/intel64" -liomp5 -lpthread -lm -ldl -std=c++0x -o $(SOURCE)

clean :
	rm -f $(SOURCE)
	rm -f *.dat
	rm -f *.txt
	rm -f ALG_Output*

run :
	make clean
	make
	./$(SOURCE) InputFile
