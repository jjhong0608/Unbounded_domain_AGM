CC = icc
AGMS=AGM.cpp
HEADERDIR=Header
HEADER=$(HEADERDIR)/AXLGEN.h $(HEADERDIR)/Basic++.h $(HEADERDIR)/Class++.h $(HEADERDIR)/Read.h $(HEADERDIR)/Stuff++.h $(HEADERDIR)/AGM.h $(HEADERDIR)/AxialData.h $(HEADERDIR)/CalcRepresenSol.h $(HEADERDIR)/ControlData.h $(HEADERDIR)/ftns.h $(HEADERDIR)/MatrixProcess.h $(HEADERDIR)/Point.h $(HEADERDIR)/util.h $(HEADERDIR)/CalcDiff.h $(HEADERDIR)/CalcNeumannpt.h
SOURCE=AGM
SERVERPATH=jhjo@165.246.142.45
TARGETPATH=$(SERVERPATH):/home/jhjo/180523
RESULTPATH=/Users/junhong/github/Unbounded_domain_AGM/results

$(SOURCE) : $(AGMS) $(HEADER)
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

move :
	scp $(AGMS) $(TARGETPATH)
	scp Makefile $(TARGETPATH)
	scp GeoInfo $(TARGETPATH)
	scp InputFile $(TARGETPATH)
	scp -r $(HEADERDIR) $(TARGETPATH)

result :
	scp $(TARGETPATH)/*.dat $(RESULTPATH)
