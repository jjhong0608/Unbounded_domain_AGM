CC = g++
SRC = AXLGEN

$(SRC): AGM.cpp
	$(CC) AGM.cpp -lm -std=c++0x -o $(SRC)

clean :
	rm -f $(SRC)
	rm -f *.dat
	rm -f *.txt
	rm -f ALG_Output*

run : clean
	make
	./$(SRC) InputFile
