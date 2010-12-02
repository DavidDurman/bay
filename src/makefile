# @project Learn Bayesian network structure from dataset
#	   then use it for prediction
# @author David Durman, 2008

EX_BAYBUILD=baybuild
EX_BAYPREDICT=baypredict
EXS=$(EX_BAYBUILD) $(EX_BAYPREDICT)

DFLAG=-DNDEBUG	
CXX=g++ 

CXX_LIBS=arrangement.cc BNGenome.cc

GA_INC_DIR= /usr/local/include
GA_LIB_DIR= /usr/local/lib

INC_DIRS= -I$(GA_INC_DIR)
LIB_DIRS= -L$(GA_LIB_DIR)


all: $(EXS)

$(EX_BAYBUILD): %: %.o
	$(CXX) $@.o -o $@ $(LIB_DIRS) -lga -lm $(CXX_LIBS)

$(EX_BAYPREDICT): %: %.o
	$(CXX) $@.o -o $@

.PHONY=clean clean-all test pack

clean:
	rm -f *~ $(EXS)
	rm -f *.o

clean-all: clean
	rm -f votes.png votes.net votes.prob votes.dot

pack: clean
	tar zcfv bay.tar.gz arrangement.cc arrangement.h \
		configuration.h BNDataset.h BNGenome.h BNGenome.cc \
		baybuild.cc baypredict.cc \
		bayname.py \
		votes.csv \
		votes.info \
		votes.names \
		makefile \
		README \
		INSTALL

test: 
	@echo "Running tests. This could take some time, depending on"
	@echo "the speed of your computer."
	./baybuild -s 333 -d votes.csv -p 100 -g 10 -m 0.2 -c 0.6 -b votes -v
	./bayname.py votes.names votes.dot >votes-names.dot
	cat votes-names.dot | dot -Tpng >votes.png
	./baypredict -n votes.net -p votes.prob -e 0=1,2=0 -c 3=1
