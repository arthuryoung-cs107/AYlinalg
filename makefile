include config/makefile.includes

# AYlinalg object files
AYOBJS:= $(addprefix $(AY_DIR), $(addsuffix .o, $(AYLINALG)))
SOLOBJS:= $(AYOBJS) $(addprefix $(SOL_DIR), $(addsuffix .o, $(SOLVERS)))

MTESTOBJS= $(AYOBJS)

STESTOBJS= $(SOLOBJS)

SCRATCHOBJS= $(AYOBJS)

JUNKOBJS=

# rules for each directory
# AYlinalg rules
$(AY_DIR)%.o: $(AY_SRC)%.c | $(AY_DIR)
	$(CC) $(IDIR) $(CFLAGS) -c $< -o $@

$(AY_DIR)%.o: $(AY_SRC)%.cc | $(AY_DIR)
	$(CXX) $(IDIR) $(CFLAGS) -c $< -o $@

$(SOL_DIR)%.o: $(SOL_SRC)%.cc | $(SOL_DIR)
	$(CXX) $(IDIR) $(CFLAGS) -c $< -o $@

all: math_test solve_test scratch_test junk_test

math_test: $(TEST_SRC)math_test.cc $(MTESTOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

solve_test: $(TEST_SRC)solve_test.cc $(STESTOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

scratch_test: $(TEST_SRC)scratch.cc $(SCRATCHOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

junk_test: $(TEST_SRC)junk.cc $(JUNKOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

$(AY_DIR) $(SOL_DIR):
	mkdir -p $@

clean_:
	rm -f *_test

clean_AY:
	rm -f $(AY_DIR)*.o

clean_SOL:
	rm -f $(SOL_DIR)*.o

clean: clean_

clean_all: clean_ clean_SOL clean_AY

clean_dat:
	rm -f *.dat
	rm -f *.aysml
	rm -f *.aydat

clean_datdir:
	rm -f ./dat_dir/*.dat
	rm -f ./dat_dir/*.aysml
	rm -f ./dat_dir/*.aydat
