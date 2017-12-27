COM_COLOR   = \033[0;34m
OBJ_COLOR   = \033[0;36m
OK_COLOR    = \033[0;32m
ERROR_COLOR = \033[0;31m
WARN_COLOR  = \033[0;33m
NO_COLOR    = \033[m

OK_STRING    = "[OK]"
ERROR_STRING = "[ERROR]"
WARN_STRING  = "[WARNING]"
COM_STRING   = "Compiling"

NAME     := heat-2d-compute-redundancy
CXX      := mpicxx
CFLAGS   := -g -Wall -O2
LDFLAGS  := -lm
STDFLAGS := -std=c++0x
DFLAGS   := -D MPI_SUPPORT

all: compute-redundancy

compute-redundancy: ft-lib
	@mkdir -p build
	@cp slurm.job build/slurm.job
	@cp torque.job build/torque.job
	@$(CXX) $(CFLAGS) -o build/$(NAME) $(NAME).cpp build/ft-lib.a $(LDFLAGS) $(STDFLAGS) $(DFLAGS); \
	RESULT=$$?; \
	if [ $$RESULT -ne 0 ]; then \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $@" "$(ERROR_COLOR)$(ERROR_STRING)$(NO_COLOR)\n"; \
	else \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $(@F)" "$(OK_COLOR)$(OK_STRING)$(NO_COLOR)\n"; \
	fi; \
	exit $$RESULT;

ft-lib: Utils Grid Redundancy Logger
	@ar cr build/$@.a build/Utils.o build/Grid.o build/Redundancy.o build/Logger.o; \
	RESULT=$$?; \
	if [ $$RESULT -ne 0 ]; then \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $@" "$(ERROR_COLOR)$(ERROR_STRING)$(NO_COLOR)\n"; \
	else \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $(@F)" "$(OK_COLOR)$(OK_STRING)$(NO_COLOR)\n"; \
	fi; \
	exit $$RESULT;

Utils:
	@mkdir -p build
	@$(CXX) $(CFLAGS) -c -o build/$@.o src/utils.cpp $(LDFLAGS) $(STDFLAGS) $(DFLAGS); \
	RESULT=$$?; \
	if [ $$RESULT -ne 0 ]; then \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $@" "$(ERROR_COLOR)$(ERROR_STRING)$(NO_COLOR)\n"; \
	else \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $(@F)" "$(OK_COLOR)$(OK_STRING)$(NO_COLOR)\n"; \
	fi; \
	exit $$RESULT;

Grid:
	@mkdir -p build
	@$(CXX) $(CFLAGS) -c -o build/$@.o src/grid-task.cpp $(LDFLAGS) $(STDFLAGS) $(DFLAGS); \
	RESULT=$$?; \
	if [ $$RESULT -ne 0 ]; then \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $@" "$(ERROR_COLOR)$(ERROR_STRING)$(NO_COLOR)\n"; \
	else \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $(@F)" "$(OK_COLOR)$(OK_STRING)$(NO_COLOR)\n"; \
	fi; \
	exit $$RESULT;

Redundancy:
	@mkdir -p build
	@$(CXX) $(CFLAGS) -c -o build/$@.o src/redundancy.cpp $(LDFLAGS) $(STDFLAGS) $(DFLAGS); \
	RESULT=$$?; \
	if [ $$RESULT -ne 0 ]; then \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $@" "$(ERROR_COLOR)$(ERROR_STRING)$(NO_COLOR)\n"; \
	else \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $(@F)" "$(OK_COLOR)$(OK_STRING)$(NO_COLOR)\n"; \
	fi; \
	exit $$RESULT;

Logger:
	@mkdir -p build
	@$(CXX) $(CFLAGS) -c -o build/$@.o src/$@.cpp $(LDFLAGS) $(STDFLAGS) $(DFLAGS); \
	RESULT=$$?; \
	if [ $$RESULT -ne 0 ]; then \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $@" "$(ERROR_COLOR)$(ERROR_STRING)$(NO_COLOR)\n"; \
	else \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $(@F)" "$(OK_COLOR)$(OK_STRING)$(NO_COLOR)\n"; \
	fi; \
	exit $$RESULT;

debug:
	@mkdir -p build
	@g++ $(CFLAGS) -o build/$@ test.cpp src/grid-task.cpp src/utils.cpp src/redundancy.cpp $(LDFLAGS) $(STDFLAGS); \
	RESULT=$$?; \
	if [ $$RESULT -ne 0 ]; then \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $@" "$(ERROR_COLOR)$(ERROR_STRING)$(NO_COLOR)\n"; \
	else  \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $(@F)" "$(OK_COLOR)$(OK_STRING)$(NO_COLOR)\n"; \
	fi

clean:
	rm -Rf build
