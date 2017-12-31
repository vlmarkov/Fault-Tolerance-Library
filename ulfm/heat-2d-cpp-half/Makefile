COM_COLOR    := \033[0;34m
OBJ_COLOR    := \033[0;36m
OK_COLOR     := \033[0;32m
ERROR_COLOR  := \033[0;31m
WARN_COLOR   := \033[0;33m
NO_COLOR     := \033[m
OK_STRING    := "[OK]"
ERROR_STRING := "[ERROR]"
WARN_STRING  := "[WARNING]"
COM_STRING   := "Compiling"
CREAT_STRING := "Creating"
CLEAN_STRING := "Cleaning"

CXX          := g++
MPICXX       := mpicxx
CFLAGS       := -g -Wall -O2
LDFLAGS      := -lm
STDFLAGS     := -std=c++0x
DMPI         := -D MPI_SUPPORT
DULFM        := -D ULFM_SUPPORT

BUILD_DIR    := build

all: heat-2d

heat-2d:
	@mkdir -p $(BUILD_DIR)
	@RESULT=$$?; \
	if [ $$RESULT -ne 0 ]; then \
		printf "%-60b%b" "$(COM_COLOR)$(CREAT_STRING)$(OBJ_COLOR) $(BUILD_DIR)" "$(ERROR_COLOR)$(ERROR_STRING)$(NO_COLOR)\n"; \
	else  \
		printf "%-60b%b" "$(COM_COLOR)$(CREAT_STRING)$(OBJ_COLOR) $(BUILD_DIR)" "$(OK_COLOR)$(OK_STRING)$(NO_COLOR)\n"; \
	fi;
	@cp slurm.job $(BUILD_DIR)/slurm.job
	@$(MPICXX) $(CFLAGS) $(LDFLAGS) $(STDFLAGS) $(DMPI) $(DULFM) \
	Utils.cpp \
	Grid.cpp \
	Task.cpp \
	-o $(BUILD_DIR)/$@ heat-2d-compute-redundancy.cpp; \
	RESULT=$$?; \
	if [ $$RESULT -ne 0 ]; then \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $@" "$(ERROR_COLOR)$(ERROR_STRING)$(NO_COLOR)\n"; \
	else  \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $(@F)" "$(OK_COLOR)$(OK_STRING)$(NO_COLOR)\n"; \
	fi;

debug:
	@mkdir -p $(BUILD_DIR)
	@RESULT=$$?; \
	if [ $$RESULT -ne 0 ]; then \
		printf "%-60b%b" "$(COM_COLOR)$(CREAT_STRING)$(OBJ_COLOR) $(BUILD_DIR)" "$(ERROR_COLOR)$(ERROR_STRING)$(NO_COLOR)\n"; \
	else  \
		printf "%-60b%b" "$(COM_COLOR)$(CREAT_STRING)$(OBJ_COLOR) $(BUILD_DIR)" "$(OK_COLOR)$(OK_STRING)$(NO_COLOR)\n"; \
	fi;
	@$(CXX) $(CFLAGS) $(LDFLAGS) $(STDFLAGS) \
	Grid.cpp \
	Task.cpp \
	-o $(BUILD_DIR)/$@ test.cpp; \
	RESULT=$$?; \
	if [ $$RESULT -ne 0 ]; then \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $@" "$(ERROR_COLOR)$(ERROR_STRING)$(NO_COLOR)\n"; \
	else  \
		printf "%-60b%b" "$(COM_COLOR)$(COM_STRING)$(OBJ_COLOR) $(@F)" "$(OK_COLOR)$(OK_STRING)$(NO_COLOR)\n"; \
	fi;

clean:
	@rm -Rf $(BUILD_DIR)
	@RESULT=$$?; \
	if [ $$RESULT -ne 0 ]; then \
		printf "%-60b%b" "$(COM_COLOR)$(CLEAN_STRING)$(OBJ_COLOR) $(BUILD_DIR)" "$(ERROR_COLOR)$(ERROR_STRING)$(NO_COLOR)\n"; \
	else  \
		printf "%-60b%b" "$(COM_COLOR)$(CLEAN_STRING)$(OBJ_COLOR) $(BUILD_DIR)" "$(OK_COLOR)$(OK_STRING)$(NO_COLOR)\n"; \
	fi;
