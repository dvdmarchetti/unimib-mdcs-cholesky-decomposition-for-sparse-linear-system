# Compiler flags
CC = g++
# CFLAGS = -Ilibs/eigen -Iinclude -lpsapi --std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0
CFLAGS = -Ilibs/eigen -Ilibs/csv/include --std=c++17 -Iinclude -lpsapi

# Project files
SRCS = main.cpp src/memory_xplatform.cpp src/csv_utils.cpp
OBJS = $(SRCS:.cpp=.o)
EXE = main

# Debug build settings
DBGDIR = build/debug
DBGEXE = $(DBGDIR)/$(EXE)
DBGOBJS = $(addprefix $(DBGDIR)/, $(OBJS))
DBGCFLAGS = -g -O0 -DDEBUG

# Release build settings
RELDIR = build/release
RELEXE = $(RELDIR)/$(EXE)
RELOBJS = $(addprefix $(RELDIR)/, $(OBJS))
RELCFLAGS = -O3 -DDEBUG

.PHONY: all clean debug prep release remake

# Default build
all: prep release

# Debug rules
debug: $(DBGEXE)

$(DBGEXE): $(DBGOBJS)
	$(CC) -o $(DBGEXE) $^ $(CFLAGS) $(DBGCFLAGS)

$(DBGDIR)/%.o: %.cpp
	$(CC) -o $@ $< -c $(CFLAGS) $(DBGCFLAGS)

# Release rules
release: $(RELEXE)

$(RELEXE): $(RELOBJS)
	$(CC) -o $(RELEXE) $^ $(CFLAGS) $(RELCFLAGS)

$(RELDIR)/%.o: %.cpp
	$(CC) -o $@ $< -c $(CFLAGS) $(RELCFLAGS)

# Other rules
prep:
	@mkdir -p $(DBGDIR) $(DBGDIR)/src $(RELDIR) $(RELDIR)/src

remake: clean all

clean:
	rm -f $(RELEXE) $(RELOBJS) $(DBGEXE) $(DBGOBJS)