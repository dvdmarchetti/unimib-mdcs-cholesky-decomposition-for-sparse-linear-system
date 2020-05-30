# Compiler flags
CC = g++
_CFLAGS = -Ilibs/eigen --std=c++17 -Iinclude
ifeq ($(OS),Windows_NT)
CFLAGS = $(_CFLAGS) -lpsapi
else
CFLAGS = $(_CFLAGS)
endif

# Project files
SRCS = main.cpp src/memory_xplatform.cpp
OBJS = $(SRCS:.cpp=.o)
EXE = main

# Debug build settings
DBGDIR = build/debug
DBGEXE = $(DBGDIR)/$(EXE)
DBGOBJS = $(addprefix $(DBGDIR)/, $(OBJS))
DBGCFLAGS = -O0 -DDEBUG -fno-omit-frame-pointer -gdwarf-2

# Release build settings
RELDIR = build/release
RELEXE = $(RELDIR)/$(EXE)
RELOBJS = $(addprefix $(RELDIR)/, $(OBJS))
RELCFLAGS = -Ofast -DDEBUG

.PHONY: all clean debug prep release remake

# Default build
all: prep release

windows:

unix:

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