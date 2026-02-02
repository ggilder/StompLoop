# Makefile for looper~ external

PD_PATH ?= /usr/local/include        # Path to Pd headers
PD_SRC  ?= looper~.c
PD_OBJ  := $(PD_SRC:.c=.o)

# Detect platform
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  EXT = pd_darwin
  CFLAGS += -fPIC -O2 -I$(PD_PATH) -Wall -W -Wshadow -Wstrict-prototypes -Wno-unused -Wno-cast-function-type
  LDFLAGS += -bundle -undefined dynamic_lookup
else
  EXT = pd_linux
  CFLAGS += -fPIC -O2 -I$(PD_PATH) -Wall -W -Wshadow -Wstrict-prototypes -Wno-unused -Wno-cast-function-type
  LDFLAGS += -shared
endif

TARGET = looper~.$(EXT)

all: $(TARGET)

$(TARGET): $(PD_SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGET) $(PD_OBJ)

.PHONY: all clean
