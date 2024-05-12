#------------------------------------------------#

PROJECT	:= qualetize
BUILD	:= build
RELEASE := release
CFLAGS	:= -O2 -Wall -Wextra -Isource -Iinclude -s
LIBS	:= -lm
RM	:= rm -rf

#------------------------------------------------#

CFILES	:= $(wildcard source/*.c)
OFILES	:= $(addprefix $(BUILD)/, $(addsuffix .o, $(CFILES)))
DFILES	:= $(OFILES:.o=.d)

OFILES_EXE := $(OFILES)
OFILES_DLL := $(filter-out $(BUILD)/source/qualetize-cli.c.o, $(OFILES))

#------------------------------------------------#

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
IS_UNIX = true
endif
ifeq ($(UNAME), Darwin)
IS_UNIX = true
endif
ifdef IS_UNIX
EXE = $(PROJECT)
DLL = lib$(PROJECT).so
else
EXE = $(PROJECT).exe
DLL = lib$(PROJECT).dll
endif

#------------------------------------------------#

all : $(RELEASE)/$(EXE) $(RELEASE)/$(DLL)

$(RELEASE)/$(EXE) : $(OFILES_EXE) | $(RELEASE)
	$(CC) $(LIBS) -o $@ $^

$(RELEASE)/$(DLL) : $(OFILES_DLL) | $(RELEASE)
	$(CC) $(LIBS) -shared -o $@ $^

$(RELEASE) :
	mkdir $@

-include $(DFILES)

#------------------------------------------------#

.PHONY: clean

clean:
	$(RM) $(RELEASE) $(BUILD)

#------------------------------------------------#

$(BUILD)/%.c.o : %.c
	mkdir -p $(dir $@)
	$(CC) -c $(CFLAGS) -DDECLSPEC="$(DDECLSPEC)" -MD -MP -MF $(BUILD)/$<.d -o $@ $<

#------------------------------------------------#
