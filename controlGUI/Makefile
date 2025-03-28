# GNU Make project makefile autogenerated by Premake

ifndef config
  config=debug_x86_64
endif

ifndef verbose
  SILENT = @
endif

.PHONY: clean prebuild prelink

ifeq ($(config),debug_x86_64)
  ifeq ($(origin CC), default)
    CC = clang
  endif
  ifeq ($(origin CXX), default)
    CXX = clang++
  endif
  ifeq ($(origin AR), default)
    AR = ar
  endif
  TARGETDIR = ../bin/Debug-linux-x86_64/controlGUI
  TARGET = $(TARGETDIR)/controlGUI
  OBJDIR = ../bin-int/Debug-linux-x86_64/controlGUI
  DEFINES += -DF_V2_API= -DF_AUX_API= -DF_CLIENTAPP -DPFM_API= -DSYS_ARCH=x86_64 -DF_OS_LINUX -DPFM_CONTROL_GUI -DAS_BUILD_APP -DX64 -DAS_DEBUG
  INCLUDES += -Iinclude -I../depend/fAux/include -I../depend/fViz2D/include -I../simulation/API
  FORCE_INCLUDE +=
  ALL_CPPFLAGS += $(CPPFLAGS) -MD -MP $(DEFINES) $(INCLUDES)
  ALL_CFLAGS += $(CFLAGS) $(ALL_CPPFLAGS) -m64 -g -Wall -Wextra
  ALL_CXXFLAGS += $(CXXFLAGS) $(ALL_CPPFLAGS) -m64 -g -Wall -Wextra -std=c++17
  ALL_RESFLAGS += $(RESFLAGS) $(DEFINES) $(INCLUDES)
  LIBS += ../bin/Debug-linux-x86_64/simulation/libsimulation.a -lfAux -lglfw -lOpenGL -lfViz2D
  LDDEPS += ../bin/Debug-linux-x86_64/simulation/libsimulation.a
  ALL_LDFLAGS += $(LDFLAGS) -L../depend/fAux/lib/Debug-linux-x86_64 -L../depend/fViz2D/lib/Debug-linux-x86_64 -L/usr/lib64 -m64
  LINKCMD = $(CXX) -o "$@" $(OBJECTS) $(RESOURCES) $(ALL_LDFLAGS) $(LIBS)
  define PREBUILDCMDS
  endef
  define PRELINKCMDS
  endef
  define POSTBUILDCMDS
	@echo Running postbuild commands
	cp -f ../res/guiCfg/imgui.ini ../bin/Debug-linux-x86_64/controlGUI
  endef
all: prebuild prelink $(TARGET)
	@:

endif

ifeq ($(config),release_x86_64)
  ifeq ($(origin CC), default)
    CC = clang
  endif
  ifeq ($(origin CXX), default)
    CXX = clang++
  endif
  ifeq ($(origin AR), default)
    AR = ar
  endif
  TARGETDIR = ../bin/Release-linux-x86_64/controlGUI
  TARGET = $(TARGETDIR)/controlGUI
  OBJDIR = ../bin-int/Release-linux-x86_64/controlGUI
  DEFINES += -DF_V2_API= -DF_AUX_API= -DF_CLIENTAPP -DPFM_API= -DSYS_ARCH=x86_64 -DF_OS_LINUX -DPFM_CONTROL_GUI -DAS_BUILD_APP -DX64 -DAS_RELEASE
  INCLUDES += -Iinclude -I../depend/fAux/include -I../depend/fViz2D/include -I../simulation/API
  FORCE_INCLUDE +=
  ALL_CPPFLAGS += $(CPPFLAGS) -MD -MP $(DEFINES) $(INCLUDES)
  ALL_CFLAGS += $(CFLAGS) $(ALL_CPPFLAGS) -m64 -O2 -Wall -Wextra
  ALL_CXXFLAGS += $(CXXFLAGS) $(ALL_CPPFLAGS) -m64 -O2 -Wall -Wextra -std=c++17
  ALL_RESFLAGS += $(RESFLAGS) $(DEFINES) $(INCLUDES)
  LIBS += ../bin/Release-linux-x86_64/simulation/libsimulation.a -lfAux -lglfw -lOpenGL -lfViz2D
  LDDEPS += ../bin/Release-linux-x86_64/simulation/libsimulation.a
  ALL_LDFLAGS += $(LDFLAGS) -L../depend/fAux/lib/Release-linux-x86_64 -L../depend/fViz2D/lib/Release-linux-x86_64 -L/usr/lib64 -m64
  LINKCMD = $(CXX) -o "$@" $(OBJECTS) $(RESOURCES) $(ALL_LDFLAGS) $(LIBS)
  define PREBUILDCMDS
  endef
  define PRELINKCMDS
  endef
  define POSTBUILDCMDS
	@echo Running postbuild commands
	cp -f ../res/guiCfg/imgui.ini ../bin/Release-linux-x86_64/controlGUI
  endef
all: prebuild prelink $(TARGET)
	@:

endif

OBJECTS := \
	$(OBJDIR)/GUIcontrolMain.o \
	$(OBJDIR)/guiMenus.o \
	$(OBJDIR)/guiTests.o \

RESOURCES := \

CUSTOMFILES := \

SHELLTYPE := posix
ifeq (.exe,$(findstring .exe,$(ComSpec)))
	SHELLTYPE := msdos
endif

$(TARGET): $(GCH) ${CUSTOMFILES} $(OBJECTS) $(LDDEPS) $(RESOURCES) | $(TARGETDIR)
	@echo Linking controlGUI
	$(SILENT) $(LINKCMD)
	$(POSTBUILDCMDS)

$(CUSTOMFILES): | $(OBJDIR)

$(TARGETDIR):
	@echo Creating $(TARGETDIR)
ifeq (posix,$(SHELLTYPE))
	$(SILENT) mkdir -p $(TARGETDIR)
else
	$(SILENT) mkdir $(subst /,\\,$(TARGETDIR))
endif

$(OBJDIR):
	@echo Creating $(OBJDIR)
ifeq (posix,$(SHELLTYPE))
	$(SILENT) mkdir -p $(OBJDIR)
else
	$(SILENT) mkdir $(subst /,\\,$(OBJDIR))
endif

clean:
	@echo Cleaning controlGUI
ifeq (posix,$(SHELLTYPE))
	$(SILENT) rm -f  $(TARGET)
	$(SILENT) rm -rf $(OBJDIR)
else
	$(SILENT) if exist $(subst /,\\,$(TARGET)) del $(subst /,\\,$(TARGET))
	$(SILENT) if exist $(subst /,\\,$(OBJDIR)) rmdir /s /q $(subst /,\\,$(OBJDIR))
endif

prebuild:
	$(PREBUILDCMDS)

prelink:
	$(PRELINKCMDS)

ifneq (,$(PCH))
$(OBJECTS): $(GCH) $(PCH) | $(OBJDIR)
$(GCH): $(PCH) | $(OBJDIR)
	@echo $(notdir $<)
	$(SILENT) $(CXX) -x c++-header $(ALL_CXXFLAGS) -o "$@" -MF "$(@:%.gch=%.d)" -c "$<"
else
$(OBJECTS): | $(OBJDIR)
endif

$(OBJDIR)/GUIcontrolMain.o: src/GUIcontrolMain.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(ALL_CXXFLAGS) $(FORCE_INCLUDE) -o "$@" -MF "$(@:%.o=%.d)" -c "$<"
$(OBJDIR)/guiMenus.o: src/guiMenus.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(ALL_CXXFLAGS) $(FORCE_INCLUDE) -o "$@" -MF "$(@:%.o=%.d)" -c "$<"
$(OBJDIR)/guiTests.o: src/guiTests.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(ALL_CXXFLAGS) $(FORCE_INCLUDE) -o "$@" -MF "$(@:%.o=%.d)" -c "$<"

-include $(OBJECTS:%.o=%.d)
ifneq (,$(PCH))
  -include $(OBJDIR)/$(notdir $(PCH)).d
endif