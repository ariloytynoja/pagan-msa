#############################################################################
# Makefile for building: papaya
# Generated by qmake (2.01a) (Qt 4.6.2) on: Thu Mar 4 22:45:10 2010
# Project:  papaya.pro
# Template: app
# Command: /usr/bin/qmake-qt4 -spec /usr/share/qt4/mkspecs/linux-g++ -unix -o Makefile papaya.pro
#############################################################################

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = -DQT_NO_DEBUG -DQT_CORE_LIB -DQT_SHARED
CFLAGS        = -pipe -O2 -Wall -W -D_REENTRANT $(DEFINES)
CXXFLAGS      = -pipe -O2 -Wall -W -D_REENTRANT $(DEFINES)
INCPATH       = -I/usr/share/qt4/mkspecs/linux-g++ -I. -I/usr/include/qt4/QtCore -I/usr/include/qt4 -I/usr/include -I.
LINK          = g++
LFLAGS        = -Wl,-O1
LIBS          = $(SUBLIBS)  -L/usr/lib -lboost_program_options-mt -lQtCore -lpthread 
AR            = ar cqs
RANLIB        = 
QMAKE         = /usr/bin/qmake-qt4
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = $(COPY)
COPY_DIR      = $(COPY) -r
STRIP         = strip
INSTALL_FILE  = install -m 644 -p
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = install -m 755 -p
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = main.cpp \
		utils/text_utils.cpp \
		utils/settings.cpp \
		utils/node.cpp \
		utils/newick_reader.cpp \
		utils/model_factory.cpp \
		utils/fasta_reader.cpp \
		utils/eigen.cpp \
		utils/db_matrix.cpp \
		main/sequence.cpp \
		main/simple_alignment.cpp \
		utils/int_matrix.cpp \
		utils/settings_handle.cpp \
		utils/xml_writer.cpp \
		utils/evol_model.cpp \
		main/reads_alignment.cpp 
OBJECTS       = main.o \
		text_utils.o \
		settings.o \
		node.o \
		newick_reader.o \
		model_factory.o \
		fasta_reader.o \
		eigen.o \
		db_matrix.o \
		sequence.o \
		simple_alignment.o \
		int_matrix.o \
		settings_handle.o \
		xml_writer.o \
		evol_model.o \
		reads_alignment.o
DIST          = /usr/share/qt4/mkspecs/common/g++.conf \
		/usr/share/qt4/mkspecs/common/unix.conf \
		/usr/share/qt4/mkspecs/common/linux.conf \
		/usr/share/qt4/mkspecs/qconfig.pri \
		/usr/share/qt4/mkspecs/features/qt_functions.prf \
		/usr/share/qt4/mkspecs/features/qt_config.prf \
		/usr/share/qt4/mkspecs/features/exclusive_builds.prf \
		/usr/share/qt4/mkspecs/features/default_pre.prf \
		/usr/share/qt4/mkspecs/features/release.prf \
		/usr/share/qt4/mkspecs/features/default_post.prf \
		/usr/share/qt4/mkspecs/features/warn_on.prf \
		/usr/share/qt4/mkspecs/features/qt.prf \
		/usr/share/qt4/mkspecs/features/unix/thread.prf \
		/usr/share/qt4/mkspecs/features/moc.prf \
		/usr/share/qt4/mkspecs/features/resources.prf \
		/usr/share/qt4/mkspecs/features/uic.prf \
		/usr/share/qt4/mkspecs/features/yacc.prf \
		/usr/share/qt4/mkspecs/features/lex.prf \
		/usr/share/qt4/mkspecs/features/include_source_dir.prf \
		papaya.pro
QMAKE_TARGET  = papaya
DESTDIR       = 
TARGET        = papaya

first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all: Makefile $(TARGET)

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

Makefile: papaya.pro  /usr/share/qt4/mkspecs/linux-g++/qmake.conf /usr/share/qt4/mkspecs/common/g++.conf \
		/usr/share/qt4/mkspecs/common/unix.conf \
		/usr/share/qt4/mkspecs/common/linux.conf \
		/usr/share/qt4/mkspecs/qconfig.pri \
		/usr/share/qt4/mkspecs/features/qt_functions.prf \
		/usr/share/qt4/mkspecs/features/qt_config.prf \
		/usr/share/qt4/mkspecs/features/exclusive_builds.prf \
		/usr/share/qt4/mkspecs/features/default_pre.prf \
		/usr/share/qt4/mkspecs/features/release.prf \
		/usr/share/qt4/mkspecs/features/default_post.prf \
		/usr/share/qt4/mkspecs/features/warn_on.prf \
		/usr/share/qt4/mkspecs/features/qt.prf \
		/usr/share/qt4/mkspecs/features/unix/thread.prf \
		/usr/share/qt4/mkspecs/features/moc.prf \
		/usr/share/qt4/mkspecs/features/resources.prf \
		/usr/share/qt4/mkspecs/features/uic.prf \
		/usr/share/qt4/mkspecs/features/yacc.prf \
		/usr/share/qt4/mkspecs/features/lex.prf \
		/usr/share/qt4/mkspecs/features/include_source_dir.prf \
		/usr/lib/libQtCore.prl
	$(QMAKE) -spec /usr/share/qt4/mkspecs/linux-g++ -unix -o Makefile papaya.pro
/usr/share/qt4/mkspecs/common/g++.conf:
/usr/share/qt4/mkspecs/common/unix.conf:
/usr/share/qt4/mkspecs/common/linux.conf:
/usr/share/qt4/mkspecs/qconfig.pri:
/usr/share/qt4/mkspecs/features/qt_functions.prf:
/usr/share/qt4/mkspecs/features/qt_config.prf:
/usr/share/qt4/mkspecs/features/exclusive_builds.prf:
/usr/share/qt4/mkspecs/features/default_pre.prf:
/usr/share/qt4/mkspecs/features/release.prf:
/usr/share/qt4/mkspecs/features/default_post.prf:
/usr/share/qt4/mkspecs/features/warn_on.prf:
/usr/share/qt4/mkspecs/features/qt.prf:
/usr/share/qt4/mkspecs/features/unix/thread.prf:
/usr/share/qt4/mkspecs/features/moc.prf:
/usr/share/qt4/mkspecs/features/resources.prf:
/usr/share/qt4/mkspecs/features/uic.prf:
/usr/share/qt4/mkspecs/features/yacc.prf:
/usr/share/qt4/mkspecs/features/lex.prf:
/usr/share/qt4/mkspecs/features/include_source_dir.prf:
/usr/lib/libQtCore.prl:
qmake:  FORCE
	@$(QMAKE) -spec /usr/share/qt4/mkspecs/linux-g++ -unix -o Makefile papaya.pro

dist: 
	@$(CHK_DIR_EXISTS) .tmp/papaya1.0.0 || $(MKDIR) .tmp/papaya1.0.0 
	$(COPY_FILE) --parents $(SOURCES) $(DIST) .tmp/papaya1.0.0/ && $(COPY_FILE) --parents utils/text_utils.h utils/settings.h utils/node.h utils/newick_reader.h utils/model_factory.h utils/fasta_reader.h utils/eigen.h utils/exceptions.h utils/db_matrix.h main/sequence.h main/simple_alignment.h utils/int_matrix.h utils/fasta_entry.h utils/settings_handle.h utils/xml_writer.h utils/evol_model.h main/reads_alignment.h .tmp/papaya1.0.0/ && $(COPY_FILE) --parents main.cpp utils/text_utils.cpp utils/settings.cpp utils/node.cpp utils/newick_reader.cpp utils/model_factory.cpp utils/fasta_reader.cpp utils/eigen.cpp utils/db_matrix.cpp main/sequence.cpp main/simple_alignment.cpp utils/int_matrix.cpp utils/settings_handle.cpp utils/xml_writer.cpp utils/evol_model.cpp main/reads_alignment.cpp .tmp/papaya1.0.0/ && (cd `dirname .tmp/papaya1.0.0` && $(TAR) papaya1.0.0.tar papaya1.0.0 && $(COMPRESS) papaya1.0.0.tar) && $(MOVE) `dirname .tmp/papaya1.0.0`/papaya1.0.0.tar.gz . && $(DEL_FILE) -r .tmp/papaya1.0.0


clean:compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) Makefile


mocclean: compiler_moc_header_clean compiler_moc_source_clean

mocables: compiler_moc_header_make_all compiler_moc_source_make_all

compiler_moc_header_make_all:
compiler_moc_header_clean:
compiler_rcc_make_all:
compiler_rcc_clean:
compiler_image_collection_make_all: qmake_image_collection.cpp
compiler_image_collection_clean:
	-$(DEL_FILE) qmake_image_collection.cpp
compiler_moc_source_make_all:
compiler_moc_source_clean:
compiler_uic_make_all:
compiler_uic_clean:
compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: 

####### Compile

main.o: main.cpp utils/settings.h \
		utils/settings_handle.h \
		utils/newick_reader.h \
		utils/exceptions.h \
		utils/node.h \
		main/sequence.h \
		utils/fasta_entry.h \
		main/simple_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/int_matrix.h \
		utils/model_factory.h \
		utils/fasta_reader.h \
		utils/xml_writer.h \
		main/reads_alignment.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp

text_utils.o: utils/text_utils.cpp utils/text_utils.h \
		utils/exceptions.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o text_utils.o utils/text_utils.cpp

settings.o: utils/settings.cpp utils/settings.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o settings.o utils/settings.cpp

node.o: utils/node.cpp utils/node.h \
		utils/exceptions.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		utils/fasta_entry.h \
		main/simple_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/int_matrix.h \
		utils/model_factory.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o node.o utils/node.cpp

newick_reader.o: utils/newick_reader.cpp utils/newick_reader.h \
		utils/exceptions.h \
		utils/node.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		utils/fasta_entry.h \
		main/simple_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/int_matrix.h \
		utils/model_factory.h \
		utils/text_utils.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o newick_reader.o utils/newick_reader.cpp

model_factory.o: utils/model_factory.cpp utils/model_factory.h \
		utils/db_matrix.h \
		utils/int_matrix.h \
		utils/evol_model.h \
		utils/settings.h \
		utils/settings_handle.h \
		utils/eigen.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o model_factory.o utils/model_factory.cpp

fasta_reader.o: utils/fasta_reader.cpp utils/fasta_reader.h \
		utils/exceptions.h \
		utils/node.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		utils/fasta_entry.h \
		main/simple_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/int_matrix.h \
		utils/model_factory.h \
		utils/text_utils.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o fasta_reader.o utils/fasta_reader.cpp

eigen.o: utils/eigen.cpp utils/eigen.h \
		utils/settings.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o eigen.o utils/eigen.cpp

db_matrix.o: utils/db_matrix.cpp utils/db_matrix.h \
		utils/settings.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o db_matrix.o utils/db_matrix.cpp

sequence.o: main/sequence.cpp utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		utils/fasta_entry.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o sequence.o main/sequence.cpp

simple_alignment.o: main/simple_alignment.cpp main/simple_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/int_matrix.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		utils/fasta_entry.h \
		utils/exceptions.h \
		utils/node.h \
		utils/model_factory.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o simple_alignment.o main/simple_alignment.cpp

int_matrix.o: utils/int_matrix.cpp utils/int_matrix.h \
		utils/settings.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o int_matrix.o utils/int_matrix.cpp

settings_handle.o: utils/settings_handle.cpp utils/settings_handle.h \
		utils/settings.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o settings_handle.o utils/settings_handle.cpp

xml_writer.o: utils/xml_writer.cpp utils/xml_writer.h \
		utils/exceptions.h \
		utils/node.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		utils/fasta_entry.h \
		main/simple_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/int_matrix.h \
		utils/model_factory.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o xml_writer.o utils/xml_writer.cpp

evol_model.o: utils/evol_model.cpp utils/evol_model.h \
		utils/db_matrix.h \
		utils/int_matrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o evol_model.o utils/evol_model.cpp

reads_alignment.o: main/reads_alignment.cpp main/reads_alignment.h \
		utils/settings.h \
		utils/settings_handle.h \
		utils/model_factory.h \
		utils/db_matrix.h \
		utils/int_matrix.h \
		utils/evol_model.h \
		utils/fasta_entry.h \
		utils/fasta_reader.h \
		utils/exceptions.h \
		utils/node.h \
		main/sequence.h \
		main/simple_alignment.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o reads_alignment.o main/reads_alignment.cpp

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:

