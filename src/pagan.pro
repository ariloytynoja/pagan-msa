# -------------------------------------------------
# Project created by QtCreator 2009-05-05T14:26:31
# -------------------------------------------------
QT -= gui
TARGET = pagan

CONFIG = debug
CONFIG += console
CONFIG -= app_bundle
TEMPLATE = app
SOURCES += main.cpp \
    utils/text_utils.cpp \
    utils/settings.cpp \
    utils/newick_reader.cpp \
    utils/model_factory.cpp \
    utils/fasta_reader.cpp \
    utils/eigen.cpp \
    utils/db_matrix.cpp \
    main/sequence.cpp \
    utils/int_matrix.cpp \
    utils/settings_handle.cpp \
    utils/xml_writer.cpp \
    utils/evol_model.cpp \
    main/node.cpp \
    main/reference_alignment.cpp \
    main/viterbi_alignment.cpp \
    main/basic_alignment.cpp \
    main/reads_aligner.cpp
HEADERS += utils/text_utils.h \
    utils/settings.h \
    utils/newick_reader.h \
    utils/model_factory.h \
    utils/fasta_reader.h \
    utils/eigen.h \
    utils/exceptions.h \
    utils/db_matrix.h \
    main/sequence.h \
    utils/int_matrix.h \
    utils/fasta_entry.h \
    utils/settings_handle.h \
    utils/xml_writer.h \
    utils/evol_model.h \
    main/node.h \
    main/reference_alignment.h \
    main/viterbi_alignment.h \
    main/basic_alignment.h \
    main/reads_aligner.h
LIBS += -lboost_program_options-mt
INCLUDEPATH += /usr/include
OTHER_FILES += missing_things.txt
