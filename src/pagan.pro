# -------------------------------------------------
# Project created by QtCreator 2009-05-05T14:26:31
# -------------------------------------------------
QT -= gui
TARGET = pagan

#CONFIG = release
CONFIG = debug

CONFIG += console
CONFIG -= app_bundle
TEMPLATE = app
SOURCES += main.cpp \
    main/node.cpp \
    main/sequence.cpp \
    main/reference_alignment.cpp \
    main/viterbi_alignment.cpp \
    main/basic_alignment.cpp \
    main/reads_aligner.cpp \
    utils/text_utils.cpp \
    utils/settings.cpp \
    utils/newick_reader.cpp \
    utils/model_factory.cpp \
    utils/fasta_reader.cpp \
    utils/eigen.cpp \
    utils/db_matrix.cpp \
    utils/int_matrix.cpp \
    utils/settings_handle.cpp \
    utils/xml_writer.cpp \
    utils/evol_model.cpp \
    utils/check_version.cpp \
    utils/exonerate_queries.cpp \
    utils/log_output.cpp \
    utils/find_anchors.cpp \
    utils/codon_translation.cpp \
    utils/input_output_parser.cpp \
    utils/mafft_alignment.cpp \
    utils/raxml_tree.cpp \
    utils/tree_node.cpp \
    utils/bppdist_tree.cpp \
    utils/bppphysamp_tree.cpp \
    utils/bppancestors.cpp
HEADERS += utils/text_utils.h \
    main/node.h \
    main/sequence.h \
    main/reference_alignment.h \
    main/viterbi_alignment.h \
    main/basic_alignment.h \
    main/reads_aligner.h \
    utils/settings.h \
    utils/newick_reader.h \
    utils/model_factory.h \
    utils/fasta_reader.h \
    utils/eigen.h \
    utils/exceptions.h \
    utils/db_matrix.h \
    utils/int_matrix.h \
    utils/fasta_entry.h \
    utils/settings_handle.h \
    utils/xml_writer.h \
    utils/evol_model.h \
    utils/check_version.h \
    utils/exonerate_queries.h \
    utils/log_output.h \
    utils/find_anchors.h \
    utils/substring_hit.h \
    utils/codon_translation.h \
    utils/input_output_parser.h \
    utils/mafft_alignment.h \
    utils/raxml_tree.h \
    utils/tree_node.h \
    utils/bppdist_tree.h \
    utils/bppphysamp_tree.h \
    utils/bppancestors.h
LIBS += -lboost_program_options -lboost_regex -lboost_thread  -lboost_system -lrt -lgomp -lcurl
INCLUDEPATH += /usr/include
OTHER_FILES +=  ../VERSION_HISTORY

QMAKE_CXXFLAGS += -fopenmp
QMAKE_CFLAGS += -fopenmp




