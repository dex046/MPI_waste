TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

SOURCES += \
    RWsgy.cpp \
    testmain.cpp \
    Partition.cpp \
    testTDFWI.cpp \
    DataTran.cpp \
    h_Coord.cpp \
    Inside.cpp \
    H_Border.cpp

INCLUDEPATH += /usr/include/mpich-x86_64/
DEPENDPATH  += /usr/include/mpich-x86_64/

LIBS += -L/lib64/mpich/lib/
include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    RWsgy.h \
    Partition.h \
    DataTran.h \
    testTDFWI.h


unix|win32: LIBS += -lfmpich

unix|win32: LIBS += -lmpich

unix|win32: LIBS += -lmpichcxx

unix|win32: LIBS += -lmpichf90

unix|win32: LIBS += -lmpl

unix|win32: LIBS += -lopa
