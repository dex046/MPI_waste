TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    RWsgy.cpp \
    testmain.cpp \
    Partition.cpp \
    testTDFWI.cpp \
    DataTran.cpp \
    h_Coord.cpp \
    Inside.cpp

INCLUDEPATH += /usr/include/mpich-x86_64/
DEPENDPATH  += /usr/include/mpich-x86_64/


include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    RWsgy.h \
    Partition.h \
    DataTran.h \
    testTDFWI.h \
    h_Coord.h \
    Inside.h

