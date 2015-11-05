TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

SOURCES += main.cpp \
    RWsgy.cpp \
    testmain.cpp \
    Partition.cpp \
    testTDFWI.cpp \
    DataTran.cpp \
    h_Coord.cpp \
    Inside.cpp \
    H_U_VW_Vp.cpp \
    H_Border.cpp

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

