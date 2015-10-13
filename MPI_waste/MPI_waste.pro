TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    RWsgy.cpp \
    TDFWI.cpp \
    testmain.cpp

INCLUDEPATH += /usr/include/mpich-x86_64/
DEPENDPATH  += /usr/include/mpich-x86_64/


include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    RWsgy.h \
    TDFWI.h

