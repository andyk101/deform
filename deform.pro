#INCLUDEPATH *= /usr/include/c++/4.6

DEV_PATH = /home/andyk/development
INCLUDEPATH *= $${DEV_PATH}/boost
DEPENDPATH  *= $${DEV_PATH}/boost
#LIBS        *= -L$${DEV_PATH}/boost/stage/lib -lboost_system
#LIBS        *= -L$${DEV_PATH}/boost/stage/lib -lboost_filesystem
LIBS        *= -L$${DEV_PATH}/boost/stage/lib -lboost_serialization
#LIBS        *= -lboost_serialization
#LIBS        *= -L$${DEV_PATH}/boost/stage/lib -lboost_wserialization

INCLUDEPATH *= $${DEV_PATH}/andyk
DEPENDPATH  *= $${DEV_PATH}/andyk
LIBS        *= -L$${DEV_PATH}/andyk/lib -landyk-core
LIBS        *= -L$${DEV_PATH}/andyk/lib -landyk-serialization

INCLUDEPATH *= $${DEV_PATH}/qwt/src
DEPENDPATH  *= $${DEV_PATH}/qwt/src
LIBS        *= -L$${DEV_PATH}/qwt/lib -lqwt

#INCLUDEPATH *= /home/andyk/development/qt/src/testlib
#DEPENDPATH  *= /home/andyk/development/qt/src/testlib

TEMPLATE = app

QT += sql \
    xml

CONFIG += warn_on \
    qtestlib \
    debug

SOURCES += deform_app.cpp \
    deform_dlg.cpp \
#    deform_model.cpp \
#    deform_plot.cpp \
#    deform_table.cpp \
#    deform_test.cpp \
    detail.cpp

HEADERS += deform_app.h \
    deform_dlg.h \
#    deform_model.h \
#    deform_plot.h \
#    deform_table.h \
    detail.h

OTHER_FILES += read_me.txt
DESTDIR = /home/andyk/Dropbox/QT/exe
TARGET = deform
