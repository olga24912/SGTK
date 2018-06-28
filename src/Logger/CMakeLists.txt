include_directories(.)

set(SOURCE_FILES logger_impl.cpp logger.hpp perfcounter.hpp
        log_writers.hpp path_helper.cpp path_helper.hpp
        stacktrace.hpp verify.hpp
        format.cc format.h
        posix.cc posix.h)

add_library(logger_lib ${SOURCE_FILES})