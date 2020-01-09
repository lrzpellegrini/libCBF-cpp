cmake_minimum_required(VERSION 3.14)
project(CBF)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(linux)

find_package(OpenSSL REQUIRED)

add_library(libCBF
        linux/libexport.h
        linux/lindef.h
        base64.cpp
        base64.h
        end.cpp
        end.h
        cbf.cpp
        cbf.h
        cbflib.h)


target_link_libraries(libCBF OpenSSL::SSL)

add_executable(appCBF test-app/test-app-cbf.cpp)
target_link_libraries(appCBF OpenSSL::SSL libCBF)
