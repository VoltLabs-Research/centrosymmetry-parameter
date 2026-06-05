cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(centrosymmetry VERSION 1.0.0 LANGUAGES CXX)

find_package(coretoolkit REQUIRED)

volt_add_plugin(centrosymmetry
    SOURCES src/centrosymmetry_engine.cpp
    MAIN_SOURCE src/plugin.cpp
)
