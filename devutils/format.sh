#!/usr/bin/bash 

find ./apps -iname "*.cpp" -o -iname "*.hpp"  -o -iname "*.h" | xargs clang-format -i --style="file:.clang-format"