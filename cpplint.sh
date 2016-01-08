#!/bin/sh

# sudo -H pip install cpplint

cpplint \
  --extensions=hpp,cpp \
  --filter=-whitespace/parens,-whitespace/braces,-whitespace/line_length,-whitespace/comments,-runtime/references,-build/include_order,-readability/todo \
  \
  src/*.cpp

cpplint \
 --filter=-whitespace/parens,-whitespace/braces,-whitespace/line_length,-whitespace/comments,-runtime/references,-build/include_order,-readability/todo,-whitespace/indent \
  --extensions=cpp,hpp \
  src/*.hpp
