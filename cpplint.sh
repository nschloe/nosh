CPPLINT=~/software/cpplint/cpplint/cpplint.py

${CPPLINT} \
  --extensions=hpp,cpp \
  --filter=-whitespace/parens,-whitespace/braces,-whitespace/line_length,-whitespace/comments,-runtime/references,-build/include_order,-readability/todo \
  \
  nosh/*.cpp

~/software/cpplint/cpplint/cpplint.py
${CPPLINT} \
 --filter=-whitespace/parens,-whitespace/braces,-whitespace/line_length,-whitespace/comments,-runtime/references,-build/include_order,-readability/todo,-whitespace/indent \
  --extensions=cpp,hpp nosh/*.hpp
