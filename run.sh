echo "extracting syms..."
python extractSym.py
echo "building..."
g++ -lfst -lboost_regex -o makebackoff -g makeBackOffLM.cpp
echo "running..."
./makebackoff
