#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
//OpenFST
#include <fst/fstlib.h>
#include <fst/fst-decl.h>
#include <fst/map.h>

//Boost
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
using namespace std;
using namespace fst;
using namespace boost;

typedef fst::LexicographicWeight<fst::TropicalWeight, fst::TropicalWeight> LexWeight;
typedef fst::LexicographicArc<fst::TropicalWeight, fst::TropicalWeight> LexArc;

int main(int argc, char *argv[]) {
    fst::VectorFst<LexArc> l_fst;
    SymbolTable *syms = fst::SymbolTable::ReadText("data/syme.sym");
    l_fst.SetInputSymbols(syms);
    l_fst.SetOutputSymbols(syms);
    int start_id = l_fst.AddState();
    l_fst.SetStart(start_id);
    int end_id = l_fst.AddState();
    //creating a lexicographic weight
    LexWeight lex_w = LexWeight(0.0, 1.0);
    //creating a lexicographic arc using the lexicographic weight
    LexArc lex_arc = LexArc(0,0, lex_w, end_id);
    //adding arc to fst
    l_fst.AddArc(start_id, lex_arc);
    
    LexWeight final_w = LexWeight(0.0, 0.0);
    l_fst.SetFinal(end_id, final_w);
    
    l_fst.Write("simple_lex.fst");
}


