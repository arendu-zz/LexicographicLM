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

typedef vector<string> string_vector;

struct arparow {
    double logprob;
    string gram;
    double bkprob;
};

typedef vector<arparow> arpa_vector;

arparow make_arparow(string_vector vec) {
    arparow arow;
    char* pEnd;
    arow.logprob = -strtod(vec[0].c_str(), &pEnd);
    arow.gram = vec[1];
    if (vec.size() == 3) {
        arow.bkprob = -strtod(vec[2].c_str(), &pEnd);
    } else {
        arow.bkprob = 0.0;
    }
    return arow;
}

void push_arparow(string_vector vec, arpa_vector *gram) {
    if (vec.size() > 1) {
        arparow arow = make_arparow(vec);
        gram->push_back(arow);
    } else {
        cout << " ignoring row" << endl;
    }
}

typedef map<string, int> Dictionary;

int main(int argc, char *argv[]) {
    //ifstream ifs("data/news.small.20.tok.arpa");
    ifstream ifs("data/lm");
    string line;
    string_vector str_vec;
    int g = 0;
    arpa_vector unigram;
    arpa_vector bigram;
    arpa_vector trigram;
    while (getline(ifs, line)) {
        trim(line);
        cout << line << endl;
        if (iequals(line, "")) {
            cout << "uigrams:" << unigram.size() << " bigrams:" << bigram.size() << " trigrams:" << trigram.size() << endl;
        } else {
            if (iequals(line.c_str(), "\\1-grams:")) {
                g = 1;
            } else if (iequals(line.c_str(), "\\2-grams:")) {
                g = 2;
            } else  if (iequals(line.c_str(), "\\3-grams:")) {
                g = 3;
            } else {
                if (g == 1) {
                    split(str_vec, line, boost::is_any_of("\t"));
                    push_arparow(str_vec, &unigram);
                } else if (g == 2) {
                    split(str_vec, line, boost::is_any_of("\t"));
                    push_arparow(str_vec, &bigram);
                } else if (g == 3) {
                    split(str_vec, line, boost::is_any_of("\t"));
                    push_arparow(str_vec, &trigram);
                } else {
                    //ignore lines
                }
            }
        }
    }
    cout << "uigrams:" << unigram.size() << " bigrams:" << bigram.size() << " trigrams:" << trigram.size() << endl;

    SymbolTable *syms = fst::SymbolTable::ReadText("data/syme.sym");
    syms->Write("data/syme.bin");
    Dictionary lm_id;
    lm_id.insert(pair<string, int>(string("_INITIAL_"), lm_id.size()));
    lm_id.insert(pair<string, int>(string("_NULL_"), lm_id.size()));
    //Make the symbol table, this iteration over the unigrams can be
    //combined with the next iteration to speed up.
    //Currently left seperate for clarity
    int k = 1;
    for (vector<arparow>::iterator it = unigram.begin(); it != unigram.end(); ++it) {
        cout << it->gram << " gram " << endl;
        //syms->AddSymbol(it->gram.c_str(), k);
        k++;
    }
    //syms->AddSymbol("<s>", k++);
    //syms->AddSymbol("</s>", k++);
    //syms->Write("data/syme.bin");
    //syms->WriteText("data/syme.txt");
    fst::StdVectorFst lm_fst;
    lm_fst.SetStart(0);
    int initial_state = lm_fst.AddState();
    int null_state = lm_fst.AddState();
    cout << "_NULL_ state id:" << null_state << endl;
    cout << "_INITIAL_ state id:" << initial_state << endl;
    lm_fst.SetInputSymbols(syms);
    lm_fst.SetOutputSymbols(syms);
    for (vector<arparow>::iterator it = unigram.begin(); it != unigram.end(); ++it) {
        cout << it->gram << " is trying to be added" << endl;
        int from_id = lm_id.find(string("_NULL_"))->second;
        int to_id = -1;
        if (lm_id.find(it->gram) == lm_id.end()) {
            lm_id.insert(std::pair<string, int>(it->gram, lm_id.size()));
            to_id  = lm_fst.AddState();
            cout << to_id << " is the new state id" << endl;
        } else {
            //to get the value associated form a map you need ->second.
            //because find returns a position pointer
            to_id = lm_id.find(it->gram)->second;
            cout << " look up to id:" << to_id << endl;
        }
        cout << "from id:" << from_id << " to id" << to_id << endl;

        cout << "finding " << it->gram << syms->Find(it->gram);
        cout << "finding " << it->gram << syms->Find(it->gram.c_str());
        lm_fst.AddArc(from_id, fst::StdArc(syms->Find(it->gram.c_str()), syms->Find(it->gram.c_str()), it->logprob, to_id));
        lm_fst.AddArc(to_id, fst::StdArc(0, 0, it->bkprob, from_id));
        if (boost::iequals(it->gram, string("</s>"))) {
            lm_fst.SetFinal(to_id, 0.0);
        }
    }

    for (vector<arparow>::iterator it = bigram.begin(); it != bigram.end(); ++it) {
        cout << it->gram << " is trying to be added" << endl;
        string_vector bigram_vec;
        split(bigram_vec, it->gram, boost::is_any_of(" "), boost::token_compress_on);
        string ng1 = bigram_vec[0];
        string ng2 = bigram_vec[1];
        string ng12 = ng1 + string(" ") + ng2;
        int to_id = -1;
        if (lm_id.find(ng12) == lm_id.end()) {
            lm_id.insert(pair<string, int>(ng12, lm_id.size()));
            to_id = lm_fst.AddState();
        } else {
            to_id = lm_id.find(ng12)->second;
        }
        int from_id = lm_id.find(ng1)->second;
        int return_id = lm_id.find(ng2)->second;
        lm_fst.AddArc(from_id, fst::StdArc(syms->Find(ng2.c_str()), syms->Find(ng2.c_str()), it->logprob, to_id));
        lm_fst.AddArc(to_id, fst::StdArc(0, 0, it->bkprob, return_id));

        if (boost::iequals(ng2, "</s>")) {
            lm_fst.SetFinal(to_id, 0.0);
        }
    }

    for (vector<arparow>::iterator it = trigram.begin(); it != trigram.end(); ++it) {
        cout << it->gram << " is trying to be added" << endl;
        string_vector trigram_vec;
        split(trigram_vec, it->gram, boost::is_any_of(" "), boost::token_compress_on);
        string ng1 = trigram_vec[0];
        string ng2 = trigram_vec[1];
        string ng3 = trigram_vec[2];
        cout << "split:" << ng1 << " " << ng2 << " " << ng3 << endl;
        string ng12 = ng1 + string(" ") + ng2;
        string ng23 = ng2 + string(" ") + ng3;
        cout << "join:" << ng12 << " " << ng23 << endl;
        int to_id = -1;
        if (lm_id.find(ng23) == lm_id.end()) {
            cout << "could not find:" << ng23 << endl;
            lm_id.insert(pair<string, int>(ng23, lm_id.size()));
            to_id = lm_fst.AddState();
            cout << "new state added:" << to_id << endl;
        } else {
            to_id  = lm_id.find(ng23)->second;
        }
        int from_id = lm_id.find(ng12)->second;

        lm_fst.AddArc(from_id, fst::StdArc(syms->Find(ng3.c_str()), syms->Find(ng3.c_str()), it->logprob, to_id));
        //Wait a min! do we have to add an arc back to null?
        int bk_to_id = lm_id.find(ng3)->second;
        lm_fst.AddArc(to_id, fst::StdArc(0, 0, 0.0, bk_to_id));
        if (boost::iequals(ng3, string("</s>"))) {
            lm_fst.SetFinal(to_id, 0.0);
        }
    }
    cout << "COMPLETED GRAMS" << endl;
    int from_id = lm_id.find("_INITIAL_")->second;
    int to_id = lm_id.find("_NULL_")->second;
    lm_fst.AddArc(from_id, fst::StdArc(0, 0, 0.0, to_id));
    lm_fst.Write("data/cpplm.fst");
}





