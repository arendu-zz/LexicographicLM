__author__ = 'arenduchintala'
import fst, re, pprint, pdb

INITIAL = '_initial_'
NULL = '_null_'


def trysplit(gramline):
    try:
        [p, ng, b] = gramline.split('\t')
        return -float(p), ng, -float(b)
    except:
        [p, ng] = gramline.split('\t')
        return -float(p), ng, None


def add_arc_pr(sym, lmfst, fid, tid, isy, osy, wt):
    lmfst.add_arc(fid, tid, isy, osy, wt)
    '''
    print 'added arc', fid, tid, sym[isy], sym[osy]
    if fid == 1961:
        count = 0
        for c, s in enumerate(lmfst.states):
            count = c
        print 'number of total states', count
    '''


if __name__ == '__main__':
    sym_e = fst.read_symbols('data/syme.bin')
    pdb.set_trace()
    lm_txt = open('data/lm', 'r').read()
    [bs, unigrams, bigrams, trigrams] = re.split('1-grams:|2-grams:|3-grams:', lm_txt)
    unigrams = re.split('\n+', unigrams)
    bigrams = re.split('\n+', bigrams)
    trigrams = re.split('\n+', trigrams)

    lm_id = {}
    lm_id[INITIAL] = len(lm_id)
    lm_fst = fst.Transducer(sym_e, sym_e)
    lm_fst.add_state()
    lm_id[NULL] = len(lm_id)

    for uni_line in unigrams:
        if uni_line.strip() != '' and len(uni_line.split('\t')) > 1:
            [p, ng, bk] = trysplit(uni_line)
            if ng not in lm_id:
                lm_fst.add_state()
            lm_id[ng] = lm_id.get(ng, len(lm_id))

            from_id = lm_id[NULL]
            to_id = lm_id[ng]
            add_arc_pr(sym_e, lm_fst, from_id, to_id, ng, ng, p)

            if bk is not None:
                add_arc_pr(sym_e, lm_fst, to_id, lm_id[NULL], sym_e.find(0), sym_e.find(0), bk)
            else:
                add_arc_pr(sym_e, lm_fst, to_id, lm_id[NULL], sym_e.find(0), sym_e.find(0), 0.0)

            if ng == '</s>':
                lm_fst[lm_id[ng]].final = True

    for bi_line in bigrams:
        if bi_line.strip() != '' and len(bi_line.split('\t')) > 1:
            [p, ng, bk] = trysplit(bi_line)
            [ng1, ng2] = ng.split()
            if (ng1, ng2) not in lm_id:
                lm_fst.add_state()
            lm_id[ng1, ng2] = lm_id.get((ng1, ng2), len(lm_id))
            from_id = lm_id[ng1]
            to_id = lm_id[ng1, ng2]

            add_arc_pr(sym_e, lm_fst, from_id, to_id, ng2, ng2, p)

            if bk is not None:
                add_arc_pr(sym_e, lm_fst, to_id, lm_id[ng2], sym_e.find(0), sym_e.find(0), bk)
            else:
                add_arc_pr(sym_e, lm_fst, to_id, lm_id[ng2], sym_e.find(0), sym_e.find(0), 0.0)

            if ng2 == '</s>':
                lm_fst[lm_id[ng1, ng2]].final = True

    for tri_line in trigrams:
        if tri_line.strip() != '' and len(tri_line.split('\t')) > 1:
            [p, ng, bk] = trysplit(tri_line)
            [ng1, ng2, ng3] = ng.split()
            if (ng2, ng3) not in lm_id:
                lm_fst.add_state()
            lm_id[ng2, ng3] = lm_id.get((ng2, ng3), len(lm_id))
            from_id = lm_id[ng1, ng2]
            to_id = lm_id[ng2, ng3]
            add_arc_pr(sym_e, lm_fst, from_id, to_id, ng3, ng3, p)
            if ng3 == '</s>':
                lm_fst[lm_id[ng2, ng3]].final = True

    pprint.pprint(lm_id)
    #connect init to start

    add_arc_pr(sym_e, lm_fst, lm_id[INITIAL], lm_id[NULL], sym_e.find(0), sym_e.find(0), 0.0)
    lm_fst.write('data/pylm.fst', sym_e, sym_e)
    """
    99 <s> -1.640621
     9 -2.411867   </s>

    for l in lm_lines:
        if l.strip() != '':
    """
