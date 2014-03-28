# coding: utf-8
__author__ = 'arenduchintala'

if __name__ == '__main__':
    lm = open('data/lm', 'r').readlines()
    sym_e = set(["<unk>","<bk>"])
    for l in lm:
        ls = l.split('\t')
        if len(ls) == 2 or len(ls) == 3:
            sym_e.update(ls[1].split())

    sym_w = open('data/syme.sym', 'w')
    sym_o = [str(sym.strip() + ' ' + str(id + 1)).strip() for id, sym in enumerate(sym_e)]
    sym_o.insert(0, '<eps> 0')
    sym_w.write('\n'.join(sym_o))
    sym_w.flush()
    sym_w.close()

