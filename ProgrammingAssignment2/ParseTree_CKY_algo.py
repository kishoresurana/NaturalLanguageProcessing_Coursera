import os
import sys
import re
import json

_RARE_ = '_RARE_'
word_cnt = dict()
q_param = dict()
q_terminal_param = dict()

def GetWordCountFromUnaryTree(node):
    if len(node) == 2:
        word = node[1]

        if word in word_cnt:
            word_cnt[word] += 1
        else:
            word_cnt[word] = 1
    else:
        GetWordCountFromUnaryTree(node[1])
        GetWordCountFromUnaryTree(node[2])
    return

def BuildWordCounts(input):
    train = open(input, "r")

    for line in train:
        tree = json.loads(line)
        GetWordCountFromUnaryTree(tree[1])
        GetWordCountFromUnaryTree(tree[2])

def GetNormalizedWord(word):
    if word not in word_cnt or word_cnt[word] < 5:
        return _RARE_
    else:
        return word

def NormalizeUnaryTree(node):
    if len(node) == 2 and type(node[1]).__name__ != 'list':
        word = node[1]
        node[1] = GetNormalizedWord(word)
    else:
        for idx, child in enumerate(node):
            if idx == 0:
                continue
            NormalizeUnaryTree(child)
    return

def NormalizeData(input, output):
    train_norm = open(output, "w")
    train = open(input, "r")

    for line in train:
        tree = json.loads(line)
        NormalizeUnaryTree(tree)
        train_norm.write("%s\n" %(json.dumps(tree)))

    train_norm.close()

def NormalizeDevData(input, output):
    train_norm = open(output, "w")
    train = open(input, "r")

    for line in train:
        words = line.split()
        for idx, word in enumerate(words):
            words[idx] = GetNormalizedWord(word)
        train_norm.write("%s\n" %(' '.join(words)))

    train_norm.close()

def Build_qparams(input):
    train_counts = open(input, "r")
    xyzrule_dict = dict()
    xrule_dict = dict()
    xterminal_dict = dict()
    xterminal_agg_dict = dict()

    #column indexes in input data
    wordcount_idx = 0
    ruletype_idx = 1
    xrule_idx = 2
    yrule_idx = 3
    zrule_idx = 4
    terminal_idx = 3

    for row in train_counts:
        cols = row.split()
        if cols[ruletype_idx] == 'NONTERMINAL':
            xrule_dict[str.format('{0}'.format(cols[xrule_idx]))] = int(cols[wordcount_idx])
        elif cols[ruletype_idx] == 'BINARYRULE':
            xyzrule_dict[str.format('{0} {1} {2}'.format(cols[xrule_idx], cols[yrule_idx], cols[zrule_idx]))] = int(cols[wordcount_idx])
        elif cols[ruletype_idx] == 'UNARYRULE':
            xterminal_dict[str.format('{0} {1}'.format(cols[xrule_idx], cols[terminal_idx]))] = int(cols[wordcount_idx])
            if (cols[xrule_idx] not in xterminal_agg_dict):
                xterminal_agg_dict[cols[xrule_idx]] = int(cols[wordcount_idx])
            else:
                xterminal_agg_dict[cols[xrule_idx]] += int(cols[wordcount_idx])

    q_param_file = open("q_param", "w")
    for xyzrule, cnt in xyzrule_dict.iteritems():
        rules = xyzrule.split()
        xrule = rules[0]
        yrule = rules[1]
        zrule = rules[2]
        yzrule = "{0} {1}".format(yrule, zrule)
        xcnt = xrule_dict[xrule]
        if (xrule not in q_param):
            q_param[xrule] = dict()
        q_param[xrule][yzrule] = (cnt * 1.0) / xcnt
        q_param_file.write("{0} {1}\n".format(xyzrule, q_param[xrule][yzrule]))

    for xterminal, cnt in xterminal_dict.iteritems():
        cols = xterminal.split()
        xrule = cols[0]
        word = cols[1]
        if word not in q_terminal_param:
            q_terminal_param[word] = dict()
        xterminal_agg_cnt = xterminal_agg_dict[xrule]
        q_terminal_param[word][xrule] = (cnt * 1.0) / xterminal_agg_cnt

def appnd2(u, v):
    return (str(u) + ' ' + str(v))

def appnd3(w, u, v):
    return (w + ' ' + u + ' ' + v)

def Build_CKYParseTree(input, input_norm, output):
    dev = open(input, "r").readlines()
    dev_norm = open(input_norm, "r").readlines()
    dev_p1_out = open(output, "w")

    #for each line in file
    for lineidx, line in enumerate(dev_norm):
        words_norm = line.split()
        words = dev[lineidx].split()
        cky_pi = dict()

        #initialize cky_pi['1 1']['N'],  cky_pi['1 1']['VP'], etc..
        for (idx, word) in enumerate(words_norm):
            key = appnd2(idx+1, idx+1)

            if key not in cky_pi:
                cky_pi[key] = dict()

            for rule, q in q_terminal_param[word].iteritems():
                cky_pi[key][rule] = dict()
                cky_pi[key][rule]['prob'] = q
                cky_pi[key][rule]['terminalword'] = words[idx]
                cky_pi[key][rule]['terminalrule'] = rule

        for wordidx in range1(1, len(words_norm)-1):
            for I in range1(1, len(words_norm)-wordidx):
                J = I + wordidx
                key = appnd2(I, J)

                if key not in cky_pi:
                    cky_pi[key] = dict()

                for xrule, yzrules in q_param.iteritems():
                    curr_cky_pi = 0.0
                    max_cky_pi = 0.0
                    maxrule = ''
                    maxs = 0

                    for yzrule, ruleprob in yzrules.iteritems():
                        rule = yzrule.split()
                        yrule = rule[0]
                        zrule = rule[1]

                        for s in range1(I, J-1):

                            if yrule in cky_pi[appnd2(I, s)] and zrule in cky_pi[appnd2(s+1, J)]:
                                curr_cky_pi = ruleprob * cky_pi[appnd2(I, s)][yrule]['prob'] * cky_pi[appnd2(s+1, J)][zrule]['prob']

                                if (curr_cky_pi >= max_cky_pi):
                                    max_cky_pi = curr_cky_pi
                                    maxrule = yzrule
                                    maxs = s

                    if max_cky_pi > 0.0:
                        cky_pi[key][xrule] = dict()
                        cky_pi[key][xrule]['prob'] = max_cky_pi
                        cky_pi[key][xrule]['yzrule'] = maxrule
                        cky_pi[key][xrule]['s'] = maxs

        parsetree = GetBestParseTree(cky_pi, len(words))
        dev_p1_out.write("{0}\n".format(parsetree))
        dev_p1_out.flush()
    dev_p1_out.close()

def GetBestParseTree(cky_pi, sentencelength):
    key = appnd2(1, sentencelength)
    max_xrule = ''
    maxprob = 0.0
    currprob = 0.0

    #Assumes multi-word sentences
    for xrule, mapping in cky_pi[key].iteritems():
        currprob = mapping['prob']
        if currprob > maxprob:
            maxprob = currprob
            max_xrule = xrule

    return '["{0}", {1}]'.format(max_xrule, GetParseTreeFragment(cky_pi, 1, sentencelength, max_xrule))
            #return GetParseTreeFragment(cky_pi, 1, sentencelength, max_xrule)

def GetParseTreeFragment(cky_pi, slicestart, sliceend, xrule):
    rulefragment = cky_pi[appnd2(slicestart, sliceend)][xrule]
    if 'terminalrule' in rulefragment:
        return '["{0}", "{1}"]'.format(rulefragment['terminalrule'], rulefragment['terminalword'])
    else:
        yzrule = rulefragment['yzrule'].split()
        yrule = yzrule[0]
        zrule = yzrule[1]
        s = int(rulefragment['s'])
        
        if slicestart == s:
            yfragment = GetParseTreeFragment(cky_pi, slicestart, s, yrule)
        else:
            yfragment = '["{0}", {1}]'.format(yrule, GetParseTreeFragment(cky_pi, slicestart, s, yrule))

        if s+1 == sliceend:
            zfragment = GetParseTreeFragment(cky_pi, s+1, sliceend, zrule)
        else:
            zfragment = '["{0}", {1}]'.format(zrule, GetParseTreeFragment(cky_pi, s+1, sliceend, zrule))
        
        return '{1}, {2}'.format(xrule, yfragment, zfragment)

range1 = lambda start, end: range(start, end+1)
train = "parse_train.dat"
#train = "tree.example"
train_norm = "parse_train.dat.norm"
train_norm_counts = "cfg.counts.norm"

BuildWordCounts(train)
NormalizeData(train, train_norm)

#Manually load "cfg.counts.norm" by running command - python count_cfg_freq.py parse_train.dat.norm > cfg.counts.norm
Build_qparams(train_norm_counts)

dev_sample = "parse_sample.dat"
dev_sample_norm = "parse_sample.dat.norm"
dev_sample_out = "parse_sample.out"
NormalizeDevData(dev_sample, dev_sample_norm)

Build_CKYParseTree(dev_sample, dev_sample_norm, dev_sample_out)

dev = "parse_dev.dat"
dev_norm = "parse_dev.dat.norm"
dev_out = "parse_dev.out"
NormalizeDevData(dev, dev_norm)

Build_CKYParseTree(dev, dev_norm, dev_out)