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
        xrule = xyzrule.split()[0]
        xcnt = xrule_dict[xrule]
        q_param[xyzrule] = (cnt * 1.0) / xcnt
        q_param_file.write("{0} {1}\n".format(xyzrule, q_param[xyzrule]))

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
        words = line.split()
        cky_pi = dict()

        #initialize cky_pi['1 1']['N'],  cky_pi['1 1']['VP'], etc..
        for (idx, word) in enumerate(words):
            key = appnd2(idx+1, idx+1)

            if key not in cky_pi:
                cky_pi[key] = dict()

            for rule, q in q_terminal_param[word].iteritems():
                cky_pi[key][rule] = q

        for l in range1(1, len(words)-1):
            for i in range1(1, len(words)-l):
                j = i + l
                key = appnd2(i, j)

                if key not in cky_pi:
                    cky_pi[key] = dict()

                curr_cky_pi = 0.0
                max_cky_pi = 0.0
                maxrule = ''
                maxs = 0

                for xyzrule, ruleprob in q_param.iteritems():
                    rule = xyzrule.split()
                    xrule = rule[0]
                    yrule = rule[1]
                    zrule = rule[2]

                    for s in range1(i, j-1):

                        if yrule in cky_pi[appnd2(i, s)] and zrule in cky_pi[appnd2(s+1, j)]:
                            curr_cky_pi = ruleprob * cky_pi[appnd2(i, s)][yrule] * cky_pi[appnd2(s+1, j)][zrule]

                            if (curr_cky_pi >= max_cky_pi):
                                max_cky_pi = curr_cky_pi
                                maxrule = xyzrule
                                maxs = s

                    cky_pi[key][xrule] = max_cky_pi


range1 = lambda start, end: range(start, end+1)
train = "parse_train.dat"
#train = "tree.example"
train_norm = "parse_train.dat.norm"
train_norm_counts = "cfg.counts.norm"

BuildWordCounts(train)
NormalizeData(train, train_norm)

#Manually load "cfg.counts.norm" by running command - python count_cfg_freq.py parse_train.dat.norm > cfg.counts.norm
Build_qparams(train_norm_counts)

dev = "parse_dev.dat"
dev_norm = "parse_dev.dat.norm"
dev_out = "parse_dev.out"
NormalizeDevData(dev, dev_norm)

Build_CKYParseTree(dev, dev_norm, dev_out)