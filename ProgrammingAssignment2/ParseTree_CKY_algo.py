import os
import sys
import re
import json

_RARE_ = '_RARE_'
word_cnt = dict()

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

train = "parse_train.dat"
#train = "tree.example"
train_norm = "parse_train.dat.norm"
train_norm_counts = "cfg.counts.norm"

BuildWordCounts(train)
NormalizeData(train, train_norm)

#Manually load "cfg.counts.norm" by running command - python count_cfg_freq.py parse_train.dat.norm > cfg.counts.norm

dev = "parse_dev.dat"
dev_norm = "parse_dev.dat.norm"
NormalizeDevData(dev, dev_norm)

#Manually load "train_counts" by running command - python count_freqs.py gene.train > gene.counts
train_counts = "gene.counts"