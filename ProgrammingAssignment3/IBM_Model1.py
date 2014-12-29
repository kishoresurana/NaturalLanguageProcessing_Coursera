import os
import sys
import re
import json
from collections import defaultdict
import pickle

t_param = defaultdict(float)
ignore_empty_lines = set()
__NULL__ = '__NULL__'

def save_t_params(tfile):
    d = dict(t_param)
    with open(tfile,'wb') as f:
        pickle.dump(d,f)

def load_t_params(tfile):
    with open(tfile,'r') as f:
        d = pickle.load(f)
        #t_param = defaultdict(lambda: defaultdict(float))
        t_param.update(d)

def Initialize_t_param(input_lines_from, input_lines_to):
    #load a from_word_dict: contains from-word, and the sentence idx list
    to_word_dict = defaultdict(set)
    all_from_words = set()
    for idx, to_line in enumerate(input_lines_to):
        to_words = to_line.split()
        if (len(to_words) == 0):
            ignore_empty_lines.add(idx)
            continue;
        for to_word in to_words:
            to_word_dict[to_word].add(idx)

    #iterate from_word_dict:
        # iterate corresponding sentences (using idx) and get unique_from_words
        # initialize t_param('to_word from_word') = 1/count(unique_from_words)
    for to_word in to_word_dict:
        unique_from_words = set()
        for idx in to_word_dict[to_word]:
            from_words = input_lines_from[idx].split()
            if (len(from_words) == 0):
                ignore_empty_lines.add(idx)
                continue;
            for from_word in from_words:
                all_from_words.add(from_word)
                unique_from_words.add(from_word)

        for unique_from_word in unique_from_words:
            t_param['{0} {1}'.format(unique_from_word, to_word)] = 1.0/len(unique_from_words)

    len_all_from_words = len(all_from_words)
    for from_word in all_from_words:
        t_param['{0} {1}'.format(from_word, __NULL__)] = 1.0/len_all_from_words

def Build_t_param(input_from, input_to, output, reuse_t_param = 1):
    input_lines_from = open(input_from, "r").readlines()
    input_lines_to = open(input_to, "r").readlines()

    t_param_file = "t_param_initial.dat"
    if not os.path.isfile(t_param_file) or reuse_t_param == 0:
        Initialize_t_param(input_lines_from, input_lines_to)
        save_t_params(t_param_file)
    else:
        load_t_params(t_param_file)

    #Test fo a sample  if totalpct is 1.0
    totalpct = 0.0
    for t in t_param:
        if (t.endswith(' possibly')):
            totalpct += t_param[t]

    totalpct = 0.0
    for t in t_param:
        if (t.endswith(__NULL__)):
            totalpct += t_param[t]

    iterations = 5
    for idx in range(iterations):
        test = 'test'
        c_to_from = defaultdict(lambda: defaultdict(float))
        #c_to_from = defaultdict(float)
        c_to = defaultdict(float)
        c_j_i_l_m = defaultdict(float)
        c_i_l_m = defaultdict(float)
        delta_k_i_j_dict = dict()

        for k, line in enumerate(input_lines_from):
            if k in ignore_empty_lines:
                continue
            words_from = line.split()
            words_to = input_lines_to[k].split()
            words_to.insert(0, __NULL__)
            l = len(words_to)
            m = len(words_from)

            delta_denominator = defaultdict(float)
            for i, word_f in enumerate(words_from):
                for j, word_t in enumerate(words_to):
                    delta_denominator[word_f] += t_param["{0} {1}".format(word_f, word_t)]

            for i, word_f in enumerate(words_from):
                for j, word_t in enumerate(words_to):
                    #delta_denominator = 0.0
                    fromto = "{0} {1}".format(word_f, word_t)
                    delta_numerator = t_param[fromto]
                    #for dj_word_t in words_to:
                    #    delta_denominator += t_param["{0} {1}".format(word_f, dj_word_t)]

                    delta_k_i_j = (delta_numerator * 1.0)/delta_denominator[word_f]
                    #delta_k_i_j_dict["{0} {1} {2}".format(k, i, j)] = delta_k_i_j

                    c_to_from[word_t][word_f] += delta_k_i_j
                    #c_to_from[fromto] += delta_k_i_j
                    c_to[word_t] += delta_k_i_j
                    c_j_i_l_m["{0} {1} {2} {3}".format(j, i, l, m)] += delta_k_i_j
                    c_i_l_m["{0} {1} {2}".format(i, l, m)] += delta_k_i_j

        for to_word in c_to_from:
            for from_word in c_to_from[to_word]:
                t_param["{0} {1}".format(from_word, to_word)] = (c_to_from[to_word][from_word] * 1.0) / c_to[to_word]

    t_param_file = "t_param_IBM_Model1.dat"
    save_t_params(t_param_file)

    output_file = open(output,'wb')
    for k, line in enumerate(input_lines_from):
        if k in ignore_empty_lines:
            continue
        words_from = line.split()
        words_to = input_lines_to[k].split()

        for i, word_f in enumerate(words_from):
            max_j = 0
            max_word_t = ''
            max_t_param = 0.0
            curr_t_param = 0.0
            
            for j, word_t in enumerate(words_to):
                curr_t_param = t_param["{0} {1}".format(word_f, word_t)]
                if curr_t_param > max_t_param:
                    max_t_param = curr_t_param
                    max_j = j
                    max_word_t = word_t
            
            if max_t_param > 0.0:
                output_file.write("{0} {1} {2}\n".format(k+1, i+1, max_j+1))
    output_file.close()
    test = 'hello'

#input_from = "corpus_sample.es"
#input_to = "corpus_sample.en"
#output = "alignment_sample.out"
#Build_t_param(input_from, input_to, output)

input_from = "corpus.es"
input_to = "corpus.en"
output = "alignment_test.p1.out"
Build_t_param(input_from, input_to, output)
test = 'hello'