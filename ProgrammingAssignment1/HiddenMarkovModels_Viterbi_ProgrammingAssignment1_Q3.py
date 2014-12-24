import os
import sys
import re

O_tag = 'O'
I_Gene_tag = 'I-GENE'
_RARE_ = '_RARE_'
_HASNUMBERS_ = '_HASNUMBERS_'
_ALLCAPITALS_ = '_ALLCAPITALS_'
_LASTCAPITAL_ = '_LASTCAPITAL_'

exv_O_dict = dict()
exv_I_Gene_dict = dict()
q_trigram = dict()

def NormalizeData(input, output, inputcnts):
    train_norm = open(output, "w")
    train = open(input, "r")
    train_counts = open(inputcnts, "r")
    cnt_dict = dict()

    #column indexes in inputcnts
    wordcount_idx = 0
    tag_idx = 2
    word_idx = 3

    #Read the records from gene.counts file
    for row in train_counts:
        cols = row.split()
        if cols[1] == 'WORDTAG':
            if cols[word_idx] in cnt_dict:
                cnt_dict[cols[word_idx]] += int(cols[wordcount_idx])
            else:
                cnt_dict[cols[word_idx]] = int(cols[wordcount_idx])

    for row in train:
        cols = row.split()
        if len(cols) == 0:
            train_norm.write('\n')
            continue
        
        word = cols[0]
        if word not in cnt_dict or cnt_dict[word] < 5:
            word = _RARE_

        if len(cols) > 1:
            train_norm.write("%s %s\n" %(word, cols[1]))
        else:
            train_norm.write("%s\n" %(word))

    train_norm.close()

def NormalizeData_Informational(input, output, inputcnts):
    train_norm = open(output, "w")
    train = open(input, "r")
    train_counts = open(inputcnts, "r")
    cnt_dict = dict()

    #column indexes in inputcnts
    wordcount_idx = 0
    tag_idx = 2
    word_idx = 3

    #Read the records from gene.counts file
    for row in train_counts:
        cols = row.split()
        if cols[1] == 'WORDTAG':
            if cols[word_idx] in cnt_dict:
                cnt_dict[cols[word_idx]] += int(cols[wordcount_idx])
            else:
                cnt_dict[cols[word_idx]] = int(cols[wordcount_idx])

    for row in train:
        cols = row.split()
        if len(cols) == 0:
            train_norm.write('\n')
            continue
        
        word = cols[0]
        if word not in cnt_dict or cnt_dict[word] < 5:
            if hasNumbers(word):
                word = _HASNUMBERS_
            elif allCaps(word):
                word = _ALLCAPITALS_
            elif lastCaps(word):
                word = _LASTCAPITAL_
            else:
                word = _RARE_

        if len(cols) > 1:
            train_norm.write("%s %s\n" %(word, cols[1]))
        else:
            train_norm.write("%s\n" %(word))

    train_norm.close()

def RegexMatches(inputString, regex):
    if re.search(regex, inputString) is not None:
        return True
    else:
        return False

def hasNumbers(inputString):
    return RegexMatches(inputString, r".*[0-9]+.*")

def allCaps(inputString):
    return RegexMatches(inputString, r'^[A-Z]+$')

def lastCaps(inputString):
    lastchar = inputString[-1]
    return allCaps(lastchar)

def Build_exv_viterbi_unigram(input):
    train_counts = open(input, "r")
    total_dict = dict()
    O_dict = dict()
    I_Gene_dict = dict()

    #column indexes in input data
    wordcount_idx = 0
    tag_idx = 2
    word_idx = 3

    for row in train_counts:
        cols = row.split()
        if cols[1] == 'WORDTAG':
                if cols[tag_idx] in total_dict:
                    total_dict[cols[tag_idx]] += int(cols[wordcount_idx])
                else:
                    total_dict[cols[tag_idx]] = int(cols[wordcount_idx])
            
                if cols[tag_idx] == O_tag:
                    O_dict[cols[word_idx]] = int(cols[wordcount_idx])
                else:
                    I_Gene_dict[cols[word_idx]] = int(cols[wordcount_idx])

    train_counts = open(input, "r")
    for row in train_counts:
        cols = row.split()
        if cols[1] == 'WORDTAG':
            O_prob = 0.0
            I_Gene_prob = 0.0

            if cols[word_idx] in O_dict:
                O_prob = O_dict[cols[word_idx]]/(total_dict[O_tag] * 1.0)

            if cols[word_idx] in I_Gene_dict:
                I_Gene_prob = I_Gene_dict[cols[word_idx]]/(total_dict[I_Gene_tag] * 1.0)

            exv_O_dict[cols[word_idx]] = O_prob
            exv_I_Gene_dict[cols[word_idx]] = I_Gene_prob

def Build_exv_viterbi_trigram(input):
    train_counts = open(input, "r")
    trigram_dict = dict()
    bigram_dict = dict()

    #column indexes in input data
    wordcount_idx = 0
    tag_idx = 2
    word_idx = 3

    for row in train_counts:
        cols = row.split()
        if cols[1] == '2-GRAM':
            bigram_dict[str.format('{0} {1}'.format(cols[2], cols[3]))] = int(cols[wordcount_idx])
        elif cols[1] == '3-GRAM':
            trigram_dict[str.format('{0} {1} {2}'.format(cols[2], cols[3], cols[4]))] = int(cols[wordcount_idx])

    train_counts = open(input, "r")
    for row in train_counts:
        cols = row.split()
        if cols[1] == '3-GRAM':
            bigram = str.format('{0} {1}'.format(cols[2], cols[3]))
            trigram = str.format('{0} {1} {2}'.format(cols[2], cols[3], cols[4]))
            
            q_trigram[trigram] = trigram_dict[trigram]/(bigram_dict[bigram] * 1.0)

def GetBestTrigram(bigram):
    q_value_dict = dict()

    for trigram, prob in q_trigram.iteritems():
        if trigram.startswith(bigram):
            q_value_dict[trigram.split()[2]] = prob

    return q_value_dict

def TestUnigramModel(input, input_norm, output):
    dev = open(input, "r").readlines()
    dev_norm = open(input_norm, "r").readlines()
    dev_p1_out = open(output, "w")

    for idx, row in enumerate(dev_norm):
        cols = row.split()

        if len(cols) == 0:
            dev_p1_out.write('\n')
            continue

        word_norm  =cols[0]
        word  =dev[idx].split()[0]

        O_prob = exv_O_dict[word_norm]
        I_Gene_prob = exv_I_Gene_dict[word_norm]

        if O_prob > I_Gene_prob:
            dev_p1_out.write("%s %s\n" %(word, O_tag))
        else:
            dev_p1_out.write("%s %s\n" %(word, I_Gene_tag))

    dev_p1_out.close()

def appnd2(u, v):
    return (u + ' ' + v)

def appnd3(w, u, v):
    return (w + ' ' + u + ' ' + v)

def TestTrigramModel_Viterbi(input, input_norm, output):
    dev = open(input, "r").readlines()
    dev_norm = open(input_norm, "r").readlines()
    dev_trigram_out = open(output, "w")
    #lines  = dev.readlines()

    k_2tag = '*'
    k_1tag = '*'
    ktag = ''
    bpk = [{}]
    Viterbi_PI = [{}]
    Viterbi_PI[0]['* *'] = 1.0
    wseq = []
    idx = 0

    validtags = [I_Gene_tag, O_tag]

    for rowidx, row in enumerate(dev_norm):
        idx = idx + 1 #Changing to index value which begins from 1 instead of default 0
        cols = row.split()

        if len(cols) == 0:

            maxu = ''
            maxv = ''
            maxVertibi_PI_For_u_v = 0.0
            currVertibi = 0.0
            for v in validtags:
                for u in validtags:
                    currVertibi = Viterbi_PI[-1][appnd2(u, v)] * q_trigram[appnd3(u, v, 'STOP')]
                    if currVertibi >= maxVertibi_PI_For_u_v:
                        maxVertibi_PI_For_u_v = currVertibi
                        maxu = u
                        maxv = v
            
            yseq = [''] * len(Viterbi_PI)
            yseq[len(Viterbi_PI)-1] = maxv
            yseq[len(Viterbi_PI)-2] = maxu

            for k in range(len(Viterbi_PI)-3, 0, -1):
                yseq[k] = bpk[k+2][appnd2(yseq[k+1], yseq[k+2])]

            for k in range(1, len(yseq)):
                word = dev[rowidx-len(yseq) + k].split()[0]
                dev_trigram_out.write("%s %s\n" %(word, yseq[k]))

            dev_trigram_out.write('\n')
            dev_trigram_out.flush()
            bpk = [{}]
            bpk.append({})
            Viterbi_PI = [{}]
            Viterbi_PI[0]['* *'] = 1.0
            wseq = []
            idx = 0
            continue

        kword_norm = cols[0]

        exv_dict = {O_tag : 0.0, I_Gene_tag : 0.0}
        if kword_norm in exv_O_dict:
            exv_dict[O_tag] = exv_O_dict[kword_norm]

        if kword_norm in exv_I_Gene_dict:
            exv_dict[I_Gene_tag] = exv_I_Gene_dict[kword_norm]

        if len(Viterbi_PI) == 1:
            Viterbi_PI.append({})
            bpk.append({})
            for v in validtags:
                Viterbi_PI[idx][appnd2('*', v)] = Viterbi_PI[idx-1]['* *'] * q_trigram['* * ' + v] * exv_dict[v]
                bpk[idx][appnd2('*', v)] = '*'
        elif len(Viterbi_PI) == 2:
            Viterbi_PI.append({})
            bpk.append({})
            for v in validtags:
                for u in validtags:
                    Viterbi_PI[idx][appnd2(u, v)] = Viterbi_PI[idx-1][appnd2('*', u)] * q_trigram[appnd3('*', u, v)] * exv_dict[v]
                    bpk[idx][appnd2(u, v)] = '*'
        else:
            Viterbi_PI.append({})
            bpk.append({})
            for v in validtags:
                for u in validtags:
                    maxw = ''
                    maxVertibi_PI_For_w = 0.0
                    currVertibi = 0.0

                    for w in validtags:
                         currVertibi = Viterbi_PI[idx-1][appnd2(w, u)] * q_trigram[appnd3(w, u, v)] * exv_dict[v]
                         if currVertibi >= maxVertibi_PI_For_w:
                             maxVertibi_PI_For_w = currVertibi
                             maxw = w

                    Viterbi_PI[idx][appnd2(u, v)] = maxVertibi_PI_For_w
                    wseq.append(maxw)
                    bpk[idx][appnd2(u, v)] = maxw

    dev_trigram_out.close()

train = "gene.train"
train_norm = "gene.train.informational.norm"
train_norm_counts = "gene.train.counts.informational.norm"

#Manually load "train_counts" by running command - python count_freqs.py gene.train > gene.counts
train_counts = "gene.counts"

NormalizeData_Informational(train, train_norm, train_counts)
#Manually load "gene.train.counts.informational.norm" by running command - python count_freqs.py gene.train.informational.norm > gene.train.counts.informational.norm

Build_exv_viterbi_unigram(train_norm_counts)
Build_exv_viterbi_trigram(train_norm_counts)

## Sample sentence for debugging purposes
#NormalizeData_Informational("sample.dev", "sample.dev.norm", train_counts)
#TestTrigramModel_Viterbi("sample.dev", "sample.dev.norm", "sample.dev.out")

dev = "gene.dev"
dev_norm = "gene.dev.informational.norm"
NormalizeData_Informational(dev, dev_norm, train_counts)
dev_unigram_output = "gene_dev.p1.informational.out"

#TestUnigramModel(dev, dev_norm, dev_unigram_output)

dev_trigram_output = "gene_trigram_dev.informational.out"
TestTrigramModel_Viterbi(dev, dev_norm, dev_trigram_output)
# After loading dev_trigram_output, evaluate F score by running command - python eval_gene_tagger.py gene.key gene_trigram_dev.informational.out
# Below is expected output of F score for dev_trigram_output
#Found 415 GENEs. Expected 642 GENEs; Correct: 222
#         precision      recall          F1-Score
#GENE:    0.534940       0.345794        0.420057

# test = "gene.test"
# test_unigram_output = "gene_test.p1.out"
# TestUnigramModel(test, test_unigram_output)