outputFile = "epron-jpron-unsupervised.wfst"
alignmentFile = "epron-jpron.alignment"
comparisonFile = "epron-jpron.data"
# EM algorithm implemented with the purpose of determining english -> japanese phenome alignments from parallel data
# Takes in:
# - Set of all japanses phenomes
# - Instance of english phenomes to translate

# Builds Derivation latice as 2 D Matrix
from parallel_data import batch_gen, listAllPossibleMappings

def em(pairings):
    jSequences = set([])
    eTokens = set([])
    eToJList = {}
    countMap = {}
    alignmentMap = {}

    eJCount = {}
    eCount = {}

    for (english, japanese) in pairings:
        alignments = listAllPossibleMappings(english, japanese)

        # Store alignments so they don't need to be recomputed
        if tuple(english) not in alignmentMap:
            alignmentMap[tuple(english)] = {}
        if tuple(japanese) not in alignmentMap[tuple(english)]:
            alignmentMap[tuple(english)][tuple(japanese)] = alignments

        for e in english:
            eTokens.add(e)
            if e not in eToJList:
                eToJList[e] = set([])

        for a in alignments:
            for i in range(len(english)):
                J = getCorrespondingJapeneseSequence(i, japanese, a)
                jSequences.add(tuple(J))
                eToJList[english[i]].add(tuple(J))

    P = {}
    for e in eTokens:
        countMap[e] = {}
        P[e] = {}
        # for J in jSequences:
        #     P[e][tuple(J)] = 1/float(len(jSequences))
        for J in eToJList[e]:
            countMap[e][J] = 0
            P[e][tuple(J)] = 1/float(len(jSequences))

    iteration = 0
    maxIterations = 5
    while iteration < maxIterations:
        for (english, japanese) in pairings:
            alignProbs = {}
            alignments = alignmentMap[tuple(english)][tuple(japanese)] # listAllPossibleMappings(english, japanese)

            # Get new probability of each alignment
            # Set alignment probabilies equal to start
            if iteration == 0:
                for a in alignments:
                    alignProbs[tuple(a)] = 1/float(len(alignments))
            else:
                for a in alignments:
                    p = 1
                    for i in range(len(english)):
                        J = getCorrespondingJapeneseSequence(i, japanese, a)
                        p *= P[english[i]][tuple(J)]
                    alignProbs[tuple(a)] = p

                # Normalize alignment probs
                total = 0.0
                for a in alignments:
                    total += alignProbs[tuple(a)]
                for a in alignments:
                    alignProbs[tuple(a)] = alignProbs[tuple(a)]/total

            # Collect fractional counts
            for a in alignments:
                for i in range(len(english)):
                    J = getCorrespondingJapeneseSequence(i, japanese, a)
                    countMap[english[i]][tuple(J)] += alignProbs[tuple(a)]

        # Normalize P(J|e)
        for e in eTokens:
            total = 0.0
            for J in eToJList[e]:
                total += countMap[e][J]
            for J in eToJList[e]:
                P[e][J] = countMap[e][J] / total

        # Clear counts map
        for e in eTokens:
            countMap[e] = {}
            for J in eToJList[e]:
                countMap[e][J] = 0

        iteration += 1

        print "ITERATION {} TOP 5 ALIGNMENTS".format(iteration)
        topAlignments = []
        for (english, japanese) in pairings[:5]:
            print english
            print japanese
            topK = getTopKAlignments(english, japanese, alignmentMap, P, 5)
            for a in topK:
                print "{} : {}".format([v + 1 for v in a[1]], a[0])
            print

    topAlignments = []
    for (english, japanese) in pairings:
        a = getTopKAlignments(english, japanese, alignmentMap, P, 1)[0]
        topAlignments.append((english, japanese, [v + 1 for v in a[1]]))

    output = {}
    for e in eTokens:
        output[e] = {}
    # for e in ["EY"]:
        for k in P[e]:
            if P[e][k] > 0.01:
                output[e][k] = P[e][k]
                # print "P[\"{}\"][\"{}\"]: {}".format(e, k, P[e][k])
    return output, topAlignments

def getTopKAlignments(english, japanese, alignmentMap, P, K):
    alignments = alignmentMap[tuple(english)][tuple(japanese)]
    k = K if K <= len(alignments) else len(alignments)
    alignmentsWithScores = []
    for a in alignments:
        score = 1
        for i in range(len(english)):
            J = getCorrespondingJapeneseSequence(i, japanese, a)
            score *= P[english[i]][tuple(J)]
        alignmentsWithScores.append((score, a))
    alignmentsWithScores.sort(key=lambda x:x[0])
    return alignmentsWithScores[::-1][:k]

def getCorrespondingJapeneseSequence(eIndex, j, alignment):
    start = alignment.index(eIndex)
    end = start
    while end < len(alignment) and alignment[end] == eIndex:
        end += 1
    J = j[start:end]
    return J

##################################################

class FSANode:
    idCounter = 0
    outgoing = {}
    incoming = {}
    ID = -1
    is_start = False
    is_end = False

    def __init__(self, start=False, end=False):
        self.ID = FSANode.idCounter 
        FSANode.idCounter = FSANode.idCounter + 1
        self.is_start = start
        self.is_end = end

        self.incoming = {}
        self.outgoing = {}
        self.startToStart = []

    def stateName(self):
        if self.is_start:
            return "START_OR_END"
        elif self.is_end:
            return "START_OR_END"
        else:
            return str(self.ID)

    def removeOutgoingTransitionToNode(self, node):
        if node.ID in self.outgoing:
            del self.outgoing[node.ID]

    def removeIncomingTransitionFromNode(self, node):
        if node.ID in self.outgoing:
            del self.incoming[node.ID]

    def getNodesUsingInputValue(self, inputValue):
        nodesUsingInputValue = []
        for _, (node, inputSymbol, _, _) in self.outgoing.iteritems():
            if inputValue == inputSymbol:
                nodesUsingInputValue.append(node)

        return nodesUsingInputValue

    def getTransitions(self):
        transitions = []
        if len(self.outgoing) == 0:
            return transitions

        # Get all transitions from current state
        # for key, (node, symbol) in self.outgoing.iteritems():
        #     transition_item = "*e*" if symbol == "*e*" else "\"" + symbol + "\""
        #     content_string = "(%s (%s %s))" % (self.stateName(), node.stateName(), transition_item)
        #     transitions.append(content_string)

        for key, (node, inputSymbol, outputSymbol, prob) in self.outgoing.iteritems():
            inputItem = "*e*" if inputSymbol == "*e*" else "\"" + inputSymbol + "\""
            outputItem = "*e*" if outputSymbol == "*e*" else "\"" + outputSymbol + "\""
            content_string = "({0} ({1} {2} {3} {4}))".format(self.stateName(), node.stateName(), inputItem, outputItem, prob)
            transitions.append(content_string)

        for (node, inputSymbol, outputSymbol, prob) in self.startToStart:
            inputItem = "*e*" if inputSymbol == "*e*" else "\"" + inputSymbol + "\""
            outputItem = "*e*" if outputSymbol == "*e*" else "\"" + outputSymbol + "\""
            content_string = "({0} ({1} {2} {3} {4}))".format(self.stateName(), node.stateName(), inputItem, outputItem, prob)
            transitions.append(content_string)

        # Get transitions for child states
        for key, (node, inputSymbol, outputSymbol, prob) in self.outgoing.iteritems():
            if node.ID != self.ID and not node.is_start:
                transitions = transitions + node.getTransitions()

        return transitions

def createLinkBetweenNodes(fromNode, toNode, inputSymbol, outputSymbol, prob):
    if fromNode == toNode:
        assert fromNode.is_start
        fromNode.startToStart += [(toNode, inputSymbol, outputSymbol, prob)]
        return
    fromNode.outgoing[toNode.ID] = (toNode, inputSymbol, outputSymbol, prob)
    toNode.incoming[fromNode.ID] = (fromNode, inputSymbol, outputSymbol, prob)

def addPath(fstStart, e, jList, prob):
    if len(jList) == 1:
        createLinkBetweenNodes(fstStart, fstStart, e, jList[0], prob)
        return

    current = FSANode()
    createLinkBetweenNodes(fstStart, current, e, jList.pop(0), prob)

    while len(jList) > 1:
        nexxt = FSANode()
        createLinkBetweenNodes(current, nexxt, "*e*", jList.pop(0), 1)
        current = nexxt

    createLinkBetweenNodes(current, fstStart, "*e*", jList.pop(0), 1)

def buildFST(toPrint):
    fstStart = FSANode(start=True, end=True)
    for e in toPrint:
        for J in toPrint[e]:
            value = toPrint[e][J]
            addPath(fstStart, e, list(J), value)

    return fstStart

def outputFST(fst):
    content = []
    content.append("START_OR_END")
    content = content + fst.getTransitions()
    with open(outputFile, 'w') as f:
        for line in content:
            f.write("%s\n" % line)

##################################################

def compareAlignments(alignmentFile, comparisonFile):
    alignmentData = []
    with open(alignmentFile, 'r') as f:
        lines = f.readlines()
        for (_, _, mapping) in batch_gen(lines, 3):
            mapping = mapping.strip().split(" ")
            alignmentData.append(mapping)

    comparisonData = []
    with open(comparisonFile, 'r') as f:
        lines = f.readlines()
        for (_, _, mapping) in batch_gen(lines, 3):
            mapping = mapping.strip().split(" ")
            comparisonData.append(mapping)

    assert len(comparisonData) == len(alignmentData)

    i = 0
    while i < len(comparisonData):
        assert len(comparisonData[i]) == len(alignmentData[i])
        correct = 0.0
        total = 0.0
        for j in range(len(comparisonData[i])):
            total += 1
            correct += 1 if comparisonData[i][j] == alignmentData[i][j] else 0
        i += 1
    print "ACCURACY: {}".format(correct/total)


##################################################

# def printAlignments(topAlignments):
#     for (e, j, a) in topAlignments[:5]:
#         print map(lambda x: "{}".format(x), e)
#         print map(lambda x: "{}".format(x), j)
#         print a

def outputAlignments(alignmentFile, topAlignments):
    with open(alignmentFile, 'w') as f:
        for (e, j, a) in topAlignments:
            f.write("%s\n" % " ".join(["\"{}\"".format(e_) for e_ in e]))
            f.write("%s\n" % " ".join(["\"{}\"".format(j_) for j_ in j]))
            f.write("%s\n" % " ".join([str(a_) for a_ in a]))


##################################################

def test():
    filename = "epron-jpron-unsupervised.data"
    pairs = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        count = 0
        for (e, j, _) in batch_gen(lines, 3):
            e = map(lambda x: x.replace("\"",""), e.strip().split(" "))
            j = map(lambda x: x.replace("\"",""), j.strip().split(" "))
            pairs.append((e, j))

    toPrint, topAlignments = em(pairs)
    outputAlignments(alignmentFile, topAlignments)

    compareAlignments(alignmentFile, comparisonFile)

    fst = buildFST(toPrint)
    outputFST(fst)