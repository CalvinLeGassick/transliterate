def getNumberOfTimesStartingElementRepeats(listOfIntegers, start):
    valueToMatch = listOfIntegers[start]
    length = 0
    while (start + length) < len(listOfIntegers) and valueToMatch == listOfIntegers[start + length]:
        length = length + 1
    return length

class PhenomePair:
    ePhenList = []
    jPhenList = []
    phenMapping = []


    def __init__(self, eng, jap, mapping):
        self.ePhenList = eng
        self.jPhenList = jap
        self.phenMapping = mapping

    def getMappingOfEnglishPhenomesToJapansePhenomeLists(self):
        assert len(self.jPhenList) == len(self.phenMapping) 

        ePhenToJPhenList = {}
        ePhenomeIndex = 0
        jPhenomeIndex = 0

        while jPhenomeIndex < len(self.jPhenList):
            length = getNumberOfTimesStartingElementRepeats(self.phenMapping, jPhenomeIndex)
            ePhenToJPhenList[self.ePhenList[ePhenomeIndex]] = self.jPhenList[jPhenomeIndex:jPhenomeIndex + length]
            jPhenomeIndex = jPhenomeIndex + length
            # print "{0}: {1}".format(self.ePhenList[ePhenomeIndex], ePhenToJPhenList[self.ePhenList[ePhenomeIndex]])
            ePhenomeIndex = ePhenomeIndex + 1

        return ePhenToJPhenList

def allMappingsHelper(e, j, eIndex, jIndex):
    toReturn = []

    # Does not work
    if (len(e) - eIndex > len(j) - jIndex):
        return (False, [])

    # Base case:
    if len(e) - eIndex == 1:
        return True, [[eIndex for x in range(len(j) - jIndex)]]

    for i in range(jIndex + 1, len(j)):
        works, mappings = allMappingsHelper(e, j, eIndex + 1, i)
        if works:
            for m in mappings:
                toReturn.append([eIndex for x in range(i - jIndex)] + m)

    return (True, toReturn)

def listAllPossibleMappings(e, j):
    assert len(e) <= len(j)
    works, mappings  = allMappingsHelper(e, j, 0, 0)
    assert works
    return mappings

def batch_gen(data, batch_size):
    for i in range(0, len(data), batch_size):
        yield data[i:i + batch_size]

def getPhenomePairsFromData(filename):
    phenomePairs = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        count = 0
        for (e, j, _) in batch_gen(lines, 3):
            if count > 5:
                break;
            e = map(lambda x: x.replace("\"",""), e.strip().split(" "))
            j = map(lambda x: x.replace("\"",""), j.strip().split(" "))
            mappings =  listAllPossibleMappings(e, j)
            print "\nWORD PAIRINGS"
            print "{}\n{}".format(e, j)
            print "----------------------\n"
            for m in mappings:
                print "{}".format(map(lambda x: x + 1, m))
            print
            count += 1

    return phenomePairs
