from Cytoscape import *
from string import strip, split
import pickle
from specialFunc import regenReaction
import sys



def read_1col(f):
    out=open(f,"r")
    data=[];
    for line in out.readlines():
        if line:
            aa=map(strip,split(line,'\t'))
            if aa[0] not in data:
                data.append(aa[0])
    return data

def drawNet(netFile,rxnFile):
    rxnList=read_1col(rxnFile)
    temp = pickle.load(open(netFile))

    reactions_info=temp[2]
    new_rInfo = {}

    for k in rxnList:
        if (k in (temp[0]['rxns'])) & (k in reactions_info.keys()):  # this is just for signaling reactions
            E = reactions_info[k]['M']
            ind = temp[0]['rxns'].index(k)
            act = temp[0]['activators'][ind]
            inh = temp[0]['inhibitors'][ind]

            if not (act == 'NA'):
                acts = act.split('|')

                for a in acts:
                    if a in temp[1].keys():
                        name = temp[1][a][0];
                    elif a in temp[0]['rashmiv']:
                        name = temp[0]['rashmik'][temp[0]['rashmiv'].index(a)]
                    else:
                        name = a;
                        print a
                    E.append(name)
            I = []
            if not (inh == 'NA'):
                acts = inh.split('|')
                for a in acts:
                    if a in temp[1].keys():
                        name = temp[1][a][0];
                    elif a in temp[0]['rashmiv']:
                        name = temp[0]['rashmik'][temp[0]['rashmiv'].index(a)]
                    else:
                        name = a;
                        print a
                    I.append(name)

            new_rInfo[k] = {'S': reactions_info[k]['S'], 'P': reactions_info[k]['P'], 'E': E, 'I': I}
    # now we want to create network on cytoscape
    # new_rInfo {'reaction_1980109': {'I': [], 'P': ['NICD1:DTX [cytosol]'], 'S': ['CNTN1:NOTCH1:DTX [plasma membrane]'],
    #'E': ['gamma-secretase complex [plasma membrane]']},
    # species_label species can be repalced
    # cofactorsID list of cofactor Ids
    cofactorsID=[]
    cys = CyNetwork()
    species_label={}
    rInfo, labels, newIds = regenReaction(new_rInfo, species_label, cofactorsID)

    node_list, edge_list = cys.createNodeEdgeCofactor(rInfo, labels, newIds)
    cys.createCysNet(node_list, edge_list)


if "__main__==__name__":
    sys.argv[1]='/Users/vikash/switchdrive2/RedoConnectMetSig/Data/Human_signaling_info.pickle'
    sys.argv[2]='/Users/vikash/switchdrive2/PyCytoscape/rxn.txt'
    drawNet(sys.argv[1],sys.argv[2])
