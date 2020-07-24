from Cytoscape import *
from libsbml import *
from specialFunc import parse_xml_reactome,regenReaction
import pickle
import scipy.io as spio

def notchOld():
    cys=CyNetwork()
    filename = '/Users/vikash/Documents/MATLAB/MetSignal/REACTOME_SBML2/Signaling_by_EGFR.xml'
    species, reactions_info = parse_xml_reactome(filename)

    cofactorsID=['AMP [cytosol]','Pi [nucleoplasm]','Ca2+ [cytosol]','H2O [cytosol]','GDP [cytosol]','ATP [cytosol]','Ub [cytosol]','ADP [nucleoplasm]','ADP [cytosol]','Pi [cytosol','ATP [nucleoplasm]','GTP [cytosol]']
    out=open('egfr.txt', 'w')
    # species dictionary
    new_species={}

    for k,v in species.items():
        new_species[v[0]]=v[0].split('[')[0]
        print v[0]
        out.write("%s\n"%v[0])
    for k,v in reactions_info.items():
        new_species[k]=''


    import pdb;
    pdb.set_trace()
    rInfo, labels, newIds=regenReaction(reactions_info, new_species, cofactorsID)
    node_list, edge_list = cys.createNodeEdgeCofactor(rInfo,labels,newIds)
    cys.createCysNet(node_list, edge_list)


def notch():
    temp=pickle.load(open('/Users/vikashpandey/Documents/MATLAB/MetSignal/TIGER/tigerSignal/ReactomeModel/HumanSignaling/Human_signaling_info.pickle'))
    filename = '/Users/vikashpandey/Documents/PyCytoscape/VisPathways/notch1.xml'
    species, reactions_info = parse_xml_reactome(filename)


    # species dictionary
    import pdb;
    pdb.set_trace()

    new_rInfo={}

    for k,vals in reactions_info.items():
        if k in temp[0]['rxns']: # this is just for signaling reactions
            E=reactions_info[k]['M']
            ind=temp[0]['rxns'].index(k)
            act=temp[0]['activators'][ind]
            inh = temp[0]['inhibitors'][ind]

            if not (act=='NA'):
                acts=act.split('|')
                for a in acts:
                    if a in temp[1].keys():
                        name=temp[1][a][0];
                    elif a in temp[0]['rashmiv']:
                        name=temp[0]['rashmik'][temp[0]['rashmiv'].index(a)]
                    else:
                        name=a;
                        print a
                    E.append(name)
            I=[]
            if not (inh=='NA'):
                acts=inh.split('|')
                for a in acts:
                    if a in temp[1].keys():
                        name=temp[1][a][0];
                    elif a in temp[0]['rashmiv']:
                        name=temp[0]['rashmik'][temp[0]['rashmiv'].index(a)]
                    else:
                        name=a;
                        print a
                    I.append(name)

            new_rInfo[k]={'S':reactions_info[k]['S'],'P':reactions_info[k]['P'],'E':E,'I':I}

    new_species = {}


    for k, v in new_rInfo.items():
        new_species[k] = ''
        for item in v['S']:
            new_species[item]=item.split('[')[0]

        for item in v['P']:
            new_species[item] = item.split('[')[0]

        for item in v['E']:
            new_species[item] = item.split('[')[0]

        for item in v['I']:
            new_species[item] = item.split('[')[0]
        # species dictionary

    pickle.dump([new_rInfo,new_species],open('/Users/vikashpandey/Documents/PyCytoscape/VisPathways/notch.pickle','w'))


def notchViz():
    Info=pickle.load(open('/Users/vikashpandey/Documents/PyCytoscape/VisPathways/notch.pickle'))
    for k in Info[0].keys():
        print k
    import pdb;pdb.set_trace()

    cofactorsID = ['AMP [cytosol]', 'Pi [nucleoplasm]', 'Ca2+ [cytosol]', 'H2O [cytosol]', 'GDP [cytosol]',
                   'ATP [cytosol]', 'Ub [cytosol]', 'ADP [nucleoplasm]', 'ADP [cytosol]', 'Pi [cytosol',
                   'ATP [nucleoplasm]', 'GTP [cytosol]']
    cys = CyNetwork()
    rInfo, labels, newIds = regenReaction(Info[0], Info[1],cofactorsID)
    node_list, edge_list = cys.createNodeEdgeCofactor(rInfo, labels, newIds)
    cys.createCysNet(node_list, edge_list)
#
# def readGraphML():
#     from pygraphml import Graph
#     from pygraphml import GraphMLParser
#     parser = GraphMLParser()
#     g = parser.parse("myGraph.graphml")

def exportTest():
    cys = CyNetwork()
    rxnInfo=cys.exportReactionInfo()
    import pdb;pdb.set_trace()
    ll=rxnInfo.keys()
    mm=rxnInfo.values()
    spio.savemat('/Users/vikashpandey/Documents/MATLAB/MetSignal/TIGER/tigerSignal/ReactomeModel/NOTCH1New/notch1.mat',rxnInfo)
    spio.savemat('/Users/vikashpandey/Documents/MATLAB/MetSignal/TIGER/tigerSignal/ReactomeModel/NOTCH1New/notch1Key.mat',mdict={'my_list': ll})

def notchViz2():
    Info=pickle.load(open('/Users/vikashpandey/Documents/PyCytoscape/VisPathways/notch.pickle'))
    for k in Info[0].keys():
        print k
    import pdb;pdb.set_trace()

    # #cofactorsID = ['AMP [cytosol]', 'Pi [nucleoplasm]', 'Ca2+ [cytosol]', 'H2O [cytosol]', 'GDP [cytosol]',
    #                'ATP [cytosol]', 'Ub [cytosol]', 'ADP [nucleoplasm]', 'ADP [cytosol]', 'Pi [cytosol',
    #                'ATP [nucleoplasm]', 'GTP [cytosol]']
    cofactorsID=[]
    cys = CyNetwork()
    rInfo, labels, newIds = regenReaction(Info[0], Info[1],cofactorsID)
    node_list, edge_list = cys.createNodeEdgeList(rInfo, labels)
    cys.createCysNet(node_list, edge_list)








if "__main__==__name__":
    notch()
    #notchViz()
    #readGraphML()
    #exportTest()

    ##
    # This is the new analysis
    #cytoscape

    ##