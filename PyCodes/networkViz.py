from Cytoscape import *
#from libsbml import *
from specialFunc import parse_xml_reactome,regenReaction
import pickle
import scipy.io as spio

from string import strip, split

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
    pdb.set_trace()\

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


def egfViz():
    # This is the xml file for the egf network. we can use this to parse any network
    fin = '/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/testModels/R-HSA-177929_EGFR.xml'

    species, reactions_info = parse_xml_reactome(fin)

    # step 2
    # we want take information  from pickle signaling model because signaling model is more complete
    import pickle
    sigf = '/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/SignalingAnalysis/PythonData/DataPickle/Signaling.pickle'
    data = pickle.load(open(sigf))
    species = data[0]
    rxnInfo = data[1]
    # find reaction and specis info of reduced model from signaling model
    rInfo = {}

    for rxn in reactions_info.keys():
        newRxn = rxn.replace('reaction_', 'R-HSA-')
        if newRxn in rxnInfo.keys():
            rInfo[newRxn] = rxnInfo[newRxn]
        else:
            print rxn

    #import pdb;pdb.set_trace()
    new_rInfo, new_species=collectNewNet(rInfo)

    # visualize with new data

    cofactorsID = ['AMP [cytosol]', 'Pi [nucleoplasm]', 'Ca2+ [cytosol]', 'H2O [cytosol]', 'GDP [cytosol]',
                   'ATP [cytosol]', 'Ub [cytosol]', 'ADP [nucleoplasm]', 'ADP [cytosol]', 'Pi [cytosol]',
                   'ATP [nucleoplasm]', 'GTP [cytosol]']

    cys = CyNetwork()
    rInfo, labels, newIds = regenReaction(new_rInfo, new_species,cofactorsID)
    f= '/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/SignalingAnalysis/PythonData/DataPickle/proteinShortName.pickle'
    proteinShortName = pickle.load(open(f))
    new_labels={}
    for k,v in labels.items():
        if k in proteinShortName.keys():
            new_labels[k]=proteinShortName[k]
        else:
            new_labels[k]=v

    node_list, edge_list = cys.createNodeEdgeCofactor(rInfo, new_labels,newIds)
    cys.createCysNet(node_list, edge_list)



def collectNewNet(reactions_info):
    # network collect information from reaction information to visualize networks

    new_rInfo={}

    for k,vals in reactions_info.items():
       new_rInfo[k]={'S':reactions_info[k]['SN'],'P':reactions_info[k]['PN'],'E':reactions_info[k]['EN'],'I':reactions_info[k]['IN']}

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
    return new_rInfo,new_species


def visFamilyShortName():
    # here we want to visualize long complexes with short name
    import pickle
    sigf = '/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/SignalingAnalysis/PythonData/DataPickle/speciesAll.pickle'
    species = pickle.load(open(sigf))
    #import pdb;pdb.set_trace()
    familyProtein={}
    for sp,item in species.items():
        if item[2] in ['DefinedSet','Complex','CandidateSet']:
            fm_protein=item[0].split(':')[0].strip()
            if len(fm_protein)>20:
                fm_protein=sp
            if fm_protein not in familyProtein.keys():
                familyProtein[fm_protein]=[]


            familyProtein[fm_protein].append([sp,item[1]])
    ll=[]
    proteinShortName={}
    proteinShortName2 = {}
    for k,item in familyProtein.items():

        for i,v in enumerate(item,1):
            proteinShortName[v[0]]=k+'*'+str(i)
            proteinShortName2[k + '*' + str(i)] = v[1]
    out=open('/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/SignalingAnalysis/PythonData/DataPickle/shortNameProtein.txt','w')
    for k,v in proteinShortName.items():
        out.write('%s\t%s\t%s\n' % (proteinShortName2[v],v,k))
    pickle.dump(proteinShortName, open('/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/SignalingAnalysis/PythonData/DataPickle/proteinShortName.pickle','w'))

    #import pdb;pdb.set_trace()


    # netfile='/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/SignalingAnalysis/PythonData/DataPickle/ExceptMetabolismNew.pickle'
    # data2 = pickle.load(open(netfile))
    # species = data2[0]
    # rxnInfo = data2[1]
    # new_species={}
    #
    # for k, v in rxnInfo.items():
    #     for i1,item in enumerate(v['S']):
    #         if item in data[0]:
    #             new_species[item]=[v['SN'][i1].split('[')[0],v['SN'][i1],data[1][data[0].index(item)]]
    #
    #     for i1, item in enumerate(v['P']):
    #         new_species[item] = [v['PN'][i1].split('[')[0], v['PN'][i1],data[1][data[0].index(item)]]
    #
    #     for i1, item in enumerate(v['E']):
    #         new_species[item] = [v['EN'][i1].split('[')[0], v['EN'][i1],data[1][data[0].index(item)]]
    #
    #     for i1, item in enumerate(v['I']):
    #         new_species[item] = [v['IN'][i1].split('[')[0], v['IN'][i1],data[1][data[0].index(item)]]
    # import pdb;pdb.set_trace('/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/SignalingAnalysis/PythonData/DataPickle/ExceptMetabolism.pickle')
    # #for c, value in enumerate(data[1]):
    # pickle.dump(new_species,open('/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/SignalingAnalysis/PythonData/DataPickle/speciesAll.pickle','w'))

    return species

def read_1col(f):
    out=open(f,"r")
    data=[];
    for line in out.readlines():
        if line:
            aa=map(strip,split(line,'\t'))
            if aa[0] not in data:
                data.append(aa[0])

    return data

def rxnWiseViz_Maria():
    # This function takes input as reactions as a text file
    f="/Users/vikash/switchdrive2/RedoConnectMetSig/Data/signaling_viz.txt"
    rxns=read_1col(f)
    import pickle
    sigf='/Users/vikash/Desktop/consign/ExceptMetabolismNew.pickle'

    # sigf = '/Users/vikash/switchdrive2/RedoConnectMetSig/Data/Signaling.pickle'
    data = pickle.load(open(sigf))
    species = data[0]
    rxnInfo = data[1]
    # find reaction and specis info of reduced model from signaling model
    rInfo = {}

    for rxn in rxns:
        if rxn in rxnInfo.keys():
            rInfo[rxn] = rxnInfo[rxn]
        else:
            print "this is not found in the siganling network=  %s"%rxn
    # import pdb;
    # pdb.set_trace()
    # import pdb;pdb.set_trace()
    new_rInfo, new_species = collectNewNet(rInfo)

    # visualize with new data

    cofactorsID = ['AMP [cytosol]', 'Pi [nucleoplasm]', 'Ca2+ [cytosol]', 'H2O [cytosol]', 'GDP [cytosol]',
                   'ATP [cytosol]', 'Ub [cytosol]', 'ADP [nucleoplasm]', 'ADP [cytosol]', 'Pi [cytosol]',
                   'ATP [nucleoplasm]', 'GTP [cytosol]']

    cys = CyNetwork()
    rInfo, labels, newIds = regenReaction(new_rInfo, new_species, cofactorsID)
    f = '/Users/vikash/switchdrive2/RedoConnectMetSig/Data/proteinShortName.pickle'
    proteinShortName = pickle.load(open(f))
    new_labels = {}
    for k, v in labels.items():
        if k in proteinShortName.keys():
            new_labels[k] = proteinShortName[k]
        else:
            new_labels[k] = v

    node_list, edge_list = cys.createNodeEdgeCofactor(rInfo, new_labels, newIds)
    cys.createCysNet(node_list, edge_list)



if "__main__==__name__":
    rxnWiseViz_Maria()
    #egfViz() #to visualize EGF signaling
    #visFamilyShortName()
    #notch()
    #notchViz()
    #readGraphML()
    #exportTest()

    ##
    # This is the new analysis
    #cytoscape

    ##
