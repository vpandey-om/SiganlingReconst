from Cytoscape import *
from libsbml import *
from specialFunc import parse_xml_reactome,regenReaction
import pickle
import scipy.io as spio


def notch1Vizualize(vizFile):
    Info=pickle.load(open(vizFile))
    # for k in Info[0].keys():
    #     print k
    # import pdb;pdb.set_trace()
    # Info[0] : rxn Info
    # Info[1]: is label
    #import pdb;pdb.set_trace()
    cofactorsID = ['AMP [cytosol]', 'Pi [nucleoplasm]', 'Ca2+ [cytosol]', 'H2O [cytosol]', 'GDP [cytosol]',
                   'ATP [cytosol]', 'Ub [cytosol]', 'ADP [nucleoplasm]', 'ADP [cytosol]', 'Pi [cytosol',
                   'ATP [nucleoplasm]', 'GTP [cytosol]']
    cys = CyNetwork()
    rInfo, labels, newIds = regenReaction(Info[0], Info[1],cofactorsID)
    node_list, edge_list = cys.createNodeEdgeCofactor(rInfo, labels, newIds)
    cys.createCysNet(node_list, edge_list)



def exportTest(Info,keyFile,speciesInfoFile):
    # Info[0] is reaction information
    # Info[1] nodes labels
    # Info[2] original reaction information
    # Info[3] original species information
    cys = CyNetwork()
    rxnInfo=cys.exportReactionInfo()
    # we want to catch Stoichiomatrix from the original


    for k in rxnInfo.keys():
        SC=[]
        k1 = k.replace('reaction_', '')
        for item in rxnInfo[k]['S']:


            if item in Info[2][k1]['S']:
                i= Info[2][k1]['S'].index(item)
                SC.append(Info[2][k1]['SC'][i])
        PC=[]
        for item in rxnInfo[k]['P']:
            k1 = k.replace('reaction_', '')
            if item in Info[2][k1]['P']:
                i = Info[2][k1]['P'].index(item)
                PC.append(Info[2][k1]['PC'][i])
        rxnInfo[k]['SC']=SC
        rxnInfo[k]['PC'] = PC
    #import pdb;pdb.set_trace()

    # ll=rxnInfo.keys()
    # mm=rxnInfo.values()
    out=open(speciesInfoFile,'w')
    for k,v in Info[3].items():
        out.write('%s\t%s\t%s\n'%(k,v[0],v[1]))
    spio.savemat(valueFile,mdict={'my_values': rxnInfo.values()})
    spio.savemat(keyFile,mdict={'my_keys': rxnInfo.keys()})



if "__main__==__name__":
    # if you want to visualize network form python to the cytoscape software
    # we want to take input as dictionary format
    # input: in form of nodes and edges
    # for a reaction network we use a dictionary:
    #### {'R1':{'I':['I1','I2'],'P':['P1','P2'],'S': ['S1','S2']},'E':['E1,E2']}
    #vizFile='/Users/vikashpandey/Documents/MATLAB/MetSignal/TIGER/tigerSignal/SignalingAnalysis/PythonData/DataPickle/Notch1VizModel.pickle'
    vizFile = '/Users/vikashpandey/Documents/MATLAB/MetSignal/TIGER/tigerSignal/SignalingAnalysis/PythonData/DataPickle/SignalingVizModel.pickle'
    #vizFile = '/Users/vikashpandey/Documents/MATLAB/MetSignal/TIGER/tigerSignal/SignalingAnalysis/PythonData/DataPickle/ExceptMetabolismVizModel.pickle'
    notch1Vizualize(vizFile)
    originalVizInfo=pickle.load(open(vizFile))
    #Info[0] is reaction information
    # Info[1] nodes labels
    # Info[2] original reaction information
    # Info[3] original species information
    #import pdb;pdb.set_trace()
    # we can export from the cytoscape
    # valueFile='/Users/vikashpandey/Documents/MATLAB/MetSignal/TIGER/tigerSignal/SignalingAnalysis/MetalabData/Data/notch1Cyto.mat'
    # keyFile='/Users/vikashpandey/Documents/MATLAB/MetSignal/TIGER/tigerSignal/SignalingAnalysis/MetalabData/Data/notch1CytoKey.mat'
    # speciesInfoFile = '/Users/vikashpandey/Documents/MATLAB/MetSignal/TIGER/tigerSignal/SignalingAnalysis/MetalabData/Data/speciesInfo.txt'

    # valueFile = '/Users/vikashpandey/Documents/MATLAB/MetSignal/TIGER/tigerSignal/SignalingAnalysis/MetalabData/Data/ExceptMetabolismCyto.mat'
    # keyFile = '/Users/vikashpandey/Documents/MATLAB/MetSignal/TIGER/tigerSignal/SignalingAnalysis/MetalabData/Data/ExceptMetabolismCytoKey.mat'
    # speciesInfoFile = '/Users/vikashpandey/Documents/MATLAB/MetSignal/TIGER/tigerSignal/SignalingAnalysis/MetalabData/Data/ExceptMetabolismspeciesInfo.txt'

    # valueFile = '/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/SignalingAnalysis/MetalabData/Data/SignalingWholeCyto.mat'
    # keyFile = '/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/SignalingAnalysis/MetalabData/Data/SignalingWholeCytoKey.mat'
    # speciesInfoFile = '/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/SignalingAnalysis/MetalabData/Data/SignalingWholespeciesInfo.txt'
    #
    # valueFile = '/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/SignalingAnalysis/MetalabData/Data/SignalingBigCompoCyto.mat'
    # keyFile = '/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/SignalingAnalysis/MetalabData/Data/SignalingBigCompoCytoKey.mat'
    # speciesInfoFile = '/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/SignalingAnalysis/MetalabData/Data/SignalingBigCompospeciesInfo.txt'
    #
    # exportTest(originalVizInfo,keyFile,speciesInfoFile) # Info is coming form cytoscape visualization

