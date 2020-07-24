import pickle
import os
import sys
import networkx as nx
import pandas as pd

code=os.getcwd()
upLevel=code.replace('scripts','') ####### we are going to upper level of code directory
sys.path.insert(0,upLevel+'/data')
sys.path.insert(1, upLevel+'/outs')
data_folder=sys.path[0]

## output folder where we will write figures and output files
out_folder=sys.path[1]

# ID conversion: we used plasmoDB to convert ID for P. Berghai

def buildGraph(data_folder,rxnList=None):
    ''' Genrate graph from list of reactions '''
    # read reactome db
    rxnInfo,species=readReactomedb(data_folder)
    if rxnList==None:
        rxnList=list(rxnInfo.keys())

    ## build directed graph
    G=nx.DiGraph()
    participations=['A','S','I', 'E']
    for rxn in rxnList:
        if rxn in rxnInfo.keys():
            for p in participations:
                for item in rxnInfo[rxn][p]:
                    G.add_edge( item, rxn )

            for item in rxnInfo[rxn]['P']:
                G.add_edge(rxn,item )
    return G


def remove_hub_modes(G,data_folder):
    remove_node_df=pd.read_csv(data_folder+'/Block_mets.txt','\t')
    nodes=[node for node in G.nodes]
    ids_all={}
    for node in nodes:
        ids_all[node.split('-')[-1]]=node

    block_nodes=remove_node_df['Id'].str.replace('species_','').to_list()
    remove_nodes=set(ids_all.keys()) & set(block_nodes)
    remove_nodes_ingraph=[]
    for item in remove_nodes:
        remove_nodes_ingraph.append(ids_all[item])
    oldG=G.copy()
    for item in remove_nodes_ingraph:
        G.remove_node(item)
    return G,oldG

def readReactomedb(data_folder):
    reactome_DB=pickle.load(open(data_folder+'/ExceptMetabolismNew.pickle','rb'))
    species = reactome_DB[0] ###  proteins
    rxnInfo = reactome_DB[1] ## reactions
    return rxnInfo,species


def getSuccessors(G, species,depth=1):
    ''' Compute successors at given depth form a source'''

    print ('Executing successors')
    rxns=[]
    visitedNodes=[]

    for i in range(1,depth+1):

        if i==1:
            # get successors

            current_rxns=[successor for successor in G.successors(species)]
            if len(current_rxns)>0:
                rxns.append(current_rxns) # these are the reaction nodes
            temp_proteins=[] # we collect all proteins form given reactions
            for rxn_item in current_rxns:
                for protein_item in [successor for successor in G.successors(rxn_item)]:
                    temp_proteins.append(protein_item)



        else:
            # if i==5:
            #     import pdb; pdb.set_trace()
            ## we want to get reactions
            if len(temp_proteins)>0:
                for rxn_item in current_rxns:
                    visitedNodes.append(rxn_item)
                temp_rxns=[]
                for protein_item in  temp_proteins:

                     current_rxns=[successor for successor in G.successors(protein_item)]
                     if len(current_rxns)>0:

                         for item in current_rxns:
                             temp_rxns.append(item)
            ### comapre visted reactions with current reactions
                current_rxns=list(set(temp_rxns)-set(visitedNodes))
                rxns.append(current_rxns) # these are the reaction nodes
                temp_proteins=[] # we collect all proteins form given reactions
                for  rxn_item in current_rxns:
                    for protein_item in  [successor for successor in G.successors(rxn_item)]:
                        temp_proteins.append(protein_item)

    return rxns


def getPredicessors(G, species,depth=1):
    ''' Compute successors at given depth form a source'''

    print ('Executing Predicessors ')
    rxns=[]
    visitedNodes=[]

    for i in range(1,depth+1):

        if i==1:
            # get successors


            current_rxns=[predecessor for predecessor in G.predecessors(species)]

            if len(current_rxns)>0:
                rxns.append(current_rxns) # these are the reaction nodes
            temp_proteins=[] # we collect all proteins form given reactions
            for rxn_item in current_rxns:
                for protein_item in [predecessor for predecessor in G.predecessors(rxn_item)]:
                    temp_proteins.append(protein_item)



        else:
            # if i==5:
            #     import pdb; pdb.set_trace()
            ## we want to get reactions
            if len(temp_proteins)>0:
                for rxn_item in current_rxns:
                    visitedNodes.append(rxn_item)
                temp_rxns=[]
                for protein_item in  temp_proteins:

                     current_rxns=[predecessor for predecessor in G.predecessors(protein_item)]
                     if len(current_rxns)>0:

                         for item in current_rxns:
                             temp_rxns.append(item)
            ### comapre visted reactions with current reactions
                current_rxns=list(set(temp_rxns)-set(visitedNodes))
                rxns.append(current_rxns) # these are the reaction nodes
                temp_proteins=[] # we collect all proteins form given reactions
                for  rxn_item in current_rxns:
                    for protein_item in  [predecessor for predecessor in G.predecessors(rxn_item)]:
                        temp_proteins.append(protein_item)

    return rxns


def getSource(G,rxnInfo,rxns):
    ''' identify source nodes from a graph'''
    ## get indegre and out degress of block_nodes
    sources=[]
    for node in G.nodes:
        if (G.in_degree(node)==0) and (G.out_degree(node)>0):
            sources.append(node)
    species_names=get_species_name(rxnInfo,rxns)

    source_ids_names={}
    for source in sources:
        if source in species_names:
            source_ids_names[source]=species_names[source]


    return source_ids_names



def get_species_name(rxn_info,rxns):
    ## we need to give reaction information from reactome db
    new_species = {}
    for rxn  in rxns:
        if rxn in rxn_info.keys():
            v=rxn_info[rxn] ### values of reaction
            for i,item in enumerate(v['SN']):
                new_species[v['S'][i]]=item
            for i,item in enumerate(v['PN']):
                new_species[v['P'][i]]=item
            for i,item in enumerate(v['AN']):
                new_species[v['A'][i]]=item
            for i,item in enumerate(v['EN']):
                new_species[v['E'][i]]=item
            for i,item in enumerate(v['IN']):
                new_species[v['I'][i]]=item


    return new_species



def stepwiseGraphAnalysis(data_folder):


    ###
    species=pickle.load(open(data_folder+'/speciesAll.pickle','rb'))
    unique_species_types=[]
    speciestype={}
    for k,v in species.items():
        speciestype[k]=v[2]
        if v[2] not in  unique_species_types:
            unique_species_types.append(v[2])


    specific_species_type=['Complex','CandidateSet','DefinedSet']

    ###

    ## read Reactome DB
    rxnInfo,species=readReactomedb(data_folder)
    ## build complete DiGraph
    G=buildGraph(data_folder,rxnList=None)
    # remove hub nodes
    G,oldG=remove_hub_modes(G,data_folder)

    il2ra='R-HSA-450046'
    rxns=getSuccessors(G, il2ra,depth=30)

    ## merge list of lists
    successors_rxns = [item for sublist in rxns for item in sublist]

    ##  make a subgraph for given reactions
    subG=buildGraph(data_folder,successors_rxns)
    subG,old_subG=remove_hub_modes(subG,data_folder)
    # get source or input node of a subgraph
    inputNodeNames=getSource(subG,rxnInfo,successors_rxns)

    ## select  specific species type
    specific_species_type=['Complex','CandidateSet','DefinedSet']

    ## get predicessors
    # il2rc_jak1='R-HSA-451905'
    list_of_predecessors_rxns=[]
    sources={}
    for  source in inputNodeNames.keys():
        if source in speciestype.keys() and speciestype[source] in specific_species_type:
            rxns=getPredicessors(G, source,depth=20)
            predicessors_rxns = [item for sublist in rxns for item in sublist]
            sources[source]=predicessors_rxns
            list_of_predecessors_rxns.append(predicessors_rxns)
    ## final all
    # all_predicessors_rxns=[item for sublist in list_of_predecessors_rxns for item in sublist]
    import pdb; pdb.set_trace()



if "__main__==__name__":
    stepwiseGraphAnalysis(data_folder)
