from string import strip, split
from libsbml import *
import json
import requests
from IPython.display import Image
#PORT_NUMBER = 1234
#BASE = 'http://localhost:' + str(PORT_NUMBER) + '/v1/'
    # Header for posting data to the server as JSON
HEADERS = {'Content-Type': 'application/json'}
import scipy.io as inp
import pandas as pd


def readXml(f):
    reader = SBMLReader()
    document = reader.readSBML(f)
    model = document.getModel()
    print model.getNumSpecies()
    #print document.getNumErrors()
    ## this code to get species infromation
    out=open('species.txt','w')
    for ind in range(model.getNumSpecies()):
        sp=model.getSpecies(ind)
        ano_string=sp.annotation_string
        notes_string=sp.notes_string
        notes=notes_string.replace('\n','')
        out.write("%s\t%s\t%s"%(sp.id,sp.name,model.getCompartment(sp.compartment).name))
        arr=ano_string.split('urn:miriam:uniprot:')
        uniprots=[]
        
        for item in arr[1:]:
            uniprots.append(item.split('"/')[0])
        if len(uniprots)>0:
            uniIds=','.join(uniprots)
        else:
            uniIds='NA'

        arr=ano_string.split('urn:miriam:kegg.compound:')
        
        
        for item in arr[1:]:
            uniprots.append(item.split('"/')[0])
        if len(uniprots)>0:
            uniIds=','.join(uniprots)
        else:
            uniIds='NA'

        arr=ano_string.split('urn:miriam:obo.chebi:')
        
        for item in arr[1:]:
            uniprots.append(item.split('"/')[0])
        if len(uniprots)>0:
            uniIds=','.join(uniprots)
        else:
            uniIds='NA'

        out.write("\t%s\t%s\n"%(uniIds,notes))
        
    ## this code to get reactions infromation
    print model.getNumReactions()

    out=open('rxns.txt','w')
    
    for ind in range(model.getNumReactions()):
        sp=model.getReaction(ind)
        #print ind
        reactants=[]
        r_stoi=[]
        products=[]
        p_stoi=[]
        modifiers=[]
        mod_ids=[]
        #import pdb;pdb.set_trace()
        ano_string=sp.annotation_string
        arr=ano_string.split('urn:miriam:reactome:')
        if len(arr)==2:
            rid=arr[1].split('"/')[0];
        else:
            rid=sp.id
        out.write("%s\t%s\t%s\t%s"%(sp.id,rid,sp.name,sp.reversible))
    
        for s in range(len(sp.getListOfReactants())):
        
            reactants.append(sp.getReactant(s).species)
            r_stoi.append(str(sp.getReactant(s).stoichiometry))
            
        if len(reactants)>0:
            uniIds=','.join(reactants)
            r_coeff=','.join(r_stoi)
        else:
            uniIds='NA'
            r_coeff='NA'
        out.write("\t%s\t%s"%(uniIds,r_coeff))

        for p in range(len(sp.getListOfProducts())):
            products.append(sp.getProduct(p).species)
            p_stoi.append(str(sp.getProduct(p).stoichiometry))

        if len(products)>0:
            uniIds=','.join(products)
            p_coeff=','.join(p_stoi)
        else:
            uniIds='NA'
            p_coeff='NA'
        out.write("\t%s\t%s"%(uniIds,p_coeff))

        for m in range(len(sp.getListOfModifiers())):
            modifiers.append(sp.getModifier(m).species)
            mod_ids.append(sp.getModifier(m).id)
            #import pdb;pdb.set_trace()
        if len(modifiers)>0:
            uniIds=','.join(modifiers)
            mm=','.join(mod_ids)
        else:
            uniIds='NA'
            mm='NA'
        out.write("\t%s\t%s\n"%(uniIds,mm))

def parse_xml_reactome(f):
    # This function is for parsing xml from reactome database
    reader = SBMLReader()
    document = reader.readSBML(f)
    model = document.getModel()
    # this script for species
    species={}
    
    for ind in range(model.getNumSpecies()):
        sp=model.getSpecies(ind)
        ano_string=sp.annotation_string
        notes_string=sp.notes_string
        notes=notes_string.replace('\n','')
        #out.write("%s\t%s\t%s"%(sp.id,sp.name,model.getCompartment(sp.compartment).name))
        arr=ano_string.split('urn:miriam:uniprot:')
        uniprots=[]
        
        for item in arr[1:]:
            uniprots.append(item.split('"/')[0])
        if len(uniprots)>0:
            uniIds=','.join(uniprots)
        else:
            uniIds='NA'

        #out.write("\t%s\t%s\n"%(uniIds,notes))
        if sp.id not in species.keys():
            species[sp.id]=[]
        species[sp.id].append(sp.name)
        species[sp.id].append(model.getCompartment(sp.compartment).name)
        species[sp.id].append(uniIds)
        species[sp.id].append(notes)

    # this script for reactions
    reactions_info={}
    for ind in range(model.getNumReactions()):
        sp=model.getReaction(ind)
        #print ind
        #out.write("%s\t%s\t%s"%(sp.id,sp.name,sp.reversible))
        reactants=[]
        r_stoi=[]
        products=[]
        p_stoi=[]
        modifiers=[]
        mod_ids=[]
        
        for s in range(len(sp.getListOfReactants())):
        
            reactants.append(sp.getReactant(s).species)
            r_stoi.append(str(sp.getReactant(s).stoichiometry))
            
        if len(reactants)>0:
            uniIds=','.join(reactants)
            r_coeff=','.join(r_stoi)
        else:
            uniIds='NA'
            r_coeff='NA'
        #out.write("\t%s\t%s"%(uniIds,r_coeff))

        for p in range(len(sp.getListOfProducts())):
            products.append(sp.getProduct(p).species)
            p_stoi.append(str(sp.getProduct(p).stoichiometry))

        if len(products)>0:
            uniIds=','.join(products)
            p_coeff=','.join(p_stoi)
        else:
            uniIds='NA'
            p_coeff='NA'
        #out.write("\t%s\t%s"%(uniIds,p_coeff))

        for m in range(len(sp.getListOfModifiers())):
            modifiers.append(sp.getModifier(m).species)
            mod_ids.append(sp.getModifier(m).id)
            #import pdb;pdb.set_trace()
        if len(modifiers)>0:
            uniIds=','.join(modifiers)
            mm=','.join(mod_ids)
        else:
            uniIds='NA'
            mm='NA'
        #out.write("\t%s\t%s\n"%(uniIds,mm))
        # now get labels for substarte
        sl=[]
        for rea in reactants:
            if rea in species.keys():
                sl.append(species[rea][0])
            else:
                sl.append('NA')
        pl=[]
        for rea in products:
            if rea in species.keys():
                pl.append(species[rea][0])
            else:
                pl.append('NA')
        ml=[]
        for rea in modifiers:
            if rea in species.keys():
                ml.append(species[rea][0])
            else:
                ml.append('NA')
                    
        if sp.id not in reactions_info.keys():
            reactions_info[sp.id]={'name': sp.name,'reversible':sp.reversible,'Sid':reactants,'Pid':products,'Mid':modifiers,'SC':r_stoi,'PC':p_stoi,'action':mm,'S':sl,'P':pl,'M':ml}

    return species,reactions_info

def  test1():
    filename='/Users/vikashpandey/Documents/Jian/Signaling/signaling_2_09_16.xml' # this is the whole signaling model
    filename='/Users/vikashpandey/Documents/Jian/Signaling/EGFR.xml'
    species,reactions_info=parse_xml_reactome(filename)
    import pdb;pdb.set_trace()


def reactomeInCytoscape():
    #filename='/Users/vikashpandey/Documents/Jian/Signaling/EGFR.xml'
    filename='/Users/vikashpandey/Documents/MATLAB/MetSignal/REACTOME_SBML/NOTCH.xml'
    #filename='/Users/vikashpandey/Documents/Jian/Signaling/EGFR.xml'
    species,reactions_info=parse_xml_reactome(filename)
    node_list,edge_list=createNodeEdgeList(reactions_info)
    
    cytoscapeDiagram(node_list,edge_list)
    import pdb;pdb.set_trace()




def createNodeEdgeList(reactions_info):
    # this function is for creating nodes and their attributes 
    # nodes attribute
    n_id=[]
    n_name=[]
    n_size=[]
    n_color=[]
    n_shape=[]
    # edge list attribute
    source=[]
    target=[]
    eid=[]
    line_type=[]
    source_arrow_shape=[]
    e_width=[]
    e_paint=[]
    ######
    for rxn,species_dict in reactions_info.items():
        if rxn not in n_id:
            n_id.append(rxn)
            #n_name.append('') #str(item[0][0])
            n_name.append(reactions_info[rxn]['name'])
            n_size.append(5)
            n_shape.append('RECTANGLE')
            n_color.append(0)
            
            for key,value in species_dict.items():
                ## substrate and product
                #import pdb;pdb.set_trace()
                if key in ['S','P','E','I','TF-A','M','TF-R']:
                    for item in value:
                        if item not in n_id:
                            n_id.append(item)
                            #n_name.append(item[0][0]) #str(metNames[i][0][0])
                            n_name.append(item)
                            n_size.append(15)
                            n_shape.append('Ellipse')
                            n_color.append(0)
                            # now we want to add different type of arrows in edge list
                        if key in ['S']:
                            source.append(item)
                            target.append(rxn)
                            eid.append(item+'_'+rxn)
                            line_type.append('SOLID')
                            e_width.append(2.0)
                            source_arrow_shape.append('ARROW')

                            e_paint.append(0)

                        elif key in ['P']:
                            source.append(rxn)
                            target.append(item)
                            eid.append(rxn+'_'+item)
                            line_type.append('SOLID')
                            e_width.append(2.0)
                            source_arrow_shape.append('ARROW')

                            e_paint.append(0)

                        elif key in ['E','M','TF-A']:
                            source.append(item)
                            target.append(rxn)
                            eid.append(item +'_'+rxn)
                            line_type.append('DOT')
                            e_width.append(2.0)
                            source_arrow_shape.append('DIAMOND')

                            e_paint.append(0)

                        elif key in ['I','TF-R']:
                            source.append(item)
                            target.append(rxn)
                            eid.append(item+'_'+rxn)
                            line_type.append('DOT')
                            e_width.append(2.0)
                            source_arrow_shape.append('T')
                            e_paint.append(0)

    ## now create node lists
    node_list=[]
    for i,item in enumerate(n_id):
        res=createNode(n_id[i],n_name[i],n_size[i],n_shape[i],n_color[i])
        node_list.append(res)
    ## create edge list
    edge_list=[]
    for i,item in enumerate(source):
        res=createEdge(source[i],target[i],line_type[i],e_width[i],source_arrow_shape[i],e_paint[i],eid[i])
        edge_list.append(res)

    return node_list,edge_list


def createNode(n_id,n_name,n_size,n_shape,n_color):
    # we will make data dictionary using all node attributes
  
    return {'data': { 'id': n_id,'n_name':n_name,'n_size':n_size,'n_shape':n_shape,'n_color':n_color }}


def createEdge(source,target,line_type,e_width,source_arrow_shape,e_paint,eid):
    # we will make data dictionary using all node attributes
  
    return {'data': { 'source': source,'target':target,'line_type':line_type,'source_arrow_shape':source_arrow_shape,'e_paint':e_paint,'e_width':e_width,'eid':eid ,'e_label':eid}}




def cytoscapeDiagram(node_list,edge_list):
    # based on node and edge list we need to plot cytoscape diagram

    PORT_NUMBER = 1234
    BASE = 'http://localhost:' + str(PORT_NUMBER) + '/v1/'

# Header for posting data to the server as JSON
    #HEADERS = {'Content-Type': 'application/json'}
# Get server status
    res = requests.get(BASE)
    print res

    empty_network = {
        'data': {
            'name': 'I\'m empty!'
        },
        'elements': {
            'nodes':[ ],
            'edges':[ ]
        }
    }

    res3 = requests.post(BASE + 'networks?collection=My%20Collection', data=json.dumps(empty_network), headers=HEADERS)
    res3_dict = json.loads(res3.content)
    net_suid = res3_dict['networkSUID']
    print('Empty network has SUID ' + str(net_suid))

    import copy

    # working with nodes

    # Create a copy of the empty network object
    small_network = copy.deepcopy(empty_network)
   
   
    
    small_network['elements']['nodes'] = node_list
    small_network['elements']['edges'] = edge_list
    small_network['data']['name'] = 'example1'
    
    res3 = requests.post(BASE + 'networks?collection=My%20Collection', data=json.dumps(small_network), headers=HEADERS)
    res3_dict = json.loads(res3.content)
    new_suid = res3_dict['networkSUID']
   
    # Apply layout
    requests.get(BASE + 'apply/layouts/force-directed/' + str(new_suid))
    res = requests.get(BASE + 'styles/default')
    d_style=json.loads(res.content)
    d_style_defalts=d_style['defaults']
    style_name = 'Changed Visual Style'
    my_style = {
      "title" : style_name,
      "defaults" : [{
    "visualProperty" : "NODE_FILL_COLOR",
    "value" : '#0099CC'
  } ],
      "mappings" : [ {
        "mappingType" : "passthrough",
        "mappingColumn" : "n_name",
        "mappingColumnType" : "String",
        "visualProperty" : "NODE_LABEL"},{
        "mappingType" : "passthrough",
        "mappingColumn" : "n_size",
        "mappingColumnType" : "Double",
        "visualProperty" : "NODE_SIZE"},## {
        ## "mappingType" : "passthrough",
        ## "mappingColumn" : "n_color",
        ## "mappingColumnType" : "Double",
        ## "visualProperty" : "NODE_FILL_COLOR"},
        {
        "mappingType" : "passthrough",
        "mappingColumn" : "n_shape",
        "mappingColumnType" : "String",
        "visualProperty" : "NODE_SHAPE"},{
        "mappingType" : "passthrough",
        "mappingColumn" : "line_type",
        "mappingColumnType" : "String",
        "visualProperty" : "EDGE_LINE_TYPE"},{
        "mappingType" : "passthrough",
        "mappingColumn" : "source_arrow_shape",
        "mappingColumnType" : "String",
        "visualProperty" : "EDGE_TARGET_ARROW_SHAPE"},## {
        ## "mappingType" : "passthrough",
        ## "mappingColumn" : "e_paint",
        ## "mappingColumnType" : "Double",
        ## "visualProperty" : "EDGE_STROKE_UNSELECTED_PAINT"},
        {
        "mappingType" : "passthrough",
        "mappingColumn" : "e_width",
        "mappingColumnType" : "Double",
        "visualProperty" : "EDGE_WIDTH"}## ,{
        ## "mappingType" : "passthrough",
        ## "mappingColumn" : "e_paint",
        ## "mappingColumnType" : "Double",
        ## "visualProperty" : "EDGE_TARGET_ARROW_UNSELECTED_PAINT"}
        ]}


    res1=json.loads(res.content)


    #Delete all style
    requests.delete(BASE + "styles")

    #Create new Visual Style
    requests.post(BASE + "styles", data=json.dumps(my_style), headers=HEADERS)

    # Apply it to current netwrok
    requests.get(BASE + 'apply/styles/' + style_name + '/' + str(net_suid))


def read_netX_reactome(species,reactions_info):
    
    # species,reactions_info=parse_xml_reactome(f)
    # res=[species,reactions_info]
    # pik = open('species_reactions_signaling.p',"w")
    # pickle.dump(res,pik)
    # pik1 = open('species_reactions_signaling.p','r')
    # res1 = pickle.load(pik1)
    # species=res1[0]
    # reactions_info=res1[1]
    # import pdb;pdb.set_trace()
    
    n_id=[]
    n_name=[]
    edges=[]
    ## codes to build nodes and edge list
    for rxn,species_dict in reactions_info.items():
        if rxn not in n_id:
            n_id.append(rxn)
            
            n_name.append(reactions_info[rxn]['name'])
            for key,value in species_dict.items():
                
                if key in ['S','P','E','I','TF-A','M','TF-R']:
                    for item in value:
                        if item not in n_id:
                            n_id.append(item)
                            n_name.append(item)
                        if key in ['S','M','E','TF-A']:
                            edges.append((item,rxn))
                        else:
                            edges.append((rxn,item))
    return n_id,n_name,edges


def structureAnalysis():
    import networkx as nx
    import copy
    import pickle
    #f='/Users/vikashpandey/Documents/Jian/Signaling/EGFR.xml'
    #f='/Users/vikashpandey/Documents/Jian/Signaling/signaling_2_09_16.xml'
    #f='/Users/vikashpandey/Documents/MATLAB/MetSignal/REACTOME_SBML/complete_siganl_transduction.xml'
    f='/Users/vikashpandey/Documents/MATLAB/MetSignal/REACTOME_SBML/NOTCH.xml'
    species,reactions_info=parse_xml_reactome(f)
    res=[species,reactions_info]
    # pik = open('/Users/vikashpandey/Documents/MATLAB/MetSignal/REACTOME_SBML/Codes_xml/species_reactions_signaling_101116.p',"w")
    # pickle.dump(res,pik)
    #pik1 = open('/Users/vikashpandey/Documents/MATLAB/MetSignal/REACTOME_SBML/Codes_xml/species_reactions_signaling_101116.p')
    #res1 = pickle.load(pik1)
    species=res[0]
    reactions_info=res[1]
    
    n_id,n_name,edges=read_netX_reactome(species,reactions_info)
    G=nx.DiGraph()
    G.add_nodes_from(n_id)
    G.add_edges_from(edges)
    import pdb;pdb.set_trace()
    ls=nx.connected_components(G)
    components = [comp for comp in nx.connected_components(G)]
    component_size = [len(comp) for comp in components]
    print G.number_of_nodes(), G.number_of_edges(), component_size
    bigCompo=copy.deepcopy(G);
    rm_nodes=[]
    for comp in components[1:]:
        for item in comp:
            if item not in rm_nodes:
                rm_nodes.append(item)
    bigCompo.remove_nodes_from(rm_nodes)
    print bigCompo.number_of_nodes()
    out=open('degree_bigComponets.txt','w')
    for k,v in nx.degree(bigCompo).items():
       out.write("%s\t%d\n"%(k,v))
    

                            
def find_Connected_Compo():
    import networkx as nx
    import copy
    import pickle
    #f='/Users/vikashpandey/Documents/Jian/Signaling/EGFR.xml'
    #f='/Users/vikashpandey/Documents/Jian/Signaling/signaling_2_09_16.xml'
    f='/Users/vikashpandey/Documents/MATLAB/MetSignal/REACTOME_SBML/complete_siganl_transduction.xml'
    pik1 = open('species_reactions_signaling.p','r')
    res1 = pickle.load(pik1)
    species=res1[0]
    reactions_info=res1[1]
    
    n_id,n_name,edges=read_netX_reactome(species,reactions_info)
    G=nx.DiGraph()
    G.add_nodes_from(n_id)
    G.add_edges_from(edges)
    ls=nx.connected_components(G)
    components = [comp for comp in nx.connected_components(G)]
    component_size = [len(comp) for comp in components]
    print G.number_of_nodes(), G.number_of_edges(), component_size
    bigCompo=copy.deepcopy(G);
    rm_nodes=[]
    for comp in components[1:]:
        for item in comp:
            if item not in rm_nodes:
                rm_nodes.append(item)
    bigCompo.remove_nodes_from(rm_nodes)
    print bigCompo.number_of_nodes()
    out=open('degree_bigComponets.txt','w')
    for k,v in nx.degree(bigCompo).items():
       out.write("%s\t%d\n"%(k,v))
    import pdb;pdb.set_trace()
def find_ExtType():

    ###################
    # This code is for fining external componets type
    #################
    import pickle
    pik1 = open('/Users/vikashpandey/Documents/Jian/Signaling/species_reactions_signaling.p','r')
    res1 = pickle.load(pik1)
    species=res1[0]
    reactions_info=res1[1]
    
    f='/Users/vikashpandey/Documents/Jian/Signaling/external_componets.txt'
    ext_comp=read_1col(f)
    
    ex_sub=[] # this is for external substrate
    ex_prod=[] # this is for external product
    ex_modi=[] # this is for external modifiers
    
    for rxn,species_dict in reactions_info.items():
        for key,value in species_dict.items():
            #if key in ['S','P','E','I','TF-A','M','TF-R']:
            if key in ['S']:
                for item in value:
                    if (item in ext_comp) and (item not in ex_sub):
                        ex_sub.append(item)
            elif key in ['P']:
                for item in value:
                    if (item in ext_comp) and (item not in ex_prod):
                        ex_prod.append(item)
            elif  key in ['M']:
                for item in value:
                    if (item in ext_comp) and (item not in ex_modi):
                        ex_modi.append(item)
            else:
                continue
    print len((set(ex_modi)& set(ex_sub))),(set(ex_modi)& set(ex_sub)), len((set(ex_modi)& set(ex_prod))),(set(ex_modi)& set(ex_prod))
    print len(set(ext_comp)),len(ex_sub),len(ex_prod),len(ex_modi),len(ex_sub)+len(ex_prod)+len(ex_modi)
    
    f='/Users/vikashpandey/Documents/Jian/Signaling/external_substrate.txt'
    write_1col(f,ex_sub)
    f='/Users/vikashpandey/Documents/Jian/Signaling/external_product.txt'
    write_1col(f,ex_prod)
    f='/Users/vikashpandey/Documents/Jian/Signaling/external_modifier.txt'
    write_1col(f,ex_modi)
    
def read_2col(f):
    out=open(f,"r")
    data={};
    for line in out.readlines():
        if line:
            aa=map(strip,split(line,'\t'))
            if aa[0] not in data.keys():
                data[aa[0]]=aa[1]

    return data

def read_1col(f):
    out=open(f,"r")
    data=[];
    for line in out.readlines():
        if line:
            aa=map(strip,split(line,'\t'))
            if aa[0] not in data:
                data.append(aa[0])

    return data

def write_1col(f,data):
    out=open(f,"w")
    
    for item in data:
        out.write("%s\n"%item)
    out.close()


def parseXml(f,o):
    reader = SBMLReader()
    document = reader.readSBML(f)
    model = document.getModel()
    print model.getNumSpecies()
    #print document.getNumErrors()
    ## this code to get species infromation
    out=open(o+'species.txt','w')
    for ind in range(model.getNumSpecies()):
        sp=model.getSpecies(ind)
        ano_string=sp.annotation_string
        notes_string=sp.notes_string
        notes=notes_string.replace('\n','')
        out.write("%s\t%s\t%s"%(sp.id,sp.name,model.getCompartment(sp.compartment).name))
        arr=ano_string.split('urn:miriam:uniprot:')
        uniprots=[]
        
        for item in arr[1:]:
            uniprots.append(item.split('"/')[0])
        if len(uniprots)>0:
            uniIds=','.join(uniprots)
        else:
            uniIds='NA'

        arr=ano_string.split('urn:miriam:kegg.compound:')
        
        
        for item in arr[1:]:
            uniprots.append(item.split('"/')[0])
        if len(uniprots)>0:
            uniIds=','.join(uniprots)
        else:
            uniIds='NA'

        arr=ano_string.split('urn:miriam:obo.chebi:')
        
        for item in arr[1:]:
            uniprots.append(item.split('"/')[0])
        if len(uniprots)>0:
            uniIds=','.join(uniprots)
        else:
            uniIds='NA'

        out.write("\t%s\t%s\n"%(uniIds,notes))
        
    ## this code to get reactions infromation
    print model.getNumReactions()

    out=open(o+'rxns.txt','w')
    
    for ind in range(model.getNumReactions()):
        sp=model.getReaction(ind)
        #print ind
        reactants=[]
        r_stoi=[]
        products=[]
        p_stoi=[]
        modifiers=[]
        mod_ids=[]
        #import pdb;pdb.set_trace()
        ano_string=sp.annotation_string
        arr=ano_string.split('urn:miriam:reactome:')
        if len(arr)==2:
            rid=arr[1].split('"/')[0];
        else:
            rid=sp.id
        out.write("%s\t%s\t%s\t%s"%(sp.id,rid,sp.name,sp.reversible))
    
        for s in range(len(sp.getListOfReactants())):
        
            reactants.append(sp.getReactant(s).species)
            r_stoi.append(str(sp.getReactant(s).stoichiometry))
            
        if len(reactants)>0:
            uniIds=','.join(reactants)
            r_coeff=','.join(r_stoi)
        else:
            uniIds='NA'
            r_coeff='NA'
        out.write("\t%s\t%s"%(uniIds,r_coeff))

        for p in range(len(sp.getListOfProducts())):
            products.append(sp.getProduct(p).species)
            p_stoi.append(str(sp.getProduct(p).stoichiometry))

        if len(products)>0:
            uniIds=','.join(products)
            p_coeff=','.join(p_stoi)
        else:
            uniIds='NA'
            p_coeff='NA'
        out.write("\t%s\t%s"%(uniIds,p_coeff))

        for m in range(len(sp.getListOfModifiers())):
            modifiers.append(sp.getModifier(m).species)
            mod_ids.append(sp.getModifier(m).id)
            #import pdb;pdb.set_trace()
        if len(modifiers)>0:
            uniIds=','.join(modifiers)
            mm=','.join(mod_ids)
        else:
            uniIds='NA'
            mm='NA'
        out.write("\t%s\t%s\n"%(uniIds,mm))

    
def pathwayWise_parse():
    
    path='/Users/vikashpandey/Documents/MATLAB/MetSignal/REACTOME_SBML/'
    path_out='/Users/vikashpandey/Documents/MATLAB/MetSignal/REACTOME_SBML/Species_rxns/'
    files=read_1col(path+'name.txt')
    
    for item in files:
        parseXml(path+item,path_out+item)
    import pdb;pdb.set_trace()

def show_InCytoscope():
    
    #filename='/Users/vikashpandey/Documents/MATLAB/MetSignal/REACTOME_SBML/Activin.xml'
    filename = '/Users/vikash/Documents/MATLAB/MetSignal/REACTOME_SBML2/notch1.xml'
    species,reactions_info=parse_xml_reactome(filename)
    import pdb;pdb.set_trace()
    node_list,edge_list=createNodeEdgeList(reactions_info)
    cytoscapeDiagram(node_list,edge_list)
   
def updateXMLmodel():
    f='/Users/vikashpandey/Documents/MATLAB/MetSignal/REACTOME_SBML/Activin.xml'
    reader = SBMLReader()
    document = reader.readSBML(f)
    model = document.getModel()
    r1 = model.createReaction()
    check(r1,                                 'create reaction')
    check(r1.setId('r1'),                     'set reaction id')
    check(r1.setReversible(False),            'set reaction reversibility flag')
    check(r1.setFast(False),                  'set reaction "fast" attribute')
 
    species_ref1 = r1.createReactant()
    check(species_ref1,                       'create reactant')
    check(species_ref1.setSpecies('s1'),      'assign reactant species')
    check(species_ref1.setConstant(True),     'set "constant" on species ref 1')
 
    species_ref2 = r1.createProduct()
    check(species_ref2,                       'create product')
    check(species_ref2.setSpecies('s2'),      'assign product species')
    check(species_ref2.setConstant(True),     'set "constant" on species ref 2')
    import pdb;pdb.set_trace()

def check(value, message):
   """If 'value' is None, prints an error message constructed using
   'message' and then exits with status code 1.  If 'value' is an integer,
   it assumes it is a libSBML return status code.  If the code value is
   LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
   prints an error message constructed using 'message' along with text from
   libSBML explaining the meaning of the code, and exits with status code 1.
   """
   if value == None:
     raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
   elif type(value) is int:
     if value == LIBSBML_OPERATION_SUCCESS:
       return
     else:
       err_msg = 'Error encountered trying to ' + message + '.' \
                 + 'LibSBML returned error code ' + str(value) + ': "' \
                 + OperationReturnValue_toString(value).strip() + '"'
       raise SystemExit(err_msg)
   else:
     return

def notch1():
    filename = '/Users/vikash/Documents/MATLAB/MetSignal/REACTOME_SBML2/notch1.xml'
    species, reactions_info = parse_xml_reactome(filename)

    # species dictionary
    new_species={}
    for k,v in species.items():
        new_species[v[0]]=v[0].split('[')[0]
        print v[0].split('[')[0]
    import pdb;
    pdb.set_trace()
    node_list, edge_list = createNodeEdgeList2(reactions_info,new_species)
    cytoscapeDiagram(node_list, edge_list)


def createNodeEdgeList2(reactions_info,species):
    # this function is for creating nodes and their attributes
    # nodes attribute
    n_id = []
    n_name = []
    n_size = []
    n_color = []
    n_shape = []
    # edge list attribute
    source = []
    target = []
    eid = []
    line_type = []
    source_arrow_shape = []
    e_width = []
    e_paint = []
    ######
    for rxn, species_dict in reactions_info.items():
        if rxn not in n_id:
            n_id.append(rxn)
            # n_name.append('') #str(item[0][0])
            n_name.append(reactions_info[rxn]['name'])
            n_size.append(5)
            n_shape.append('Ellipse')
            n_color.append(0)

            for key, value in species_dict.items():
                ## substrate and product
                # import pdb;pdb.set_trace()
                if key in ['S', 'P', 'E', 'I', 'TF-A', 'M', 'TF-R']:
                    for item in value:
                        if item not in n_id:
                            n_id.append(item)
                            # n_name.append(item[0][0]) #str(metNames[i][0][0])
                            if item in species.keys():
                                n_name.append(species[item])
                            else:
                                n_name.append(item)
                            n_size.append(20)
                            n_shape.append('RECTANGLE')
                            n_color.append(0)
                            # now we want to add different type of arrows in edge list
                        if key in ['S']:
                            source.append(item)
                            target.append(rxn)
                            eid.append(item + '_' + rxn)
                            line_type.append('SOLID')
                            e_width.append(2.0)
                            source_arrow_shape.append('ARROW')

                            e_paint.append(0)

                        elif key in ['P']:
                            source.append(rxn)
                            target.append(item)
                            eid.append(rxn + '_' + item)
                            line_type.append('SOLID')
                            e_width.append(2.0)
                            source_arrow_shape.append('ARROW')

                            e_paint.append(0)

                        elif key in ['E', 'M', 'TF-A']:
                            source.append(item)
                            target.append(rxn)
                            eid.append(item + '_' + rxn)
                            line_type.append('DOT')
                            e_width.append(2.0)
                            source_arrow_shape.append('DIAMOND')

                            e_paint.append(0)

                        elif key in ['I', 'TF-R']:
                            source.append(item)
                            target.append(rxn)
                            eid.append(item + '_' + rxn)
                            line_type.append('DOT')
                            e_width.append(2.0)
                            source_arrow_shape.append('T')
                            e_paint.append(0)

    ## now create node lists
    node_list = []
    for i, item in enumerate(n_id):
        res = createNode(n_id[i], n_name[i], n_size[i], n_shape[i], n_color[i])
        node_list.append(res)
    ## create edge list
    edge_list = []
    for i, item in enumerate(source):
        res = createEdge(source[i], target[i], line_type[i], e_width[i], source_arrow_shape[i], e_paint[i], eid[i])
        edge_list.append(res)

    return node_list, edge_list


def parse_reaction():
    fin='/Users/vikash/Documents/MATLAB/MetSignal/homo_sapiens.2.xml'
    fout='/Users/vikash/Documents/MATLAB/MetSignal/homo_sapiens.2_out.txt'
    parseXml(fin, fout)

if "__main__==__name__":
    #updateXMLmodel()
    #show_InCytoscope()
    #pathwayWise_parse()
    #filename='/Users/vikashpandey/Documents/Jian/Signaling/signaling_2_09_16.xml'
    #readXml(filename)
    #test1()
    #reactomeInCytoscape()
    #structureAnalysis()
    #find_Connected_Compo()
    #find_ExtType()
    #notch1()

    parse_reaction()