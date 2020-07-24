import json
import requests
#from IPython.display import Image
PORT_NUMBER = 1234
BASE = 'http://localhost:' + str(PORT_NUMBER) + '/v1/'
    # Header for posting data to the server as JSON
HEADERS = {'Content-Type': 'application/json'}
import scipy.io as inp
#import pandas as pd
from string import strip,split


# Utulity to POST object
def create(param, dict_data):
    res=requests.post(BASE + param, data=json.dumps(dict_data), headers=HEADERS)
    return res

def update(param, dict_data):
    res=requests.put(BASE + param, data=json.dumps(dict_data), headers=HEADERS)
    return res

def get_basic_mapping(map_type, column, column_type, vp):
    new_mapping = {
         'mappingType': map_type,
        'mappingColumn': column,
        'mappingColumnType': column_type,
        'visualProperty': vp,
    }
    return new_mapping
    
def get_discrete_mapping(column, column_type, vp):
    mapping = get_basic_mapping('discrete', column, column_type, vp)
    mapping['map'] = []
    return mapping

def get_continuous_mapping(column, column_type, vp):
    mapping = get_basic_mapping('continuous', column, column_type, vp)
    mapping['points'] = []
    return mapping

def get_passthrough_mapping(column, column_type, vp):
    return get_basic_mapping('passthrough', column, column_type, vp)


def readmatFile(matfile,varName):
    # this is a model structure file
    #varName='toymodel'
    model_dict=inp.loadmat(matfile,struct_as_record=True)# this is a model dictionary
    model_info=model_dict[varName]
    metNames=model_dict[varName]['metNames'][0][0]
    rxns=model_dict[varName]['rxns'][0][0]
    rxnNames=model_dict[varName]['rxnNames'][0][0]
    mets=model_dict[varName]['mets'][0][0]
    genes=model_dict[varName]['genes'][0][0]
    flux_data=model_dict[varName]['flux_data'][0][0][0].tolist()
    conc_data=model_dict[varName]['conc_data'][0][0][0].tolist()
    gene_data=model_dict[varName]['gene_data'][0][0][0].tolist()
    S=model_dict[varName]['S'][0][0]
    G=model_dict[varName]['G'][0][0]
    # now we want to collect all nodes from metabolic network which are simply : metaboliets genes and enzymes or genes
    # nodes attribute (){'data': {'id': 'a','name':'vikash','n_size':15,'n_color':'#0099CC','n_shape':'ROUND_RECTANGLE'}
    n_id=[]
    n_name=[]
    n_size=[]
    n_color=[]
    n_shape=[]
    #import pdb;pdb.set_trace()
    for i,item in enumerate(mets):
       
        n_id.append(str(item[0][0]))
        #n_name.append(item[0][0]) #str(metNames[i][0][0])
        n_name.append(str(metNames[i][0][0]))
        n_size.append(15)
        n_shape.append('Ellipse')
        n_color.append(conc_data[i])
    for i,item in enumerate(rxns):
        n_id.append(str(item[0][0]))
        n_name.append('') #str(item[0][0])
        #n_name.append(str(rxnNames[i][0][0]))
        n_size.append(5)
        n_shape.append('RECTANGLE')
        n_color.append(0)
    
    for i,item in enumerate(genes):
        n_id.append(str(item[0][0]))
        n_name.append(str(item[0][0]))
        n_size.append(15)
        n_shape.append('Triangle')
        n_color.append(gene_data[i])
 
    # function definition to make node list
       
    node_list=[]
    for i,item in enumerate(n_id):
        res=createNode(n_id[i],n_name[i],n_size[i],n_shape[i],n_color[i])
        node_list.append(res)

    # prepare to build edge list {'data': {  'source': 'a', 'target': 'b','line_type':'SOLID','source_arrow_shape':'ARROW','e_width':5,'e_paint':'#FF0000' }
    source=[]
    target=[]
    eid=[]
    line_type=[]
    source_arrow_shape=[]
    e_width=[]
    e_paint=[]
   
    for i,item in enumerate(rxns):
        reactants=(S[:,i]<0).nonzero()[0]
        for j in range(reactants.shape[0]):
            source.append(str(mets[reactants[j]][0][0]))
            target.append(str(rxns[i][0][0]))
            eid.append((str(mets[reactants[j]][0][0]))+'_'+str(rxns[i][0][0]))
            line_type.append('SOLID')
            e_width.append(2.0)
            source_arrow_shape.append('ARROW')
            #import pdb;pdb.set_trace()
            e_paint.append(flux_data[i])
          
      
        reactants=(S[:,i]>0).nonzero()[0]
        for j in range(reactants.shape[0]):
            target.append(str(mets[reactants[j]][0][0]))
            source.append(str(rxns[i][0][0]))
            eid.append((str(rxns[i][0][0])+'_'+str(mets[reactants[j]][0][0])))        
            line_type.append('SOLID')
            e_width.append(2.0)
            source_arrow_shape.append('ARROW')
            e_paint.append(flux_data[i])
          
           
        reactants=(G[i,:]>0).nonzero()[0]
        for j in range(reactants.shape[0]):
            source.append(str(genes[reactants[j]][0][0]))
            target.append(str(rxns[i][0][0]))
            eid.append((str(genes[reactants[j]][0][0]))+'_'+str(rxns[i][0][0]))
            line_type.append('DOT')
            e_width.append(2.0)
            source_arrow_shape.append('DIAMOND')
            e_paint.append(flux_data[i])
            
    
    edge_list=[]
    for i,item in enumerate(source):
        res=createEdge(source[i],target[i],line_type[i],e_width[i],source_arrow_shape[i],e_paint[i],eid[i])
        edge_list.append(res)
    return node_list,edge_list



def createNet_fromMat():

    ## this will create an empty network in cytoscape  and map information from mat file in    
    # Basic Setup
    #filename='/Users/vikashpandey/Documents/MATLAB/Visualize/SO_visualize.mat'
    #filename='/Users/vikashpandey/Documents/MATLAB/matfiles/SO_visualize.mat'
    #filename='/home/vikash/Programming/PyScripts/cyREST_copy/toy_visualize.mat'
    filename='/Users/vikashpandey/Documents/MATLAB/Visualize/Ceramide_visualize_human.mat'
    filename='/Users/vikashpandey/Documents/MATLAB/Visualize/phospho_visualize_human.mat'
    #filename='/Users/vikashpandey/Documents/MATLAB/Visualize/phospho_aj_visualize_human.mat'
    filename='/Users/vikashpandey/Documents/MATLAB/Visualize/phospho_aj_visualize_human.mat'
    filename='/Users/vikashpandey/Documents/MATLAB/Visualize/ps_hs_aj_visualize_human.mat'
    filename='/Users/vikashpandey/Documents/MATLAB/Visualize/ps_hs_pwd_visualize_human.mat'
    filename='/Users/vikashpandey/Documents/MATLAB/Visualize/ps_hs_pwd_aj_visualize_human.mat'
    # structure save as 'toymodel'
    #node_list,edge_list=readmatFile(filename,'toymodel') # this function is for building newtwork from .mat file
    #node_list,edge_list=readmatFile(filename,'SO')
    #node_list,edge_list=readmatFile(filename,'CR')
    ###
    # for puja data
    ###
    filename='/Users/vikashpandey/Documents/MATLAB/Puja/MiNEAFBA/Vis/tag.mat'
    filename='/Users/vikash/Documents/PyCytoscape/toy_visualize.mat'
    filename='/Users/vikashpandey/Documents/MATLAB/Nousin/A290517/Res190717/GEMVisualize/GEMecoli.mat'

    node_list,edge_list=readmatFile(filename,'toymodel')
    
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
        
   


def createNode(n_id,n_name,n_size,n_shape,n_color):
    # we will make data dictionary using all node attributes
  
    return {'data': { 'id': n_id,'n_name':n_name,'n_size':n_size,'n_shape':n_shape,'n_color':n_color }}

def createEdge(source,target,line_type,e_width,source_arrow_shape,e_paint,eid):
    # we will make data dictionary using all node attributes
  
    return {'data': { 'source': source,'target':target,'line_type':line_type,'source_arrow_shape':source_arrow_shape,'e_paint':e_paint,'e_width':e_width,'eid':eid ,'e_label':eid}}


def make_Ids(matfile,varName):
     # this is a model structure file
    node_list,edge_list=readmatFile(matfile,varName) # this function is for building newtwork from .mat file
    node_dict={}
    edge_dict={}
    for item in node_list:
        if item['data']['id'] not in node_dict.keys():
            node_dict[item['data']['id']]=item['data']
    for item in edge_list:
        edge=item['data']['source']+'_'+item['data']['target']
        if edge not in edge_dict.keys():
            edge_dict[edge]=item['data']
    return node_dict,edge_dict


def apply_data_to_Nodes(node_suids,BASE,net_suid,node_dict):
    label_list=[]
    color_list=[]
    size_list=[]
    import random
    
    for n_suid in node_suids:
        res = requests.get(BASE + 'networks/' + str(net_suid) + '/nodes/'+str(n_suid))
        res_dict = json.loads(res.content)
        if res_dict['data']['id'] in node_dict.keys():
           
            node_att=node_dict[res_dict['data']['id']]
            n_size=int(node_att['n_size']*random.random())
            label_entry={ 'key':str(n_suid),'value':node_att['n_name']}
            label_list.append(label_entry)
            #color_entry={ 'key':str(n_suid),'value':node_att['n_color']}
            #color_list.append(color_entry)
            size_entry={ 'key':str(n_suid),'value':100}
            size_list.append(size_entry)
    
    # color_mapping = get_discrete_mapping('SUID', 'Long', 'NODE_FILL_COLOR')
    # color_mapping = get_passthrough_mapping('SUID', 'Long', 'NODE_FILL_COLOR')
    # color_mapping = get_continuous_mapping('SUID', 'Long', 'NODE_FILL_COLOR')
    label_mapping=get_discrete_mapping('SUID', 'Long', 'NODE_SIZE')
    label_mapping['map'] = label_list
    size_mapping=get_discrete_mapping('SUID', 'Integer', 'NODE_SIZE')
    size_mapping['map'] = size_list
    #color_mapping=get_continuous_mapping('SUID', 'Double', 'NODE_FILL_COLOR')
    #color_mapping['points'] = color_list
    return [size_mapping]
    #return [label_mapping,size_mapping,color_mapping]
            
def add_defaultnode_table(net_suid,node_dict):
    import pdb;
    pdb.set_trace()
    res = requests.get(BASE + 'networks/' + str(net_suid) + '/nodes')
    node_suids = res.json()
    change_nodes=[]
    for n_suid in node_suids:
        res = requests.get(BASE + 'networks/' + str(net_suid) + '/nodes/'+str(n_suid))
        res_dict = json.loads(res.content)
        if res_dict['data']['id'] in node_dict.keys():
            entry={'id':node_dict[res_dict['data']['id']]['id'],'n_name':node_dict[res_dict['data']['id']]['n_name'],'new_color': node_dict[res_dict['data']['id']]['n_color']}
            change_nodes.append(entry)

    new_table_data = {
    "key": "id",
    "dataKey": "n_name",
    "data" : change_nodes
    }
    update_table_url =  BASE + "networks/" + str(net_suid) + "/tables/defaultnode"
    res=requests.put(update_table_url, data=json.dumps(new_table_data), headers=HEADERS)
    print res


def add_defaultedge_table(net_suid,node_dict): 
    res = requests.get(BASE + 'networks/' + str(net_suid) + '/edges')
    node_suids = res.json()
    change_nodes=[]
    for n_suid in node_suids:
        res = requests.get(BASE + 'networks/' + str(net_suid) + '/edges/'+str(n_suid))
        res_dict = json.loads(res.content)
        if res_dict['data']['eid'] in node_dict.keys():
            entry={'eid':node_dict[res_dict['data']['eid']]['eid'],'e_label':node_dict[res_dict['data']['eid']]['e_label'],'new_paint': node_dict[res_dict['data']['eid']]['e_paint']}
            change_nodes.append(entry)
            
    new_table_data = {
    "key": "eid",
    "dataKey": "e_label",
    "data" : change_nodes
    }
    update_table_url =  BASE + "networks/" + str(net_suid) + "/tables/defaultedge"
    res=requests.put(update_table_url, data=json.dumps(new_table_data), headers=HEADERS)
    print res

def load_and_update_cys():
    ## this function load network from cys file and then update what you have in mat file
    # first get information from cytoscape
    # then bulid python dictionary betwwen IDS and node id
    # parse your mat structure similarlily like building network
    # then map everytjing
    #then show network
    print "hi"
    # -------------------load network cytoscape-------------------------
     
   
    
    res = requests.get(BASE + 'networks')
    res3_dict = json.loads(res.content)

    net_suid = res3_dict[0]
    print net_suid

     # make ids for nodes and edges from mat files
    matfile='/Users/vikashpandey/Documents/MATLAB/matfiles/toy_visualize.mat'
    matfile='/Users/vikash/Documents/PyCytoscape/toy_visualize.mat'
    varName='toymodel'
    node_dict,edge_dict=make_Ids(matfile,varName)
    node_dict={'RASHMI': {'n_color': 0, 'n_name': 'BIBI', 'n_size': 5, 'id': 'RASHMI', 'n_shape': 'RECTANGLE'}}
    import pdb;pdb.set_trace()
    add_defaultnode_table(net_suid,node_dict)
   
    #add_defaultedge_table(net_suid,edge_dict)

def add_update():
    res = requests.get(BASE + 'networks')
    res3_dict = json.loads(res.content)
    net_suid = res3_dict[0]
    print net_suid
    res= requests.get(BASE + 'networks/'+str(net_suid)+'/nodes')
    res3_list = json.loads(res.content)
    change_nodes=[]
    res = requests.get(BASE + 'networks/' + str(net_suid) + '/nodes/'+str(res3_list[3]))
    res_dict = json.loads(res.content)
    import pdb;pdb.set_trace()
    entry={'id':5127,'n_name':'bibi','n_color':'#00FFCC'}
    new_table_data = { "key": "id", "dataKey": "n_name", "data" : change_nodes }
    update_table_url =  BASE + "networks/" + str(net_suid) + "/tables/defaultnode"
    res=requests.put(update_table_url, data=json.dumps(new_table_data), headers=HEADERS)
    
def load_add_node_edges():
    res = requests.get(BASE + 'networks')
    res3_dict = json.loads(res.content)
    import pdb; pdb.set_trace()
    net_suid = res3_dict[0]
    print net_suid
    nodes=['RASHMI']
    {'data': {u'name': u'vikash1', u'SUID': 5127, u'selected': False, u'id': u'a', u'n_size': 15.0, u'shared_name': u'vikash', u'n_color': u'#0099CC', u'n_shape': u'ROUND_RECTANGLE'}}

   
    res= requests.post(BASE + 'networks/'+str(net_suid)+'/nodes',data=json.dumps(nodes), headers=HEADERS)
    res= requests.get(BASE + 'networks/'+str(net_suid)+'/nodes')
    
    res3_list = json.loads(res.content)
    edges=[{'source': res3_list[0],'target':res3_list[3]}]
    res= requests.post(BASE + 'networks/'+str(net_suid)+'/edges',data=json.dumps(edges), headers=HEADERS)

def getContentNodes(net_suid):
    # get nodes and edge content
    res= requests.get(BASE + 'networks/'+str(net_suid))
    res_dict = json.loads(res.content)
    node_list=res_dict['elements']['nodes']
    name_suid={}
    
    for node in node_list:
        name_suid[node['data']['name']]=node['data']['SUID']
        
    return name_suid

def getContentEdges(net_suid):
  
    # for edges
    res= requests.get(BASE + 'networks/'+str(net_suid))
    res_dict = json.loads(res.content)
    edge_list=res_dict['elements']['edges']
    
    edge_suid={}
    for edge in edge_list:
        edge_suid[(edge['data']['source'],edge['data']['target'])]=edge['data']['SUID']
    return edge_suid


def add_subSys():
    filename='/Users/vikashpandey/Documents/MATLAB/Visualize/SO_visualize_human.mat'
    node_list,edge_list=readmatFile(filename,'SO')
    #import pdb;pdb.set_trace()
    ## get existing network
    
    res = requests.get(BASE + 'networks')
    res3_dict = json.loads(res.content)
    net_suid = res3_dict[0]
    print net_suid
    
    node_suid=getContentNodes(net_suid) ## store nodes and edges of existing network
    
    edge_suid=getContentEdges(net_suid)
    # compare how many are new
    
    add_nodes=[]
    up_nodes=[] ## update all
    for item in node_list:
        if item['data']['id'] not in node_suid.keys():
            add_nodes.append(item['data']['id'])
            up_nodes.append(item['data'])
        else:
            up_nodes.append(item['data'])
            
    res= requests.post(BASE + 'networks/'+str(net_suid)+'/nodes',data=json.dumps(add_nodes), headers=HEADERS)
    node_suid=getContentNodes(net_suid) ## get All suid
    
    add_edges=[] # [{"source": SOURCE_NODE_SUID, "target": TARGET_NODE_SUID}]
    up_edges=[] ## update all
    
    for item in edge_list:
        if (item['data']['source'],item['data']['target']) not in edge_suid.keys():
            new_edge={"source": node_suid[item['data']['source']], "target": node_suid[item['data']['target']],"interaction":item['data']['source']+'_'+item['data']['target']}
            add_edges.append(new_edge)
            up_edges.append(item['data'])
        else:
            up_edges.append(item['data'])
            
    res= requests.post(BASE + 'networks/'+str(net_suid)+'/edges',data=json.dumps(add_edges), headers=HEADERS)
    ## update shape color etc for add edge and add nodes
    import pdb;pdb.set_trace()
    new_table_data = { "key": "name", "dataKey": "id", "data" : up_nodes }
    update_table_url =  BASE + "networks/" + str(net_suid) + "/tables/defaultnode"
    res=requests.put(update_table_url, data=json.dumps(new_table_data), headers=HEADERS)
    # for edge
    new_table_data = { "key": "interaction", "dataKey": "eid", "data" : up_edges }
    update_table_url =  BASE + "networks/" + str(net_suid) + "/tables/defaultedge"
    res=requests.put(update_table_url, data=json.dumps(new_table_data), headers=HEADERS)
    print res

def del_subSys():
    filename='/Users/vikashpandey/Documents/MATLAB/Visualize/mouse_sub.mat'
    node_list,edge_list=readmatFile(filename,'human_sub')

def read_2col(f):
    out=open(f,"r")
    data={};
    for line in out.readlines():
        if line:
            aa=map(strip,split(line,'\t'))
            if aa[0] not in data.keys():
                data[aa[0]]=aa[1]

    return data
    

def update_attribute():
    #f='/Users/vikashpandey/Documents/MATLAB/Visualize/hfr.csv'
    f='/Users/vikashpandey/Documents/MATLAB/Puja/MiNEAFBA/Vis/tagTxt/101.txt'
    data=read_2col(f)
    res = requests.get(BASE + 'networks')
    res3_dict = json.loads(res.content)
    net_suid = res3_dict[0]
    print net_suid
    lab=[]
    node_suid=getContentNodes(net_suid)
    
    for k,v in node_suid.items():
        if k in data.keys():
            entry={"name":k,"labelCol":int(data[k])}
        else:
            entry={"name":k,"labelCol":int(0)}
            
        lab.append(entry)
    import pdb;pdb.set_trace()
    
    
    new_table_data = { "key": "name", "dataKey": "name", "data" :lab }
    
    
    update_table_url =  BASE + "networks/" + str(net_suid) + "/tables/defaultnode"
    res=requests.put(update_table_url, data=json.dumps(new_table_data), headers=HEADERS)




if "__main__==__name__":
    print 'hi'

    createNet_fromMat()
   
    #load_and_update_cys()
    
    #load_add_node_edges()
    # after add load and update
   # add_update()
    #add_subSys()
    #update_attribute()
 
