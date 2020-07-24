from string import strip,split
from Cytoscape import *
#from libsbml import *


def notch1():
    cys=CyNetwork()
    filename = '/Users/vikash/Documents/MATLAB/MetSignal/REACTOME_SBML2/notch1.xml'
    species, reactions_info = parse_xml_reactome(filename)

    # species dictionary
    new_species={}
    for k,v in species.items():
        new_species[v[0]]=v[0].split('[')[0]
        print v[0].split('[')[0]
    import pdb;
    pdb.set_trace()
    node_list, edge_list = cys.createNodeEdgeList(reactions_info,new_species)
    cys.cytoscapeDiagram(node_list, edge_list)

def notch():
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


def parse_xml_reactome(f):
    # This function is for parsing xml from reactome database
    reader = SBMLReader()
    document = reader.readSBML(f)
    model = document.getModel()
    # this script for species
    species = {}

    for ind in range(model.getNumSpecies()):
        sp = model.getSpecies(ind)
        ano_string = sp.annotation_string
        notes_string = sp.notes_string
        notes = notes_string.replace('\n', '')
        # out.write("%s\t%s\t%s"%(sp.id,sp.name,model.getCompartment(sp.compartment).name))
        arr = ano_string.split('urn:miriam:uniprot:')
        uniprots = []

        for item in arr[1:]:
            uniprots.append(item.split('"/')[0])
        if len(uniprots) > 0:
            uniIds = ','.join(uniprots)
        else:
            uniIds = 'NA'

        # out.write("\t%s\t%s\n"%(uniIds,notes))
        if sp.id not in species.keys():
            species[sp.id] = []
        species[sp.id].append(sp.name)
        species[sp.id].append(model.getCompartment(sp.compartment).name)
        species[sp.id].append(uniIds)
        species[sp.id].append(notes)

    # this script for reactions
    reactions_info = {}
    for ind in range(model.getNumReactions()):
        sp = model.getReaction(ind)
        # print ind
        # out.write("%s\t%s\t%s"%(sp.id,sp.name,sp.reversible))
        reactants = []
        r_stoi = []
        products = []
        p_stoi = []
        modifiers = []
        mod_ids = []

        for s in range(len(sp.getListOfReactants())):
            reactants.append(sp.getReactant(s).species)
            r_stoi.append(str(sp.getReactant(s).stoichiometry))

        if len(reactants) > 0:
            uniIds = ','.join(reactants)
            r_coeff = ','.join(r_stoi)
        else:
            uniIds = 'NA'
            r_coeff = 'NA'
        # out.write("\t%s\t%s"%(uniIds,r_coeff))

        for p in range(len(sp.getListOfProducts())):
            products.append(sp.getProduct(p).species)
            p_stoi.append(str(sp.getProduct(p).stoichiometry))

        if len(products) > 0:
            uniIds = ','.join(products)
            p_coeff = ','.join(p_stoi)
        else:
            uniIds = 'NA'
            p_coeff = 'NA'
        # out.write("\t%s\t%s"%(uniIds,p_coeff))

        for m in range(len(sp.getListOfModifiers())):
            modifiers.append(sp.getModifier(m).species)
            mod_ids.append(sp.getModifier(m).id)
            # import pdb;pdb.set_trace()
        if len(modifiers) > 0:
            uniIds = ','.join(modifiers)
            mm = ','.join(mod_ids)
        else:
            uniIds = 'NA'
            mm = 'NA'
        # out.write("\t%s\t%s\n"%(uniIds,mm))
        # now get labels for substarte
        sl = []
        for rea in reactants:
            if rea in species.keys():
                sl.append(species[rea][0])
            else:
                sl.append('NA')
        pl = []
        for rea in products:
            if rea in species.keys():
                pl.append(species[rea][0])
            else:
                pl.append('NA')
        ml = []
        for rea in modifiers:
            if rea in species.keys():
                ml.append(species[rea][0])
            else:
                ml.append('NA')

        if sp.id not in reactions_info.keys():
            reactions_info[sp.id] = {'name': sp.name, 'reversible': sp.reversible, 'Sid': reactants, 'Pid': products,
                                     'Mid': modifiers, 'SC': r_stoi, 'PC': p_stoi, 'action': mm, 'S': sl, 'P': pl,
                                     'M': ml}

    return species, reactions_info

def regenReaction(rDict,species,cofactorsID):
    # rDict is reaction dictionaries
    # node_label :species
    # cofactors is going to be add for each reactions
    rInfo={}
    new_species={}
    newIds=[];

    for r,vals in rDict.items():
        sub = []
        prod = []
        I = []
        E = []
        for v in vals['S']:
            if v in species.keys():
                tmp=species[v]
            else:
                tmp=v;
            if v in cofactorsID:
                sub.append(v+r)
                new_species[v+r]=tmp;
                newIds.append(v+r)
            else:
                sub.append(v)
                new_species[v]=tmp;
        for v in vals['P']:
            if v in species.keys():
                tmp=species[v]
            else:
                tmp=v;
            if v in cofactorsID:
                prod.append(v+r)
                new_species[v+r]=tmp;
                newIds.append(v + r)
            else:
                prod.append(v)
                new_species[v]=tmp;
        if 'E' in vals.keys():
            for v in vals['E']:
                if v in species.keys():
                    tmp=species[v]
                else:
                    tmp=v
                if v in cofactorsID:
                    E.append(v+r)
                    new_species[v+r]=tmp;
                    newIds.append(v + r)
                else:
                    E.append(v)
                    new_species[v]=tmp;
        if 'I' in vals.keys():
            for v in vals['I']:
                if v in species.keys():
                    tmp=species[v]
                else:
                    tmp=v;
                if v in cofactorsID:
                    I.append(v+r)
                    new_species[v+r]=tmp;
                    newIds.append(v + r)
                else:
                    I.append(v)
                    new_species[v]=tmp;
        if r in species.keys():
            new_species[r]=species[r]
        else:
            new_species[r]=r
        rInfo[r]={'S':sub,'P':prod,'E':E,'I':I}

    return rInfo,new_species,newIds


def save_update_example():
    rInfo={}
    rInfo['R1']={'S':['A','B'],'P':['C']}
    new_species={}
    cofactorsID=[]
    cys = CyNetwork()
    rInfo, labels, newIds = regenReaction(rInfo, new_species, cofactorsID)
    node_list, edge_list = cys.createNodeEdgeCofactor(rInfo, labels, newIds)
    cys.cytoscapeDiagram(node_list, edge_list)


def addNodeExample():
    cys = CyNetwork()
    nodes = ['F']
    cys.addNewNode(nodes)

def addEdgeExample():
    cys = CyNetwork()
    edge=['F','A']
    cys.addNewEdge(edge)


def updateNodeExample():
    cys = CyNetwork()
    node_style = {'id': 'F', 'n_name': 'F', 'n_size': 20, 'n_color': 0, 'n_shape': 'Ellipse', 'n_font': 12}
    cys.updateNode(node_style)


def updateEdgeExample():
    cys = CyNetwork()
    edge_style = {'source': 'A', 'target': 'R1', 'line_type': 'DOT',
                  'source_arrow_shape': 'ARROW', 'e_paint': 0, 'e_width': 2.0, 'eid': 'A (-) R1',
                  'e_label': 'A (-) R1'}
    cys.updateEdge(edge_style)


def add_update():
    cys = CyNetwork()
    BASE=cys.BASE
    HEADERS=cys.HEADERS
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




if "__main__==__name__":
    print 'Executing ........'
    #notch()
    #addNodeExample()
    #save_update_example()
    #add_update()
    #addEdgeExample()
    #updateNodeExample()
    #updateEdgeExample()
