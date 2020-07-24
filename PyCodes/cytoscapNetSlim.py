
import json
import requests


# PORT_NUMBER = 1234
# BASE = 'http://localhost:' + str(PORT_NUMBER) + '/v1/'
# Header for posting data to the server as JSON
HEADERS = {'Content-Type': 'application/json'}
def createNodeEdgeList(reactions_info):
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
            n_shape.append('RECTANGLE')
            n_color.append(0)

            for key, value in species_dict.items():
                ## substrate and product
                # import pdb;pdb.set_trace()
                if key in ['S', 'P', 'E', 'I', 'TF-A', 'M', 'TF-R']:
                    for item in value:
                        if item not in n_id:
                            n_id.append(item)
                            # n_name.append(item[0][0]) #str(metNames[i][0][0])
                            n_name.append(item)
                            n_size.append(15)
                            n_shape.append('Ellipse')
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


def createNode(n_id, n_name, n_size, n_shape, n_color):
    # we will make data dictionary using all node attributes

    return {'data': {'id': n_id, 'n_name': n_name, 'n_size': n_size, 'n_shape': n_shape, 'n_color': n_color}}


def createEdge(source, target, line_type, e_width, source_arrow_shape, e_paint, eid):
    # we will make data dictionary using all node attributes

    return {
        'data': {'source': source, 'target': target, 'line_type': line_type, 'source_arrow_shape': source_arrow_shape,
                 'e_paint': e_paint, 'e_width': e_width, 'eid': eid, 'e_label': eid}}


def cytoscapeDiagram(node_list, edge_list):
    # based on node and edge list we need to plot cytoscape diagram

    PORT_NUMBER = 1234
    BASE = 'http://localhost:' + str(PORT_NUMBER) + '/v1/'

    # Header for posting data to the server as JSON
    # HEADERS = {'Content-Type': 'application/json'}
    # Get server status
    res = requests.get(BASE)
    print
    res

    empty_network = {
        'data': {
            'name': 'I\'m empty!'
        },
        'elements': {
            'nodes': [],
            'edges': []
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
    d_style = json.loads(res.content)
    d_style_defalts = d_style['defaults']
    style_name = 'Changed Visual Style'
    my_style = {
        "title": style_name,
        "defaults": [{
            "visualProperty": "NODE_FILL_COLOR",
            "value": '#0099CC'
        }],
        "mappings": [{
            "mappingType": "passthrough",
            "mappingColumn": "n_name",
            "mappingColumnType": "String",
            "visualProperty": "NODE_LABEL"}, {
            "mappingType": "passthrough",
            "mappingColumn": "n_size",
            "mappingColumnType": "Double",
            "visualProperty": "NODE_SIZE"},  ## {
            ## "mappingType" : "passthrough",
            ## "mappingColumn" : "n_color",
            ## "mappingColumnType" : "Double",
            ## "visualProperty" : "NODE_FILL_COLOR"},
            {
                "mappingType": "passthrough",
                "mappingColumn": "n_shape",
                "mappingColumnType": "String",
                "visualProperty": "NODE_SHAPE"}, {
                "mappingType": "passthrough",
                "mappingColumn": "line_type",
                "mappingColumnType": "String",
                "visualProperty": "EDGE_LINE_TYPE"}, {
                "mappingType": "passthrough",
                "mappingColumn": "source_arrow_shape",
                "mappingColumnType": "String",
                "visualProperty": "EDGE_TARGET_ARROW_SHAPE"},  ## {
            ## "mappingType" : "passthrough",
            ## "mappingColumn" : "e_paint",
            ## "mappingColumnType" : "Double",
            ## "visualProperty" : "EDGE_STROKE_UNSELECTED_PAINT"},
            {
                "mappingType": "passthrough",
                "mappingColumn": "e_width",
                "mappingColumnType": "Double",
                "visualProperty": "EDGE_WIDTH"}  ## ,{
            ## "mappingType" : "passthrough",
            ## "mappingColumn" : "e_paint",
            ## "mappingColumnType" : "Double",
            ## "visualProperty" : "EDGE_TARGET_ARROW_UNSELECTED_PAINT"}
        ]}

    res1 = json.loads(res.content)

    # Delete all style
    requests.delete(BASE + "styles")

    # Create new Visual Style
    requests.post(BASE + "styles", data=json.dumps(my_style), headers=HEADERS)

    # Apply it to current netwrok
    requests.get(BASE + 'apply/styles/' + style_name + '/' + str(net_suid))
