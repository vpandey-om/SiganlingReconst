import json
import requests

class CyNetwork:
    def __init__(self):
        PORT_NUMBER = 1234
        self.BASE = 'http://localhost:' + str(PORT_NUMBER) + '/v1/'
        self.HEADERS = {'Content-Type': 'application/json'}

    def createBasedOnNodeEdge(self,node_list,edge_list):
        #node list and edge_list are lists and each element of the list is an another list
        # node list contains node_attributes
        # edge list contains edges attributes
        ## now create node lists
        n_id=node_list[0][0]
        n_name=node_list[0][1]
        n_size=node_list[0][2]
        n_shape=node_list[0][3]
        n_color=node_list[0][4]
        n_font=node_list[0][5]

        # edge arribute
        source=edge_list[0][0]
        target=edge_list[0][1]
        line_type=edge_list[0][2]
        e_width=edge_list[0][3]
        source_arrow_shape=edge_list[0][4]
        e_paint=edge_list[0][5]
        eid = edge_list[0][6]

        node_list = []
        for i, item in enumerate(n_id):
            res = self.createNode1(n_id[i], n_name[i], n_size[i], n_shape[i], n_color[i], n_font[i])
            node_list.append(res)
        ## create edge list
        edge_list = []
        for i, item in enumerate(source):
            res = self.createEdge(source[i], target[i], line_type[i], e_width[i], source_arrow_shape[i], e_paint[i],eid[i])
            edge_list.append(res)

        return node_list, edge_list
    # we want to create cytoscape based on nodes and edges
    def createNodeEdgeCofactor(self,reactions_info,node_label,newIds):
        # this function is for creating nodes and their attributes
        # node lable is labels for node
        n_id = []
        n_name = []
        n_size = []
        n_color = []
        n_shape = []
        n_font=[]
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
                #import pdb;pdb.set_trace()
                if rxn in node_label.keys():
                    n_name.append(node_label[rxn])
                else:
                    n_name.append(rxn)

                #n_name.append(reactions_info[rxn]['name'])
                n_size.append(15) ## node size 5
                n_font.append(6)
                n_shape.append('Ellipse')
                n_color.append(1)

                for key, value in species_dict.items():
                    ## substrate and product
                    # import pdb;pdb.set_trace()
                    if key in ['S', 'P', 'E', 'I', 'TF-A', 'M', 'TF-R']:
                        for item in value:
                            if item not in n_id:


                                n_id.append(item)
                                # n_name.append(item[0][0]) #str(metNames[i][0][0])
                                if item in node_label.keys():
                                    n_name.append(node_label[item])
                                else:
                                    n_name.append(item)
                                if item in newIds:
                                    n_size.append(40)  ## node size 5
                                    n_font.append(8)
                                else:
                                    n_size.append(80)  ## node size 50
                                    n_font.append(12)
                                n_shape.append('RECTANGLE')
                                n_color.append(0)
                                # now we want to add different type of arrows in edge list
                            if key in ['S']:
                                source.append(item)
                                target.append(rxn)
                                eid.append(item + ' (-) ' + rxn)
                                line_type.append('SOLID')
                                e_width.append(2.0)
                                source_arrow_shape.append('ARROW')

                                e_paint.append(0)

                            elif key in ['P']:
                                source.append(rxn)
                                target.append(item)
                                eid.append(rxn + ' (-) ' + item)
                                line_type.append('SOLID')
                                e_width.append(2.0)
                                source_arrow_shape.append('ARROW')

                                e_paint.append(0)

                            elif key in ['E', 'M', 'TF-A']:
                                source.append(item)
                                target.append(rxn)
                                eid.append(item + ' (-) ' + rxn)
                                line_type.append('DOT')
                                e_width.append(2.0)
                                source_arrow_shape.append('DIAMOND')

                                e_paint.append(0)

                            elif key in ['I', 'TF-R']:
                                source.append(item)
                                target.append(rxn)
                                eid.append(item + ' (-) ' + rxn)
                                line_type.append('DOT')
                                e_width.append(2.0)
                                source_arrow_shape.append('T')
                                e_paint.append(0)

        ## now create node lists
        node_list = []
        for i, item in enumerate(n_id):
            res = self.createNode1(n_id[i], n_name[i], n_size[i], n_shape[i], n_color[i],n_font[i])
            node_list.append(res)
        ## create edge list
        edge_list = []
        for i, item in enumerate(source):
            res = self.createEdge(source[i], target[i], line_type[i], e_width[i], source_arrow_shape[i], e_paint[i], eid[i])
            edge_list.append(res)

        return node_list, edge_list


    def createNodeEdgeCofactorSize(self,reactions_info,node_label,newIds):
        # this function is for creating nodes and their attributes
        # node lable is labels for node
        n_id = []
        n_name = []
        n_size = []
        n_color = []
        n_shape = []
        n_font=[]
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
                #import pdb;pdb.set_trace()
                if rxn in node_label.keys():
                    n_name.append(node_label[rxn])
                else:
                    n_name.append(rxn)

                #n_name.append(reactions_info[rxn]['name'])
                n_size.append(5) ## node size 5
                n_font.append(6)
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
                                if item in node_label.keys():
                                    n_name.append(node_label[item])
                                else:
                                    n_name.append(item)
                                if item in newIds:
                                    n_size.append(40)  ## node size 5
                                    n_font.append(8)
                                else:
                                    n_size.append(80)  ## node size 50
                                    n_font.append(12)
                                n_shape.append('RECTANGLE')
                                n_color.append(0)
                                # now we want to add different type of arrows in edge list
                            if key in ['S']:
                                source.append(item)
                                target.append(rxn)
                                eid.append(item + ' (-) ' + rxn)
                                line_type.append('SOLID')
                                e_width.append(2.0)
                                source_arrow_shape.append('ARROW')

                                e_paint.append(0)

                            elif key in ['P']:
                                source.append(rxn)
                                target.append(item)
                                eid.append(rxn + ' (-) ' + item)
                                line_type.append('SOLID')
                                e_width.append(2.0)
                                source_arrow_shape.append('ARROW')

                                e_paint.append(0)

                            elif key in ['E', 'M', 'TF-A']:
                                source.append(item)
                                target.append(rxn)
                                eid.append(item + ' (-) ' + rxn)
                                line_type.append('DOT')
                                e_width.append(2.0)
                                source_arrow_shape.append('DIAMOND')

                                e_paint.append(0)

                            elif key in ['I', 'TF-R']:
                                source.append(item)
                                target.append(rxn)
                                eid.append(item + ' (-) ' + rxn)
                                line_type.append('DOT')
                                e_width.append(2.0)
                                source_arrow_shape.append('T')
                                e_paint.append(0)

        ## now create node lists
        node_list = []
        for i, item in enumerate(n_id):
            res = self.createNode1(n_id[i], n_name[i], n_size[i], n_shape[i], n_color[i],n_font[i])
            node_list.append(res)
        ## create edge list
        edge_list = []
        for i, item in enumerate(source):
            res = self.createEdge(source[i], target[i], line_type[i], e_width[i], source_arrow_shape[i], e_paint[i], eid[i])
            edge_list.append(res)

        return node_list, edge_list







    def createNodeEdgeList(self,reactions_info,node_label):
        # this function is for creating nodes and their attributes
        # node lable is labels for node
        n_id = []
        n_name = []
        n_size = []
        n_color = []
        n_shape = []
        # edge list attribute
        #n_font=[]
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
                #import pdb;pdb.set_trace()
                if rxn in node_label.keys():
                    n_name.append(node_label[rxn])
                else:
                    n_name.append(rxn)

                #n_name.append(reactions_info[rxn]['name'])
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
                                if item in node_label.keys():
                                    n_name.append(node_label[item])

                                else:
                                    n_name.append(item)
                                n_size.append(15)
                                n_shape.append('Ellipse')
                                n_color.append(0)
                                # now we want to add different type of arrows in edge list
                            if key in ['S']:
                                source.append(item)
                                target.append(rxn)
                                eid.append(item + ' (-) ' + rxn)
                                line_type.append('SOLID')
                                e_width.append(2.0)
                                source_arrow_shape.append('ARROW')

                                e_paint.append(0)

                            elif key in ['P']:
                                source.append(rxn)
                                target.append(item)
                                eid.append(rxn + ' (-) ' + item)
                                line_type.append('SOLID')
                                e_width.append(2.0)
                                source_arrow_shape.append('ARROW')

                                e_paint.append(0)

                            elif key in ['E', 'M', 'TF-A']:
                                source.append(item)
                                target.append(rxn)
                                eid.append(item + ' (-) ' + rxn)
                                line_type.append('DOT')
                                e_width.append(2.0)
                                source_arrow_shape.append('DIAMOND')

                                e_paint.append(0)

                            elif key in ['I', 'TF-R']:
                                source.append(item)
                                target.append(rxn)
                                eid.append(item + ' (-) ' + rxn)
                                line_type.append('DOT')
                                e_width.append(2.0)
                                source_arrow_shape.append('T')
                                e_paint.append(0)

        ## now create node lists
        node_list = []
        for i, item in enumerate(n_id):
            res = self.createNode(n_id[i], n_name[i], n_size[i], n_shape[i], n_color[i])
            node_list.append(res)
        ## create edge list
        edge_list = []
        for i, item in enumerate(source):
            res = self.createEdge(source[i], target[i], line_type[i], e_width[i], source_arrow_shape[i], e_paint[i], eid[i])
            edge_list.append(res)

        return node_list, edge_list

    def createNode(self,n_id, n_name, n_size, n_shape, n_color):
        # we will make data dictionary using all node attributes

        return {'data': {'id': n_id, 'n_name': n_name, 'n_size': n_size, 'n_shape': n_shape, 'n_color': n_color}}

    def createNode1(self,n_id, n_name, n_size, n_shape, n_color,n_font):
        # we will make data dictionary using all node attributes

        return {'data': {'id': n_id, 'n_name': n_name, 'n_size': n_size, 'n_shape': n_shape, 'n_color': n_color,'n_font':n_font}}

    def createEdge(self,source, target, line_type, e_width, source_arrow_shape, e_paint, eid):
        # we will make data dictionary using all node attributes

        return {'data': {'source': source, 'target': target, 'line_type': line_type,
                         'source_arrow_shape': source_arrow_shape, 'e_paint': e_paint, 'e_width': e_width, 'eid': eid,
                         'e_label': eid}}

    def addNewNode(self,node,node_style=None):

        # function will add new node in th exiting network the network
        # input: list of nodes ids ['D','E','F']
        if node_style==None:
            node_style={'id':node[0],'n_name':node[0],'n_size':20,'n_color':0,'n_shape':'RECTANGLE','n_font':12}

        res = requests.get(self.BASE + 'networks')
        res3_dict = json.loads(res.content)
        net_suid = res3_dict[0]
        names=self.getContentNodes(net_suid) # this is for getting content of nodes in network
        # if name matches then not add or add
        if node not in names.keys():
            res2 = requests.post(self.BASE + 'networks/' + str(net_suid) + '/nodes', data=json.dumps(node), headers=self.HEADERS)

            new_table_data = {
                "key": "name",
                "dataKey": "id",
                "data": [node_style]
            }
            update_table_url = self.BASE + "networks/" + str(net_suid) + "/tables/defaultnode"
            res = requests.put(update_table_url, data=json.dumps(new_table_data), headers=self.HEADERS)
            print res
            style_name = 'Changed Visual Style'
            my_style = self.getVisualStyle(style_name)



            # Delete all style
            requests.delete(self.BASE + "styles")

            # Create new Visual Style
            requests.post(self.BASE + "styles", data=json.dumps(my_style), headers=self.HEADERS)

            # Apply it to current netwrok
            requests.get(self.BASE + 'apply/styles/' + style_name + '/' + str(net_suid))

    def updateNode(self,node_style):

        # function will add new node in th exiting network the network
        # input: list of nodes ids ['D','E','F']
        res = requests.get(self.BASE + 'networks')
        res3_dict = json.loads(res.content)
        net_suid = res3_dict[0]
        names=self.getContentNodes(net_suid) # this is for getting content of nodes in network
        # if name matches then not add or add
        if node_style['id'] in names.keys():

            new_table_data = {
                "key": "name",
                "dataKey": "id",
                "data": [node_style]
            }
            update_table_url = self.BASE + "networks/" + str(net_suid) + "/tables/defaultnode"
            res = requests.put(update_table_url, data=json.dumps(new_table_data), headers=self.HEADERS)
            print res
            style_name = 'Changed Visual Style'
            my_style = self.getVisualStyle(style_name)



            # Delete all style
            requests.delete(self.BASE + "styles")

            # Create new Visual Style
            requests.post(self.BASE + "styles", data=json.dumps(my_style), headers=self.HEADERS)

            # Apply it to current netwrok
            requests.get(self.BASE + 'apply/styles/' + style_name + '/' + str(net_suid))
        else:
            print '%s is not in existing network'%node_style['id']





    def updateEdge(self,edge_style):

        res = requests.get(self.BASE + 'networks')
        res3_dict = json.loads(res.content)
        net_suid = res3_dict[0]

        # get exixting edges
        names,edge_suid=self.getContentEdges(net_suid)

        if (edge_style['eid']  in names):


            new_table_data = {
                "key": "eid",
                "dataKey": "eid",
                "data": [edge_style]
            }
            update_table_url = self.BASE + "networks/" + str(net_suid) + "/tables/defaultedge"
            res = requests.put(update_table_url, data=json.dumps(new_table_data), headers=self.HEADERS)
            print res
            style_name = 'Changed Visual Style'
            my_style = self.getVisualStyle(style_name)



            # Delete all style
            requests.delete(self.BASE + "styles")

            # Create new Visual Style
            requests.post(self.BASE + "styles", data=json.dumps(my_style), headers=self.HEADERS)
        else:
            print '%s is not exist in network' %edge_style['eid']

    def addNewEdge(self,edge,edge_style=None):



        # function will add new node in th exiting network the network
        # input: list edge
        eid=edge[0]+' (-) '+edge[1]
        # get suids for souece and target




        # get existing network
        res = requests.get(self.BASE + 'networks')
        res3_dict = json.loads(res.content)
        net_suid = res3_dict[0]

        # get exixting edges
        names,edge_suid=self.getContentEdges(net_suid)
        nodes_suids = self.getContentNodes(net_suid)
        # identify source and target
        fg=1;
        if edge[0] in nodes_suids.keys():
            source=nodes_suids[edge[0]];
        else:
            fg=0;
            print edge[0]
            print 'is not in the network'

        if edge[1] in nodes_suids.keys():
            target=nodes_suids[edge[1]];
        else:
            fg=0;
            print edge[1]
            print 'is not in the network'


        if fg==1 and (eid not in names):
            edges = [{'source': source, 'target': target}]
            res = requests.post(self.BASE + 'networks/' + str(net_suid) + '/edges', data=json.dumps(edges), headers=self.HEADERS)
            print res

            if edge_style == None:
                edge_style = {'source': source, 'target': target, 'line_type': 'SOLID',
                              'source_arrow_shape': 'ARROW', 'e_paint': 0, 'e_width': 2.0, 'eid': eid,
                              'e_label': eid}


            new_table_data = {
                "key": "name",
                "dataKey": "eid",
                "data": [edge_style]
            }
            update_table_url = self.BASE + "networks/" + str(net_suid) + "/tables/defaultedge"
            res = requests.put(update_table_url, data=json.dumps(new_table_data), headers=self.HEADERS)
            print res
            style_name = 'Changed Visual Style'
            my_style = self.getVisualStyle(style_name)



            # Delete all style
            requests.delete(self.BASE + "styles")

            # Create new Visual Style
            requests.post(self.BASE + "styles", data=json.dumps(my_style), headers=self.HEADERS)

    def cytoscapeDiagram(self,node_list, edge_list):
        # based on node and edge list we need to plot cytoscape diagram

        BASE=self.BASE
        HEADERS=self.HEADERS
        # Header for posting data to the server as JSON
        # HEADERS = {'Content-Type': 'application/json'}
        # Get server status
        res = requests.get(BASE)
        print res

        empty_network = {
            'data': {
                'name': 'I\'m empty!'
            },
            'elements': {
                'nodes': [],
                'edges': []
            }
        }

        res3 = requests.post(BASE + 'networks?collection=My%20Collection', data=json.dumps(empty_network),
                             headers=HEADERS)
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

        res3 = requests.post(BASE + 'networks?collection=My%20Collection', data=json.dumps(small_network),
                             headers=HEADERS)
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

    def createCysNet(self,node_list, edge_list):
        # based on node and edge list we need to plot cytoscape diagram

        BASE=self.BASE
        HEADERS=self.HEADERS
        # Header for posting data to the server as JSON
        # HEADERS = {'Content-Type': 'application/json'}
        # Get server status
        res = requests.get(BASE)
        print res

        empty_network = {
            'data': {
                'name': 'I\'m empty!'
            },
            'elements': {
                'nodes': [],
                'edges': []
            }
        }

        res3 = requests.post(BASE + 'networks?collection=My%20Collection', data=json.dumps(empty_network),
                             headers=HEADERS)
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

        res3 = requests.post(BASE + 'networks?collection=My%20Collection', data=json.dumps(small_network),
                             headers=HEADERS)
        res3_dict = json.loads(res3.content)
        new_suid = res3_dict['networkSUID']

        # Apply layout
        requests.get(BASE + 'apply/layouts/force-directed/' + str(new_suid))
        res = requests.get(BASE + 'styles/default')
        d_style = json.loads(res.content)
        d_style_defalts = d_style['defaults']
        style_name = 'Changed Visual Style'
        my_style=self.getVisualStyle(style_name)

        res1 = json.loads(res.content)

        # Delete all style
        requests.delete(BASE + "styles")

        # Create new Visual Style
        requests.post(BASE + "styles", data=json.dumps(my_style), headers=HEADERS)

        # Apply it to current netwrok
        requests.get(BASE + 'apply/styles/' + style_name + '/' + str(net_suid))


    def exportReactionInfo(self):
        # get nodes and edges
        res = requests.get(self.BASE + 'networks')
        res3_dict = json.loads(res.content)
        net_suid = res3_dict[0]

        # get nodes
        res = requests.get(self.BASE + 'networks/' + str(net_suid))
        res_dict = json.loads(res.content)
        node_list = res_dict['elements']['nodes']
        # build node_dictionary

        node_dict={}
        for node in node_list:
            node_dict[node['data']['id_original']]=node['data']

        res = requests.get(self.BASE + 'networks/' + str(net_suid))
        res_dict = json.loads(res.content)
        edge_list = res_dict['elements']['edges']
        edge_dict = {}
        for edge in edge_list:
            edge_dict[edge['data']['eid']] = edge['data']

        # find  reaction nodes
        r_nodes = filter(lambda s: s.startswith('reaction_'), node_dict.keys())
        # find source and targets
        edges = map(lambda s: s.split(' (-) '), edge_dict.keys())
        source=[]
        target=[]
        for edge in edges:
            if len(edge)==2:
                source.append(edge[0])
                target.append(edge[1])
            else:
                print edge

        # find products
        edgeKeys=edge_dict.keys()
        rxnInfo={}
        for item in r_nodes:
            #prods = filter(lambda s: s.startswith(item), source)


            P=[]
            pinds = [i for i, x in enumerate(source) if x == item]
            P=[target[i] for i in pinds]
            PN=[node_dict[target[i]]['n_name'] for i in pinds]
            # substrate,enzyme, and inhibitors
            SEIinds = [i for i, x in enumerate(target) if x == item]
            SEI = [source[i] for i in SEIinds]
            S=[]
            E=[]
            I=[]
            SN=[]
            EN=[]
            IN=[]
            for ind in SEIinds:

                if (edge_dict[edgeKeys[ind]]['line_type']=='SOLID') and (edge_dict[edgeKeys[ind]]['source_arrow_shape']=='ARROW'):
                    S.append(source[ind])

                    SN.append(node_dict[source[ind]]['n_name'])
                elif (edge_dict[edgeKeys[ind]]['line_type']=='DOT') and (edge_dict[edgeKeys[ind]]['source_arrow_shape']=='DIAMOND'):
                    E.append(source[ind])
                    EN.append(node_dict[source[ind]]['n_name'])
                elif (edge_dict[edgeKeys[ind]]['line_type']=='DOT') and (edge_dict[edgeKeys[ind]]['source_arrow_shape']=='T'):
                    I.append(source[ind])
                    IN.append(node_dict[source[ind]]['n_name'])
                else:
                    print 'cannot recognize'
            rxnInfo[item]={'S':S,'E':E,'P':P,'I':I,'SN':SN,'EN':EN,'PN':PN,'IN':IN}

        return rxnInfo


    def getContentNodes(self,net_suid):
        # get nodes and edge content
        res = requests.get(self.BASE + 'networks/' + str(net_suid))
        res_dict = json.loads(res.content)
        node_list = res_dict['elements']['nodes']
        name_suid = {}

        for node in node_list:
            name_suid[node['data']['name']] = node['data']['SUID']

        return name_suid

    def getContentEdges(self,net_suid):

        # for edges
        res = requests.get(self.BASE + 'networks/' + str(net_suid))
        res_dict = json.loads(res.content)
        edge_list = res_dict['elements']['edges']
        eids=[]
        edge_suid = {}
        suids=[];
        for edge in edge_list:
            edge_suid[(edge['data']['source'], edge['data']['target'])] = edge['data']['SUID']
            eids.append(edge['data']['eid'])
            suids.append(edge['data']['SUID'])
        return eids,suids

    def getVisualStyle(self,style_name):


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


        return my_style





