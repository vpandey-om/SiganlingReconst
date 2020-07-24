from Cytoscape import *
from string import strip, split
import pickle


def highConservedNet():
    # slist : vertex  vertex_label  edges edge_arrows
    slist=pickle.load(open('/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/ModuleAnalysis/STE_modules_withweight.pickle'))
    # tarnsfer picjle file

    nodes = slist[0]
    nlables = slist[1]
    edges = slist[2]
    arrow = slist[3]
    terminal=slist[4]
    # now we want to add node attributes
    n_id=[]
    n_name=[]
    n_size=[]
    n_shape=[]
    n_color=[]
    n_font=[]

    # edge arribute
    source=[]
    target=[]
    line_type=[]
    e_width=[]
    source_arrow_shape=[]
    e_paint=[]
    e_id=[]
    nodes_list=[]
    edges_list=[]

    for i, item in enumerate(nodes):
        n_id.append(item)
        if item in terminal:
            n_color.append(1)
        else:
            n_color.append(0)

        #n_name.append(nlables[i])
        if nodes[i]== nlables[i]:
            n_shape.append('Ellipse')

            n_size.append(8)
            n_font.append(6)
            n_name.append('')
        else:
            n_name.append(nlables[i])
            n_shape.append('RECTANGLE')

            n_size.append(40)
            n_font.append(20)
        nodes_list.append([n_id,n_name,n_size,n_shape,n_color,n_font])
    for i,item in enumerate(edges):
        source.append(item[0])
        target.append(item[1])
        e_width.append(2.0)
        e_paint.append(0)
        e_id.append(item[0] + ' (-) ' + item[1])
        if arrow[i]==1:
            source_arrow_shape.append('DIAMOND')
            line_type.append('DOT')
        elif arrow[i]==-1:
            source_arrow_shape.append('T')
            line_type.append('DOT')
        else:
            source_arrow_shape.append('None')
            #source_arrow_shape.append('ARROW')
            line_type.append('SOLID')
        edges_list.append([source,target,line_type,e_width,source_arrow_shape,e_paint,e_id])
    return nodes_list,edges_list




if "__main__==__name__":
    nodes_list, edges_list=highConservedNet()
    fout=open('/Users/vikashpandey/Documents/MATLAB/reDoMetSignal/ModuleAnalysis/STE_modules_withweight_nodes_edge_list_cytoscape.pickle','w')
    pickle.dump([nodes_list,edges_list],fout)
    import pdb;

    pdb.set_trace()
    cys = CyNetwork()
    node_list, edge_list = cys.createBasedOnNodeEdge(nodes_list,edges_list)
    cys.createCysNet(node_list, edge_list)