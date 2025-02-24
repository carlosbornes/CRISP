import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import plotly.graph_objects as go
import plotly.io as pio


pio.renderers.default = 'svg'
pio.renderers.default = 'notebook'


def aggregate_data(data, index_map, N):
    # Aggregates hydrogen bond data
    node_frequency = {node: 0 for node in index_map.keys()}
    edge_weight = {}
    all_pairs = []

    for frame, group in data.groupby('Frame'):
        pairs = group[['Donor', 'Acceptor']].values
        all_pairs.extend(pairs)
        
        for donor, acceptor in pairs:
            if donor in index_map and acceptor in index_map:
                node_frequency[donor] += 1
                node_frequency[acceptor] += 1
                edge = tuple(sorted([donor, acceptor]))
                edge_weight[edge] = edge_weight.get(edge, 0) + 1

    G = nx.Graph()

    for node, freq in node_frequency.items():
        G.add_node(node, size=freq)
    
    for (donor, acceptor), weight in edge_weight.items():
        G.add_edge(donor, acceptor, weight=weight)

    corr_matrix = np.zeros((N, N), dtype=int)
    for (donor, acceptor), weight in edge_weight.items():
        donor_idx = index_map[donor]
        acceptor_idx = index_map[acceptor]
        corr_matrix[donor_idx, acceptor_idx] = weight
        corr_matrix[acceptor_idx, donor_idx] = weight

    return corr_matrix, G, all_pairs

def process_frame(data, water_indices, frame_index, index_map, N):
    # Processes data for a specific frame
    frame_data = data[data['Frame'] == frame_index]
    pairs = frame_data[['Donor', 'Acceptor']].values
    corr_matrix = np.zeros((N, N), dtype=int)

    for donor, acceptor in pairs:
        if donor in index_map and acceptor in index_map:
            donor_idx = index_map[donor]
            acceptor_idx = index_map[acceptor]
            corr_matrix[donor_idx, acceptor_idx] += 1
            corr_matrix[acceptor_idx, donor_idx] += 1

    G = nx.Graph()
    for i, donor_idx in enumerate(range(len(water_indices))):
        for j, acceptor_idx in enumerate(range(len(water_indices))):
            if corr_matrix[i, j] > 0:
                G.add_edge(water_indices[i], water_indices[j], weight=corr_matrix[i, j])

    return corr_matrix, G, pairs

def visualize_hydrogen_bonds_matrix(corr_matrix, water_indices=None, frame_index=None, average=False):
    # Visualizes the hydrogen bond correlation matrix using Matplotlib
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr_matrix, annot=True, fmt="d", cmap='viridis', xticklabels=water_indices, yticklabels=water_indices)
    
    if average:
        plt.title("Average Hydrogen Bond Correlation Matrix Across All Frames")
    else:
        plt.title(f"Hydrogen Bond Correlation Matrix for Frame {frame_index}")
    
    plt.xlabel("Water Index")
    plt.ylabel("Water Index")
    plt.show()

def visualize_hydrogen_bonds_plotly(G, water_indices=None, frame_index=None, average=False):
    # Visualizes the hydrogen bond network using Plotly
    seed = 42
    pos = nx.spring_layout(G, seed=seed)
    
    node_size = [G.nodes[node].get('size', 20) for node in G.nodes()]
    node_color = [G.degree(node) for node in G.nodes()]
    edge_width = [G[u][v]['weight'] * 0.5 for u, v in G.edges()]

    x_nodes = [pos[node][0] for node in G.nodes()]
    y_nodes = [pos[node][1] for node in G.nodes()]

    edge_trace = []
    for (u, v), width in zip(G.edges(), edge_width):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_trace.append(go.Scatter(x=[x0, x1], y=[y0, y1], mode='lines', line=dict(width=width, color='Magenta'), name=f'Bond {u}-{v}'))
    
    node_trace = go.Scatter(x=x_nodes, y=y_nodes, mode='markers', 
                            marker=dict(size=node_size, color=node_color, colorscale='Viridis',
                                        colorbar=dict(thickness=15, title='Node Connections',
                                                      xanchor='left', titleside='right')),
                            text=[str(node) for node in G.nodes()],
                            textposition='top center', name='Water Molecules')
    
    fig = go.Figure(data=edge_trace + [node_trace],
                    layout=go.Layout(title=f'<br>{"Average Hydrogen Bond Network Across All Frames" if average else f"Hydrogen Bond Network for Frame {frame_index}"}',
                                     titlefont_size=16, showlegend=False, hovermode='closest',
                                     margin=dict(b=20, l=5, r=5, t=40),
                                     annotations=[dict(text="Python code: <a href='https://plotly.com/python/network-graphs/'> https://plotly.com/python/network-graphs/</a>",
                                                       showarrow=False, xref="paper", yref="paper", x=0.005, y=-0.002)],
                                     xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                     yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
    
    # Save the figure to an HTML file
    fig.write_html('hydrogen_bonds_network.html')
    print("Figure saved as 'hydrogen_bonds_network.html'")


def visualize_hydrogen_bonds_matplotlib(corr_matrix, water_indices, pairs, frame_index, index_map, average=False):
    # Plot hydrogen bonds matrix and network graph using Matplotlib
    G = nx.Graph()

    for donor, acceptor in pairs:
        if donor in index_map and acceptor in index_map:
            donor_idx = index_map[donor]
            acceptor_idx = index_map[acceptor]
            
            if not G.has_node(donor):
                G.add_node(donor, size=0)
            if not G.has_node(acceptor):
                G.add_node(acceptor, size=0)
            
            if G.has_edge(donor, acceptor):
                G[donor][acceptor]['weight'] += 1
            else:
                G.add_edge(donor, acceptor, weight=1)

    node_size = {node: len([1 for donor, acceptor in pairs if donor == node or acceptor == node]) * 100 for node in G.nodes()}
    edge_width = [G[u][v]['weight'] for u, v in G.edges()]

    seed = 42
    pos = nx.spring_layout(G, seed=seed)
    
    plt.figure(figsize=(12, 12))
    nx.draw_networkx_edges(G, pos, alpha=0.5, width=edge_width, edge_color="m")
    node_sizes = [node_size[n] for n in G.nodes()]
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color="#210070", alpha=0.9)
    nx.draw_networkx_labels(G, pos, font_size=10, font_color="white", verticalalignment='center')

    plt.title(f"Hydrogen Bond Network for Frame {frame_index}" if not average else "Average Hydrogen Bond Network Across All Frames")
    plt.axis("off")
    plt.show()

def visualize_hydrogen_bonds(csv_file, wat_array_path, frame_index=None, average=False):
    # Visualizes hydrogen bonds from a CSV file containing donor-acceptor pairs
    data = pd.read_csv(csv_file)
    water_indices = np.load(wat_array_path)    
    index_map = {idx: i for i, idx in enumerate(water_indices)}
    N = len(water_indices)
    
    if average:
        corr_matrix, G, pairs = aggregate_data(data, index_map, N)
        visualize_hydrogen_bonds_matrix(corr_matrix, water_indices=water_indices, average=True)
        visualize_hydrogen_bonds_plotly(G, average=True)
        visualize_hydrogen_bonds_matplotlib(corr_matrix, water_indices, pairs, frame_index=0, index_map=index_map, average=True)
    else:
        if frame_index is None:
            raise ValueError("frame_index must be provided when average=False")
        
        corr_matrix, G, pairs = process_frame(data, water_indices, frame_index, index_map, N)
        visualize_hydrogen_bonds_matrix(corr_matrix, water_indices=water_indices, frame_index=frame_index)
        visualize_hydrogen_bonds_plotly(G, frame_index=frame_index)
        visualize_hydrogen_bonds_matplotlib(corr_matrix, water_indices, pairs, frame_index, index_map)

