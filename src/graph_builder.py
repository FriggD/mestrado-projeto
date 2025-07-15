import polars as pl
import numpy as np
import torch
from torch_geometric.data import Data
from pyvis.network import Network
from sklearn.preprocessing import LabelEncoder
from tqdm import tqdm
import os


class LDGraphBuilder:
    def __init__(self, r2_threshold=0.8):
        self.r2_threshold = r2_threshold
        self.node_encoder = LabelEncoder()
        
    def load_ld_data(self, file_path, batch_size=10000):
        """Load LD data from PLINK output file in batches"""
        file_size = os.path.getsize(file_path)
        chunks = []
        
        with open(file_path, 'r') as f:
            with tqdm(total=file_size, desc="Loading LD data", unit='B', unit_scale=True) as pbar:
                while True:
                    lines = []
                    for _ in range(batch_size):
                        line = f.readline()
                        if not line:
                            break
                        lines.append(line)
                        pbar.update(len(line.encode('utf-8')))
                    
                    if not lines:
                        break
                    
                    # Process batch
                    batch_data = []
                    for line in lines:
                        parts = line.strip().split()
                        if len(parts) >= 7 and parts[6] != 'R2':  # Skip header
                            batch_data.append(parts[:7])
                    
                    if batch_data:
                        batch_df = pl.DataFrame(
                            batch_data,
                            schema=['CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B', 'SNP_B', 'R2']
                        )
                        chunks.append(batch_df)
        
        return pl.concat(chunks) if chunks else pl.DataFrame()
    
    def create_nodes_from_ids(self, node_ids):
        """Create node DataFrame from node IDs"""
        nodes_data = []
        for node_id in node_ids:
            chr_str, bp_str = node_id.split('_')
            nodes_data.append([chr_str, bp_str, node_id])
        
        return pl.DataFrame(
            nodes_data,
            schema=['chr', 'bp', 'node_id']
        )
    
    def build_pytorch_geometric_graph(self, file_path, batch_size=10000):
        """Build PyTorch Geometric graph from LD data"""
        print("Loading LD data...")
        df = self.load_ld_data(file_path, batch_size)
        
        print("Filtering by R² threshold...")
        filtered_chunks = []
        for batch_df in tqdm(df.iter_slices(n_rows=batch_size), desc="Filtering batches"):
            filtered_batch = batch_df.filter(
                pl.col('R2').str.contains(r'^\d+\.?\d*$').fill_null(False) &
                (pl.col('R2').cast(pl.Float64) >= self.r2_threshold)
            )
            if len(filtered_batch) > 0:
                filtered_chunks.append(filtered_batch)
        
        print("Concatenating filtered data...")
        df = pl.concat(filtered_chunks) if filtered_chunks else pl.DataFrame()
        
        print("Creating unique nodes...")
        # Create nodes in batches to avoid memory issues
        unique_node_ids = set()
        for batch_df in tqdm(df.iter_slices(n_rows=batch_size), desc="Collecting nodes"):
            batch_nodes_a = batch_df.select([
                pl.concat_str([pl.col('CHR_A'), pl.col('BP_A')], separator='_').alias('node_id')
            ])['node_id'].to_list()
            
            batch_nodes_b = batch_df.select([
                pl.concat_str([pl.col('CHR_B'), pl.col('BP_B')], separator='_').alias('node_id')
            ])['node_id'].to_list()
            
            unique_node_ids.update(batch_nodes_a + batch_nodes_b)
        
        print(f"Found {len(unique_node_ids)} unique nodes")
        node_ids = list(unique_node_ids)
        self.node_encoder.fit(node_ids)
        
        # Create edges in batches
        print("Creating edges...")
        edge_indices = []
        edge_attrs = []
        
        for batch_df in tqdm(df.iter_slices(n_rows=batch_size), desc="Processing edges"):
            source_ids = batch_df.with_columns(
                pl.concat_str([pl.col('CHR_A'), pl.col('BP_A')], separator='_').alias('source')
            )['source'].to_list()
            
            target_ids = batch_df.with_columns(
                pl.concat_str([pl.col('CHR_B'), pl.col('BP_B')], separator='_').alias('target')
            )['target'].to_list()
            
            batch_edge_index = torch.tensor([
                self.node_encoder.transform(source_ids),
                self.node_encoder.transform(target_ids)
            ], dtype=torch.long)
            
            batch_edge_attr = torch.tensor(batch_df['R2'].cast(pl.Float64).to_numpy(), dtype=torch.float)
            
            edge_indices.append(batch_edge_index)
            edge_attrs.append(batch_edge_attr)
        
        edge_index = torch.cat(edge_indices, dim=1) if edge_indices else torch.empty((2, 0), dtype=torch.long)
        edge_attr = torch.cat(edge_attrs) if edge_attrs else torch.empty(0, dtype=torch.float)
        
        # Node features (chromosome and position)
        print("Creating node features...")
        nodes = self.create_nodes_from_ids(node_ids)
        x = torch.tensor(nodes.select(['chr', 'bp']).cast(pl.Float64).to_numpy(), dtype=torch.float)
        
        return Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
    
    def visualize_graph(self, file_path, output_file='ld_graph.html', batch_size=10000):
        """Create interactive visualization using Pyvis"""
        print("Loading data for visualization...")
        df = self.load_ld_data(file_path, batch_size)
        df = df.filter(
            pl.col('R2').str.contains(r'^\d+\.?\d*$').fill_null(False) &
            (pl.col('R2').cast(pl.Float64) >= self.r2_threshold)
        )
        
        net = Network(height='600px', width='100%', bgcolor='#222222', font_color='white')
        
        # Add nodes in batches
        unique_node_ids = set()
        for batch_df in df.iter_slices(n_rows=batch_size):
            batch_nodes_a = batch_df.select([
                pl.concat_str([pl.col('CHR_A'), pl.col('BP_A')], separator='_').alias('node_id')
            ])['node_id'].to_list()
            batch_nodes_b = batch_df.select([
                pl.concat_str([pl.col('CHR_B'), pl.col('BP_B')], separator='_').alias('node_id')
            ])['node_id'].to_list()
            unique_node_ids.update(batch_nodes_a + batch_nodes_b)
        
        nodes = self.create_nodes_from_ids(list(unique_node_ids))
        print("Adding nodes to visualization...")
        for batch_nodes in tqdm(nodes.iter_slices(n_rows=batch_size), desc="Adding nodes"):
            for row in batch_nodes.iter_rows(named=True):
                net.add_node(
                    row['node_id'],
                    label=f"Chr{row['chr']}:{row['bp']}",
                    color='#97c2fc'
                )
        
        # Add edges in batches
        print("Adding edges to visualization...")
        for batch_df in tqdm(df.iter_slices(n_rows=batch_size), desc="Adding edges"):
            for row in batch_df.iter_rows(named=True):
                source = f"{row['CHR_A']}_{row['BP_A']}"
                target = f"{row['CHR_B']}_{row['BP_B']}"
                net.add_edge(source, target, weight=float(row['R2']), title=f"R²={float(row['R2']):.3f}")
        
        net.save_graph(output_file)
        return output_file


def main():
    builder = LDGraphBuilder(r2_threshold=0.8)
    
    # Example usage
    graph = builder.build_pytorch_geometric_graph('../data/chr1_ld.ld')
    builder.visualize_graph('../data/chr1_ld.ld', '../output/ld_network.html')


if __name__ == "__main__":
    main()