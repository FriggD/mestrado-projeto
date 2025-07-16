# Máscara baseada em MAF (Minor Allele Frequency)
#     Como funciona: mascara preferencialmente SNPs com baixa frequência alélica
#     Justificativa: esses SNPs são mais difíceis de imputar e mais suscetíveis a erros em sequenciamento real
#     Aplicável em: avaliação do modelo em variantes raras
# Máscara por cluster (bloco de haplótipo)
#     Como funciona: remove blocos inteiros de SNPs que pertencem a um mesmo haplótipo
#     Aplicação: simula falhas por regiões inteiras em painéis SNP arrays
#     Pode usar: segmentação por LD ou mapas de haplótipos reais

import polars as pl
import numpy as np
from tqdm import tqdm
import os
from collections import defaultdict
from sklearn.cluster import AgglomerativeClustering


class MAFClusterMasking:
    def __init__(self, maf_threshold=0.05, missing_rate=0.1, r2_threshold=0.8, cluster_size_min=3):
        self.maf_threshold = maf_threshold
        self.missing_rate = missing_rate
        self.r2_threshold = r2_threshold
        self.cluster_size_min = cluster_size_min
    
    def simulate_maf(self, snp_positions):
        """Simulate MAF values using Markov model based on genomic position"""
        np.random.seed(42)
        n_snps = len(snp_positions)
        
        # Markov transition probabilities for MAF categories
        # States: 0=rare (MAF<0.05), 1=common (MAF>=0.05)
        transition_matrix = np.array([
            [0.7, 0.3],  # rare -> [rare, common]
            [0.2, 0.8]   # common -> [rare, common]
        ])
        
        # Initialize first SNP state
        states = [np.random.choice([0, 1], p=[0.3, 0.7])]
        
        # Generate states using Markov chain
        for i in range(1, n_snps):
            current_state = states[-1]
            next_state = np.random.choice([0, 1], p=transition_matrix[current_state])
            states.append(next_state)
        
        # Convert states to MAF values
        maf_values = []
        for state in states:
            if state == 0:  # rare
                maf = np.random.beta(1, 20)  # skewed towards low MAF
            else:  # common
                maf = np.random.beta(2, 2) * 0.5  # uniform between 0-0.5
            maf_values.append(min(maf, 0.5))
        
        return np.array(maf_values)
    
    def create_ld_clusters(self, ld_data, snp_positions):
        """Create haplotype clusters based on LD structure"""
        print("Creating LD-based clusters...")
        
        # Build LD matrix
        snp_to_idx = {snp: i for i, snp in enumerate(snp_positions)}
        n_snps = len(snp_positions)
        ld_matrix = np.eye(n_snps)  # Initialize with identity
        
        for _, row in ld_data.iter_rows(named=True):
            snp_a = f"{row['CHR_A']}_{row['BP_A']}"
            snp_b = f"{row['CHR_B']}_{row['BP_B']}"
            
            if snp_a in snp_to_idx and snp_b in snp_to_idx:
                idx_a, idx_b = snp_to_idx[snp_a], snp_to_idx[snp_b]
                r2_val = float(row['R2'])
                ld_matrix[idx_a, idx_b] = r2_val
                ld_matrix[idx_b, idx_a] = r2_val
        
        # Convert to distance matrix for clustering
        distance_matrix = 1 - ld_matrix
        
        # Hierarchical clustering
        n_clusters = max(2, n_snps // 10)  # Adaptive number of clusters
        clustering = AgglomerativeClustering(
            n_clusters=n_clusters,
            metric='precomputed',
            linkage='average'
        )
        
        cluster_labels = clustering.fit_predict(distance_matrix)
        
        # Group SNPs by cluster
        clusters = defaultdict(list)
        for i, label in enumerate(cluster_labels):
            clusters[label].append(snp_positions[i])
        
        # Filter clusters by minimum size
        valid_clusters = {k: v for k, v in clusters.items() if len(v) >= self.cluster_size_min}
        
        print(f"Created {len(valid_clusters)} valid clusters")
        return valid_clusters
    
    def create_mask_from_ld(self, ld_file_path, output_path, batch_size=10000):
        """Create combined MAF and cluster-based mask"""
        print("Loading LD data...")
        
        # Load LD data
        ld_chunks = []
        unique_snps = set()
        file_size = os.path.getsize(ld_file_path)
        
        with open(ld_file_path, 'r') as f:
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
                    
                    batch_data = []
                    for line in lines:
                        parts = line.strip().split()
                        if len(parts) >= 7 and parts[6] != 'R2':
                            try:
                                r2_val = float(parts[6])
                                if r2_val >= self.r2_threshold:
                                    snp_a = f"{parts[0]}_{parts[1]}"
                                    snp_b = f"{parts[3]}_{parts[4]}"
                                    batch_data.append(parts[:7])
                                    unique_snps.update([snp_a, snp_b])
                            except ValueError:
                                continue
                    
                    if batch_data:
                        batch_df = pl.DataFrame(
                            batch_data,
                            schema=['CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B', 'SNP_B', 'R2']
                        )
                        ld_chunks.append(batch_df)
        
        ld_data = pl.concat(ld_chunks) if ld_chunks else pl.DataFrame()
        snp_positions = sorted(list(unique_snps))
        
        print(f"Found {len(snp_positions)} unique SNPs")
        
        # Simulate MAF values using Markov model
        print("Simulating MAF values with Markov model...")
        maf_values = self.simulate_maf(snp_positions)
        
        # Create LD-based clusters
        clusters = self.create_ld_clusters(ld_data, snp_positions)
        
        # Select SNPs to mask based on MAF and clusters
        print("Selecting SNPs to mask...")
        masked_snps = set()
        
        # 1. MAF-based masking (prioritize rare variants)
        rare_snps = [snp for i, snp in enumerate(snp_positions) 
                    if maf_values[i] < self.maf_threshold]
        
        n_maf_mask = int(len(rare_snps) * self.missing_rate * 2)  # Higher rate for rare SNPs
        if rare_snps:
            np.random.seed(42)
            maf_masked = np.random.choice(rare_snps, size=min(n_maf_mask, len(rare_snps)), replace=False)
            masked_snps.update(maf_masked)
        
        # 2. Cluster-based masking (remove entire haplotype blocks)
        cluster_list = list(clusters.values())
        n_clusters_to_mask = max(1, int(len(cluster_list) * self.missing_rate))
        
        np.random.seed(43)
        selected_clusters = np.random.choice(len(cluster_list), size=n_clusters_to_mask, replace=False)
        
        for cluster_idx in selected_clusters:
            cluster_snps = cluster_list[cluster_idx]
            masked_snps.update(cluster_snps)
        
        print(f"Masking {len(masked_snps)} SNPs total")
        
        # Create detailed mask data
        mask_data = []
        snp_to_cluster = {}
        for cluster_id, snps in clusters.items():
            for snp in snps:
                snp_to_cluster[snp] = cluster_id
        
        for i, snp in enumerate(tqdm(snp_positions, desc="Creating mask")):
            chr_pos = snp.split('_')
            is_masked = snp in masked_snps
            maf_val = maf_values[i]
            is_rare = maf_val < self.maf_threshold
            cluster_id = snp_to_cluster.get(snp, -1)
            
            mask_data.append([
                chr_pos[0], chr_pos[1], snp, is_masked, maf_val, is_rare, cluster_id
            ])
        
        # Save mask
        mask_df = pl.DataFrame(
            mask_data,
            schema=['CHR', 'POS', 'SNP_ID', 'MASKED', 'MAF', 'IS_RARE', 'CLUSTER_ID']
        )
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        mask_df.write_csv(output_path)
        
        print(f"Mask saved to: {output_path}")
        print(f"Total SNPs: {len(snp_positions)}")
        print(f"Masked SNPs: {len(masked_snps)}")
        print(f"Rare SNPs (MAF<{self.maf_threshold}): {sum(maf_values < self.maf_threshold)}")
        print(f"Clusters created: {len(clusters)}")
        
        return mask_df


def main():
    masker = MAFClusterMasking(
        maf_threshold=0.05,
        missing_rate=0.1,
        r2_threshold=0.8,
        cluster_size_min=3
    )
    
    mask_df = masker.create_mask_from_ld(
        '../data/chr1_ld.ld',
        '../output/maf_cluster_mask.csv'
    )


if __name__ == "__main__":
    main()