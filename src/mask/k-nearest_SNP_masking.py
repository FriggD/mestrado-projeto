# Máscara baseada em distância (k-nearest SNP masking)
#     Como funciona: mascara SNPs com base na proximidade no genoma ou no grafo de LD
#     Por quê: SNPs com alta correlação são mais prováveis de serem ausentes juntos (ex: falhas locais no sequenciamento)
#     Implementação: seleciona um SNP e remove também seus k-vizinhos

import polars as pl
import numpy as np
from tqdm import tqdm
import os
from collections import defaultdict


class KNearestSNPMasking:
    def __init__(self, k=5, missing_rate=0.1, r2_threshold=0.8):
        self.k = k
        self.missing_rate = missing_rate
        self.r2_threshold = r2_threshold
    
    def create_mask_from_ld(self, ld_file_path, output_path, batch_size=10000):
        """Create k-nearest SNP mask from LD data"""
        print("Loading LD data and building neighbor graph...")
        
        # Build SNP neighbor graph based on LD
        snp_neighbors = defaultdict(list)
        unique_snps = set()
        file_size = os.path.getsize(ld_file_path)
        
        with open(ld_file_path, 'r') as f:
            with tqdm(total=file_size, desc="Building LD graph", unit='B', unit_scale=True) as pbar:
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
                    
                    for line in lines:
                        parts = line.strip().split()
                        if len(parts) >= 7 and parts[6] != 'R2':
                            try:
                                r2_value = float(parts[6])
                                if r2_value >= self.r2_threshold:
                                    snp_a = f"{parts[0]}_{parts[1]}"
                                    snp_b = f"{parts[3]}_{parts[4]}"
                                    
                                    snp_neighbors[snp_a].append((snp_b, r2_value))
                                    snp_neighbors[snp_b].append((snp_a, r2_value))
                                    unique_snps.update([snp_a, snp_b])
                            except ValueError:
                                continue
        
        # Sort neighbors by R² value (highest first)
        for snp in snp_neighbors:
            snp_neighbors[snp].sort(key=lambda x: x[1], reverse=True)
        
        snp_list = sorted(list(unique_snps))
        print(f"Found {len(snp_list)} unique SNPs")
        print(f"Built neighbor graph with {len(snp_neighbors)} connected SNPs")
        
        # Select seed SNPs to mask
        n_seeds = int(len(snp_list) * self.missing_rate)
        np.random.seed(42)
        seed_snps = np.random.choice(snp_list, size=n_seeds, replace=False)
        
        # Find k-nearest neighbors for each seed SNP
        masked_snps = set()
        for seed_snp in tqdm(seed_snps, desc="Finding k-nearest neighbors"):
            masked_snps.add(seed_snp)
            
            # Get k-nearest neighbors based on LD
            neighbors = snp_neighbors.get(seed_snp, [])
            k_nearest = neighbors[:self.k]
            
            for neighbor_snp, r2_val in k_nearest:
                masked_snps.add(neighbor_snp)
        
        print(f"Masking {len(masked_snps)} SNPs total (seeds: {len(seed_snps)}, k={self.k})")
        
        # Create mask data
        mask_data = []
        for snp in tqdm(snp_list, desc="Creating mask"):
            chr_pos = snp.split('_')
            is_masked = snp in masked_snps
            is_seed = snp in seed_snps
            
            # Count neighbors for this SNP
            n_neighbors = len(snp_neighbors.get(snp, []))
            
            mask_data.append([
                chr_pos[0], chr_pos[1], snp, is_masked, is_seed, n_neighbors
            ])
        
        # Save mask
        mask_df = pl.DataFrame(
            mask_data,
            schema=['CHR', 'POS', 'SNP_ID', 'MASKED', 'IS_SEED', 'N_NEIGHBORS']
        )
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        mask_df.write_csv(output_path)
        
        print(f"Mask saved to: {output_path}")
        print(f"Total SNPs: {len(snp_list)}")
        print(f"Seed SNPs: {len(seed_snps)}")
        print(f"Masked SNPs: {len(masked_snps)}")
        print(f"Masking rate: {len(masked_snps)/len(snp_list)*100:.1f}%")
        
        return mask_df


def main():
    masker = KNearestSNPMasking(k=5, missing_rate=0.1, r2_threshold=0.8)
    mask_df = masker.create_mask_from_ld(
        '../data/chr1_ld.ld',
        '../output/k_nearest_mask.csv'
    )


if __name__ == "__main__":
    main()