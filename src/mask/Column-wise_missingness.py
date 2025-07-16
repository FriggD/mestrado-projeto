# Máscara por SNP (Column-wise Missingness)
#     Como funciona: remove uma fração de genótipos de um SNP específico
#     Útil para: avaliar a capacidade de imputar SNPs inteiros com base em vizinhos (LD)
#     Aplicável em: seleção de tag SNPs, avaliação de robustez

import polars as pl
import numpy as np
from tqdm import tqdm
import os


class ColumnwiseMissingness:
    def __init__(self, missing_rate=0.1):
        self.missing_rate = missing_rate
    
    def create_mask_from_ld(self, ld_file_path, output_path, batch_size=10000):
        """Create column-wise missingness mask from LD data"""
        print("Loading LD data...")
        
        # Get unique SNP positions
        unique_snps = set()
        file_size = os.path.getsize(ld_file_path)
        
        with open(ld_file_path, 'r') as f:
            with tqdm(total=file_size, desc="Collecting SNPs", unit='B', unit_scale=True) as pbar:
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
                            snp_a = f"{parts[0]}_{parts[1]}"
                            snp_b = f"{parts[3]}_{parts[4]}"
                            unique_snps.update([snp_a, snp_b])
        
        snp_list = sorted(list(unique_snps))
        print(f"Found {len(snp_list)} unique SNPs")
        
        # Select SNPs to mask
        n_mask = int(len(snp_list) * self.missing_rate)
        np.random.seed(42)
        masked_snps = np.random.choice(snp_list, size=n_mask, replace=False)
        
        print(f"Masking {len(masked_snps)} SNPs ({self.missing_rate*100:.1f}%)")
        
        # Create mask data
        mask_data = []
        for snp in tqdm(snp_list, desc="Creating mask"):
            chr_pos = snp.split('_')
            is_masked = snp in masked_snps
            mask_data.append([chr_pos[0], chr_pos[1], snp, is_masked])
        
        # Save mask
        mask_df = pl.DataFrame(
            mask_data,
            schema=['CHR', 'POS', 'SNP_ID', 'MASKED']
        )
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        mask_df.write_csv(output_path)
        
        print(f"Mask saved to: {output_path}")
        print(f"Total SNPs: {len(snp_list)}")
        print(f"Masked SNPs: {len(masked_snps)}")
        
        return mask_df


def main():
    masker = ColumnwiseMissingness(missing_rate=0.1)
    mask_df = masker.create_mask_from_ld(
        '../data/chr1_ld.ld',
        '../output/columnwise_mask.csv'
    )


if __name__ == "__main__":
    main()