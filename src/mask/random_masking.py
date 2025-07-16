# Random masking
#     Como funciona: aplica máscaras aleatórias aos dados para obscurecer informações sensíveis
#     Aplicação: reduz correlação entre dados sensíveis e informações observáveis
#     Implementação: gera números aleatórios e os aplica sistematicamente aos dados

import polars as pl
import numpy as np
from tqdm import tqdm
import os


class RandomMasking:
    def __init__(self, missing_rate=0.1, seed=42):
        self.missing_rate = missing_rate
        self.seed = seed
    
    def generate_random_masks(self, n_snps):
        """Generate random masks using deterministic RNG"""
        np.random.seed(self.seed)
        
        # Generate multiple sets of random masks for enhanced security
        mask_sets = []
        for i in range(3):  # Multiple stages of masking
            mask_set = np.random.random(n_snps)
            mask_sets.append(mask_set)
        
        return mask_sets
    
    def apply_structured_masking(self, snp_list, mask_sets):
        """Apply masks in structured manner with multiple stages"""
        n_snps = len(snp_list)
        n_to_mask = int(n_snps * self.missing_rate)
        
        # Stage 1: Primary random selection
        primary_mask = mask_sets[0]
        primary_indices = np.argsort(primary_mask)[:n_to_mask]
        
        # Stage 2: Secondary randomization
        secondary_mask = mask_sets[1]
        secondary_selection = secondary_mask[primary_indices]
        
        # Stage 3: Final mask application with different order
        final_mask = mask_sets[2]
        
        # Combine masks using structured approach
        masked_indices = set()
        
        # Apply primary selection
        for idx in primary_indices:
            if secondary_selection[np.where(primary_indices == idx)[0][0]] > 0.3:
                masked_indices.add(idx)
        
        # Apply additional random selection if needed
        remaining_needed = n_to_mask - len(masked_indices)
        if remaining_needed > 0:
            remaining_candidates = [i for i in range(n_snps) if i not in masked_indices]
            if remaining_candidates:
                additional_indices = np.random.choice(
                    remaining_candidates, 
                    size=min(remaining_needed, len(remaining_candidates)), 
                    replace=False
                )
                masked_indices.update(additional_indices)
        
        return masked_indices
    
    def quantify_masking_strength(self, original_snps, masked_indices):
        """Quantify masking strength (QMS) - resistance against correlation attacks"""
        n_total = len(original_snps)
        n_masked = len(masked_indices)
        
        # Calculate distribution uniformity
        positions = [int(snp.split('_')[1]) for snp in original_snps]
        masked_positions = [positions[i] for i in masked_indices]
        
        # QMS based on spatial distribution and randomness
        if len(masked_positions) < 2:
            return 0.0
        
        # Calculate spacing uniformity
        sorted_masked = sorted(masked_positions)
        spacings = [sorted_masked[i+1] - sorted_masked[i] for i in range(len(sorted_masked)-1)]
        
        if len(spacings) == 0:
            uniformity = 1.0
        else:
            mean_spacing = np.mean(spacings)
            spacing_variance = np.var(spacings) if mean_spacing > 0 else 0
            uniformity = 1.0 / (1.0 + spacing_variance / (mean_spacing**2 + 1e-10))
        
        # QMS combines masking rate and spatial uniformity
        qms = (n_masked / n_total) * uniformity
        
        return qms
    
    def create_mask_from_ld(self, ld_file_path, output_path, batch_size=10000):
        """Create random mask from LD data"""
        print("Loading LD data for random masking...")
        
        # Extract unique SNPs
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
        
        # Generate random masks
        print("Generating random masks...")
        mask_sets = self.generate_random_masks(len(snp_list))
        
        # Apply structured masking
        print("Applying structured masking...")
        masked_indices = self.apply_structured_masking(snp_list, mask_sets)
        
        # Quantify masking strength
        qms = self.quantify_masking_strength(snp_list, masked_indices)
        
        print(f"Masking {len(masked_indices)} SNPs ({len(masked_indices)/len(snp_list)*100:.1f}%)")
        print(f"Quantitative Masking Strength (QMS): {qms:.4f}")
        
        # Create mask data
        mask_data = []
        for i, snp in enumerate(tqdm(snp_list, desc="Creating mask")):
            chr_pos = snp.split('_')
            is_masked = i in masked_indices
            
            # Include mask stage information
            stage_1_val = mask_sets[0][i] if i < len(mask_sets[0]) else 0
            stage_2_val = mask_sets[1][i] if i < len(mask_sets[1]) else 0
            stage_3_val = mask_sets[2][i] if i < len(mask_sets[2]) else 0
            
            mask_data.append([
                chr_pos[0], chr_pos[1], snp, is_masked, 
                stage_1_val, stage_2_val, stage_3_val, qms
            ])
        
        # Save mask
        mask_df = pl.DataFrame(
            mask_data,
            schema=['CHR', 'POS', 'SNP_ID', 'MASKED', 
                   'STAGE1_VAL', 'STAGE2_VAL', 'STAGE3_VAL', 'QMS']
        )
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        mask_df.write_csv(output_path)
        
        print(f"Random mask saved to: {output_path}")
        print(f"Total SNPs: {len(snp_list)}")
        print(f"Masked SNPs: {len(masked_indices)}")
        print(f"Masking rate: {len(masked_indices)/len(snp_list)*100:.1f}%")
        print(f"QMS Score: {qms:.4f}")
        
        return mask_df


def main():
    masker = RandomMasking(missing_rate=0.1, seed=42)
    mask_df = masker.create_mask_from_ld(
        '../data/chr1_ld.ld',
        '../output/random_mask.csv'
    )


if __name__ == "__main__":
    main()