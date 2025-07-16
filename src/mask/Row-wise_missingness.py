# Row-wise missingness
#     Como funciona: determina a proporção de valores ausentes dentro de cada linha do dataset
#     Aplicação: avalia qualidade dos dados e impacto na análise estatística
#     Implementação: conta entradas ausentes por linha e calcula percentual de missingness

import polars as pl
import numpy as np
from tqdm import tqdm
import os
from collections import defaultdict


class RowwiseMissingness:
    def __init__(self, missing_rate=0.1, high_missingness_threshold=0.3):
        self.missing_rate = missing_rate
        self.high_missingness_threshold = high_missingness_threshold
    
    def simulate_sample_data(self, snp_list, n_samples=1000):
        """Simulate genotype data for samples to demonstrate row-wise missingness"""
        np.random.seed(42)
        
        # Create sample IDs
        sample_ids = [f"SAMPLE_{i:04d}" for i in range(n_samples)]
        
        # Simulate genotype data (0, 1, 2 for AA, AB, BB)
        genotype_data = {}
        for sample_id in sample_ids:
            genotypes = np.random.choice([0, 1, 2], size=len(snp_list), p=[0.25, 0.5, 0.25])
            genotype_data[sample_id] = genotypes
        
        return sample_ids, genotype_data
    
    def apply_rowwise_missing_pattern(self, sample_ids, genotype_data, snp_list):
        """Apply row-wise missing patterns to simulate realistic missing data"""
        np.random.seed(42)
        
        # Classify samples into different missingness categories
        n_samples = len(sample_ids)
        
        # Low missingness samples (0-10% missing)
        n_low = int(n_samples * 0.7)
        # Medium missingness samples (10-30% missing)  
        n_medium = int(n_samples * 0.2)
        # High missingness samples (30%+ missing)
        n_high = n_samples - n_low - n_medium
        
        sample_categories = (['low'] * n_low + 
                           ['medium'] * n_medium + 
                           ['high'] * n_high)
        np.random.shuffle(sample_categories)
        
        # Apply missing data patterns
        missing_data = {}
        sample_missingness = {}
        
        for i, sample_id in enumerate(tqdm(sample_ids, desc="Applying row-wise missingness")):
            category = sample_categories[i]
            original_genotypes = genotype_data[sample_id].copy()
            
            # Determine missing rate for this sample
            if category == 'low':
                sample_missing_rate = np.random.uniform(0.0, 0.1)
            elif category == 'medium':
                sample_missing_rate = np.random.uniform(0.1, 0.3)
            else:  # high
                sample_missing_rate = np.random.uniform(0.3, 0.6)
            
            # Apply missing data
            n_to_mask = int(len(snp_list) * sample_missing_rate)
            if n_to_mask > 0:
                missing_indices = np.random.choice(len(snp_list), size=n_to_mask, replace=False)
                original_genotypes[missing_indices] = -9  # Missing value code
            
            missing_data[sample_id] = original_genotypes
            
            # Calculate actual missingness percentage
            n_missing = np.sum(original_genotypes == -9)
            missingness_pct = (n_missing / len(snp_list)) * 100
            sample_missingness[sample_id] = {
                'n_total': len(snp_list),
                'n_missing': n_missing,
                'missingness_pct': missingness_pct,
                'category': category,
                'high_missingness': missingness_pct >= (self.high_missingness_threshold * 100)
            }
        
        return missing_data, sample_missingness
    
    def assess_missing_patterns(self, sample_missingness):
        """Assess patterns of missingness (MAR vs MNAR)"""
        missingness_values = [info['missingness_pct'] for info in sample_missingness.values()]
        
        # Basic statistics
        mean_missingness = np.mean(missingness_values)
        std_missingness = np.std(missingness_values)
        
        # Categorize samples
        low_miss_samples = sum(1 for info in sample_missingness.values() 
                              if info['missingness_pct'] < 10)
        medium_miss_samples = sum(1 for info in sample_missingness.values() 
                                 if 10 <= info['missingness_pct'] < 30)
        high_miss_samples = sum(1 for info in sample_missingness.values() 
                               if info['missingness_pct'] >= 30)
        
        # Assess if pattern suggests MAR or MNAR
        # High variance might suggest MNAR (systematic reasons for missingness)
        cv = std_missingness / mean_missingness if mean_missingness > 0 else 0
        pattern_type = "MNAR" if cv > 1.0 else "MAR"
        
        return {
            'mean_missingness': mean_missingness,
            'std_missingness': std_missingness,
            'cv': cv,
            'pattern_type': pattern_type,
            'low_miss_samples': low_miss_samples,
            'medium_miss_samples': medium_miss_samples,
            'high_miss_samples': high_miss_samples
        }
    
    def create_mask_from_ld(self, ld_file_path, output_path, batch_size=10000):
        """Create row-wise missingness mask from LD data"""
        print("Loading LD data for row-wise missingness analysis...")
        
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
        
        # Simulate sample data
        print("Simulating sample genotype data...")
        sample_ids, genotype_data = self.simulate_sample_data(snp_list, n_samples=1000)
        
        # Apply row-wise missing patterns
        print("Applying row-wise missing patterns...")
        missing_data, sample_missingness = self.apply_rowwise_missing_pattern(
            sample_ids, genotype_data, snp_list
        )
        
        # Assess missing patterns
        pattern_stats = self.assess_missing_patterns(sample_missingness)
        
        print(f"Row-wise missingness analysis complete:")
        print(f"  Mean missingness: {pattern_stats['mean_missingness']:.2f}%")
        print(f"  Pattern type: {pattern_stats['pattern_type']}")
        print(f"  Low missingness samples: {pattern_stats['low_miss_samples']}")
        print(f"  Medium missingness samples: {pattern_stats['medium_miss_samples']}")
        print(f"  High missingness samples: {pattern_stats['high_miss_samples']}")
        
        # Create mask data
        mask_data = []
        for sample_id in tqdm(sample_ids, desc="Creating row-wise mask"):
            info = sample_missingness[sample_id]
            
            mask_data.append([
                sample_id,
                info['n_total'],
                info['n_missing'],
                info['missingness_pct'],
                info['category'],
                info['high_missingness'],
                pattern_stats['pattern_type']
            ])
        
        # Save mask
        mask_df = pl.DataFrame(
            mask_data,
            schema=['SAMPLE_ID', 'N_TOTAL_SNPS', 'N_MISSING_SNPS', 
                   'MISSINGNESS_PCT', 'CATEGORY', 'HIGH_MISSINGNESS', 'PATTERN_TYPE']
        )
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        mask_df.write_csv(output_path)
        
        print(f"Row-wise missingness mask saved to: {output_path}")
        print(f"Total samples: {len(sample_ids)}")
        print(f"Samples with high missingness (>{self.high_missingness_threshold*100}%): {pattern_stats['high_miss_samples']}")
        
        return mask_df, pattern_stats


def main():
    masker = RowwiseMissingness(missing_rate=0.1, high_missingness_threshold=0.3)
    mask_df, stats = masker.create_mask_from_ld(
        '../data/chr1_ld.ld',
        '../output/rowwise_mask.csv'
    )


if __name__ == "__main__":
    main()