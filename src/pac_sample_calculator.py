"""
PAC Sample Size Calculator
Calcula quantidade ideal de amostras de treinamento baseado na teoria PAC
"""

import numpy as np
import math


class PACCalculator:
    def __init__(self):
        pass
    
    def finite_hypothesis_class(self, H_size, epsilon=0.1, delta=0.05):
        """
        Calcula amostras para classe de hipóteses finita
        m >= (1/ε) * ln(|H|/δ)
        """
        m = (1/epsilon) * math.log(H_size/delta)
        return math.ceil(m)
    
    def vc_dimension_bound(self, vc_dim, epsilon=0.1, delta=0.05):
        """
        Calcula amostras baseado na dimensão VC
        m >= O(d/ε * log(d/ε) + 1/ε * log(1/δ))
        """
        term1 = (vc_dim/epsilon) * math.log(vc_dim/epsilon)
        term2 = (1/epsilon) * math.log(1/delta)
        m = 8 * (term1 + term2)  # Constante empírica
        return math.ceil(m)
    
    def agnostic_learning_bound(self, vc_dim, epsilon=0.1, delta=0.05):
        """
        Bound para aprendizado agnóstico
        m >= O(d/ε² * log(1/ε) + 1/ε² * log(1/δ))
        """
        term1 = (vc_dim/(epsilon**2)) * math.log(1/epsilon)
        term2 = (1/(epsilon**2)) * math.log(1/delta)
        m = 32 * (term1 + term2)  # Constante empírica
        return math.ceil(m)
    
    def rademacher_complexity_bound(self, complexity, epsilon=0.1, delta=0.05):
        """
        Bound baseado na complexidade de Rademacher
        m >= O(complexity²/ε² + log(1/δ)/ε²)
        """
        term1 = (complexity**2) / (epsilon**2)
        term2 = math.log(1/delta) / (epsilon**2)
        m = 2 * (term1 + term2)
        return math.ceil(m)
    
    def calculate_optimal_training_size(self, total_data, 
                                      hypothesis_type="linear",
                                      n_features=None,
                                      n_classes=2,
                                      epsilon=0.1, 
                                      delta=0.05,
                                      validation_split=0.2):
        """
        Calcula tamanho ótimo de treinamento dado total de dados
        """
        # Reservar dados para validação
        available_for_training = int(total_data * (1 - validation_split))
        
        # Calcular bounds baseado no tipo de hipótese
        bounds = {}
        
        if hypothesis_type == "linear":
            if n_features:
                # Dimensão VC para classificadores lineares = n_features + 1
                vc_dim = n_features + 1
                bounds['vc_bound'] = self.vc_dimension_bound(vc_dim, epsilon, delta)
                bounds['agnostic_bound'] = self.agnostic_learning_bound(vc_dim, epsilon, delta)
        
        elif hypothesis_type == "polynomial":
            if n_features:
                # Dimensão VC aproximada para polinômios de grau d
                degree = 2  # Assumindo grau 2
                vc_dim = math.comb(n_features + degree, degree)
                bounds['vc_bound'] = self.vc_dimension_bound(vc_dim, epsilon, delta)
        
        elif hypothesis_type == "decision_tree":
            if n_features:
                # Estimativa conservadora para árvores de decisão
                vc_dim = n_features * math.log2(n_features)
                bounds['vc_bound'] = self.vc_dimension_bound(int(vc_dim), epsilon, delta)
        
        elif hypothesis_type == "neural_network":
            if n_features:
                # Estimativa para redes neurais (muito conservadora)
                n_params = n_features * 10 + 10 * n_classes  # Rede simples
                bounds['param_bound'] = n_params * math.log(n_params) / epsilon
        
        # Bound genérico baseado no número de classes
        if n_classes:
            bounds['finite_bound'] = self.finite_hypothesis_class(
                2**n_features if n_features else 1000, epsilon, delta
            )
        
        # Selecionar bound mais conservador
        if bounds:
            recommended_size = min(max(bounds.values()), available_for_training)
        else:
            # Heurística padrão: 80% dos dados disponíveis
            recommended_size = int(available_for_training * 0.8)
        
        # Garantir mínimo razoável
        min_samples = max(100, 10 * (n_features or 10))
        recommended_size = max(recommended_size, min_samples)
        
        return {
            'total_data': total_data,
            'available_for_training': available_for_training,
            'recommended_training_size': min(recommended_size, available_for_training),
            'validation_size': total_data - available_for_training,
            'bounds_calculated': bounds,
            'parameters': {
                'epsilon': epsilon,
                'delta': delta,
                'hypothesis_type': hypothesis_type,
                'n_features': n_features,
                'n_classes': n_classes
            }
        }
    
    def learning_curve_analysis(self, total_data, n_features=None):
        """
        Análise de curva de aprendizado para diferentes tamanhos
        """
        sizes = np.logspace(2, np.log10(total_data), 10, dtype=int)
        sizes = np.unique(sizes)
        
        results = []
        for size in sizes:
            result = self.calculate_optimal_training_size(
                size, n_features=n_features
            )
            results.append({
                'total_size': size,
                'recommended': result['recommended_training_size'],
                'efficiency': result['recommended_training_size'] / size
            })
        
        return results


def main():
    calculator = PACCalculator()
    
    # Exemplo de uso
    total_samples = 111442046
    n_features = 250
    
    result = calculator.calculate_optimal_training_size(
        total_data=total_samples,
        hypothesis_type="linear",
        n_features=n_features,
        epsilon=0.1,
        delta=0.05
    )
    
    print("PAC Sample Size Calculator Results:")
    print(f"Total data available: {result['total_data']}")
    print(f"Recommended training size: {result['recommended_training_size']}")
    print(f"Validation size: {result['validation_size']}")
    print(f"Training efficiency: {result['recommended_training_size']/result['total_data']*100:.1f}%")
    print(f"Bounds calculated: {result['bounds_calculated']}")


if __name__ == "__main__":
    main()