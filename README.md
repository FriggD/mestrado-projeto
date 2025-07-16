Saída do PLINK: Essa saída mostra pares de variantes localizadas no mesmo cromossomo e o grau de correlação entre elas, medido principalmente pelo R².

| Coluna  | Descrição                                                  |
| ------- | ---------------------------------------------------------- |
| `CHR_A` | Cromossomo da primeira variante                            |
| `BP_A`  | Posição da primeira variante (em pares de base)            |
| `SNP_A` | Identificador da primeira variante (aqui aparece como `.`) |
| `CHR_B` | Cromossomo da segunda variante                             |
| `BP_B`  | Posição da segunda variante                                |
| `SNP_B` | Identificador da segunda variante (aqui também `.`)        |
| `R2`    | Valor de **R²**, medida de correlação entre os dois SNPs   |

Exemplo:
1        10505    .      1        10506    .            1 

A variante A está no cromossomo 1, posição 10505
A variante B está no cromossomo 1, posição 10506
O valor de R² = 1 → desequilíbrio de ligação perfeito (totalmente correlacionadas)

Sobre o R²:
    É uma medida que varia de 0 a 1, onde:
        0 = sem correlação entre as variantes (aleatórias)
        1 = variantes totalmente correlacionadas (herdadas juntas)
    Um valor acima de 0.8 ou 0.9 costuma indicar forte LD em estudos genômicos.

Variações no LD:
    Os valores intermediários (ex: 0.333067, 0.4998) indicam que há correlação parcial — essas variantes podem estar ligadas, mas com recombinação ocasional.

Por que calcular LD?
    Usado para:
        Identificar blocos de haplótipos
        Reduzir variáveis redundantes em GWAS
        Escolher SNPs representativos (tag SNPs)
        Construção de grafos em imputação genômica