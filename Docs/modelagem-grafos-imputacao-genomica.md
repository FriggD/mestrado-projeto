# Modelagem de Grafos para Imputação de Dados Genômicos

Este documento apresenta diferentes abordagens para modelar grafos genômicos no contexto do sistema de imputação de dados genômicos usando GNNs e VAEs.

## 1. Grafo baseado em Desequilíbrio de Ligação (LD)

```mermaid
graph TD
    subgraph "Cromossomo 1"
        SNP1((rs123))
        SNP2((rs456))
        SNP3((rs789))
        SNP4((rs012))
        SNP5((rs345))
        
        SNP1 -- "r²=0.85" --- SNP2
        SNP2 -- "r²=0.72" --- SNP3
        SNP1 -- "r²=0.65" --- SNP3
        SNP3 -- "r²=0.91" --- SNP4
        SNP4 -- "r²=0.45" --- SNP5
        
        style SNP3 fill:#f66,stroke:#333,stroke-width:2px
        style SNP5 fill:#f66,stroke:#333,stroke-width:2px
    end
    
    classDef known fill:#6b6,stroke:#333,stroke-width:1px
    classDef unknown fill:#f66,stroke:#333,stroke-width:2px
    
    class SNP1,SNP2,SNP4 known
    class SNP3,SNP5 unknown
```

*Neste modelo, os nós representam SNPs (polimorfismos de nucleotídeo único) e as arestas representam o desequilíbrio de ligação (r²) entre eles. Os nós vermelhos são SNPs desconhecidos que precisam ser imputados, enquanto os nós verdes são SNPs conhecidos.*

## 2. Grafo Hierárquico com Anotações Funcionais

```mermaid
graph TD
    subgraph "Região Genômica"
        SNP1((rs123))
        SNP2((rs456))
        SNP3((rs789))
        SNP4((rs012))
        
        Gene1[Gene BRCA1]
        Enhancer1[Enhancer E1]
        Promoter1[Promotor P1]
        
        SNP1 -- "Localizado em" --> Gene1
        SNP2 -- "Localizado em" --> Enhancer1
        SNP3 -- "Localizado em" --> Promoter1
        SNP4 -- "Localizado em" --> Gene1
        
        SNP1 -- "r²=0.85" --- SNP2
        SNP2 -- "r²=0.72" --- SNP3
        SNP1 -- "r²=0.65" --- SNP3
        SNP3 -- "r²=0.91" --- SNP4
        
        style SNP3 fill:#f66,stroke:#333,stroke-width:2px
    end
    
    classDef known fill:#6b6,stroke:#333,stroke-width:1px
    classDef unknown fill:#f66,stroke:#333,stroke-width:2px
    classDef feature fill:#69f,stroke:#333,stroke-width:1px
    
    class SNP1,SNP2,SNP4 known
    class SNP3 unknown
    class Gene1,Enhancer1,Promoter1 feature
```

*Este modelo incorpora informações funcionais, onde os SNPs estão conectados tanto por LD quanto por elementos genômicos funcionais (genes, enhancers, promotores). Esta estrutura permite que o GNN aprenda relações funcionais além das estatísticas.*

## 3. Grafo Multi-populacional

```mermaid
graph TD
    subgraph "População Europeia"
        E_SNP1((rs123))
        E_SNP2((rs456))
        E_SNP3((rs789))
        
        E_SNP1 -- "r²=0.85" --- E_SNP2
        E_SNP2 -- "r²=0.72" --- E_SNP3
        E_SNP1 -- "r²=0.65" --- E_SNP3
        
        style E_SNP3 fill:#f66,stroke:#333,stroke-width:2px
    end
    
    subgraph "População Africana"
        A_SNP1((rs123))
        A_SNP2((rs456))
        A_SNP3((rs789))
        
        A_SNP1 -- "r²=0.45" --- A_SNP2
        A_SNP2 -- "r²=0.38" --- A_SNP3
        A_SNP1 -- "r²=0.22" --- A_SNP3
        
        style A_SNP3 fill:#f66,stroke:#333,stroke-width:2px
    end
    
    E_SNP1 -- "Mesmo SNP" --- A_SNP1
    E_SNP2 -- "Mesmo SNP" --- A_SNP2
    E_SNP3 -- "Mesmo SNP" --- A_SNP3
    
    classDef known fill:#6b6,stroke:#333,stroke-width:1px
    classDef unknown fill:#f66,stroke:#333,stroke-width:2px
    
    class E_SNP1,E_SNP2,A_SNP1,A_SNP2 known
    class E_SNP3,A_SNP3 unknown
```

*Este modelo representa a estrutura de LD em diferentes populações, permitindo que o modelo aprenda padrões específicos de cada população e transfira conhecimento entre elas para melhorar a imputação.*

## 4. Grafo de Haplótipos

```mermaid
graph LR
    subgraph "Haplótipos"
        H1[Haplótipo 1: ACTG]
        H2[Haplótipo 2: ACTA]
        H3[Haplótipo 3: GCTG]
        H4[Haplótipo 4: GCTA]
        
        H1 -- "Distância=1" --- H2
        H1 -- "Distância=2" --- H3
        H2 -- "Distância=2" --- H3
        H3 -- "Distância=1" --- H4
        H2 -- "Distância=1" --- H4
        H1 -- "Distância=3" --- H4
    end
    
    subgraph "SNPs"
        SNP1((Pos 1: A/G))
        SNP2((Pos 2: C))
        SNP3((Pos 3: T))
        SNP4((Pos 4: A/G))
        
        style SNP1 fill:#6b6,stroke:#333,stroke-width:1px
        style SNP2 fill:#6b6,stroke:#333,stroke-width:1px
        style SNP3 fill:#6b6,stroke:#333,stroke-width:1px
        style SNP4 fill:#f66,stroke:#333,stroke-width:2px
    end
    
    H1 -- "Contém" --> SNP1
    H1 -- "Contém" --> SNP2
    H1 -- "Contém" --> SNP3
    H1 -- "Contém" --> SNP4
    
    H2 -- "Contém" --> SNP1
    H2 -- "Contém" --> SNP2
    H2 -- "Contém" --> SNP3
    H2 -- "Contém" --> SNP4
    
    H3 -- "Contém" --> SNP1
    H3 -- "Contém" --> SNP2
    H3 -- "Contém" --> SNP3
    H3 -- "Contém" --> SNP4
    
    H4 -- "Contém" --> SNP1
    H4 -- "Contém" --> SNP2
    H4 -- "Contém" --> SNP3
    H4 -- "Contém" --> SNP4
```

*Este modelo representa haplótipos (sequências de alelos em um cromossomo) como nós, com arestas representando a similaridade entre haplótipos. Os SNPs são conectados aos haplótipos que os contêm, permitindo a imputação baseada em estruturas de haplótipos completos.*

## 5. Grafo Espaço-Temporal para Recombinação

```mermaid
graph TD
    subgraph "Geração Atual"
        G0_SNP1((rs123))
        G0_SNP2((rs456))
        G0_SNP3((rs789))
        G0_SNP4((rs012))
        
        G0_SNP1 -- "r²=0.85" --- G0_SNP2
        G0_SNP2 -- "r²=0.72" --- G0_SNP3
        G0_SNP1 -- "r²=0.65" --- G0_SNP3
        G0_SNP3 -- "r²=0.91" --- G0_SNP4
        
        style G0_SNP3 fill:#f66,stroke:#333,stroke-width:2px
    end
    
    subgraph "Geração Ancestral"
        G1_SNP1((rs123))
        G1_SNP2((rs456))
        G1_SNP3((rs789))
        G1_SNP4((rs012))
        
        G1_SNP1 -- "r²=0.92" --- G1_SNP2
        G1_SNP2 -- "r²=0.88" --- G1_SNP3
        G1_SNP1 -- "r²=0.81" --- G1_SNP3
        G1_SNP3 -- "r²=0.95" --- G1_SNP4
    end
    
    G0_SNP1 -- "Evolução" --> G1_SNP1
    G0_SNP2 -- "Evolução" --> G1_SNP2
    G0_SNP3 -- "Evolução" --> G1_SNP3
    G0_SNP4 -- "Evolução" --> G1_SNP4
    
    classDef known fill:#6b6,stroke:#333,stroke-width:1px
    classDef unknown fill:#f66,stroke:#333,stroke-width:2px
    classDef ancestral fill:#aaf,stroke:#333,stroke-width:1px
    
    class G0_SNP1,G0_SNP2,G0_SNP4 known
    class G0_SNP3 unknown
    class G1_SNP1,G1_SNP2,G1_SNP3,G1_SNP4 ancestral
```

*Este modelo incorpora informações evolutivas, conectando SNPs entre gerações para capturar padrões de recombinação e deriva genética ao longo do tempo. Isso permite que o modelo aprenda como os padrões de LD mudam através das gerações.*

## 6. Grafo Heterogêneo com Múltiplos Tipos de Nós

```mermaid
graph TD
    subgraph "Estrutura Genômica"
        SNP1((rs123))
        SNP2((rs456))
        SNP3((rs789))
        SNP4((rs012))
        
        Gene1[Gene BRCA1]
        Pathway1{Via de Sinalização}
        Disease1[Doença: Câncer de Mama]
        
        SNP1 -- "Localizado em" --> Gene1
        SNP2 -- "Localizado em" --> Gene1
        SNP3 -- "Localizado em" --> Gene1
        SNP4 -- "Próximo a" --> Gene1
        
        Gene1 -- "Participa" --> Pathway1
        Pathway1 -- "Associado a" --> Disease1
        
        SNP1 -- "r²=0.85" --- SNP2
        SNP2 -- "r²=0.72" --- SNP3
        SNP1 -- "r²=0.65" --- SNP3
        SNP3 -- "r²=0.91" --- SNP4
        
        style SNP3 fill:#f66,stroke:#333,stroke-width:2px
    end
    
    classDef snp fill:#6b6,stroke:#333,stroke-width:1px
    classDef unknown fill:#f66,stroke:#333,stroke-width:2px
    classDef gene fill:#69f,stroke:#333,stroke-width:1px
    classDef pathway fill:#f96,stroke:#333,stroke-width:1px
    classDef disease fill:#c9f,stroke:#333,stroke-width:1px
    
    class SNP1,SNP2,SNP4 snp
    class SNP3 unknown
    class Gene1 gene
    class Pathway1 pathway
    class Disease1 disease
```

*Este modelo heterogêneo integra múltiplos tipos de informação biológica, incluindo SNPs, genes, vias de sinalização e doenças. Esta abordagem permite que o GNN capture relações complexas entre diferentes entidades biológicas para melhorar a imputação.*

## 7. Grafo de Janela Deslizante para Processamento de Genoma Completo

```mermaid
graph LR
    subgraph "Janela 1 (Chr1:1000-2000)"
        W1_SNP1((rs123))
        W1_SNP2((rs456))
        W1_SNP3((rs789))
        
        W1_SNP1 -- "r²=0.85" --- W1_SNP2
        W1_SNP2 -- "r²=0.72" --- W1_SNP3
        W1_SNP1 -- "r²=0.65" --- W1_SNP3
        
        style W1_SNP3 fill:#f66,stroke:#333,stroke-width:2px
    end
    
    subgraph "Janela 2 (Chr1:1500-2500)"
        W2_SNP3((rs789))
        W2_SNP4((rs012))
        W2_SNP5((rs345))
        
        W2_SNP3 -- "r²=0.91" --- W2_SNP4
        W2_SNP4 -- "r²=0.45" --- W2_SNP5
        W2_SNP3 -- "r²=0.38" --- W2_SNP5
        
        style W2_SNP3 fill:#f66,stroke:#333,stroke-width:2px
    end
    
    W1_SNP3 -- "Mesmo SNP" --- W2_SNP3
    
    classDef known fill:#6b6,stroke:#333,stroke-width:1px
    classDef unknown fill:#f66,stroke:#333,stroke-width:2px
    
    class W1_SNP1,W1_SNP2,W2_SNP4,W2_SNP5 known
    class W1_SNP3,W2_SNP3 unknown
```

*Este modelo representa uma abordagem de janela deslizante para processar genomas completos, onde janelas sobrepostas são modeladas como subgrafos conectados. Os SNPs que aparecem em múltiplas janelas servem como pontos de conexão, permitindo a propagação de informação entre janelas.*