# ceRNA Pipeline Workflow Diagram

```
                    DUAL-INPUT ceRNA DISCOVERY PIPELINE
                                      
┌─────────────────┬─────────────────┐
│   RNA-seq       │   miRNA-seq     │  INPUT FILES
│  (mRNA/lncRNA)  │   (miRNAs)      │  
└─────────┬───────┴─────────┬───────┘
          │                 │
          └─────────┬───────┘
                    │
          ┌─────────▼──────────┐
          │  QC & Normalization │  STAGE 1: Data preprocessing
          │  • Sample alignment │  • CPM normalization  
          │  • Gene ID mapping  │  • Low expression filtering
          │  • Quality control  │  • Log2 transformation
          └─────────┬──────────┘
                    │
          ┌─────────▼──────────┐
          │ Database Download   │  STAGE 2: Knowledge bases
          │  • miRTarBase       │  • miRNA-mRNA interactions
          │  • LncBase         │  • miRNA-lncRNA interactions  
          │  • Gene annotation  │  • Biotype information
          └─────────┬──────────┘
                    │
          ┌─────────▼──────────┐
          │ Feature Engineering │  STAGE 3: Candidate generation
          │  • Triplet creation │  • Expression correlations
          │  • SPONGE scores    │  • Partial correlations
          │  • Quality filters  │  • Sequence features
          └─────────┬──────────┘
                    │
          ┌─────────▼──────────┐
          │   ML Training      │  STAGE 4: Predictive modeling
          │  • XGBoost model   │  • Feature selection
          │  • Label generation │  • Cross-validation
          │  • Model validation │  • Hyperparameter tuning
          └─────────┬──────────┘
                    │
          ┌─────────▼──────────┐
          │ Triplet Prediction │  STAGE 5: Apply ML model
          │  • Probability     │  • Confidence scoring
          │    scoring         │  • Ranking predictions
          │  • High-confidence │  • Initial filtering
          │    filtering       │
          └─────────┬──────────┘
                    │
          ┌─────────▼──────────┐
          │Statistical Validation│  STAGE 6: Mediation analysis
          │  • Sobel test      │  • Indirect effect testing
          │  • Mediation types │  • Multiple testing correction
          │  • FDR correction  │  • Effect size calculation
          └─────────┬──────────┘
                    │
          ┌─────────▼──────────┐
          │ Network Analysis   │  STAGE 7: Network construction  
          │  • Graph building  │  • Centrality analysis
          │  • Hub detection   │  • Topology metrics
          │  • Community finding│  • Export formats
          └─────────┬──────────┘
                    │
          ┌─────────▼──────────┐
          │ Report Generation  │  STAGE 8: Results summary
          │  • Interactive HTML│  • Statistical plots
          │  • Network plots   │  • Top discoveries
          │  • Summary tables  │  • Quality metrics
          └─────────┬──────────┘
                    │
┌─────────────────────▼────────────────────┐
│              FINAL OUTPUTS               │
│ • validated_triplets_ranked.csv         │
│ • cerna_analysis_report.html            │  
│ • cerna_network.graphml                 │
│ • centrality_scores.csv                 │
└──────────────────────────────────────────┘

                    ceRNA BIOLOGICAL MODEL

    lncRNA ◄──┐                    ┌──► mRNA
              │                    │
              ▼                    ▲
            miRNA ◄────────────────┘
              │                     
              └─── Competition ─────┘
              
    • lncRNA "sponges" miRNA away from mRNA
    • Reduces miRNA-mediated repression  
    • Results in positive lncRNA-mRNA correlation
    • miRNA acts as regulatory mediator
```