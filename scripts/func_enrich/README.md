# Functional enrichment
The scripts included in this folder were used to perform the network enrichment analysis tests ([NEAT](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1203-6)) described in our study. These enrichments used the[FunCoup 5.0](https://funcoup5.scilifelab.se/search/) *H. sapiens* as the reference network. The network file, which can be downloaded [here](https://funcoup5.scilifelab.se/downloads/download.action?type=network&instanceID=24480085&fileName=FC5.0_H.sapiens_full.gz), must be placed into `data/funCoup/`. The following enrichment analyses were performed:
## 1. Functional enrichment of the predictor genes from our transcriptomic clock
Genes with non-zero coefficients obtained in the second round of fitting, after running successfully `fitModel_second_round.sh` in the [`model_training`](../../scripts/model_training) module, were used as input for NEAT. The gene sets used were the Gene Ontology terms, as they appear in `org.Hs.eg.db` R package (v3.17.0). For performing this enrichment run:
```bash
sbatch launch_NEAT_enrich_modGenes.sh
```
Result will be saved into `results/enrich_NEAT/mod_second_round/`. 
## 2. Functional enrichment of the DEGs induced by the three compound treatment in mice
We used DEGs obtained by comparing cortex transcriptional profiles of mice before and after the treatment with azacitidine, tranylcypromine and JNK-IN-8. Different enrichments were done with these DEGs.
### 2.1. On Gene Ontology terms
GO terms are taken from the `org.Hs.eg.db` R package (v3.17.0).
#### 2.1.1. With all the genes
Run:
```bash
sbatch launch_NEAT_enrich_mice_DEGs_all_GO.sh
```
#### 2.1.2. With upregulated genes
Run:
```bash
sbatch launch_NEAT_enrich_mice_DEGs_up_GO.sh
```
#### 2.1.3. With downregulated genes
Run:
```bash
sbatch launch_NEAT_enrich_mice_DEGs_down_GO.sh
```
### 2.2. On MSigDB gene sets
These analyses were performed using MSigDB gene sets.
#### 2.2.1. With all the genes
Run:
```bash
sbatch launch_NEAT_enrich_mice_DEGs_all.sh
```
#### 2.2.2. With upregulated genes
Run:
```bash
sbatch launch_NEAT_enrich_mice_DEGs_up.sh
```
#### 2.3.2. With downregulated genes
Run:
```bash
sbatch launch_NEAT_enrich_mice_DEGs_down.sh
```