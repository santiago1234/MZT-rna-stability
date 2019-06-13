## Expression Data *Owens et al., 2016, Cell Reports*

get the data

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65785/suppl/GSE65785_clutchApolyA_absolute_TPE_gene_isoform.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65785/suppl/GSE65785_clutchA_rdRNA_absolute_TPE_gene_isoform.txt.gz
gunzip *gz
```


## sequence data

the sequence data was obtained from biomart website

```bash
fasta_to_tab -in martquery_0613173836_403.txt >xenopus_3utr.txt
fasta_to_tab -in martquery_0613173535_586.txt >xenopus_coding.txt
```


