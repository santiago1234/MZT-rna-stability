# Arielome Data Sets

In this folder I have a copy of the Arielome data.

This data set was supplied by Ariel Bazzini. Arielome in zebra fish.

```bash
# A copied the data from another folder
cp -r  ../../../rna-stability-model/experiments/18-05-30-ArielomeAsUTRome/data/Arielome/ Arielome
```

Description of data folder: **Arielome |**

This folder contains the file/folder(s):

+ **2h |** 2hrs arielome reads
+ **8h |** 8 hrs arielome reads
+ **8hMO |** 8 hrs with Morpholine arielome reads, one replicate is split for some
unknow reason to me.

for each of the file above we have 3 replicates

The files *Arielome/metadata.csv* contains the metadata information for each file.

# Procesing the libraries

```bash
python ./processdata.py
```
