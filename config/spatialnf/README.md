10X Visium spatial transcriptomics data was processed with SpatialNF (revision: 8db7eaadc7):

https://github.com/aertslab/SpatialNF


```
nextflow pull https://github.com/aertslab/SpatialNF

nextflow -C spatialnf_multi_sample.config run SpatialNF \
           -entry multi_sample -resume -with-report
```
