# Converting MeSH to NXOntology

Key links:

- [`HHS/meshrdf` GitHub Repository](https://github.com/HHS/meshrdf).
  Issues below:
    - [Difference between mappedTo and preferredMappedTo predicates?](https://github.com/HHS/meshrdf/issues/155)
    - [Comment by Dan Davis on creating a TopicalDescriptor Ontology](https://github.com/HHS/meshrdf/issues/156#issuecomment-752226217)
    - [Pagination / offset of results for the SPARQL endpoint is faulty](https://github.com/HHS/meshrdf/issues/150)
- [MeSH RDF Technical Documentation](https://hhs.github.io/meshrdf/)
- [MeSH Homepage](https://www.nlm.nih.gov/mesh/meshhome.html)
- [MeSH SPARQL Endpoint](https://id.nlm.nih.gov/mesh/query)

## References

Works related to converting MeSH to a directed acyclic graph:

1. **Desiderata for an authoritative Representation of MeSH in RDF**   
Rainer Winnenburg, Olivier Bodenreider  
*AMIA* (2014-11-14) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4419968/>   
PMID: [25954433](https://www.ncbi.nlm.nih.gov/pubmed/25954433) · PMCID: [PMC4419968](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4419968)

2. **Transforming the Medical Subject Headings into Linked Data: Creating the Authorized Version of MeSH in RDF**   
Barbara Bushman, David Anderson, Gang Fu  
*Journal of library metadata* (2015) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4749162/>   
DOI: [10.1080/19386389.2015.1099967](https://doi.org/10.1080/19386389.2015.1099967) · PMID: [26877832](https://www.ncbi.nlm.nih.gov/pubmed/26877832) · PMCID: [PMC4749162](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4749162)

3. **MeSHHeading2vec: a new method for representing MeSH headings as vectors based on graph embedding algorithm**   
Zhen-Hao Guo, Zhu-Hong You, De-Shuang Huang, Hai-Cheng Yi, Kai Zheng, Zhan-Heng Chen, Yan-Bin Wang  
*Briefings in Bioinformatics* (2020-03-31) <https://doi.org/ghnr5n>   
DOI: [10.1093/bib/bbaa037](https://doi.org/10.1093/bib/bbaa037) · PMID: [32232320](https://www.ncbi.nlm.nih.gov/pubmed/32232320)

4. [Semantic Measures Library & ToolKit Medical Subject Headings (MeSH)](https://www.semantic-measures-library.org/sml/index.php?q=doc&page=mesh)

5. **MeSHSim: An R/Bioconductor package for measuring semantic similarity over MeSH headings and MEDLINE documents**   
Zhou Jing, Shui Yuxuan, Peng Shengwen, Li Xuhui, Mamitsuka Hiroshi, Zhu Shanfeng  
*Institute of Electrical and Electronics Engineers (IEEE)* (2015-07) <https://doi.org/ghnr5x>   
DOI: [10.1109/chicc.2015.7260989](https://doi.org/10.1109/chicc.2015.7260989)

6. **pyMeSHSim: an integrative python package for biomedical named entity recognition, normalization, and comparison of MeSH terms**   
Zhi-Hui Luo, Meng-Wei Shi, Zhuang Yang, Hong-Yu Zhang, Zhen-Xia Chen  
*BMC Bioinformatics* (2020-06-18) <https://doi.org/ghnvj9>   
DOI: [10.1186/s12859-020-03583-6](https://doi.org/10.1186/s12859-020-03583-6) · PMID: [32552728](https://www.ncbi.nlm.nih.gov/pubmed/32552728) · PMCID: [PMC7301509](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7301509)

### Key Quotations

Note that MeSH violates a core assumption of ontologies.
From the "Contextual hierarchical structure" section of Winnenburg & Bodenreider, 2014:

> MeSH also uses a non-standard hierarchical organization.
The hierarchy among MeSH descriptors is indicated through "tree numbers" assigned to descriptors.
Tree number inclusion reflects that the descriptor with the longer tree number is narrower than that with the shorter tree number.
For example, the tree number for *Liver* [A03.620] has an additional node (.620) compared to that of *Digestive System* [A03],
indicating the narrower relation between the two.
Note that tree numbers are not the unique identifiers of descriptors.
Descriptors often have multiple tree numbers reflecting particular aspects of the descriptors,
each aspect being assigned specific broader and narrower descriptors.
For example, the descriptor *Eye* (`D005123`) has two tree numbers, `A01.456.505.420` and `A09.371`.
In the `A01` tree, *Eye* is narrower than *Head* [`A01.456`] and broader than *Eyebrows* [`A01.456.505.420.338`] and *Eyelids* [`A01.456.505.420.504`],
whereas, in the `A09` tree, *Eye* is narrower than *Sense Organs* [`A09`] and broader than *Eyelids* [`A09.371.337`], *Retina* [`A09.371.729`], *Uvea* [`A09.371.894`] and nine other descriptors.
Note that, although *Head* is broader than *Eye* (`A01` tree),
some descendants of *Eye* in the `A09` tree (e.g., *Retina*) do not have *Head* as their ancestor.
For all practical purposes, **the broader/narrower relationship among MeSH descriptors is not transitive**.
