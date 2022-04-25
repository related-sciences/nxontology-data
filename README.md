# NXOntology data: making ontologies accessible as simple JSON files

[![GitHub Actions CI Build Status](https://img.shields.io/github/workflow/status/related-sciences/nxontology/Build/main?label=actions&style=for-the-badge&logo=github&logoColor=white)](https://github.com/related-sciences/nxontology/actions)  
[![Software License](https://img.shields.io/github/license/related-sciences/nxontology?style=for-the-badge&logo=Apache&logoColor=white)](https://github.com/related-sciences/nxontology/blob/main/LICENSE)  
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?style=for-the-badge&logo=Python&logoColor=white)](https://github.com/psf/black)  

This repository imports public ontologies/taxonomies into Python [NXOntology](https://github.com/related-sciences/nxontology) objects
and writes the ontologies in the JSON-based [node-link](https://networkx.org/documentation/stable/reference/readwrite/generated/networkx.readwrite.json_graph.node_link_graph.html) data format.
The goal is to standardize and simplify data access to ontologies.

For ontologies that have been imported into NXOntology and exported to JSON,
see the `output/*` branches on GitHub,
for example [`output/pubchem`](https://github.com/related-sciences/nxontology-data/tree/output/pubchem).

Once you find the ontology you'd like to read,
you can read in Python
(after installing any dependenies like `pip install nxontology`):

```py
# URL to the exported dataset.
# Here we read the ChEMBL protein/target classification hierarchy.
url = "https://github.com/related-sciences/nxontology-data/raw/output/pubchem/087_chembl_target_tree.json"
# Versioning with the commit hash is a good idea, since we might change the branch structure where data is stored.
url = "https://github.com/related-sciences/nxontology-data/raw/71cf538dc5c258ada880d58663b0205b7b7f8561/087_chembl_target_tree.json"

# To read as an NXOntology object,
# which encapsulates the networkx graph.
# Will also work for the gzip compressed files.
from nxontology import NXOntology
nxo = NXOntology.read_node_link_json(url)

# To read as a networkx.DiGraph
import requests
from networkx.readwrite.json_graph import node_link_graph
digraph = node_link_graph(requests.get(url).json())
```

## Sources

The data sources that are currently imported are listed below.
Please open an issue if you are interested in contributing support for additional sources.

### MeSH

MeSH (Medical Subject Headings) is created by the National Library Medicine and integrated into many projects including PubMed.
See [nxontology_data/mesh](nxontology_data/mesh) for a detailed README.

### PubChem

We import ontologies from the PubChem Classifications service
(see [browser](https://pubchem.ncbi.nlm.nih.gov/classification/) & [docs](https://pubchem.ncbi.nlm.nih.gov/classification/docs/classification_help.html)).
Most ontologies indexed by service do not originate with PubChem,
but PubChem provides convenient and standardized bulk access.

## Development

```shell
# Install the environment
poetry install --no-root

# Update the lock file
poetry update

# Run tests
pytest

# Set up the git pre-commit hooks.
# `git commit` will now trigger automatic checks including linting.
pre-commit install

# Run all pre-commit checks (CI will also run this).
pre-commit run --all
```

## License

This source code in this repository is released under an Apache License 2.0 License
(see [LICENSE.md](LICENSE.md)).
Source code refers to the contents of the `main` branch and any other development branches containing code and documentation.

The output branches contain data from external ontologies.
Please refer to each respective ontology for its data license.
If available, we include license information in the graph metadata for each ontology,
but often license information is not supplied in the ontology data we ingest.
Please attribute the source ontology when reusing data obtained from this project,
and as best practice mention that the data was obtained via [NXOntology data](https://github.com/related-sciences/nxontology-data).

Any original data produced by this repository is released under a [CC0 1.0 Universal (CC0 1.0) Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/).
As noted above, the underlying ontology data is not original to this repository and upstream licenses should be consulted.
