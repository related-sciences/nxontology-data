# Import ontologies to NXOntology

This repository imports public ontologies/taxonomies into Python NXOntology objects
and writes the ontologies in the JSON-based [node-link](https://networkx.org/documentation/stable/reference/readwrite/generated/networkx.readwrite.json_graph.node_link_graph.html) data format.
The goal is to standardize and simplify data access to ontologies.

For ontologies that have been imported into NXOntology and exported to JSON,
see the `output/*` branches on GitHub,
for example [`output/pubchem`](https://github.com/related-sciences/nxontology-imports/tree/output/pubchem).

Once you find the ontology you'd like to read,
you can read in Python
(after installing any dependenies like `pip install nxontology`):

```py
# URL to the exported dataset.
# Here we read the ChEMBL protein/target classification hierarchy.
url = "https://github.com/related-sciences/nxontology-imports/raw/output/pubchem/087_chembl_target_tree.json"
# Versioning with the commit hash is a good idea, since we might change the branch structure where data is stored.
url = "https://github.com/related-sciences/nxontology-imports/raw/71cf538dc5c258ada880d58663b0205b7b7f8561/087_chembl_target_tree.json"

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

This repository is released under an Apache License 2.0 License (see [LICENSE.md](LICENSE.md)).
Please see the source licenses of ontologies imported into NXOntology objects.
