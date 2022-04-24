import functools
import logging
import pathlib
from typing import Iterator

import fsspec
import networkx as nx
import nxontology
import pandas as pd
import rdflib
from nxontology import NXOntology
from rdflib.term import URIRef

from nxontology_data.utils import sparql_results_to_df

logger = logging.getLogger(__name__)

"""
Regex pattern for valid MeSH identifiers.
References:
- https://registry.identifiers.org/registry/mesh
- https://www.wikidata.org/wiki/Property:P486.
- https://id.nlm.nih.gov/mesh/describe?uri=http%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2Fvocab%23identifier
"""
mesh_id_pattern: str = r"^[CD][0-9]{6}([0-9]{3}|)$"


class MeshLoader:
    @staticmethod
    @functools.lru_cache
    def get_mesh_rdf(
        directory: str, nt_filename: str = "mesh2020.nt.gz"
    ) -> rdflib.Graph:
        """
        directory: local or remote directory with raw MeSH RDF files.
        For example, <ftp://ftp.nlm.nih.gov/online/mesh/rdf/2020>.
        """
        logging.info(f"Loading MeSH into rdflib from {directory}")
        rdf = rdflib.Graph()
        rdf.namespace_manager.bind("meshv", "http://id.nlm.nih.gov/mesh/vocab#")
        # load MeSH vocabulary (takes ~2 seconds)
        with fsspec.open(f"{directory}/vocabulary_1.0.0.ttl", "rt") as src:
            # https://github.com/HHS/meshrdf/issues/153
            rdf.parse(source=src, format="n3")
        # load MeSH triples (takes ~30 minutes)
        logging.info(f"Loading triples from {nt_filename}")
        with fsspec.open(
            f"{directory}/{nt_filename}", "rb", compression="infer"
        ) as src:
            # read in binary mode https://github.com/RDFLib/rdflib/issues/1144
            rdf.parse(source=src, format="nt")
        return rdf

    @staticmethod
    def get_query(name: str) -> str:
        """
        Read SPARQL query from text file.
        """
        return pathlib.Path(__file__).parent.joinpath(f"queries/{name}.rq").read_text()

    @classmethod
    def run_query(
        cls, rdf: rdflib.Graph, name: str, cache: bool = False
    ) -> pd.DataFrame:
        """
        Run SPARQL query on an rdflib.Graph instance of th4e MeSH RDF,
        returning the results as a pandas.DataFrame.
        Enable cache to cache results by the query text (not name/path).
        """
        if cache and not hasattr(rdf, "cached_query"):
            rdf.cached_query = functools.lru_cache(rdf.query)
        query = cls.get_query(name)
        results = (rdf.cached_query if cache else rdf.query)(query)
        return sparql_results_to_df(results)

    @staticmethod
    def mesh_uri_to_id(uri: URIRef) -> str:
        if not uri.startswith("http://id.nlm.nih.gov/mesh/"):
            raise ValueError(f"{uri} does not look like a MeSH identifier")
        return str(uri).rsplit("/", 1)[1]

    @staticmethod
    def get_relationship_triples(
        rdf: rdflib.Graph,
    ) -> Iterator[tuple[URIRef, URIRef, URIRef]]:
        for rel_type in "broaderDescriptor", "preferredMappedTo", "mappedTo":
            # expand CURIE to URI https://github.com/RDFLib/rdflib/issues/626
            predicate = rdflib.util.from_n3(
                f"meshv:{rel_type}", nsm=rdf.namespace_manager
            )
            yield from rdf.triples((None, predicate, None))

    @classmethod
    def get_id_to_tree_numbers(cls, rdf: rdflib.Graph) -> dict[str, str]:
        tree_number_df = cls.run_query(rdf, "tree-numbers")
        return {
            mesh_id: tns.to_list()
            for mesh_id, tns in tree_number_df.groupby("mesh_id").tree_number
        }

    _node_classes = {
        "CheckTag",  # 2 disconnected terms: male and female
        "GeographicalDescriptor",
        "PublicationType",
        "SCR_Chemical",
        "SCR_Disease",
        "SCR_Organism",
        "SCR_Protocol",
        "TopicalDescriptor",
    }
    """Create ontology nodes for these MeSH classes."""

    _node_attrs = [
        "mesh_id",
        "mesh_class",
        "mesh_uri",
        "mesh_label",
    ]

    @classmethod
    def create_nxo(cls, rdf: rdflib.Graph, year_yyyy: str) -> NXOntology[str]:
        nxo: NXOntology[str] = NXOntology()
        nxo.graph.graph["name"] = "MeSH"
        nxo.graph.graph["year"] = str(year_yyyy)
        nxo.set_graph_attributes(
            node_name_attribute="mesh_label",
            node_identifier_attribute="mesh_id",
            node_url_attribute="mesh_uri",
        )
        # add nodes
        id_df = cls.run_query(rdf, "identifiers")
        id_to_tns = cls.get_id_to_tree_numbers(rdf)
        for row in (
            id_df[cls._node_attrs]
            .query("mesh_class in @cls._node_classes")
            .to_dict(orient="records")
        ):
            mesh_id = row["mesh_id"]
            if tns := id_to_tns.get(mesh_id):
                row["tree_numbers"] = tns
            nxo.add_node(mesh_id, **row)
        # add edges
        for s, p, o in cls.get_relationship_triples(rdf):
            try:
                nxo.add_edge(
                    u_of_edge=cls.mesh_uri_to_id(o),
                    v_of_edge=cls.mesh_uri_to_id(s),
                    predicate=rdf.namespace_manager.qname_strict(p),
                )
            except nxontology.exceptions.NodeNotFound:
                # meshv:AllowedDescriptorQualifierPair nodes like D014199Q000031 aren't included as nodes
                pass
        return nxo

    @classmethod
    def create_vocab_digraph(cls, rdf: rdflib.Graph) -> nx.DiGraph:
        subclass_df = cls.run_query(rdf, "vocab-subclasses")
        nx_subclass = nx.DiGraph()
        for row in subclass_df.itertuples():
            nx_subclass.add_edge(row.object_suffix, row.subject_suffix)
        return nx_subclass

    @classmethod
    def create_top_level_map_df(cls, nxo: NXOntology[str]) -> pd.DataFrame:
        """
        Create a table of mesh_id-top_mesh_id pairs.
        Used to associate mesh terms with their top-level categories.
        Includes a `top_is_disease` column based on whether the top-level category
        is a recognized disease category. Filter for `top_is_disease` to get
        assignments of just diseases to their therapeutic areas.
        """
        nxo.freeze()
        rows = list()
        for node in nxo.graph:
            info = nxo.node_info(node)
            for root in info.roots:
                root_info = nxo.node_info(root).data
                (tree_number,) = root_info["tree_numbers"]
                rows.append(
                    {
                        "mesh_id": info.data["mesh_id"],
                        "mesh_label": info.data["mesh_label"],
                        "mesh_class": info.data["mesh_class"],
                        "top_mesh_id": root,
                        "top_tree_number": tree_number,
                        "top_mesh_label": root_info["mesh_label"],
                        "top_is_disease": cls.is_disease(tree_number),
                        "depth": nx.shortest_path_length(nxo.graph, root, node),
                    }
                )
        return pd.DataFrame(rows).sort_values(
            ["top_tree_number", "depth", "mesh_class", "mesh_id"],
            ascending=[True, True, False, True],
        )

    @staticmethod
    def is_disease(tree_number: str) -> bool:
        """
        Whether a MeSH tree number corresponds to a top-level disease category.
        RS internal issue: 370
        """
        import re

        # patterns to include
        include = [
            # Category C is for diseases
            # https://www.nlm.nih.gov/bsd/indexing/training/CATC_010.html
            r"C[0-9]{2}",
            # F03 mental disorders
            # https://www.nlm.nih.gov/bsd/indexing/training/CATF_010.html
            "F03",
        ]
        # patterns to exclude
        # exclude takes precedent over include
        exclude = [
            "C26",  # Wounds and Injuries
        ]
        for pattern in exclude:
            if re.match(pattern, tree_number):
                return False
        for pattern in include:
            if re.match(pattern, tree_number):
                return True
        return False
