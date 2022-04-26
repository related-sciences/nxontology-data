from __future__ import annotations

import functools
import logging
import pathlib
import re
import tempfile
from typing import Iterator
from urllib.request import urlretrieve

import bioversions
import fsspec
import networkx as nx
import nxontology
import pandas as pd
import rdflib
from nxontology import NXOntology
from rdflib.term import URIRef

from nxontology_data.utils import (
    get_output_dir,
    sparql_results_to_df,
    write_dataframe,
    write_ontology,
)

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

    MESH_RDF_ROOT = "https://nlmpubs.nlm.nih.gov/projects/mesh/rdf"

    @classmethod
    def get_mesh_rdf(cls, year_yyyy: str) -> rdflib.Graph:
        """
        Read MeSH into rdflib from the MeSH RDF FPT site.
        https://www.nlm.nih.gov/databases/download/mesh.html
        """
        # The .nt.gz file is around 115 MB.
        # Reading from HTTPS/FTP to rdflib was causing timeout errors,
        # so using a local temporary directory instead.
        temp_dir = tempfile.mkdtemp(suffix="_nxontology_data_mesh")
        logging.info(f"Copying mesh rdf files to {temp_dir}")
        nt_filename = f"mesh{year_yyyy}.nt.gz"
        for filename in "vocabulary_1.0.0.ttl", nt_filename:
            url = f"{cls.MESH_RDF_ROOT}/{year_yyyy}/{filename}"
            urlretrieve(url, f"{temp_dir}/{filename}")
        return cls._read_mesh_rdf(directory=temp_dir, nt_filename=nt_filename)

    @staticmethod
    @functools.cache
    def _read_mesh_rdf(directory: str, nt_filename: str) -> rdflib.Graph:
        """
        directory: local directory with raw MeSH RDF files.
        """
        logger.info(f"Loading MeSH into rdflib from {directory}")
        rdf = rdflib.Graph()
        rdf.namespace_manager.bind("meshv", "http://id.nlm.nih.gov/mesh/vocab#")
        # load MeSH vocabulary (takes ~2 seconds)
        with fsspec.open(f"{directory}/vocabulary_1.0.0.ttl", "rt") as src:
            # https://github.com/HHS/meshrdf/issues/153
            rdf.parse(source=src, format="n3")
        # load MeSH triples (takes ~30 minutes)
        logger.info(f"Loading triples from {nt_filename}")
        with fsspec.open(
            f"{directory}/{nt_filename}", mode="rb", compression="infer"
        ) as src:
            # read in binary mode https://github.com/RDFLib/rdflib/issues/1144
            rdf.parse(source=src, format="nt")
        # When directory is an HTTPS or FTP URL, we encountered several issues:
        # - Hit error with fsspec/aiohttp: Can not decode content-encoding: gzip
        # https://github.com/fsspec/filesystem_spec/issues/389
        # - gzip transfer-encoding is specified such that compression is automatically decoded.
        # https://github.com/HHS/meshrdf/issues/193
        # - requests.get with `stream=True` and `response.raw.decode_content=True` choked
        # rdflib.exceptions.ParserError: Invalid line: <http://id.nlm.nih.gov/mesh/2021/D00895
        logger.info("Reading rdflib.Graph is complete.")
        return rdf

    @staticmethod
    def _get_query(name: str) -> str:
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
            rdf.cached_query = functools.cache(rdf.query)
        query = cls._get_query(name)
        results = (rdf.cached_query if cache else rdf.query)(query)
        return sparql_results_to_df(results)

    @staticmethod
    def _mesh_uri_to_id(uri: URIRef) -> str:
        if not uri.startswith("http://id.nlm.nih.gov/mesh/"):
            raise ValueError(f"{uri} does not look like a MeSH identifier")
        return str(uri).rsplit("/", 1)[1]

    @staticmethod
    def _get_relationship_triples(
        rdf: rdflib.Graph,
    ) -> Iterator[tuple[URIRef, URIRef, URIRef]]:
        for rel_type in "broaderDescriptor", "preferredMappedTo", "mappedTo":
            # expand CURIE to URI https://github.com/RDFLib/rdflib/issues/626
            predicate = rdflib.util.from_n3(
                f"meshv:{rel_type}", nsm=rdf.namespace_manager
            )
            yield from rdf.triples((None, predicate, None))

    @classmethod
    def _get_id_to_tree_numbers(cls, rdf: rdflib.Graph) -> dict[str, str]:
        tree_number_df = cls.run_query(rdf, "tree-numbers")
        return {
            mesh_id: tns.to_list()
            for mesh_id, tns in tree_number_df.groupby("mesh_id").tree_number
        }

    @classmethod
    def get_identifier_df(cls, rdf: rdflib.Graph) -> pd.DataFrame:
        id_df = cls.run_query(rdf, "identifiers")
        id_df["tree_numbers"] = id_df["mesh_id"].map(cls._get_id_to_tree_numbers(rdf))
        return id_df

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
        "tree_numbers",
    ]

    @classmethod
    def create_nxo(
        cls, rdf: rdflib.Graph, year_yyyy: str
    ) -> tuple[NXOntology[str], pd.DataFrame]:
        nxo: NXOntology[str] = NXOntology()
        nxo.graph.graph["name"] = "mesh"
        nxo.graph.graph["description"] = "Medical Subject Headings"
        nxo.graph.graph["mesh_year"] = str(year_yyyy)
        nxo.set_graph_attributes(
            node_name_attribute="mesh_label",
            node_identifier_attribute="mesh_id",
            node_url_attribute="mesh_uri",
        )
        # add nodes
        id_df = cls.get_identifier_df(rdf)
        for row in (
            id_df[cls._node_attrs]
            .query("mesh_class in @cls._node_classes")
            .to_dict(orient="records")
        ):
            mesh_id = row["mesh_id"]
            nxo.add_node(mesh_id, **row)
        # add edges
        for s, p, o in cls._get_relationship_triples(rdf):
            try:
                nxo.add_edge(
                    u_of_edge=cls._mesh_uri_to_id(o),
                    v_of_edge=cls._mesh_uri_to_id(s),
                    predicate=rdf.namespace_manager.qname_strict(p),
                )
            except nxontology.exceptions.NodeNotFound:
                # meshv:AllowedDescriptorQualifierPair nodes like D014199Q000031 aren't included as nodes
                pass
        # TODO: should we remove disconnected nodes
        # nxo.graph = nxo.graph.remove_nodes_from(nxo.isolates)
        id_df["in_nxo"] = id_df.mesh_id.isin(set(nxo.graph))
        return nxo, id_df

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
                if root_info["mesh_class"] != "TopicalDescriptor":
                    continue
                (tree_number,) = root_info["tree_numbers"]
                rows.append(
                    {
                        "mesh_id": info.data["mesh_id"],
                        "mesh_label": info.data["mesh_label"],
                        "mesh_class": info.data["mesh_class"],
                        "top_mesh_id": root,
                        "top_tree_number": tree_number,
                        "top_mesh_label": root_info["mesh_label"],
                        "top_is_disease": cls._is_disease(tree_number),
                        "depth": nx.shortest_path_length(nxo.graph, root, node),
                    }
                )
        return pd.DataFrame(rows).sort_values(
            ["top_tree_number", "depth", "mesh_class", "mesh_id"],
            ascending=[True, True, False, True],
        )

    @staticmethod
    def _is_disease(tree_number: str) -> bool:
        """
        Whether a MeSH tree number corresponds to a top-level disease category.
        RS internal issue: 370
        """
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

    @classmethod
    def export_mesh_outputs(cls, year_yyyy: str | None = None) -> None:
        if year_yyyy is None:
            year_yyyy = bioversions.get_version("mesh")
        year_yyyy = str(year_yyyy)  # protect against fire
        output_dir = get_output_dir().joinpath("mesh")
        output_dir.mkdir(parents=True, exist_ok=True)
        rdf = cls.get_mesh_rdf(year_yyyy)
        logger.info(f"Creating NXOntology for mesh {year_yyyy}.")
        nxo, id_df = cls.create_nxo(rdf=rdf, year_yyyy=year_yyyy)
        nxo_path = write_ontology(nxo=nxo, output_dir=output_dir)
        logger.info(f"Wrote mesh nxontology to {nxo_path}.")
        write_dataframe(df=id_df, path=output_dir.joinpath("mesh_identifiers.json.gz"))
        logger.info(f"Creating top-level term mapping for mesh {year_yyyy}.")
        top_map_df = cls.create_top_level_map_df(nxo)
        write_dataframe(
            df=top_map_df, path=output_dir.joinpath("mesh_top_level_map.json.gz")
        )
