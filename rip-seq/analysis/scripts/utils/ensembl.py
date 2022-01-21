from collections import defaultdict
from functools import lru_cache
from typing import Tuple, Set

import pandas as pd

from .paths import RESOURCES

BIOMART_MAPPING = RESOURCES.joinpath("ENSEMBL-gene-transcripts-info.tsv")


@lru_cache(maxsize=1)
def _load_id_to_name():
    df = pd.read_csv(BIOMART_MAPPING, sep="\t")
    id_to_name = {k: v for k, v in zip(df['Gene stable ID'].values, df["Gene name"].values)}
    id_to_name.update({k: v for k, v in zip(df['Gene stable ID version'].values, df["Gene name"].values)})
    return id_to_name


def id_to_name(id: str) -> str:
    id_to_name = _load_id_to_name()
    if id in id_to_name:
        return id_to_name[id]
    raise ValueError(f"Unknown gene ID: {id}")


@lru_cache(maxsize=1)
def _load_synonym_to_name() -> Tuple[dict, dict, dict]:
    df = pd.read_csv(BIOMART_MAPPING, sep="\t")
    synonym_to_name = {str(k): str(v) for k, v in zip(df['Gene Synonym'].values, df['Gene name'].values)}
    synonym_to_name.update({str(k): str(k) for k in df['Gene name'].values})
    synonym_to_name.update({
        str(k): str(v) for k, v in zip(df['NCBI gene (formerly Entrezgene) accession'].values, df['Gene name'].values)
    })

    lower_synonym_to_name = {k.lower(): v for k, v in synonym_to_name.items()}
    upper_synonym_to_name = {k.upper(): v for k, v in synonym_to_name.items()}
    return synonym_to_name, lower_synonym_to_name, upper_synonym_to_name


def synonym_to_name(gene: str) -> str:
    synonym_to_name, lower_synonym_to_name, upper_synonym_to_name = _load_synonym_to_name()
    if gene in synonym_to_name:
        return synonym_to_name[gene]
    lgene = gene.lower()
    if lgene in lower_synonym_to_name:
        return lower_synonym_to_name[lgene]
    ugene = gene.upper()
    if ugene in upper_synonym_to_name:
        return upper_synonym_to_name[ugene]
    raise ValueError(f"Unknown gene synonym: {gene}")


@lru_cache(maxsize=1)
def _load_name_to_id():
    df = pd.read_csv(BIOMART_MAPPING, sep="\t")
    name_to_id = {k: v for k, v in zip(df['Gene name'].values, df["Gene stable ID"].values)}
    return name_to_id


def name_to_id(gene: str) -> str:
    gene = synonym_to_name(gene)

    name_to_id = _load_name_to_id()
    if gene in name_to_id:
        return name_to_id[gene]
    raise ValueError(f"Unknown gene name: {gene}")


@lru_cache(maxsize=1)
def _load_transcript_to_id():
    df = pd.read_csv(BIOMART_MAPPING, sep="\t")
    transcript_to_id = {k: v for k, v in zip(df['Transcript stable ID'].values, df["Gene stable ID"].values)}
    transcript_to_id.update(
        {k: v for k, v in zip(df['Transcript stable ID version'].values, df['Gene stable ID'].values)}
    )
    return transcript_to_id


def transcript_to_gene_name(transcriptid: str) -> str:
    transcript_to_id = _load_transcript_to_id()
    if transcriptid not in transcript_to_id:
        transcriptid = transcriptid.split(".")[0]
        if transcriptid not in transcript_to_id:
            raise ValueError(f"Unknown transcript ID: {transcriptid}")

    geneid = transcript_to_id[transcriptid]
    return id_to_name(geneid)


def transcript_to_gene_id(transcriptid: str) -> str:
    transcript_to_id = _load_transcript_to_id()
    if transcriptid not in transcript_to_id:
        transcriptid = transcriptid.split(".")[0]

    if transcriptid not in transcript_to_id:
        raise ValueError(f"Unknown transcript ID: {transcriptid}")
    return transcript_to_id[transcriptid]


@lru_cache(maxsize=1)
def _load_gene_id_to_type():
    df = pd.read_csv(BIOMART_MAPPING, sep="\t")
    geneid_to_type = {k: v for k, v in zip(df['Gene stable ID'].values, df["Gene type"].values)}
    geneid_to_type.update(
        {k: v for k, v in zip(df["Gene stable ID version"].values, df["Gene type"].values)}
    )
    return geneid_to_type


def gene_id_to_type(geneid: str) -> str:
    geneid_to_type = _load_gene_id_to_type()
    if geneid not in geneid_to_type:
        geneid = geneid.split(".")[0]

    if geneid not in geneid_to_type:
        raise ValueError(f"Unknown gene ID: {geneid}")
    return geneid_to_type[geneid]


@lru_cache(maxsize=1)
def _load_transcript_id_to_type():
    df = pd.read_csv(BIOMART_MAPPING, sep="\t")
    transcript_id_to_type = {k: v for k, v in zip(df['Transcript stable ID'].values, df["Transcript type"].values)}
    transcript_id_to_type.update(
        {k: v for k, v in zip(df["Transcript stable ID version"].values, df["Transcript type"].values)}
    )
    return transcript_id_to_type


def transcript_id_to_type(tid: str) -> str:
    transcript_id_to_type = _load_transcript_id_to_type()
    if tid not in transcript_id_to_type:
        tid = tid.split(".")[0]

    if tid not in transcript_id_to_type:
        raise ValueError(f"Unknown Transcript ID: {tid}")
    return transcript_id_to_type[tid]


@lru_cache(maxsize=1)
def _load_gene_id_to_transcripts():
    df = pd.read_csv(BIOMART_MAPPING, sep="\t")

    mapping = defaultdict(set)
    for gid, tid in zip(df['Gene stable ID'], df['Transcript stable ID']):
        mapping[gid].add(tid)
    return mapping


def gene_id_to_transcripts(tid: str) -> Set[str]:
    mapping = _load_gene_id_to_transcripts()
    return mapping.get(tid, [])
