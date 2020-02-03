import pysam

from idbd_bio_utils import NcbiTaxonomy


def get_class(taxid, ncbi=None):
    if ncbi is None:
        ncbi = NcbiTaxonomy()
    human_taxid = 9606
    bacteria_taxid = 2
    virus_taxid = 10239
    eukaryota_taxid = 2759
    parasite_taxids = {
        7563,
        188941,
        6029,
        5653,
        6935,
        6178,
        5794,
        6308,
        31277,
        119088,
        6199,
        85819,
        33083,
        33084,
        75966,
        41165,
        7509,
        6236,
        198624,
        33634,
        5988,
        6249,
        5738,
        1489900,
        740972,
        1485168,
        37104,
        10232
    }
    fungus_taxid = 4751
    archaea_taxid = 2157

    taxid_path = ncbi.get_path(taxid)
    if human_taxid in taxid_path:
        org_class = 'human'
    elif fungus_taxid in taxid_path:
        org_class = 'fungus'
    elif parasite_taxids.intersection(set(taxid_path)) != set():
        org_class = 'parasite'
    elif bacteria_taxid in taxid_path:
        org_class = 'bacteria'
    elif virus_taxid in taxid_path:
        org_class = 'virus'
    elif archaea_taxid in taxid_path:
        org_class = 'archaea'
    else:
        org_class = 'unclassified'
    return org_class


def get_genome_size(genome_path):
    with pysam.FastxFile(genome_path) as fh:
        genome_size = 0
        for entry in fh:
            genome_size += len(entry.sequence)
    return genome_size