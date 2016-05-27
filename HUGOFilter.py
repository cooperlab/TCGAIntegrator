import json
import requests


def HUGOFilter(Symbols, SearchStr="http://rest.genenames.org/search/status/"):
    """Performs filtering of a list of gene symbols against the HGNC HUGO
    database to remove non-protein coding genes and unapproved gene symbols.

    Parameters
    ----------
    Symbols : list
        A list of gene symbol strings.
    SearchStr : string
        Base string used to search the genenames.org server. Default value =
        "http://rest.genenames.org/search/status/".

    Returns
    -------
    Approved : list
        A list of the approved HUGO gene symbols representing protein-coding
        genes from input list 'Symbols'.
    Indices : array_like
        The corresponding indices of the approved symbols representing their
        position from the input list 'Symbols'. Used for additional filtering
        outside this function.
    """

    # string for finding approved genes with protein products
    AppStr = "Approved+AND+locus_type:%22gene%20with%20protein%20product%22"

    # build set of official protein-coding gene symbols
    HGNC = requests.get(SearchStr + AppStr,
                        headers={'Accept': 'application/json'})
    if HGNC.status_code != requests.codes.ok:
        print "error"
    else:
        HUGO = json.loads(HGNC.content)
        HUGO = [str(gene['symbol']) for gene in HUGO['response']['docs']]

    # filter amplification and deletion symbols
    Approved = [Gene for Gene in Symbols if Gene in HUGO]
    Indices = [Symbols.index(Gene) for Gene in Symbols if Gene in HUGO]

    # return outputs
    return Approved, Indices
