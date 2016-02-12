class Attribute(object):
    """
    Stores attributes from the gencode attribute file.
    """

    __slots__ = ("gene_id", "gene_name", "gene_type", "transcript_id", "transcript_type")

    def __init__(self, gene_id, gene_name, gene_type, transcript_id, transcript_type):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_type = gene_type
        self.transcript_id = transcript_id
        self.transcript_type = transcript_type


def get_transcript_attribute_dict(attribute_file):
    """
    Returns a dictionary mapping the transcript ID to an Attribute object.
    This stores all of the relevant information from the gencode attributes file.
    """
    attribute_dict = {}
    with open(attribute_file) as f:
        for line in f:
            line = line.split()
            gene_id, gene_name, gene_type, transcript_id, transcript_type = line
            attribute_dict[transcript_id] = Attribute(gene_id, gene_name, gene_type, transcript_id, transcript_type)
    return attribute_dict
