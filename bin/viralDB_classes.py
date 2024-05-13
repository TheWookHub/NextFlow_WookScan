
#========================== Organism Class ==========================#
# a class that stores organism info via genbank format file (.gb)

class OrganismClass:
    def __init__(choppa, locus,version,organism,full_seq,accession):
        choppa.locus = str(locus)
        choppa.accession = accession
        choppa.version = str(version)
        choppa.organism = str(organism)
        choppa.full_seq = full_seq
        choppa.taxon_id = None
        choppa.host = None
        choppa.cds_list = []
    def __str__(choppa):
        a       = f"Accession:          {choppa.accession}\n"
        l       = f"Locus:              {choppa.locus}\n"
        v       = f"Version No.:        {choppa.version}\n"
        o       = f"Organism:           {choppa.organism}\n"
        h       = f"Host:               {choppa.host}\n"
        t       = f"Taxon ID:           {choppa.taxon_id}\n"
        cdslist = f"No. of CDS region:  {len(choppa.cds_list)}"
        return(a+l+v+o+h+t+cdslist)
    def addTaxonID(choppa,taxon_id):
        choppa.taxon_id = str(taxon_id)
    def removeTaxonID(choppa):
        choppa.taxon_id = None
    def addHost(choppa,host):
        choppa.host = str(host)
    def removeHost(choppa):
        choppa.host = None
    def appendCDS(choppa,ProteinObject):
        choppa.cds_list.append(ProteinObject)
    def clearCDS(choppa):
        choppa.cds_list = []

#========================== Protein Class ==========================#
# a class that stores protein information

class ProteinClass:
    def __init__(nico,location):
        nico.location = location
        nico.protein = None
        nico.protein_id = None
        nico.protein_seq = None
        # nico.coding_seq = None
    def __str__(nico):
        p       = f"Protein:        {nico.protein}\n"
        pid     = f"Protein ID:     {nico.protein_id}\n"
        cds     = f"CDS info:       {nico.location}\n"
        seq     = f"Protein Seq:    {nico.protein_seq}\n"
        # cds_seq = f"Coding Seq:     {nico.coding_seq}"
        # return(p+pid+cds+seq+cds_seq)
        return(p+pid+cds+seq)
    def addProteinName(nico,protein):
        nico.protein = str(protein)
    def removeProteinName(nico):
        nico.protein = None
    def addProteinSeq(nico,protein_seq):
        nico.protein_seq = str(protein_seq)
    def removeProteinSeq(nico):
        nico.protein_seq = None
    def addProteinID(nico,protein_id):
        nico.protein_id = str(protein_id)
    def removeProteinID(nico):
        nico.protein_id = None
    # def addCodingSeq(nico,coding_seq):
    #     nico.coding_seq = str(coding_seq)
    # def removeCodingSeq(nico):
    #     nico.coding_seq = None
        

#========================== Peptide Comparison Class ==========================#
# a class that stores peptide comparison information

class PeptideCompare:
    def __init__(franky, virus_a, virus_b):
        franky.virus_a = str(virus_a)
        franky.virus_b = str(virus_b)
        franky.intersection = None
        franky.intersection_size = 0
        franky.exclusiveA = None
        franky.exclusiveA_size = 0
        franky.exclusiveB = None
        franky.exclusiveB_size = 0
        franky.total_pep_count = 1
    def __str__(franky):        
        va_name     = f"Virus 1:                    {franky.virus_a}\n"
        va_count    = f"No. of exclusive pepides:   {franky.exclusiveA_size}\n"
        vb_name     = f"Virus 2:                    {franky.virus_b}\n"
        vb_count    = f"No. of exclusive peptides:  {franky.exclusiveB_size}\n"
        inter       = f"No. of shared peptides:     {franky.intersection_size}\n"
        total_pep   = f"No. of peptides in VirLib:  {franky.total_pep_count}\n"
        return(va_name + va_count + vb_name + vb_count + inter + total_pep)
    def add_intersection_set(franky, intersection):
        franky.intersection = intersection
        franky.intersection_size = len(intersection)
    def add_exclusiveA_set(franky, exclusiveA):
        franky.exclusiveA = exclusiveA
        franky.exclusiveA_size = len(exclusiveA)
    def add_exclusiveB_set(franky, exclusiveB):
        franky.exclusiveB = exclusiveB
        franky.exclusiveB_size = len(exclusiveB)
    def update_total_pep_count(franky, pep_count):
        franky.total_pep_count = int(pep_count)
    
        
#========================"My Input File Error Class"=======================#

#Error classes for handling input errors
class InputFileError(Exception):
    """
        Input File Error class to notify when something is wrong
        with the input files given.    
    """
    """Exception for errors in the input.
        Attributes:
        expr -- the input that caused the error
        msg  -- message giving details of the error
    """
    def __init__(nami, expr, msg):
        nami.expr = expr
        nami.msg = msg    
    def __str__(nami):
        return repr(str(nami.msg)+str(nami.expr))
