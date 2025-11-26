###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

from dataclasses import dataclass
from typing import List


@dataclass
class MetadataNCBI:
    """NCBI metadata for a genome."""
    relation_to_type_material: str
    refseq_category: str
    excluded_from_refseq: List[str]
    assem_level: str


# Criteria for excluding or retaining genomes based on NCBI metadata
NCBI_EXCLUSION_FILTERING_CRITERIA = set(['partial'])
NCBI_COMPLETE_GENOME = 'Complete Genome'
NCBI_ENV_GENOME = set(['derived from metagenome', 'derived from single cell'])

# Criteria for establishing species representative genomes 
NCBI_TYPE_SPECIES = set(['assembly from type material',
                             'assembly from neotype material',
                             'assembly designated as neotype',
                             'assembly designated as reftype'])
NCBI_PROXYTYPE = set(['assembly from proxytype material'])
NCBI_TYPE_SUBSP = set(['assembly from synonym type material', 'assembly from heterotypic synonym type material'])

NCBI_EXCLUDED_AS_TYPE = set([
    'not used as type'
])

REFSEQ_CATEGORIES = set([
    'representative genome',
    'reference genome',
])


def parse_gid_to_ncbi_sp(ncbi_taxonomy_file: str) -> Dict[str, str]:
    """Parse NCBI taxonomy file to determine species classification of each genome."""

    gid_to_sp = {}
    with open(ncbi_taxonomy_file) as f:
        for line in f:
            tokens = line.strip().split('\t')

            gid = tokens[0]
            if gid.startswith('GCF_'):
                # GTDB Fungi DB is build strictly from genomes in GenBank
                continue

            sp = tokens[1].split(';')[-1]
            assert sp.startswith('s__')
            gid_to_sp[gid] = sp

    return gid_to_sp