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
from typing import Dict, List, Set


YEAR_UNKNOWN = 10000
NOT_IN_MYCOBANK = 'Taxon not in MycoBank'
UNDETERMINED_NAMING_PRIORITY = 'undetermined'

NAME_STATUS = set([
    'Legitimate',
    'Invalid',
    'Orthographic variant',
    'Illegitimate',
    'Deleted',
    'Unavailable',
    'Uncertain'
])

VALID_NAME_STATUS = set(['Legitimate', 'Orthographic variant'])


@dataclass
class MycoBankMetadata:
    rank: str
    year_publication: int
    name_status: str
    is_current_name: bool
    current_name: str
    basionym: str
    taxonomic_synonyms: List[str]
    obligate_synonyms: List[str]


def parse_mycobank_synonym_data(synonym_data, rank):
    """Parse synonym data from MycoBank.
    
    *** THIS IS INCOMPLETE - I BELIEVE THIS WAS PREVIOUSLY DONE BY PIERRE ***
    """

    SECTION_TOKEN = '<NEW_SECTION>'
    SYNONYM_SECTIONS = ['Current name:', 'Basionym:', 'Taxonomic synonyms:', 'Obligate synonyms:']

    # mark synonym string with the section token at the start of each section
    for synonym_section in SYNONYM_SECTIONS:
        synonym_data = synonym_data.replace(synonym_section, f'{SECTION_TOKEN}{synonym_section}')

    # split synonym string on section token and then process each sections data accordingly
    current_name = None
    basionym = None
    taxonomic_synonyms = set()
    obligate_synonyms = set()
    for section_data in synonym_data.split(SECTION_TOKEN):
        if section_data.startswith('Current name:'):
            current_name_tokens = [t.strip() for t in section_data.replace('Current name:', '').split()]
            if rank == 'sp.':
                current_name = f"s__{' '.join(current_name_tokens[0:2])}"
            elif rank == 'gen.':
                current_name = f"g__{' '.join(current_name_tokens[0:1])}"
            else:
                print('[Error] Invalid rank.')
                sys.exit(1)
        elif section_data.startswith('Basionym:'):
            pass
        elif section_data.startswith('Taxonomic synonyms:'):
            pass
        elif section_data.startswith('Obligate synonyms:'):
            pass

    return current_name, basionym, taxonomic_synonyms, obligate_synonyms


def parse_mycobank_data(mycobank_file: str) -> Dict[str, MYCOBANK_METADATA]:
    """Parse Mycobank nomenclature information."""

    mycobank_data = {}
    with open(mycobank_file, errors='replace') as f:
        header = f.readline().strip().split('\t')

        taxon_name_idx = header.index('Taxon_name')
        rank_idx = header.index('Rank')
        year_idx = header.index('Year_of_effective_publication')
        name_status_idx = header.index('Name_status')
        synonymy_idx = header.index('Synonymy')

        for line in f:
            tokens = line.strip().split('\t')

            taxon_name = tokens[taxon_name_idx]

            rank = tokens[rank_idx]
            if rank == 'gen.':
                taxon_name = f'g__{taxon_name}'
            elif rank == 'sp.':
                taxon_name = f's__{taxon_name}'
            else:
                continue

            current_name, basionym, taxonomic_synonyms, obligate_synonyms = parse_mycobank_synonym_data(tokens[synonymy_idx], rank)

            year = tokens[year_idx]
            if '/' in year: # some data have the format dd/mm/yyyy
                year = year.split('/')[-1]

            mycobank_data[taxon_name] = MycoBankMetadata(
                rank,
                int(year) if year != '?' else YEAR_UNKNOWN,
                tokens[name_status_idx],
                taxon_name == current_name,
                current_name,
                basionym,
                taxonomic_synonyms,
                obligate_synonyms
            )

    return mycobank_data

def read_sanctioned_names(sanctioned_names_file: str) -> Set[str]:
    """Parse file with sanctioned names."""

    with open(sanctioned_names_file, errors='replace') as f:
        header = f.readline().strip().split('\t')

        taxon_idx = header.index('Taxon name')

        for line in f:
            tokens = line.strip().split('\t')

            taxon = tokens[taxon_idx]
            if len(taxon.split()) == 2:
                taxon = f's__{taxon}'
            else:
                taxon = f'g__{taxon}'

            sanctioned_names.add(taxon)