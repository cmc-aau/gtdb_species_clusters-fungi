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

import re
import gzip
import logging
from dataclasses import dataclass
from collections import defaultdict
from typing import Dict, List, Set, Optional
from unidecode import unidecode


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
    current_name: Optional[str]
    basionym: Optional[str]
    taxonomic_synonyms: List[str]
    obligate_synonyms: List[str]


def strip_nonalnum(text: str) -> str:
    """Remove leading and trailing non-alphanumeric characters."""
    return re.sub(r"^\W+|\W+$", "", text)


def extract_binomial_name(text, rank, taxon_name):
    """
    Extract the appropriate taxonomic name based on rank.

    Args:
        text: The text containing the taxonomic name
        rank: The taxonomic rank (e.g., 'sp.', 'var.')
        taxon_name: The original taxon name

    Returns:
        Cleaned binomial or generic name
    """
    words = text.strip().split()

    # For species or variety ranks, extract binomial names
    if rank == 'sp.' or rank == 'var.' or len(taxon_name.strip().split()) == 2:
        # Include variety designation if present
        if len(words) > 2 and words[2] == 'var.':
            name = ' '.join(words[:4])
        else:
            name = ' '.join(words[:2])

        name = strip_nonalnum(name)

        # check if name is valid by making sure it consists of 
        # only alphabetical characters, spaces, periods, or hyphens
        #  Examples:
        #    Lecidea alboatra var. margaritacea
        #    Diacanthodes novo-guineensis
        if not all(c.isalpha() or c.isspace() or c in ['.', '-'] for c in name):
            name = None
    else:
        # For higher ranks, extract only genus name
        name = words[0] if words else ''
        name = strip_nonalnum(name)

    return name


def parse_mycobank_synonym_data(synonym_data: str, rank: str, taxon_name: str) -> Dict[str, List[str]]:
    """
    Parse synonym information from GTDB synonymy string.

    Args:
        synonymy: String containing synonym information
        rank: Taxonomic rank
        taxon_name: Name of the taxon

    Returns:
        Dictionary mapping synonym types to lists of names
    """

    synonym_data = unidecode(synonym_data)

    result = defaultdict(list)

    # Define fields to search for
    fields = [
        ('Current name:', 'current_name'),
        ('Basionym:', 'basionym'),
        ('Taxonomic synonyms:', 'taxonomic_synonyms'),
        ('Taxonomic synonym:', 'taxonomic_synonyms'),
        ('Obligate synonyms:', 'obligate_synonyms'),
        ('Obligate synonym:', 'obligate_synonyms')
    ]

    # Find all field positions
    field_positions = {}
    for field_text, field_name in fields:
        index = synonym_data.find(field_text)
        if index != -1:
            field_positions[index] = {
                'name': field_name,
                'length': len(field_text)
            }

    if not field_positions:
        return result

    # Process each field in order
    sorted_positions = sorted(field_positions.keys())

    for idx, pos in enumerate(sorted_positions):
        # Extract text between this field and the next (or end of string)
        start = pos + field_positions[pos]['length']
        end = sorted_positions[idx + 1] if idx < len(sorted_positions) - 1 else len(synonym_data)
        text_segment = synonym_data[start:end]

        # Split by MycoBank reference pattern [MB#######]
        name_segments = re.split(r'\[MB#\d+\]', text_segment)

        field_name = field_positions[pos]['name']

        for segment in name_segments:
            segment = strip_nonalnum(segment)
            if segment:
                name = extract_binomial_name(segment, rank, taxon_name)
                if name:
                    result[field_name].append(name)

    return result


def parse_mycobank_data(mycobank_file: str) -> Dict[str, MycoBankMetadata]:
    """Parse Mycobank nomenclature information."""

    mycobank_data = {}
    with gzip.open(mycobank_file, 'rt', errors='replace') as f:
        header = f.readline().strip().split('\t')

        taxon_name_idx = header.index('Taxon name')
        rank_idx = header.index('Rank.Rank name')
        year_idx = header.index('Year of effective publication')
        name_status_idx = header.index('Name status')
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

            synonym_data = parse_mycobank_synonym_data(tokens[synonymy_idx], rank, taxon_name)

            if len(synonym_data['current_name']) == 0:
                current_name = None
            elif len(synonym_data['current_name']) == 1:
                current_name = synonym_data['current_name'][0]
            else:
                self.log = logging.getLogger('rich')
                self.log.error(f'Mycobank data indicated multiple current names: {synonym_data['current_name']}')
                sys.exit(1)

            basionym = None
            if len(synonym_data['basionym']) == 0:
                pass
            elif len(synonym_data['basionym']) == 1:
                basionym = synonym_data['basionym'][0]
            else:
                self.log = logging.getLogger('rich')
                self.log.error(f'Mycobank data indicated multiple basionyms: {synonym_data['basionym']}')
                sys.exit(1)
            
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
                synonym_data['taxonomic_synonyms'],
                synonym_data['obligate_synonyms']
            )

    return mycobank_data


def read_sanctioned_names(sanctioned_names_file: str) -> Set[str]:
    """Parse file with sanctioned names."""

    sanctioned_names = set()
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

    return sanctioned_names