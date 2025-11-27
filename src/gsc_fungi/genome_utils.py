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

import os
from collections import namedtuple
from typing import Dict, Set, Optional

from gtdblib.util.bio.accession import canonical_gid


def read_genome_path(genome_path_file: str, gids_of_interest: Optional[Set[str]] = None) -> Dict[str, str]:
    """Determine path to genomic FASTA file for each genome."""

    genome_files = {}
    for line in open(genome_path_file):
        line_split = line.strip().split('\t')

        gid = line_split[0]
        if gids_of_interest is not None and gid not in gids_of_interest:
            continue

        genome_path = line_split[1]
        accession = os.path.basename(os.path.normpath(genome_path))

        genome_files[gid] = os.path.join(
            genome_path, accession + '_genomic.fna.gz')

    return genome_files