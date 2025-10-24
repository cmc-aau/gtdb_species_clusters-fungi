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
from typine import List

from gsc_fungi.ncbi import MetadataNCBI


@dataclass
class GenomeMetadata:
    """Metadata for a genome including quality metrics."""
    completness: float
    contamination: float
    genome_quality: float
    genome_size: int
    contig_count: int
    n50_contigs: int
    assemblye_score: float
    ncbi: MetadataNCBI


def assembly_score(ncbi_assem_level: str,
                    excluded_from_refseq: List[str], 
                    completeness: float, 
                    contamination: float,
                    contig_count: int,
                    ambiguous_base_perc: float) -> float:
        """Calculate assembly score for genome."""

        score = 0
        if ncbi_assem_level == NCBI_COMPLETE_GENOME:
            score += 100

        score += completeness - 5*contamination

        for exclude in excluded_from_refseq:
            if exclude in NCBI_ENV_GENOME:
                score -= 100

        score -= 5 * contig_count/100.0
        score -= 5 * ambiguous_base_perc

        return score