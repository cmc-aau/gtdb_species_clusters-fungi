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

from gsc_fungi.ncbi import MetadataNCBI, NCBI_COMPLETE_GENOME, NCBI_ENV_GENOME


@dataclass
class GenomeMetadata:
    """Metadata for a genome including quality metrics."""
    phylo_marker_perc: float
    lineage_dup_marker_perc: float
    genome_size: int
    contig_count: int
    n50_contigs: int
    ambiguous_base_perc: float
    assembly_score: float
    ncbi: MetadataNCBI


def assembly_score(ncbi_assem_level: str,
                    excluded_from_refseq: List[str],
                    phylo_marker_perc: float,
                    lineage_dup_marker_perc: float,
                    contig_count: int,
                    ambiguous_base_perc: float) -> float:
        """Calculate assembly score for genome."""

        score = 0
        if ncbi_assem_level == NCBI_COMPLETE_GENOME:
            score += 100

        score += phylo_marker_perc - 5*lineage_dup_marker_perc

        for exclude in excluded_from_refseq:
            if exclude in NCBI_ENV_GENOME:
                score -= 100

        score -= 5 * contig_count/1000.0
        score -= 5 * ambiguous_base_perc

        return score