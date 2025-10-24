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
import gzip
import logging
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Set, Tuple

from gtdblib.util.shell.execute import check_dependencies

from gsc_fungi.genome import GenomeMetadata
from gsc_fungi.genome_utils import read_genome_path
from gsc_fungi.ncbi import parse_gid_to_ncbi_sp, MetadataNCBI, NCBI_TYPE_SPECIES, NCBI_PROXYTYPE, NCBI_TYPE_SUBSP, NCBI_EXCLUDED_AS_TYPE, REFSEQ_CATEGORIES
from gsc_fungi.skani import skani as Skani


# Species with multiple NCBI representative / reference genomes
# were manually inspected. In most cases, this was found to
# be the results of genomes being representatives of hybrids.
NCBI_REP_GENOMES_TO_EXCLUDE = {
    'GCA_006992865.1', # Cryptococcus neoformans AD hybrid
    'GCA_019381815.1', # Ogataea polymorpha x Ogataea parapolymorpha
    'GCA_009666385.1', # Saccharomyces eubayanus x Saccharomyces uvarum
    'GCA_001936135.1', # Saccharomycopsis fibuligera x Saccharomycopsis cf. fibuligera
    'GCA_017656815.1', # Fusarium meridionale x Fusarium asiaticum
    'GCA_020280125.1', # Kurtzmaniella quercitrusa var. comoensis (nom. nud.)
    'GCA_020280135.1', # Kurtzmaniella quercitrusa var. filamentosus (nom. nud.)
    'GCA_030578135.1', # Saccharomyces cerevisiae x Saccharomyces kudriavzevii
    'GCA_013180185.1', # Saccharomyces cerevisiae x Saccharomyces uvarum
    'GCA_013180725.1', # Saccharomyces cerevisiae x Saccharomyces kudriavzevii
    'GCA_009666275.1', # Saccharomyces cerevisiae x Saccharomyces eubayanus
    'GCA_009665555.1', # Saccharomyces cerevisiae x Saccharomyces eubayanus
    'GCA_009666465.1', # Saccharomyces cerevisiae x Saccharomyces kudriavzevii
    'GCA_009667055.1', # Saccharomyces cerevisiae x Saccharomyces eubayanus
}


@dataclass
class RepresentativeSelectionStats:
    sp_has_type_or_rep_genome: int = 0
    sp_type_only: int = 0
    sp_rep_only: int = 0
    sp_no_type_genome: int = 0

    multi_sp_type_genomes: int = 0
    multi_sp_type_ncbi_rep_resolved: int = 0
    multi_sp_type_similar_genomes: int = 0

    multi_subsp_type_genomes: int = 0
    multi_subsp_type_ncbi_rep_resolved: int = 0
    multi_subsp_type_similar_genomes: int = 0

    ani99_resolved: int = 0

    manual_resolution: int = 0
    single_sp_type_genome: int = 0
    single_subsp_type_genome: int = 0
    sp_no_type_single: int = 0


class SelectSpeciesRepresentative():
    """Select representative genome for each NCBI fungal species."""

    def __init__(self, out_dir: str):
        """Initialization.

        :param out_dir: Output directory.
        """

        self.log = logging.getLogger('rich')
        self.out_dir = out_dir

        check_dependencies(['skani'])
        self.log.info('Using skani v{}.'.format(Skani.get_version()))

    def check_pairwise_ani_99(self, gids: Set[str], genome_files: Dict[str, str]) -> Tuple[bool, Dict[str, List[float]]]:
        """Check if all genomes have an ANI >= 99%"""

        similar_genomes = True
        anis = defaultdict(list)
        for gid1, gid2 in combinations(gids, 2):
            result = Skani.dist(gid1, gid2, genome_files[gid1], genome_files[gid2])

            anis[gid1].append(result.ani)
            anis[gid2].append(result.ani)

            if result.ani < 99.0:
                similar_genomes = False
            
        return similar_genomes, anis

    def metadata_pass_qc(self, qc_pass_file: str) -> Dict[str, GenomeMetadata]:
        """Get metadata for genomes passing QC."""

        genome_metadata = {}
        with open(qc_pass_file) as f:
            header = f.readline().strip().split('\t')

            assem_score_idx = header.index('assembly_score')
            comp_idx = header.index('completeness')
            cont_idx = header.index('contamination')
            qual_idx = header.index('quality')
            contig_count_idx = header.index('contig_count')
            n50_contigs_idx = header.index('n50_contigs')
            ambig_perc_idx = header.index('ambiguous_base_perc')
            ncbi_type_material_idx = header.index('ncbi_relation_to_type_material')
            ncbi_refseq_cat_idx = header.index('ncbi_refseq_category')
            ncbi_exclude_idx = header.index('ncbi_excluded_from_refseq')
            ncbi_assem_level_idx = header.index('ncbi_assem_level')

            for line in f:
                tokens = line.strip().split('\t')

                gid = tokens[0]
                
                genome_metadata[gid] = GenomeMetadata(
                    float(tokens[assem_score_idx]),
                    float(tokens[comp_idx]),
                    float(tokens[cont_idx]),
                    float(tokens[qual_idx]),
                    int(tokens[contig_count_idx]),
                    int(tokens[n50_contigs_idx]),
                    float(tokens[ambig_perc_idx]),
                    MetadataNCBI(
                        tokens[ncbi_type_material_idx],
                        tokens[ncbi_refseq_cat_idx],
                        [t.strip() for t in tokens[ncbi_exclude_idx].split(';')],
                        tokens[ncbi_assem_level_idx]
                    )
                )

        return genome_metadata

    def determine_type_material(self,
                                ncbi_assembly_file: str, 
                                genome_metadata: Dict[str, GenomeMetadata], 
                                gid_to_ncbi_sp: Dict[str, str]) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], Dict[str, List[str]]]:
        """Determine genomes that are type strain of species, type strain of subspecies, or NCBI species representative."""

        sp_type_genomes = defaultdict(list)
        subsp_type_genomes = defaultdict(list)
        sp_rep_genomes = defaultdict(list)
        gid_to_type_material = {}
        gid_to_rs_category = {}

        with open(ncbi_assembly_file) as f:
            f.readline()
            header = f.readline().strip().split('\t')

            org_name_idx = header.index('organism_name')

            for line in f:
                tokens = line.strip().split('\t')

                gid = tokens[0]
                if gid not in genome_metadata:
                    # only process genomes passing QC
                    continue

                sp = gid_to_ncbi_sp[gid]
                if sp == 's__':
                    # can't be type material of an unspecified species
                    continue

                type_material = genome_metadata[gid].ncbi.relation_to_type_material
                rs_category = genome_metadata[gid].ncbi.refseq_category

                gid_to_type_material[gid] = type_material
                gid_to_rs_category[gid] = rs_category

                if rs_category in REFSEQ_CATEGORIES and gid not in NCBI_REP_GENOMES_TO_EXCLUDE:
                    sp_rep_genomes[sp].append(gid)

                for exclude in genome_metadata[gid].ncbi.excluded_from_refseq:
                    if exclude in NCBI_EXCLUDED_AS_TYPE:
                        continue

                if type_material in NCBI_TYPE_SPECIES:
                    org_name = tokens[org_name_idx]
                    name_tokens = org_name.split()
                    if len(name_tokens) == 4 and name_tokens[2] in ['var.', 'f.', 'pv.', 'subsp.']:
                        if name_tokens[1] == name_tokens[3]:
                            # e.g. Schwanniomyces occidentalis var. occidentalis
                            sp_type_genomes[sp].append(gid)
                        else:
                            # e.g. Schwanniomyces occidentalis var. persoonii
                            subsp_type_genomes[sp].append(gid)
                    else:
                        sp_type_genomes[sp].append(gid)
                elif type_material in NCBI_TYPE_SUBSP:
                    subsp_type_genomes[sp].append(gid)
                    
        return sp_type_genomes, subsp_type_genomes, sp_rep_genomes

    def determine_best_assembly(gids: List[str], genome_metadata: Dict[str, GenomeMetadata]) -> str:
        """Determine genome with highest assembly score."""
        
        gtdb_rid = None
        highest_score = -1e9
        for gid in gids:
            if metadata[gid].assembly_score > highest_score:
                highest_score = metadata[gid].assembly_score
                gtdb_rid = gid

        return gtdb_rid

    def sel_sp_rep_from_type(self, 
                                genome_metadata: Dict[str, GenomeMetadata], 
                                sp_type_genomes: Dict[str, List[str]],
                                subsp_type_genomes: Dict[str, List[str]],
                                sp_rep_genomes: Dict[str, List[str]],
                                genome_files: Dict[str, str],
                                stats: RepresentativeSelectionStats) -> str:
        """Select representative genome for species with type strain or NCBI representative genomes."""

        gtdb_rid = None

        stats.sp_has_type_or_rep_genome += 1

        if sp not in sp_rep_genomes:
            stats.sp_type_only += 1

        if len(sp_type_genomes[sp]) == 1:
            # select the only type strain of species genome as the GTDB representative
            stats.single_sp_type_genome += 1
            gtdb_rid = sp_type_genomes[sp][0]
        elif len(sp_type_genomes[sp]) > 1:
            stats.multi_sp_type_genomes += 1

            rep_inter = set(sp_type_genomes.get(sp, [])).intersection(sp_rep_genomes.get(sp, []))
            if len(rep_inter) > 0:
                # select the type strain genome that is also the NCBI representative as the GTDB representative
                assert len(rep_inter) == 1
                gtdb_rid = list(rep_inter)[0]
                stats.multi_sp_type_ncbi_rep_resolved += 1
            else:
                # this is the most complicated case since we have multiple genomes that are assembled
                # from the type strain of the species and none of these genomes are the NCBI representative genome
                self.log.info(f' - {sp} has {len(sp_type_genomes[sp]):,} type species genomes and cannot be resolved with NCBI representative genome information')

                # determine ANI between genomes and select highest quality genome assembly
                # if all genomes are highly similar to each other
                similar_genomes, anis = self.check_pairwise_ani_99(sp_type_genomes[sp], genome_files)
                if similar_genomes:
                    stats.multi_sp_type_similar_genomes += 1
                    gtdb_rid = self.determine_best_assembly(sp_type_genomes[sp], genome_metadata)

        elif len(subsp_type_genomes[sp]) == 1:
            # select the only type strain of subspecies genome as the GTDB representative
            stats.single_subsp_type_genome += 1
            gtdb_rid = subsp_type_genomes[sp][0]
        elif len(subsp_type_genomes[sp]) > 1:
            stats.multi_subsp_type_genomes += 1

            rep_inter = set(subsp_type_genomes.get(sp, [])).intersection(sp_rep_genomes.get(sp, []))
            if len(rep_inter) > 0:
                # select the type strain of subspecies genome that is also the NCBI representative as the GTDB representative
                assert len(rep_inter) == 1
                gtdb_rid = list(rep_inter)[0]
                stats.multi_subsp_type_ncbi_rep_resolved += 1
            else:
                # this is the most complicated case since we have multiple genomes that are assembled
                # from the type strain of the subspecies and none of these genomes are the NCBI representative genome
                self.log.info(f' - {sp} has {len(subsp_type_genomes[sp]):,} type subspecies genomes and cannot be resolved with NCBI representative genome information')

                # determine ANI between genomes and select highest quality genome assembly
                # if all genomes are highly similar to each other
                similar_genomes, anis = check_pairwise_ani_99(subsp_type_genomes[sp], genome_files)
                if similar_genomes:
                    stats.multi_subsp_type_similar_genomes += 1
                    gtdb_rid = self.determine_best_assembly(subsp_type_genomes[sp], genome_metadata)
        elif len(sp_rep_genomes[sp]) >= 1:
            # select the NCBI representative genome as the GTDB representative
            assert len(sp_rep_genomes[sp]) == 1
            stats.sp_rep_only += 1
            gtdb_rid = sp_rep_genomes[sp][0]
        else:
            self.log.error('This should never occur as the above cases should be exhaustive.')
            sys.exit(1)

        return gtdb_rid

    def select_sp_representative_genome(self,
                                        genome_metadata: Dict[str, GenomeMetadata], 
                                        sp_type_genomes: Dict[str, List[str]],
                                        subsp_type_genomes: Dict[str, List[str]],
                                        sp_rep_genomes: Dict[str, List[str]],
                                        genome_files: Dict[str, str]):
        """Select representative genome for each NCBI species."""

        sp_rep_out_file = os.path.join(self.out_dir, 'sp_reps.tsv')
        fout = open(sp_rep_out_file, 'w')
        fout.write('ncbi_sp\tgid\tis_gtdb_rep\tncbi_is_type_strain_of_species\tncbi_relation_to_type_material\tncbi_refseq_category\tassembly_score\tncbi_excluded_from_refseq\n')

        manual_rep_out_file = os.path.join(self.out_dir, 'manual_reps.tsv')
        fout_man = open(manual_rep_out_file, 'w')
        fout_man.write('ncbi_sp\tgid\tncbi_is_type_strain_of_species\tis_highest_scoring\tncbi_relation_to_type_material\tncbi_refseq_category\tassembly_score\tintraspecific_ani\tncbi_excluded_from_refseq\n')

        stats = RepresentativeSelectionStats()

        for sp, gids in ncbi_sp_to_gids.items():
            if sp == 's__':
                continue

            if sp in sp_type_genomes or sp in subsp_type_genomes or sp in sp_rep_genomes:
                gtdb_rid = self.sel_sp_rep_from_type(genome_metadata, 
                                                        sp_type_genomes,
                                                        subsp_type_genomes,
                                                        sp_rep_genomes,
                                                        genome_files,
                                                        stats)

                union_gids = set(sp_type_genomes.get(sp, [])).union(subsp_type_genomes.get(sp, [])).union(sp_rep_genomes.get(sp, []))
                if gtdb_rid is not None:
                    union_gids.add(gtdb_rid)

                    for gid in union_gids:
                        fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\n'.format(
                            sp,
                            gid,
                            gid == gtdb_rid,
                            gid in sp_type_genomes[sp],
                            gid_to_type_material[gid],
                            gid_to_rs_category[gid],
                            genome_metadata[gid].assembly_score,
                            '; '.join(genome_metadata[gid].ncbi.excluded_from_refseq),
                        ))
                else:
                    # manual resolution of GTDB representative is required
                    stats.manual_resolution += 1

                    gtdb_rid = self.determine_best_assembly(union_gids, genome_metadata)

                    for gid in union_gids:
                        fout_man.write('{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\n'.format(
                            sp,
                            gid,
                            gid == gtdb_rid,
                            gid in sp_type_genomes[sp],
                            gid_to_type_material[gid],
                            gid_to_rs_category[gid],
                            metadata[gid].assembly_score,
                            dict(anis),
                            '; '.join(metadata[gid].ncbi.excluded_from_refseq),
                        ))

                        fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\n'.format(
                            sp,
                            gid,
                            gid == gtdb_rid,
                            gid in sp_type_genomes[sp],
                            gid_to_type_material[gid],
                            gid_to_rs_category[gid],
                            metadata[gid].assembly_score,
                            '; '.join(metadata[gid].ncbi.excluded_from_refseq),
                        ))
            else:
                # select genome with highest quality assembly to be the
                # representative genome as there is no type material or
                # NCBI representative genomes for this species
                if len(gids) == 1:
                    sp_no_type_single += 1

                sp_no_type_genome += 1

                gtdb_rid = self.determine_best_assembly(gids, genome_metadata)

                for gid in gids:
                    fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\n'.format(
                        sp,
                        gid,
                        gid == gtdb_rid,
                        gid in sp_type_genomes[sp],
                        'na',
                        'na',
                        metadata[gid].assembly_score,
                        '; '.join(metadata[gid].ncbi.excluded_from_refseq),
                    ))

        fout.close()
        fout_man.close()

        return stats

    def run(self,
            qc_pass_file: str, 
            ncbi_assembly_file: str, 
            ncbi_taxonomy_file: str,
            genome_path_file: str) -> None:
        """Select representative genome for each NCBI fungal species."""

        # read path to genomic FASTA files
        self.log.info('Reading path to genomic FASTA files:')
        genome_files = read_genome_path(genome_path_file)
        self.log.info(f' - read path for {len(genome_files):,} genomes')

        # get metadata for genomes passing QC criteria
        self.log.info('Reading metadata for genomes passing QC:')
        genome_metadata = metadata_pass_qc(qc_pass_file)
        self.log.info(f' - read metadata for {len(metadata):,} genomes')

        # get fungal species under NCBI classification
        self.log.info('Parsing NCBI species assignments:')
        gid_to_ncbi_sp = parse_gid_to_ncbi_sp(ncbi_taxonomy_file)

        ncbi_sp_to_gids = defaultdict(set)
        for gid, sp in gid_to_sp.items():
            ncbi_sp_to_gids[sp].add(gid)

        self.log.info(f' - identified {len(ncbi_sp_to_gids):,} NCBI species across {len(gid_to_ncbi_sp):,} genomes')

        # determine genomes assembled from type material for named NCBI species
        self.log.info('Determining genomes assembled from type material for NCBI species:')
        sp_type_genomes, subsp_type_genomes, sp_rep_genomes = self.determine_type_material(ncbi_assembly_file, 
                                                                                            genome_metadata, 
                                                                                            gid_to_ncbi_sp)
        self.log.info(f' - identified {len(sp_type_genomes):,} species with one or more type genomes')
        self.log.info(f' - identified {len(sp_rep_genomes):,} species with one or more NCBI representative genomes')

        subsp_only_type_genomes = set(subsp_type_genomes) - set(sp_type_genomes)
        self.log.info(f' - identified {len(subsp_only_type_genomes):,} species with only subspecies type genomes')

        self.log.info('The following species have multiple type strain of species genomes:')
        for sp, gids in sp_type_genomes.items():
            if len(gids) > 1:
                self.log.info(f' - {sp}: {gids}')

        self.log.info('The following species have multiple NCBI representative genomes:')
        for sp, gids in sp_rep_genomes.items():
            if len(gids) > 1:
                self.log.info(f' - {sp}: {gids}')

        assert len(set(sp_type_genomes) - set(ncbi_sp_to_gids)) == 0

        # select and write out type material information
        self.log.info('Selecting representative genome for each NCBI species:')
        selection_stats = self.select_sp_representative_genome(genome_metadata,
                                                                ncbi_sp_to_gids, 
                                                                sp_type_genomes, 
                                                                subsp_type_genomes, 
                                                                sp_rep_genomes,
                                                                genome_files)

        self.log.info(' - sp_has_type_or_rep_genome', selection_stats.sp_has_type_or_rep_genome)
        self.log.info('   - sp_type_only', selection_stats.sp_type_only)
        self.log.info('   - sp_rep_only', selection_stats.sp_rep_only)
        self.log.info('   - single_sp_type_genome', selection_stats.single_sp_type_genome)
        self.log.info('   - single_subsp_type_genome', selection_stats.single_subsp_type_genome)
        self.log.info('   - multi_sp_type_genomes', selection_stats.multi_sp_type_genomes)
        self.log.info('     - multi_sp_type_similar_genomes', selection_stats.multi_sp_type_similar_genomes)
        self.log.info('     - multi_sp_type_ncbi_rep_resolved', selection_stats.multi_sp_type_ncbi_rep_resolved)
        self.log.info('   - multi_subsp_type_genomes', selection_stats.multi_subsp_type_genomes)
        self.log.info('     - multi_subsp_type_similar_genomes', selection_stats.multi_subsp_type_similar_genomes)
        self.log.info('     - multi_subsp_type_ncbi_rep_resolved', selection_stats.multi_subsp_type_ncbi_rep_resolved)
        self.log.info(' - manual_resolution', selection_stats.manual_resolution)
        self.log.info(' - sp_no_type_genome', selection_stats.sp_no_type_genome)
        self.log.info('   - sp_no_type_single', selection_stats.sp_no_type_single)
                        