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
import logging
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import Dict, Set, Tuple

from gtdblib.util.bio.accession import canonical_gid

from gsc_fungi.genome import GenomeMetadata, assembly_quality, assembly_score
from gsc_fungi.ncbi import parse_gid_to_ncbi_sp, MetadataNCBI, NCBI_EXCLUSION_FILTERING_CRITERIA, NCBI_COMPLETE_GENOME, NCBI_ENV_GENOME


@dataclass
class QcCriteria:
    """Criteria for passing QC.

    :param min_comp: Minimum estimated genome completeness to pass QC.
    :param max_cont: Maximum estimated genome contamination to pass QC.
    :param min_quality: Minimum estimated genome quality to pass QC, defined as completeness - 5*contamination.
    :param max_contigs: Maximum number of contigs to pass QC.
    :param min_N50: Minimum N50 to pass QC.
    :param max_ambiguous_perc: Maximum percentage of ambiguous bases to pass QC.
    """

    min_comp: float
    max_cont: float
    min_quality: float
    max_contigs: int
    min_N50: int
    max_ambiguous_perc: float


@dataclass
class BuscoStats:
    comp: float
    cont: float
    quality: float


class QcGenomes():
    """Quality check all potential GTDB genomes."""

    def __init__(self, out_dir: str):
        """Initialization.

        :param out_dir: Output directory.
        """

        self.log = logging.getLogger('rich')
        self.out_dir = out_dir

    def parse_ncbi_assembly_summary_file(self, ncbi_assembly_file: str) -> Tuple[Dict[str, MetadataNCBI], Dict[str, int]]:
        """Parse NCBI assembly summary file."""

        metadata_ncbi = {}
        included_gids = {}
        excluded_from_refseq_count = defaultdict(int)

        with open(ncbi_assembly_file) as f:
            f.readline()
            header = f.readline().strip().split('\t')

            type_material_idx = header.index('relation_to_type_material')
            refseq_category_idx = header.index('refseq_category')
            excluded_from_refseq_idx = header.index('excluded_from_refseq')
            assem_level_idx = header.index('assembly_level')

            for line in f:
                tokens = line.strip().split('\t')

                gid = tokens[0]

                cur_canonical_gid = canonical_gid(gid)
                if cur_canonical_gid in included_gids:
                    prev_gid = included_gids[cur_canonical_gid]
                    self.log.warning(f'Skipping {gid} as {prev_gid} already included.')
                    continue

                excluded_from_refseq = [t.strip() for t in tokens[excluded_from_refseq_idx].split(';') if t != 'na']
                for reason in excluded_from_refseq:
                    excluded_from_refseq_count[reason] += 1

                included_gids[cur_canonical_gid] = gid

                type_material = tokens[type_material_idx]
                rs_category = tokens[refseq_category_idx]
                assem_level = tokens[assem_level_idx]
                       
                metadata_ncbi[gid] = MetadataNCBI(type_material,
                                                    rs_category,
                                                    excluded_from_refseq,
                                                    assem_level)

        return metadata_ncbi, excluded_from_refseq_count

    def parse_busco_results(self, busco_file: str) -> Dict[str, BuscoStats]:
        """Parse Busco results and determine genome quality statistics."""

        busco_stats = {}
        with open(busco_file) as f:
            f.readline()

            for line in f:
                tokens = line.strip().split('\t')

                gid = tokens[0]

                marker_results = Counter(tokens[1:])

                num_complete = marker_results['C']
                num_missing = marker_results['M']
                num_duplicate = marker_results['D']
                num_fragmented = marker_results['F']

                total_markers = sum(marker_results.values())
                assert total_markers == num_complete + num_missing + num_duplicate + num_fragmented

                comp = 100 * (num_complete + num_duplicate) / total_markers
                cont = 100 * num_duplicate / total_markers

                busco_stats[gid] = BuscoStats(
                    comp,
                    cont,
                    assembly_quality(comp, cont)
                )

        return busco_stats

    def parse_genome_metadata(self, 
                                genome_metadata_file: str, 
                                busco_results: Dict[str, BuscoStats], 
                                metadata_ncbi: Dict[str, MetadataNCBI]):
        """Read metadata for genome assemblies and combined with BUSCO and NCBI metadata."""

        genome_metadata = {}
        skip_count = 0
        with open(genome_metadata_file) as f:
            header = f.readline().strip().split('\t')

            genome_idx = header.index('accession')
            genome_size_idx = header.index('genome_size')
            contig_count_idx = header.index('contig_count')
            n50_idx = header.index('n50')
            ambiguous_bases_idx = header.index('ambiguous_bases')
                
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = line_split[genome_idx]
                if gid not in metadata_ncbi or gid not in busco_results:
                    skip_count += 1
                    continue

                genome_size = int(line_split[genome_size_idx])
                contig_count = int(line_split[contig_count_idx])
                n50 = int(line_split[n50_idx])
                ambiguous_base_perc = (100.0 * int(line_split[ambiguous_bases_idx])) / genome_size

                assem_score = assembly_score(metadata_ncbi[gid].assem_level,
                                                 metadata_ncbi[gid].excluded_from_refseq,
                                                 busco_results[gid].comp,
                                                 busco_results[gid].cont,
                                                 contig_count,
                                                 ambiguous_base_perc)
                
                genome_metadata[gid] = GenomeMetadata(busco_results[gid].comp,
                                                        busco_results[gid].cont,
                                                        busco_results[gid].quality,
                                                        genome_size,
                                                        contig_count,
                                                        n50,
                                                        ambiguous_base_perc,
                                                        assem_score,
                                                        metadata_ncbi[gid])

        if skip_count > 0:
            self.log.warning(f'Skipped {skip_count:,} genomes without NCBI or BUSCO metadata.')

        return genome_metadata

    def qc_genomes(self, 
                    genome_metadata: Dict[str, GenomeMetadata], 
                    gid_to_sp: Dict[str, str],
                    qc_criteria: QcCriteria) -> Tuple[Set[str], Set[str], Dict[str, int]]:
        """Determine genomes passing and failing QC."""

        fout_passed = open(os.path.join(self.out_dir, 'qc_passed.tsv'), 'w')
        fout_failed = open(os.path.join(self.out_dir, 'qc_failed.tsv'), 'w')

        header = 'gid\tncbi_species'
        header += '\tassembly_score\tcompleteness\tcontamination\tquality'
        header += '\tgenome_size\tcontig_count\tn50_contigs\tambiguous_base_perc'
        header += '\tncbi_relation_to_type_material\tncbi_refseq_category'
        header += '\tncbi_excluded_from_refseq\tncbi_assem_level'

        fout_passed.write(header + '\n')
        fout_failed.write(header)
        fout_failed.write(
            '\tfailed_completeness\tfailed_contamination\tfailed_quality')
        fout_failed.write(
            '\tfailed_contig_count\tfailed_contig_n50\tfailed_ambiguous_bases_perc\tfailed_refseq_exclude\n')

        qc_pass = set()
        qc_fail = set()
        failed_gids = defaultdict(set)
        for gid, m in genome_metadata.items():
            # check if genome should be filtered based on the RefSeq excluded assignments at NCBI
            refseq_exclusion_tokens = m.ncbi.excluded_from_refseq
            refseq_exclude = False
            for t in refseq_exclusion_tokens:
                if t in NCBI_EXCLUSION_FILTERING_CRITERIA:
                    refseq_exclude = True

            if (m.completness >= qc_criteria.min_comp
                and m.contamination <= qc_criteria.max_cont
                and m.genome_quality >= qc_criteria.min_quality
                and m.contig_count <= qc_criteria.max_contigs
                and m.n50_contigs >= qc_criteria.min_N50
                and m.ambiguous_base_perc <= qc_criteria.max_ambiguous_perc
                and not refseq_exclude):
                qc_pass.add(gid)

                fout_passed.write('{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    gid,
                    gid_to_sp[gid],
                    m.assembly_score,
                    m.completness,
                    m.contamination,
                    m.genome_quality,
                    m.genome_size,
                    m.contig_count,
                    m.n50_contigs,
                    m.ambiguous_base_perc,
                    m.ncbi.relation_to_type_material,
                    m.ncbi.refseq_category,
                    m.ncbi.excluded_from_refseq,
                    m.ncbi.assem_level,
                ))
            else:
                qc_fail.add(gid)

                fout_failed.write('{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    gid,
                    gid_to_sp[gid],
                    m.assembly_score,
                    m.completness,
                    m.contamination,
                    m.genome_quality,
                    m.genome_size,
                    m.contig_count,
                    m.n50_contigs,
                    m.ambiguous_base_perc,
                    m.ncbi.relation_to_type_material,
                    m.ncbi.refseq_category,
                    m.ncbi.excluded_from_refseq,
                    m.ncbi.assem_level,
                ))

                if m.completness < qc_criteria.min_comp:
                    failed_gids['completeness'].add(gid)
                    fout_failed.write('\tTrue')
                else:
                    fout_failed.write('\tFalse')

                if m.contamination > qc_criteria.max_cont:
                    failed_gids['contamination'].add(gid)
                    fout_failed.write('\tTrue')
                else:
                    fout_failed.write('\tFalse')

                if m.genome_quality < qc_criteria.min_quality:
                    failed_gids['genome_quality'].add(gid)
                    fout_failed.write('\tTrue')
                else:
                    fout_failed.write('\tFalse')

                if m.contig_count > qc_criteria.max_contigs:
                    failed_gids['contig_count'].add(gid)
                    fout_failed.write('\tTrue')
                else:
                    fout_failed.write('\tFalse')

                if m.n50_contigs < qc_criteria.min_N50:
                    failed_gids['n50_contigs'].add(gid)
                    fout_failed.write('\tTrue')
                else:
                    fout_failed.write('\tFalse')

                if m.ambiguous_base_perc > qc_criteria.max_ambiguous_perc:
                    failed_gids['ambiguous_base_perc'].add(gid)
                    fout_failed.write('\tTrue')
                else:
                    fout_failed.write('\tFalse')

                if refseq_exclude:
                    failed_gids['refseq_exclude'].add(gid)
                    fout_failed.write('\tTrue')
                else:
                    fout_failed.write('\tFalse')
                    
                fout_failed.write('\n')

        fout_passed.close()
        fout_failed.close()

        return qc_pass, qc_fail, failed_gids

    def species_qc_status(self, 
                            genome_metadata: Dict[str, GenomeMetadata], 
                            sp_to_gids: Dict[str, Set[str]], 
                            qc_pass: Set[str], 
                            qc_fail: Set[str],
                            failed_gids: Dict[str, Set[str]],
                            qc_criteria: QcCriteria) -> Dict[str, int]:
        """Determine species lost after QC."""

        fout_failed_sp = open(os.path.join(self.out_dir, 'species_fail_qc.tsv'), 'w')

        header = 'ncbi_species\tgid'
        header += '\tassembly_score\tcompleteness\tcontamination\tquality'
        header += '\tgenome_size\tcontig_count\tn50_contigs\tambiguous_base_perc'
        header += '\tncbi_relation_to_type_material\tncbi_refseq_category'
        header += '\tncbi_excluded_from_refseq\tncbi_assem_level'

        fout_failed_sp.write(header)
        fout_failed_sp.write(
            '\tfailed_completeness\tfailed_contamination\tfailed_quality')
        fout_failed_sp.write(
            '\tfailed_contig_count\tfailed_contig_n50\tfailed_ambiguous_bases_perc\n')

        lost_sp = set()
        lost_sp_by_reason = defaultdict(int)
        for sp, gids in sp_to_gids.items():
            if gids - qc_fail == set():
                lost_sp.add(sp)

                if gids - failed_gids['completeness'] == set():
                    lost_sp_by_reason['completeness'] += 1
                if gids - failed_gids['contamination'] == set():
                    lost_sp_by_reason['contamination'] += 1
                if gids - failed_gids['genome_quality'] == set():
                    lost_sp_by_reason['genome_quality'] += 1
                if gids - failed_gids['contig_count'] == set():
                    lost_sp_by_reason['contig_count'] += 1
                if gids - failed_gids['n50_contigs'] == set():
                    lost_sp_by_reason['n50_contigs'] += 1
                if gids - failed_gids['ambiguous_base_perc'] == set():
                    lost_sp_by_reason['ambiguous_base_perc'] += 1

                # write out information for each genome in failed species
                for gid in qc_fail.intersection(gids):
                    m = genome_metadata[gid]

                    fout_failed_sp.write('{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        sp,
                        gid,
                        m.assembly_score,
                        m.completness,
                        m.contamination,
                        m.genome_quality,
                        m.genome_size,
                        m.contig_count,
                        m.n50_contigs,
                        m.ambiguous_base_perc,
                        m.ncbi.relation_to_type_material,
                        m.ncbi.refseq_category,
                        m.ncbi.excluded_from_refseq,
                        m.ncbi.assem_level,
                    ))

                    fout_failed_sp.write(f'\t{m.completness < qc_criteria.min_comp}')
                    fout_failed_sp.write(f'\t{m.contamination > qc_criteria.max_cont}')
                    fout_failed_sp.write(f'\t{m.genome_quality < qc_criteria.min_quality}')
                    fout_failed_sp.write(f'\t{m.contig_count > qc_criteria.max_contigs}')
                    fout_failed_sp.write(f'\t{m.n50_contigs < qc_criteria.min_N50}')
                    fout_failed_sp.write(f'\t{m.ambiguous_base_perc > qc_criteria.max_ambiguous_perc}')
                    fout_failed_sp.write(f'\n')

        fout_failed_sp.close()

        return lost_sp, lost_sp_by_reason
         
    def run(self, 
            genome_metadata_file: str,
            busco_file: str, 
            ncbi_assembly_file: str, 
            ncbi_taxonomy_file: str,
            qc_criteria: QcCriteria) -> None:
        """Plot of common genomic statistics."""

        # parse NCBI assembly summary file
        self.log.info('Parsing NCBI assembly summary file:')
        metadata_ncbi, excluded_from_refseq_count = self.parse_ncbi_assembly_summary_file(ncbi_assembly_file)
        self.log.info(f' - read NCBI metadata for {len(metadata_ncbi):,} genomes')

        self.log.info('Excluded from RefSeq counts:')
        for reason, count in sorted(excluded_from_refseq_count.items(), key=lambda kv: kv[1], reverse=True):
            self.log.info(f' - {reason}: {count}')

        # parse BUSCO resutls file
        self.log.info('Paring BUSCO results file:')
        busco_results = self.parse_busco_results(busco_file)
        self.log.info(f' - read BUSCO results for {len(busco_results):,} genomes')

        # get genome metadata
        self.log.info('Parsing genome assembly metadata file:')
        genome_metadata = self.parse_genome_metadata(genome_metadata_file, busco_results, metadata_ncbi)
        self.log.info(f' - read metadata for {len(genome_metadata):,} genomes')

        # get genomes assigned to each species
        self.log.info('Parsing NCBI species assignments:')
        gid_to_sp = parse_gid_to_ncbi_sp(ncbi_taxonomy_file)

        sp_to_gids = defaultdict(set)
        for gid, sp in gid_to_sp.items():
            if gid in genome_metadata:
                sp_to_gids[sp].add(gid)
            else:
                self.log.warning(f'Skipping genome {gid} for species {sp} as no genome metadata available.')

        self.log.info(f' - identified {len(sp_to_gids):,} NCBI species across {len(gid_to_sp):,} genomes')

        # determine genomes passing and failing QC
        self.log.info('Determining genomes passing and failing QC.')
        qc_pass, qc_fail, failed_gids = self.qc_genomes(genome_metadata, gid_to_sp, qc_criteria)
        self.log.info(f' - identified {len(qc_pass):,} genomes passing and {len(qc_fail):,} genomes failing QC')

        # determine species with no remaining genomes after QC
        self.log.info('Determining species without genomes after QC:')
        lost_sp, lost_sp_by_reason = self.species_qc_status(genome_metadata, 
                                                            sp_to_gids, 
                                                            qc_pass, 
                                                            qc_fail, 
                                                            failed_gids, 
                                                            qc_criteria)
        self.log.info(f' - identified {len(lost_sp):,} lost species')

        self.log.info('Genomes and species failing each QC criterion:')
        for reason, gids in sorted(failed_gids.items(), key=lambda kv: len(kv[1]), reverse=True):
            count = len(gids)
            sp_count = lost_sp_by_reason[reason]
            self.log.info(f' - {reason} | {count:,} ({100.0*count/len(genome_metadata):.1f}%) | {sp_count:,} ({100.0*sp_count/len(sp_to_gids):.1f}%)')
