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
import sys
import gzip
import logging
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import Dict, Set, Tuple

from gtdblib.util.bio.accession import canonical_gid

from gsc_fungi.genome import assembly_score
from gsc_fungi.ncbi import (parse_ncbi_taxonomy_file,
                            is_hybrid_species,
                            MetadataNCBI,
                            NCBI_EXCLUSION_FILTERING_CRITERIA,
                            NCBI_COMPLETE_GENOME,
                            NCBI_ENV_GENOME)


@dataclass
class QcCriteria:
    """Criteria for passing QC.

    :param max_perc_dup_markers: maximum percentage of duplicate lineage-specific, single-copy genes.
    :param min_perc_phylo_markers: minimum percentage of phylogenetic marker genes to pass QC.
    :param remove_ncbi_exclusions: Filter genomes based on NCBI metadata.
    :param remove_hybrids_on_list: Filter genomes from species on the hybrid species list.
    :param keep_ncbi_hybrid_species: Do not filter genomes from NCBI species with hybrid names (e.g. 'A x B').
    """

    max_perc_dup_markers: float
    min_perc_phylo_markers: float
    remove_ncbi_exclusions: bool
    remove_hybrids_on_list: bool
    keep_ncbi_hybrid_species: bool


@dataclass
class QcLedgerEntry:
    """Entry in the QC exception ledger.

    :param species: species name given in the ledger.
    :param include: retain (True) or remove (False) genome regardless of QC criteria.
    :param route: route by which the exception was established.
    :param reason: reason for the exception.
    """

    species: str
    include: bool
    route: str
    reason: str


@dataclass
class GenomeQualityMetadata:
    phylo_marker_perc: float
    lineage_dup_marker_perc: float
    genome_size: int
    contig_count: int
    n50_contigs: int
    ambiguous_base_perc: float
    assembly_score: float
    ncbi: MetadataNCBI


@dataclass
class GenomeQualityStats:
    phylo_marker_perc: float
    lineage_dup_marker_perc: float


class QcGenomes():
    """Quality check all potential GTDB genomes."""

    def __init__(self, out_dir: str):
        """Initialization.

        :param out_dir: Output directory.
        """

        self.log = logging.getLogger('rich')
        self.out_dir = out_dir

    def parse_ncbi_assembly_summary_file(self, ncbi_assembly_file: str) -> Dict[str, MetadataNCBI]:
        """Parse NCBI assembly summary file."""

        metadata_ncbi = {}
        included_gids = {}

        open_fn = gzip.open if str(ncbi_assembly_file).endswith('.gz') else open
        with open_fn(ncbi_assembly_file, 'rt') as f:
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

                included_gids[cur_canonical_gid] = gid

                type_material = tokens[type_material_idx]
                rs_category = tokens[refseq_category_idx]
                assem_level = tokens[assem_level_idx]
                       
                metadata_ncbi[gid] = MetadataNCBI(type_material,
                                                    rs_category,
                                                    excluded_from_refseq,
                                                    assem_level)

        return metadata_ncbi

    def parse_qc_exception_list(self, qc_exception_list: str) -> Dict[str, QcLedgerEntry]:
        """Parse QC exception list and return dictionary of genomes to retain or remove."""

        qc_ledger = {}

        open_fn = gzip.open if str(qc_exception_list).endswith('.gz') else open
        with open_fn(qc_exception_list, 'rt') as f:
            header = f.readline().rstrip('\r\n').split('\t')

            genome_id_idx = header.index('Genome ID')
            species_idx = header.index('Species')
            include_idx = header.index('Include')
            route_idx = header.index('Route')
            reason_idx = header.index('Reason for exception')

            for line in f:
                tokens = line.rstrip('\r\n').split('\t')
                if not tokens[genome_id_idx].strip():
                    continue

                gid = tokens[genome_id_idx].strip()

                include = tokens[include_idx].strip().upper()
                if include not in ('TRUE', 'FALSE'):
                    self.log.error(f'Invalid Include value for {gid} in QC exception list: {tokens[include_idx]}')
                    sys.exit(1)

                if gid in qc_ledger:
                    self.log.error(f'Genome {gid} appears multiple times in QC exception list.')
                    sys.exit(1)

                qc_ledger[gid] = QcLedgerEntry(tokens[species_idx].strip(),
                                               include == 'TRUE',
                                               tokens[route_idx].strip(),
                                               tokens[reason_idx].strip())

        return qc_ledger

    def parse_hybrid_species_list(self, hybrid_species_list: str) -> Dict[str, str]:
        """Parse hybrid species list and return dictionary of genomes to retain."""

        hybrid_species = {}
        hybrid_class_counts = defaultdict(int)

        open_fn = gzip.open if str(hybrid_species_list).endswith('.gz') else open
        with open_fn(hybrid_species_list, 'rt') as f:
            header = f.readline().strip().split('\t')

            sp_idx = header.index('species')
            hybrid_class_idx = header.index('hybrid_class')

            for line in f:
                tokens = line.strip().split('\t')

                sp = tokens[sp_idx]
                if not sp.startswith('s__'):
                    sp = 's__' + sp

                hybrid_class = tokens[hybrid_class_idx]
                hybrid_species[sp] = hybrid_class
                hybrid_class_counts[hybrid_class] += 1

        return hybrid_species, hybrid_class_counts

    def parse_marker_stats_file(self, marker_set_stats_file: str) -> Dict[str, GenomeQualityStats]:
        """Parse marker set stats file."""

        marker_stats = {}

        open_fn = gzip.open if str(marker_set_stats_file).endswith('.gz') else open
        with open_fn(marker_set_stats_file, 'rt') as f:
            header = f.readline().strip().split('\t')

            phylo_comp_idx = header.index('Fun77_SingleCopy_pct_AAI100')
            lineage_dup_marker_idx = header.index('SpecSet_Duplicated_pct_AAI99')

            for line in f:
                tokens = line.strip().split('\t')

                gid = tokens[0]

                phylo_marker_perc = float(tokens[phylo_comp_idx])
                lineage_dup_marker_perc = float(tokens[lineage_dup_marker_idx])

                marker_stats[gid] = GenomeQualityStats(
                    phylo_marker_perc,
                    lineage_dup_marker_perc,
                )

        return marker_stats

    def parse_genome_metadata(self, 
                                genome_metadata_file: str, 
                                marker_stat_results: Dict[str, GenomeQualityStats],
                                metadata_ncbi: Dict[str, MetadataNCBI]):
        """Read metadata for genome assemblies and combined with BUSCO and NCBI metadata."""

        genome_qual_metadata = {}
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
                if gid not in metadata_ncbi:
                    skip_count += 1
                    continue

                genome_size = int(line_split[genome_size_idx])
                contig_count = int(line_split[contig_count_idx])
                n50 = int(line_split[n50_idx])
                ambiguous_base_perc = (100.0 * int(line_split[ambiguous_bases_idx])) / genome_size

                assem_score = assembly_score(metadata_ncbi[gid].assem_level,
                                                metadata_ncbi[gid].excluded_from_refseq,
                                                marker_stat_results[gid].phylo_marker_perc,
                                                marker_stat_results[gid].lineage_dup_marker_perc,
                                                contig_count,
                                                ambiguous_base_perc)

                genome_qual_metadata[gid] = GenomeQualityMetadata(
                                                            marker_stat_results[gid].phylo_marker_perc,
                                                            marker_stat_results[gid].lineage_dup_marker_perc,
                                                            genome_size,
                                                            contig_count,
                                                            n50,
                                                            ambiguous_base_perc,
                                                            assem_score,
                                                            metadata_ncbi[gid])

        if skip_count > 0:
            self.log.warning(f'Skipped {skip_count:,} genomes without NCBI or BUSCO metadata.')

        return genome_qual_metadata

    def qc_genomes(self, 
                    genome_qual_metadata: Dict[str, GenomeQualityMetadata], 
                    gid_to_ncbi_sp: Dict[str, str],
                    gid_to_ncbi_taxonomy: Dict[str, str],
                    hybrid_list_gids: Set[str],
                    ncbi_hybrid_gids: Set[str],
                    qc_ledger: Dict[str, QcLedgerEntry],
                    qc_criteria: QcCriteria) -> Tuple[Set[str], Set[str], Dict[str, int]]:
        """Determine genomes passing and failing QC."""

        fout_passed = open(os.path.join(self.out_dir, 'qc_passed.tsv'), 'w')
        fout_failed = open(os.path.join(self.out_dir, 'qc_failed.tsv'), 'w')

        header = 'gid\tncbi_species\tncbi_taxonomy\tqc_ledger_include'
        header += '\tassembly_score\tphylo_marker_perc\tlineage_marker_dup_perc'
        header += '\tgenome_size\tcontig_count\tn50_contigs\tambiguous_base_perc'
        header += '\tncbi_relation_to_type_material\tncbi_refseq_category'
        header += '\tncbi_excluded_from_refseq\tncbi_assem_level'

        fout_passed.write(header + '\n')
        fout_failed.write(header)
        fout_failed.write('\tfailed_phylo_marker_perc\tfailed_lineage_dup_marker_perc')
        fout_failed.write('\tfailed_refseq_exclude\tfailed_hybrid_on_list\tfailed_ncbi_hybrid_species')
        fout_failed.write('\tfailed_ledger_exclusion\n')

        qc_pass = set()
        qc_fail = set()
        failed_gids = defaultdict(set)
        for gid, m in genome_qual_metadata.items():
            # check if genome should be failed as being from a species on the hybrid list
            hybrid_on_list_exclude = qc_criteria.remove_hybrids_on_list and gid in hybrid_list_gids

            # check if genome should be failed as being from an NCBI species with a hybrid name
            ncbi_hybrid_exclude = not qc_criteria.keep_ncbi_hybrid_species and gid in ncbi_hybrid_gids

            # check if genome should be filtered based on the RefSeq excluded assignments at NCBI
            refseq_exclusion_tokens = m.ncbi.excluded_from_refseq
            refseq_exclude = False
            if qc_criteria.remove_ncbi_exclusions:
                for t in refseq_exclusion_tokens:
                    if t in NCBI_EXCLUSION_FILTERING_CRITERIA:
                        refseq_exclude = True

            pass_qc = (
                m.phylo_marker_perc >= qc_criteria.min_perc_phylo_markers
                and m.lineage_dup_marker_perc <= qc_criteria.max_perc_dup_markers
                and not refseq_exclude
                and not hybrid_on_list_exclude
                and not ncbi_hybrid_exclude)

            # genomes on the QC exception ledger are retained or removed regardless of the QC criteria
            ledger_include = ''
            ledger_exclude = False
            if gid in qc_ledger:
                ledger_include = 'TRUE' if qc_ledger[gid].include else 'FALSE'
                ledger_exclude = not qc_ledger[gid].include

            if (pass_qc and not ledger_exclude) or ledger_include == 'TRUE':
                qc_pass.add(gid)

                fout_passed.write('{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    gid,
                    gid_to_ncbi_sp[gid],
                    gid_to_ncbi_taxonomy[gid],
                    ledger_include,
                    m.assembly_score,
                    m.phylo_marker_perc,
                    m.lineage_dup_marker_perc,
                    m.genome_size,
                    m.contig_count,
                    m.n50_contigs,
                    m.ambiguous_base_perc,
                    m.ncbi.relation_to_type_material,
                    m.ncbi.refseq_category,
                    "; ".join(m.ncbi.excluded_from_refseq),
                    m.ncbi.assem_level,
                ))
            else:
                qc_fail.add(gid)

                fout_failed.write('{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    gid,
                    gid_to_ncbi_sp[gid],
                    gid_to_ncbi_taxonomy[gid],
                    ledger_include,
                    m.assembly_score,
                    m.phylo_marker_perc,
                    m.lineage_dup_marker_perc,
                    m.genome_size,
                    m.contig_count,
                    m.n50_contigs,
                    m.ambiguous_base_perc,
                    m.ncbi.relation_to_type_material,
                    m.ncbi.refseq_category,
                    "; ".join(m.ncbi.excluded_from_refseq),
                    m.ncbi.assem_level,
                ))

                if m.phylo_marker_perc < qc_criteria.min_perc_phylo_markers:
                    failed_gids['phylo_marker_perc'].add(gid)
                    fout_failed.write('\tTrue')
                else:
                    fout_failed.write('\tFalse')

                if m.lineage_dup_marker_perc > qc_criteria.max_perc_dup_markers:
                    failed_gids['lineage_dup_marker_perc'].add(gid)
                    fout_failed.write('\tTrue')
                else:
                    fout_failed.write('\tFalse')

                if refseq_exclude:
                    failed_gids['refseq_exclude'].add(gid)
                    fout_failed.write('\tTrue')
                else:
                    fout_failed.write('\tFalse')

                if hybrid_on_list_exclude:
                    failed_gids['hybrid_on_list'].add(gid)
                    fout_failed.write('\tTrue')
                else:
                    fout_failed.write('\tFalse')

                if ncbi_hybrid_exclude:
                    failed_gids['ncbi_hybrid_species'].add(gid)
                    fout_failed.write('\tTrue')
                else:
                    fout_failed.write('\tFalse')

                if ledger_exclude:
                    failed_gids['ledger_exclusion'].add(gid)
                    fout_failed.write('\tTrue')
                else:
                    fout_failed.write('\tFalse')

                fout_failed.write('\n')

        fout_passed.close()
        fout_failed.close()

        return qc_pass, qc_fail, failed_gids

    def species_qc_status(self, 
                            genome_qual_metadata: Dict[str, GenomeQualityMetadata], 
                            sp_to_gids: Dict[str, Set[str]], 
                            qc_pass: Set[str], 
                            qc_fail: Set[str],
                            failed_gids: Dict[str, Set[str]],
                            qc_criteria: QcCriteria) -> Dict[str, int]:
        """Determine species lost after QC."""

        fout_failed_sp = open(os.path.join(self.out_dir, 'species_fail_qc.tsv'), 'w')

        header = 'ncbi_species\tgid'
        header += '\tassembly_score\tphylo_marker_perc\tlineage_marker_dup_perc'
        header += '\tgenome_size\tcontig_count\tn50_contigs\tambiguous_base_perc'
        header += '\tncbi_relation_to_type_material\tncbi_refseq_category'
        header += '\tncbi_excluded_from_refseq\tncbi_assem_level'

        fout_failed_sp.write(header)
        fout_failed_sp.write('\tfailed_phylo_marker_perc\tfailed_lineage_dup_marker_perc')
        fout_failed_sp.write('\tfailed_refseq_exclude\tfailed_hybrid_on_list\tfailed_ncbi_hybrid_species')
        fout_failed_sp.write('\tfailed_ledger_exclusion\n')

        lost_sp = set()
        lost_sp_by_reason = defaultdict(int)
        for sp, gids in sp_to_gids.items():
            if gids - qc_fail == set():
                lost_sp.add(sp)

                if gids - failed_gids['phylo_marker_perc'] == set():
                    lost_sp_by_reason['phylo_marker_perc'] += 1
                if gids - failed_gids['lineage_dup_marker_perc'] == set():
                    lost_sp_by_reason['lineage_dup_marker_perc'] += 1
                if gids - failed_gids['refseq_exclude'] == set():
                    lost_sp_by_reason['refseq_exclude'] += 1
                if gids - failed_gids['hybrid_on_list'] == set():
                    lost_sp_by_reason['hybrid_on_list'] += 1
                if gids - failed_gids['ncbi_hybrid_species'] == set():
                    lost_sp_by_reason['ncbi_hybrid_species'] += 1
                if gids - failed_gids['ledger_exclusion'] == set():
                    lost_sp_by_reason['ledger_exclusion'] += 1

                # write out information for each genome in failed species
                for gid in qc_fail.intersection(gids):
                    m = genome_qual_metadata[gid]

                    fout_failed_sp.write('{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        sp,
                        gid,
                        m.assembly_score,
                        m.phylo_marker_perc,
                        m.lineage_dup_marker_perc,
                        m.genome_size,
                        m.contig_count,
                        m.n50_contigs,
                        m.ambiguous_base_perc,
                        m.ncbi.relation_to_type_material,
                        m.ncbi.refseq_category,
                        '; '.join(m.ncbi.excluded_from_refseq),
                        m.ncbi.assem_level,
                    ))

                    fout_failed_sp.write(f'\t{gid in failed_gids["phylo_marker_perc"]}')
                    fout_failed_sp.write(f'\t{gid in failed_gids["lineage_dup_marker_perc"]}')
                    fout_failed_sp.write(f'\t{gid in failed_gids["refseq_exclude"]}')
                    fout_failed_sp.write(f'\t{gid in failed_gids["hybrid_on_list"]}')
                    fout_failed_sp.write(f'\t{gid in failed_gids["ncbi_hybrid_species"]}')
                    fout_failed_sp.write(f'\t{gid in failed_gids["ledger_exclusion"]}\n')

        fout_failed_sp.close()

        return lost_sp, lost_sp_by_reason
         
    def run(self, 
            genome_metadata_file: str,
            marker_set_stats_file: str,
            ncbi_assembly_file: str, 
            ncbi_taxonomy_file: str,
            qc_exception_list: str,
            hybrid_species_list: str,
            qc_criteria: QcCriteria) -> None:
        """Plot of common genomic statistics."""

        # parse NCBI assembly summary file
        self.log.info('Parsing NCBI assembly summary file:')
        metadata_ncbi = self.parse_ncbi_assembly_summary_file(ncbi_assembly_file)
        self.log.info(f' - read NCBI metadata for {len(metadata_ncbi):,} genomes')

        # parse QC exception file
        self.log.info('Parsing QC exception list:')
        qc_ledger = self.parse_qc_exception_list(qc_exception_list)
        include_count = sum(1 for e in qc_ledger.values() if e.include)
        self.log.info(f' - identified {include_count:,} genomes to retain and {len(qc_ledger) - include_count:,} genomes to remove')

        route_count = defaultdict(lambda: defaultdict(int))
        for e in qc_ledger.values():
            route_count[e.route][e.include] += 1

        self.log.info('QC exception routes (genomes to retain | genomes to remove):')
        for route, counts in sorted(route_count.items(), key=lambda kv: sum(kv[1].values()), reverse=True):
            self.log.info(f' - {route} | {counts[True]:,} | {counts[False]:,}')

        # parse hybrid species list
        self.log.info('Parsing hybrid species list:')
        hybrid_species, hybrid_class_counts = self.parse_hybrid_species_list(hybrid_species_list)
        self.log.info(f' - identified {len(hybrid_species):,} hybrid species')

        self.log.info('Hybrid classes:')
        for hybrid_class, count in sorted(hybrid_class_counts.items(), key=lambda kv: kv[1], reverse=True):
            self.log.info(f' - {hybrid_class}: {count}')


        # parse marker set stats file
        self.log.info('Parsing marker set stats file:')
        marker_set_stats = self.parse_marker_stats_file(marker_set_stats_file)
        self.log.info(f' - read marker set stats for {len(marker_set_stats):,} genomes')

        # get genome metadata
        self.log.info('Parsing genome assembly metadata file:')
        genome_qual_metadata = self.parse_genome_metadata(genome_metadata_file, marker_set_stats, metadata_ncbi)
        self.log.info(f' - read metadata for {len(genome_qual_metadata):,} genomes')

        # get genomes assigned to each species
        self.log.info('Parsing NCBI species assignments:')
        gid_to_ncbi_taxonomy, gid_to_ncbi_sp = parse_ncbi_taxonomy_file(ncbi_taxonomy_file)

        # sanity check that genomes in the QC exception list are assigned to the expected
        # NCBI species and that genomes to retain aren't from a hybrid species
        for gid, e in qc_ledger.items():
            if gid not in gid_to_ncbi_sp:
                self.log.warning(f'Genome {gid} in QC exception list has no NCBI species assignment.')
                continue

            ncbi_sp = gid_to_ncbi_sp[gid]
            if ncbi_sp != 's__' and ncbi_sp.replace('s__', '') != e.species:
                self.log.warning(f'Genome {gid} in QC exception list is from {e.species}, but assigned to NCBI species {ncbi_sp}.')

            if e.include and ncbi_sp in hybrid_species:
                self.log.warning(f'Genome {gid} to retain in QC exception list is from a hybrid species: {ncbi_sp}')
                sys.exit(0)

        # determine hybrid genomes
        self.log.info('Identifying genomes from hybrid species:')
        hybrid_list_gids = set()
        hybrid_list_sp = set()
        ncbi_hybrid_gids = set()
        ncbi_hybrid_sp = set()
        for gid, sp in gid_to_ncbi_sp.items():
            if sp in hybrid_species:
                hybrid_list_gids.add(gid)
                hybrid_list_sp.add(sp)
            if is_hybrid_species(sp):
                ncbi_hybrid_gids.add(gid)
                ncbi_hybrid_sp.add(sp)
        self.log.info(f' - identified {len(hybrid_list_gids):,} genomes from {len(hybrid_list_sp):,} species on the hybrid species list')
        self.log.info(f' - identified {len(ncbi_hybrid_gids):,} genomes from {len(ncbi_hybrid_sp):,} NCBI species with hybrid names')

        # warn about species on the hybrid list that will not be filtered as they match no NCBI species,
        # and report NCBI species with the same specific name as these may indicate a genus transfer
        if qc_criteria.remove_hybrids_on_list:
            ncbi_sp_by_epithet = defaultdict(set)
            for sp in set(gid_to_ncbi_sp.values()):
                ncbi_sp_by_epithet[sp.split()[-1]].add(sp)

            unmatched_sp = sorted(set(hybrid_species) - hybrid_list_sp)
            if unmatched_sp:
                self.log.warning(f'Identified {len(unmatched_sp):,} species on the hybrid species list that match no NCBI species:')
            for sp in unmatched_sp:
                same_epithet_sp = sorted(ncbi_sp_by_epithet.get(sp.split()[-1], set()))
                if same_epithet_sp:
                    self.log.warning(f' - {sp}: NCBI species with same specific name (possible genus transfer): {", ".join(same_epithet_sp)}')
                else:
                    self.log.warning(f' - {sp}: no NCBI species with same specific name')

        # determine genomes in each species
        self.log.info('Determining genomes assigned to each NCBI species:')
        sp_to_gids = defaultdict(set)
        for gid, sp in gid_to_ncbi_sp.items():
            if gid in genome_qual_metadata:
                sp_to_gids[sp].add(gid)
            else:
                self.log.warning(f'Skipping genome {gid} for species {sp} as no genome metadata available.')

        self.log.info(f' - identified {len(sp_to_gids):,} NCBI species across {len(gid_to_ncbi_sp):,} genomes')

        # determine genomes passing and failing QC
        self.log.info('Determining genomes passing and failing QC:')
        qc_pass, qc_fail, failed_gids = self.qc_genomes(genome_qual_metadata, 
                                                        gid_to_ncbi_sp,
                                                        gid_to_ncbi_taxonomy,
                                                        hybrid_list_gids,
                                                        ncbi_hybrid_gids,
                                                        qc_ledger,
                                                        qc_criteria)
        self.log.info(f' - identified {len(qc_pass):,} genomes passing and {len(qc_fail):,} genomes failing QC')

        # report genomes with each excluded from RefSeq annotation and how many failed QC for any reason
        excluded_from_refseq_count = defaultdict(int)
        excluded_from_refseq_filtered = defaultdict(int)
        for gid, m in genome_qual_metadata.items():
            for reason in set(m.ncbi.excluded_from_refseq):
                excluded_from_refseq_count[reason] += 1
                if gid in qc_fail:
                    excluded_from_refseq_filtered[reason] += 1

        self.log.info('Excluded from RefSeq annotations (genomes with annotation | genomes filtered):')
        for reason, count in sorted(excluded_from_refseq_count.items(), key=lambda kv: kv[1], reverse=True):
            self.log.info(f' - {reason} | {count:,} | {excluded_from_refseq_filtered[reason]:,}')

        # determine species with no remaining genomes after QC
        self.log.info('Determining species without genomes after QC:')
        lost_sp, lost_sp_by_reason = self.species_qc_status(genome_qual_metadata, 
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
            self.log.info(f' - {reason} | {count:,} ({100.0*count/len(genome_qual_metadata):.1f}%) | {sp_count:,} ({100.0*sp_count/len(sp_to_gids):.1f}%)')
