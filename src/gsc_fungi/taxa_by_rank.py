#! /usr/bin/env python3

import os
import logging
from collections import defaultdict
from typing import Dict, Set, Tuple

from gtdblib.util.bio.accession import canonical_gid


RANKS = {
    'd': 'domain',
    'sd': 'subdomain',
    'p': 'phylum',
    'sp': 'subphylum',
    'c': 'class',
    'sc': 'subclass',
    'o': 'order',
    'so': 'suborder',
    'f': 'family',
    'sf': 'subfamily',
    'g': 'genus',
    'sg': 'subgenus',
    's': 'species'
}


class TaxaByRank():
    """Get number of genomes and taxa at each rank."""

    def __init__(self, out_dir: str) -> None:
        """Initialization."""

        self.log = logging.getLogger('rich')
        self.out_dir = out_dir

    def identify_env_genomes(self, assembly_summary_file: str) -> set:
        """Identify MAGs and SAGs from assembly summary file."""

        env_genomes = set()
        with open(assembly_summary_file) as f:
            f.readline()
            header = f.readline().strip().split('\t')

            refseq_idx = header.index('excluded_from_refseq')

            for line in f:
                tokens = line.strip().split('\t')

                gid = canonical_gid(tokens[0])

                rs_status = tokens[refseq_idx]
                if ('derived from metagenome' in rs_status
                    or 'derived from single cell' in rs_status):
                    env_genomes.add(gid)

        return env_genomes

    def determine_missing_classifications(self, 
                                            ncbi_taxonomy_file: str, 
                                            env_genomes: set) -> Tuple[Dict[str, int], Dict[str, int]]:
        """Determine number of genomes and environmental genomes with missing classifications at each rank."""

        missing = defaultdict(int)
        env_missing = defaultdict(int)
        with open(ncbi_taxonomy_file) as f:
            for line in f:
                tokens = line.strip().split('\t')

                gid = canonical_gid(tokens[0])
                taxa = tokens[1].split(';')

                for taxon in taxa:
                    if len(taxon) == 3:
                        rank = taxon.split('__')[0]
                        missing[rank] += 1
                        if gid in env_genomes:
                            env_missing[rank] += 1
                    
        return missing, env_missing

    def parse_ncbi_taxonomy(self, ncbi_taxonomy_file: str) -> Tuple[Dict[str, Set[str]], Dict[str, int], Dict[str, int], int]:
        """Parse NCBI taxonomy file."""

        taxa_in_rank = defaultdict(set)
        taxon_count = defaultdict(int)
        genome_count_by_rank = defaultdict(int)
        total_genomes = 0
        ranks = set()
        with open(ncbi_taxonomy_file) as f:
            for line in f:
                tokens = line.strip().split('\t')

                taxa = tokens[1].split(';')
                total_genomes += 1

                for cur_taxon in taxa:
                    rank, taxon = cur_taxon.split('__')
                    taxa_in_rank[rank].add(cur_taxon)
                    ranks.add(rank)

                    if taxon:
                        genome_count_by_rank[rank] += 1
                        taxon_count[cur_taxon] += 1
                    
        assert ranks == set(RANKS.keys()), f'Ranks in taxonomy file do not match expected ranks: {ranks}'

        return taxa_in_rank, taxon_count, genome_count_by_rank, total_genomes

    def write_missing_classification_table(self, missing, env_missing, total_genomes, total_env_genomes) -> None:
        """Write missing classification table."""

        fout = open(os.path.join(self.out_dir, 'missing_classifications_by_rank.tsv'), 'w')
        fout.write('rank\tnum_genomes_missing\tgenomes_missing_perc\tnum_env_genomes_missing\tenv_genomes_missing_perc\n')

        for rank_prefix, rank_label in RANKS.items():
            # skip subranks as only the standard 7 ranks have an
            # expectation of being defined across all genomes
            if rank_label.startswith('sub'):
                continue

            num_missing = missing.get(rank_prefix, 0)
            num_env_missing = env_missing.get(rank_prefix, 0)

            fout.write('{}\t{}\t{}\t{}\t{}\n'.format(
                rank_label,
                num_missing,
                f'{(num_missing / total_genomes * 100):.2f}%',
                num_env_missing,
                f'{(num_env_missing / total_env_genomes * 100):.2f}%'
            ))

        fout.close()

    def write_taxa_by_rank_table(self, 
                                    taxa_in_rank: Dict[str, Set[str]], 
                                    taxon_count: Dict[str, int], 
                                    genome_count_by_rank: Dict[str, int], 
                                    total_genomes: int) -> None:
        """Write taxa by rank table."""

        fout = open(os.path.join(self.out_dir, 'taxa_by_rank.tsv'), 'w')
        fout.write('rank\tnum_genomes\tgenomes_perc\tnum_taxa\ttaxa\n')

        for rank_prefix, rank_label in RANKS.items():
            taxon_str = {}
            for taxon in taxa_in_rank[rank_prefix]:
                taxon_str[taxon] = taxon_count[taxon]

            taxon_str = ', '.join([f'{k}:{v}' for k, v in sorted(taxon_str.items(), key = lambda kv: kv[1], reverse=True)])

            fout.write('{}\t{}\t{}\t{}\t{}\n'.format(
                rank_label,
                genome_count_by_rank[rank_prefix],
                f'{(genome_count_by_rank[rank_prefix] / total_genomes * 100):.2f}%',
                len(taxa_in_rank[rank_prefix]),
                taxon_str
            ))

        fout.close()

    def run(self, ncbi_taxonomy_file: str, assembly_summary_file: str) -> None:
        """Get number of genomes and taxa at each rank."""

        # determine environmental genomes
        self.log.info('Identifying environmental genomes (MAGs and SAGs):')
        env_genomes = self.identify_env_genomes(assembly_summary_file)
        self.log.info(f' - identified {len(env_genomes):,} environmental genomes')

        # determine missing classification
        self.log.info('Determining missing classifications at each rank:')
        missing, env_missing = self.determine_missing_classifications(ncbi_taxonomy_file, env_genomes)
        total_missing = sum(missing.values())
        total_env_missing = sum(env_missing.values())
        self.log.info(f' - determined {total_missing:,} missing classifications ({total_env_missing:,} from environmental genomes)')

        # determine taxa and genome counts at each rank
        self.log.info('Parsing taxonomic classification for genomes:')
        taxa_in_rank, taxon_count, genome_count_by_rank, total_genomes = self.parse_ncbi_taxonomy(ncbi_taxonomy_file)
        self.log.info(f' - parsed taxonomy for {total_genomes:,} genomes')

        # write results tables
        self.write_missing_classification_table(missing, env_missing, total_genomes, len(env_genomes))
        self.write_taxa_by_rank_table(taxa_in_rank, taxon_count, genome_count_by_rank, total_genomes)


