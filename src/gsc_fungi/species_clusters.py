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
import pickle
from collections import defaultdict

from gtdblib.util.shell.execute import check_dependencies

from gsc_fungi import defaults as Defaults
from gsc_fungi.genome_utils import read_genome_path
from gsc_fungi.skani import Skani
from gsc_fungi import mycobank as Mycobank
from gsc_fungi.ncbi import parse_gid_to_ncbi_sp


class SpeciesClusters():
    """Create ANI-based species clusters."""

    def __init__(self, out_dir: str):
        """Initialization.

        :param out_dir: Output directory.
        """

        self.log = logging.getLogger('rich')
        self.out_dir = out_dir

        check_dependencies(['skani'])

    def parse_pass_qc_assembly_score(self, qc_pass_file: str) -> Set[str]:
        """Get assembly score of genomes passing QC."""

        gids_pass_qc = set()
        with open(qc_pass_file) as f:
            header = f.readline().strip().split('\t')
            assem_score_idx = header.index('assembly_score')

            for line in f:
                tokens = line.strip().split('\t')
                gids_pass_qc[tokens[0]] = float(tokens[assem_score_idx])

        return gids_pass_qc

    def calculate_pairwise_ani(self, gids: Set[str], genome_files: Dict[str, str]):
        """Calculate pairwise ANI between genomes."""

        if True: # *** DEBUGGING
            genomi_files_qc = {gid: gf for gid, gf in genome_files.items() if gid in gids}
            ani_af = self.skani.triangle(genomi_files_qc, 
                                            self.out_dir, 
                                            preset = Defaults.SKANI_PRESET,
                                            min_af = Defaults.AF_SP,
                                            min_sketch_ani = Defaults.SKANI_PREFILTER_THRESHOLD)
            pickle.dump(ani_af, open(os.path.join(
                self.out_dir, 'ani_af.pkl'), 'wb'))
        else:
            ani_af = pickle.load(
                open(os.path.join(self.out_dir, 'ani_af.pkl'), 'rb'))

        return ani_af

    def normalize_species_name(self, sp_name: str) -> str:
        """Determine normalized species name.
        
        Currently this only removed brackets ([, ]) that
        can surround the generic name in NCBI species names. 
        This indicates NCBI disagrees with the genus assignment,
        but prevents matching with data from MycoBank.
        """

        return sp_name.replace('[', '').replace(']', '')

    def resolve_naming_priority(sp1: str, 
                                year1: int, 
                                name_status1: str, 
                                type_strain1: bool, 
                                sanctioned_name1: bool,
                                sp2: str, 
                                year2: int, 
                                name_status2: str, 
                                type_strain2: bool, 
                                sanctioned_name2: bool) -> Tuple[str, str]:
        """Resolve naming priority between two species.
        
        Priority order:
        - name is a legitimate fungal name according to MycoBank  
        - genomes assembled from type strain have priority over strains that are not type material 
        - sanctioned names that have priority over earlier homonyms and competing synonyms
        - year of effective publication (early = high naming priority)
        """

        # being a valid name has the highest priority
        if name_status1 in Mycobank.VALID_NAME_STATUS and name_status2 not in Mycobank.VALID_NAME_STATUS:
            return sp1, 'legitimate or orthographic variant has priority; orthographic variant should be corrected'
        elif name_status1 not in Mycobank.VALID_NAME_STATUS and name_status2 in Mycobank.VALID_NAME_STATUS:
            return sp2, 'legitimate or orthographic variant has priority; orthographic variant should be corrected'

        # being the type strain of the species has the next highest naming priority as
        # other strains have no standing in nomenclature
        if type_strain1 and not type_strain2:
            return sp1, 'type strain of species has priority'
        elif not type_strain1 and type_strain2:
            return sp2, 'type strain of species has priority'

        # sanctioned names have priority over unsanctioned names
        if sanctioned_name1 and not sanctioned_name2:
            return sp1, 'sanctioned name has priority'
        elif not sanctioned_name1 and sanctioned_name2:
            return sp2, 'sanctioned name has priority'
        
        # earlier published names have priority if they are a legitimate name
        if year1 < year2 and name_status1 in Mycobank.VALID_NAME_STATUS:
            return sp1, 'earlier effective publication date'
        elif year2 < year1 and name_status2 in Mycobank.VALID_NAME_STATUS:
            return sp2, 'earlier effective publication date'
        
        return Mycobank.UNDETERMINED_NAMING_PRIORITY, '<priority must be manually resolved>'

    def merge_sp_reps(self,
                        ani_af: Dict[str, Dict[str, List[float]]], 
                        rid_to_sp: Dict[str, str], 
                        sp_to_rid: Dict[str, str], 
                        sp_type_strains: Dict[str, Set[str]], 
                        ani_threshold: float, 
                        af_threshold: float) -> List[str]:
        """Resolve naming priority for species representatives merged at a given ANI and AF threshold."""

        # determine representatives that will be merged at ANI and AF thresholds
        merged_spp = set()
        merged_pair = defaultdict(lambda: {})
        for gid1, gid2 in combinations(rid_to_sp, 2):
            ani, af = Skani.symmetric_ani_af(ani_af, gid1, gid2)

            if ani >= ani_threshold and af >= af_threshold:
                sp1 = rid_to_sp[gid1]
                sp2 = rid_to_sp[gid2]

                merged_spp.add(sp1)
                merged_spp.add(sp2)

                # save in alphabetical order so identical pairs
                # are reported only once
                if sp1 > sp2:
                    sp1, sp2 = sp2, sp1

                merged_pair[sp1][sp2] = (ani, af)

        perc_merged = 100.0 * len(merged_spp) / len(rid_to_sp)
        self.log.info(' - merging stats (no. species / perc):', len(merged_spp), f'{perc_merged:.1f}')

        # get genus level counts for merged species
        genus_merged_spp = defaultdict(set)
        for sp in merged_spp:
            genus = sp.split()[0].replace('s__', '')
            genus_merged_spp[genus].add(sp)

        fout = open(os.path.join(SP_REP_DIR, f'species_merged-genus_counts-ani{ani_threshold}_af{af_threshold}.tsv'), 'w')
        fout.write('ncbi_genus\tnum_species_merged\tmerged_species_perc\tncbi_species\n')
        for genus, count in sorted(genus_merged_spp.items(), key=lambda kv: len(kv[1]), reverse=True):
            fout.write('{}\t{}\t{:.2f}\t{}\n'.format(
                genus,
                len(genus_merged_spp[genus]),
                100.0*len(genus_merged_spp[genus]) / len(merged_spp),
                ','.join(sorted(genus_merged_spp[genus]))
            ))
        fout.close()

        # write out species that are subject to merging at given ANI and AF criteria
        fout = open(os.path.join(SP_REP_DIR, f'species_merged-naming_priority-ani{ani_threshold}_af{af_threshold}.tsv'), 'w')
        fout.write('species1\tgid1\tyear_publication1\tname_status1\tis_type_strain_of_species1\tsanctioned_name1')
        fout.write('\tspecies2\tgid2\tyear_publication2\tname_status2\tis_type_strain_of_species2\tsanctioned_name2')
        fout.write('\tpriority_species\treason\tani\taf\n')

        priority_rids = []
        undetermined_priority_count = 0
        for sp1 in sorted(merged_pair):
            for sp2 in sorted(merged_pair[sp1]):
                ani, af = merged_pair[sp1][sp2]

                sp1_norm = normalize_species_name(sp1)
                year1 = mycobank_data[sp1_norm].year_publication if sp1_norm in mycobank_data else YEAR_UNKNOWN
                name_status1 = mycobank_data[sp1_norm].name_status if sp1_norm in mycobank_data else NOT_IN_MYCOBANK
                type_strain1 = sp1 in sp_type_strains
                sanctioned_name1 = sp1 in sanctioned_names

                sp2_norm = normalize_species_name(sp2)
                year2 = mycobank_data[sp2_norm].year_publication if sp2_norm in mycobank_data else YEAR_UNKNOWN
                name_status2 = mycobank_data[sp2_norm].name_status if sp2_norm in mycobank_data else NOT_IN_MYCOBANK
                type_strain2 = sp2 in sp_type_strains
                sanctioned_name2 = sp2 in sanctioned_names

                priority_sp, reason = resolve_naming_priority(
                    sp1, year1, name_status1, type_strain1, sanctioned_name1,
                    sp2, year2, name_status2, type_strain2, sanctioned_name2)
                
                if priority_sp == sp1:
                    priority_rids.append((sp_to_rid[sp1], sp_to_rid[sp2]))
                elif priority_sp == sp2:
                    priority_rids.append((sp_to_rid[sp2], sp_to_rid[sp1]))
                elif priority_sp == UNDETERMINED_NAMING_PRIORITY:
                    undetermined_priority_count += 1

                fout.write('{}\t{}\t{}\t{}\t{}\t{}'.format(
                    sp1,
                    sp_to_rid[sp1],
                    year1,
                    name_status1,
                    type_strain1,
                    sanctioned_name1
                ))

                fout.write('\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    sp2,
                    sp_to_rid[sp2],
                    year2,
                    name_status2,
                    type_strain2,
                    sanctioned_name2
                ))

                fout.write('\t{}\t{}\t{:.2f}\t{:.2f}\n'.format(
                    priority_sp,
                    reason,
                    ani,
                    af
                ))

        fout.close()

        self.log.info(' - undetermined_priority_count', undetermined_priority_count)

        # generate list of species representative genomes sorted such
        # that genomes with naming priority allows come before any
        # genome they would be merged with. Note that this is NOT
        # a fully sorted list by naming priority. The only guarantee
        # in ordering is between species (representative genomes) 
        # that will be merged.
        rids_by_priority = deque()
        for rid1, rid2 in priority_rids:
            if rid1 not in rids_by_priority:
                rids_by_priority.appendleft(rid1)
            
            if rid2 not in rids_by_priority:
                rids_by_priority.append(rid2)

            if rid1 in rids_by_priority and rid2 in rids_by_priority:
                rid1_idx = rids_by_priority.index(rid1)
                rid2_idx = rids_by_priority.index(rid2)

                # check if elements need to be swapped to ensure correct ordering
                if rid1_idx > rid2_idx:
                    rids_by_priority[rid1_idx], rids_by_priority[rid2_idx] = rids_by_priority[rid2_idx], rids_by_priority[rid1_idx]

            # ensure expected sorting gaurantee is obtained
            assert rids_by_priority.index(rid1) < rids_by_priority.index(rid2)

        # sanity check that sorting gaurantee is obtained
        for rid1, rid2 in priority_rids:
            rids_by_priority.index(rid1) < rids_by_priority.index(rid2)

        return list(rids_by_priority)

    def create_sp_cluters(self, 
                            ani_af: Dict[str, Dict[str, List[float]]], 
                            merged_rids_by_priority: List[str], 
                            sp_rids: Dict[str, str], 
                            gids_pass_qc: Dict[str, float], 
                            ani_threshold: float) -> Dict[str, List[str]]:
        """Create GTDB species clusters."""

        # sort all genomes by their assembly score, except we want
        # to process all NCBI species representatives first so we
        # add a large value to their score and need to ensure we
        # select the representative genomes from the species
        # with naming priority for species that will be merged
        # at the ANI and AF criteria
        scores = {}
        for gid, score in gids_pass_qc.items():
            scores[gid] = score
            if gid in sp_rids:
                scores[gid] += 1_000_000

        gids_sorted_by_score = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)

        # put merged representative species first as these are all species
        # representatives and are in the order required to preserve naming priority
        gids_sorted = merged_rids_by_priority
        for gid, _score in gids_sorted_by_score:
            if gid not in merged_rids_by_priority:
                gids_sorted.append(gid)

        # create species cluters
        clusters = {}
        for gid in gids_sorted:
            # check if genomes can be assigned to an existing cluster
            is_clustered = False
            for rid in clusters:
                ani, af = Skani.symmetric_ani_af(ani_af, rid, gid)

                if af >= Defaults.AF_SP and ani >= ani_threshold:
                    is_clustered = True
                    break

            if not is_clustered:
                # create new species cluster
                clusters[gid] = [gid]

        # assign genomes to closest species cluter
        for gid in gids_sorted:
            if gid in clusters:
                continue

            # find closest cluster
            closest_ani = 0
            closest_rid = None
            for rid in clusters:
                ani = anis[rid].get(gid, 0)
                af = max(afs[rid].get(gid, 0), afs[gid].get(rid, 0))

                if af >= AF_SP_THRESHOLD and ani >= closest_ani:
                    closest_ani = ani
                    closest_rid = rid

            if closest_ani >= ani_threshold:
                clusters[closest_rid].append(gid)
            else:
                self.log.error(f'Failed to assign genome to species cluster: {closest_ani}')
                sys.exit(1)

        return clusters

    def parse_species_representatives(self, sp_rep_file: str) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, Set[str]]]:
        """Parse metadata for species representative genomes."""

        rid_to_sp = {}
        sp_to_rid = {}
        sp_type_strains = set()
        with open(sp_rep_file) as f:
            header = f.readline().strip().split('\t')

            sp_idx = header.index('ncbi_sp')
            gid_idx = header.index('gid')
            is_rep_idx = header.index('is_gtdb_rep')
            assem_score_idx = header.index('assembly_score')
            is_type_strain_idx = header.index('ncbi_is_type_strain_of_species')

            for line in f:
                tokens = line.strip().split('\t')

                sp = tokens[sp_idx]

                is_rep = tokens[is_rep_idx] == 'True'
                if is_rep:
                    gid = tokens[gid_idx]
                    assert gid in gids_pass_qc
                    assert gids_pass_qc[gid] == float(tokens[assem_score_idx])
                    rid_to_sp[gid] = sp
                    sp_to_rid[sp] = gid

                    if tokens[is_type_strain_idx].lower().startswith('t'):
                        sp_type_strains[sp].add(gid)

        return rid_to_sp, sp_to_rid, sp_type_strains

    def write_sp_clusters(self, 
                            sp_clusters: Dict[str, List[str]], 
                            gid_to_sp: Dict[str, str], 
                            rid_to_sp: Dict[str, str], 
                            out_file: str) -> None:
        """Write out species clusters to file."""

        fout = open(out_file, 'w')
        fout.write('representative_id\tncbi_species\tnum_clustered\tclustered_gids')
        fout.write('\tnum_merged_species\tsp_reps_in_cluster\tncbi_sp_classification\n')

        for rid, cids in sp_clusters.items():
            # generate string indicating number of times each NCBI species
            # classification occurs in the cluter
            ncbi_sp_counts = defaultdict(int)
            for gid in cids:
                sp = gid_to_sp[gid]
                ncbi_sp_counts[sp] += 1

            ncbi_sp_counts_sorted = sorted(ncbi_sp_counts.items(), key=lambda kv: kv[1], reverse=True)
            ncbi_sp_str = []
            for sp, count in ncbi_sp_counts_sorted:
                ncbi_sp_str.append(f'{sp}: {count}')
            ncbi_sp_str = ', '.join(ncbi_sp_str)

            # generate string indicating species representatives in cluster
            cur_sp_rids = set(cids).intersection(rid_to_sp)
            rep_sp = [rid_to_sp[rid] for rid in cur_sp_rids]

            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                rid,
                rid_to_sp.get(rid, "<de novo>"),
                len(cids),
                ','.join(cids),
                len(rep_sp)-1,
                ','.join(rep_sp),
                ncbi_sp_str
            ))

            fout.close()

    def run(self,
            qc_pass_file: str,
            sp_rep_file: str,
            ncbi_taxonomy_file: str,
            genome_path_file: str,
            mycobank_file: str,
            mycobank_sanctioned_names_persoon_file: str,
            mycobank_sanctioned_names_ef_file: str) -> None:
        """Create ANI-based species clusters."""

        # read path to genomic FASTA files
        self.log.info('Reading path to genomic FASTA files:')
        genome_files = read_genome_path(genome_path_file)
        self.log.info(f' - read path for {len(genome_files):,} genomes')

        # read genomes passing QC
        self.log.info('Reading genomes passing QC:')
        gids_pass_qc = self.parse_pass_qc_assembly_score(qc_pass_file)
        self.log.info(f' - identified {len(gids_pass_qc):,} genomes passing QC')

        # calculate pairwise ANI values between all genomes
        ani_af = self.calculate_pairwise_ani(gids_qc, genome_files)

        # read list of sanctioned names that always have naming priority
        self.log.info('Reading sanctioned names:')
        sanctioned_names = Mycobank.read_sanctioned_names(mycobank_sanctioned_names_persoon_file)
        sanctioned_names.update(Mycobank.read_sanctioned_names(mycobank_sanctioned_names_ef_file))
        self.log.info(f' - read {len(sanctioned_names):,} sanctioned names')

        # read MycoBank nomenclature information used to 
        # resolve naming priority
        self.log.info('Reading MycoBank data:')
        mycobank_data = Mycobank.parse_mycobank_data(mycobank_file)
        self.log.info(f' - read data for {len(mycobank_data):,} genera and species')

        # read GTDB representatives for named NCBI fungal species
        self.log.info('Reading GTDB species representatives:')
        rid_to_sp, sp_to_rid, sp_type_strains = parse_species_representatives(sp_rep_file)
        self.log.info(f' - identified {len(rid_to_sp):,} representative genomes')

        # get NCBI species classification for genomes
        self.log.info('Reading NCBI species assignment for all genomes:')
        gid_to_ncbi_sp = parse_gid_to_ncbi_sp(ncbi_taxonomy_file)
        self.log.info(f' - identified NCBI species for {len(gid_to_ncbi_sp):,} genomes')

        # create species cluters
        self.log.info('Creating species cluters at different ANI thresholds:')
        for ani_threshold in [95]: #*** range(92, 100):
            self.log.info(f' - {ani_threshold}')

            # determine representatives that will be merged at ANI and AF thresholds
            self.logger.info('- determining species representatives to merged')
            merged_rids_by_priority = self.merge_sp_reps(ani_af,
                                                            rid_to_sp, 
                                                            sp_to_rid, 
                                                            sp_type_strains, 
                                                            ani_threshold, 
                                                            Defaults.AF_SP)

            # cluster genomes making sure to select species representative genomes with naming priority
            self.logger.info(' - determining species clusters')
            sp_clusters = self.create_sp_cluters(ani_af, merged_rids_by_priority, rid_to_sp, gids_pass_qc, ani_threshold)

            # write out species clusters
            self.logger.info(' - writting species clusters to file')
            out_file = os.path.join(self.out_dir, f'sp_clusters-ani{ani_threshold}.tsv')
            self.write_sp_clusters(sp_clusters, gid_to_sp, rid_to_sp, out_file)
