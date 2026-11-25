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
import re
import logging
import subprocess
import shutil
from collections import defaultdict
from dataclasses import dataclass
from typing import Tuple, List, Dict, TypeAlias

from rich import get_console
from rich.progress import (
    Progress,
    TextColumn,
    BarColumn,
    MofNCompleteColumn,
    TaskProgressColumn,
    TimeElapsedColumn,
    TimeRemainingColumn)

from gtdblib.util.shell.execute import check_dependencies

# progress message written to stderr by skani search after every 100 queries
SKANI_QUERIES_PROCESSED_RE = re.compile(r'(\d+) query sequences processed')


@dataclass
class SkaniResult:
    """Metadata for a genome including quality metrics."""

    ani: float
    af_ref: float
    af_query: float


SkaniResults = Dict[str, Dict[str, List]]


class ANIError(Exception):
    """Exception for failed execution of ANI/AF."""
    pass


class Skani():
    """Calculate ANI between genomes using Skani."""

    def __init__(self, cpus: int) -> None:
        """Initialization."""

        check_dependencies(['skani'])

        self.cpus = cpus
        self.log = logging.getLogger('rich')

        self.log.info('Using skani v{}.'.format(Skani.get_version()))

    @staticmethod
    def symmetric_ani_af(ani_af, gid1: str, gid2: str) -> Tuple[float, float]:
        """Calculate symmetric ANI and AF statistics between genomes."""

        if gid1 == gid2:
            return 100.0, 100.0

        if gid1 in ani_af and gid2 in ani_af[gid1]:
            ani, af_r, af_q = ani_af[gid1][gid2]
            return ani, max(af_r, af_q)
        elif gid2 in ani_af and gid1 in ani_af[gid2]:
            ani, af_r, af_q = ani_af[gid2][gid1]
            return ani, max(af_r, af_q)

        return 0.0, 0.0

    def _run_with_query_progress(self, cmd: List[str], num_queries: int) -> Tuple[int, str]:
        """Run skani search, reporting progress from its stderr.

        skani search writes 'N query sequences processed.' to stderr after every
        100 queries, so progress is reported by reading stderr as it is written
        rather than after skani exits. This adds no work for skani. A progress
        bar is shown when writing to a terminal, and progress is logged roughly
        every 10% so it also reaches the log file.

        Returns the skani return code and its stderr.
        """

        stderr_lines = []
        next_log_pct = 10

        progress = Progress(
            TextColumn('   '),
            BarColumn(),
            MofNCompleteColumn(),
            TaskProgressColumn(),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            disable=not get_console().is_terminal,
            transient=True)

        with progress:
            task = progress.add_task('skani search', total=num_queries)

            proc = subprocess.Popen(cmd,
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.PIPE,
                                    encoding='utf-8')
            for line in proc.stderr:
                stderr_lines.append(line)

                m = SKANI_QUERIES_PROCESSED_RE.search(line)
                if not m:
                    continue

                processed = int(m.group(1))
                progress.update(task, completed=processed)

                pct = 100 * processed / num_queries
                if pct >= next_log_pct:
                    self.log.info(f'   - processed {processed:,} of {num_queries:,} query genomes ({pct:.0f}%)')
                    next_log_pct = 10 * (int(pct) // 10 + 1)

            proc.wait()
            progress.update(task, completed=num_queries)

        return proc.returncode, ''.join(stderr_lines)

    @staticmethod
    def get_version() -> str:
        """Returns the version of skani on the system path.

        Returns
        -------
        str
            The string containing the skani version.
        """
        try:
            proc = subprocess.Popen(['skani', '-V'], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, encoding='utf-8')
            stdout, stderr = proc.communicate()
            version = re.search(r'skani (.+)', stdout)
            return version.group(1)
        except Exception as e:
            print(e)
            return 'unknown'

    @staticmethod
    def dist(qid: str, rid: str, q_gf: str, r_gf: str, preset: str = None) -> SkaniResult:
        """CalculateANI between a pair of genomes."""

        # run skani and write results to stdout
        try:
            cmd = ['skani', 'dist',
                    '-q', q_gf,
                    '-r', r_gf,
                    '-o', '/dev/stdout']

            if preset is not None:
                cmd += [preset]

            proc = subprocess.Popen(cmd,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        encoding='utf-8')
            stdout, stderr = proc.communicate()

            if proc.returncode != 0:  # skani returned an error code
                print(stderr)
                raise ANIError(f"skani exited with code {proc.returncode}.")
        except Exception as e:
            print(e)
            raise

        result_lines = stdout.splitlines()
        if len(result_lines) == 2:
            tokens = result_lines[1].split('\t')
            ani = float(tokens[2])
            af_r = float(tokens[3])
            af_q = float(tokens[4])
            result = SkaniResult(ani, af_r, af_q)
        else:
            # genomes too divergent to determine ANI and AF
            # with skani so default to zeros
            result = SkaniResult(0.0, 0.0, 0.0)

        return result

    def triangle(self, genome_paths: Dict[str, str], 
                        output_dir: str, 
                        preset: str,
                        min_af: float = 50,
                        min_sketch_ani:float = 85) -> SkaniResults:
        """Calculate ANI/AF between genome pairs in lower triangle.
        
        This method is fast, but has large memory requirements if large
        numbers of genomes need to be compared.
        """

        os.makedirs(output_dir, exist_ok=True)

        self.log.info(f'Calculating ANI between {len(genome_paths):,} genomes in lower triangle ({preset}; min_af = {min_af}; s = {min_sketch_ani}):')

        try:
            # create file with path to genomes
            self.log.info(' - creating file with path to genomes')
            genome_path_file = os.path.join(output_dir, 'skani_genome_paths.tsv')
            fout = open(genome_path_file, 'w')
            for gp in genome_paths.values():
                fout.write(f'{gp}\n')
            fout.close()
            
            # compute ANI/AF between genome pairs in lower triangle
            self.log.info(' - calculating ANI and AF between genomes')
            results_file = os.path.join(output_dir, 'skani_dist.tsv')
            cmd = ['skani', 'triangle',
                    '-t', str(self.cpus), 
                    '-l', genome_path_file,
                    '-o', results_file,
                    '--sparse',
                    '--min-af', str(min_af),
                    '-s', str(min_sketch_ani)]

            if preset is not None:
                cmd += [preset]

            proc = subprocess.Popen(cmd,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    encoding='utf-8')
            _stdout, stderr = proc.communicate()
            if proc.returncode != 0:  # skani returned an error code
                print(stderr)
                raise ANIError(f"skani triangle exited with code {proc.returncode}.")

        except Exception as e:
            print(e)
            raise

        # get map from genome path to genome ID
        gp_to_gid_map = {}
        for gid, gp in genome_paths.items():
            gp_to_gid_map[gp] = gid

        # parse skani results
        self.log.info(' - parsing results')
        ani_af = defaultdict(dict)
        with open(results_file) as f:
            header = f.readline().strip().split('\t')

            ref_file_idx = header.index('Ref_file')
            query_file_idx = header.index('Query_file')
            ani_idx = header.index('ANI')
            af_r_idx = header.index('Align_fraction_ref')
            af_q_idx = header.index('Align_fraction_query')

            for line in f:
                tokens = line.strip().split('\t')

                rid = gp_to_gid_map[tokens[ref_file_idx]]
                qid = gp_to_gid_map[tokens[query_file_idx]]
                ani = float(tokens[ani_idx])
                af_r = float(tokens[af_r_idx])
                af_q = float(tokens[af_q_idx])

                ani_af[qid][rid] = (ani, af_r, af_q)

        return ani_af

    def search(self, 
                query_paths: Dict[str, str], 
                ref_paths: Dict[str, str],
                output_dir: str, 
                preset: str,
                min_af: float = 50,
                min_sketch_ani: float = 85,
                remove_sketches: bool = False) -> SkaniResults:
        """Calculate skani between query and reference genomes.
        
        Since skani is symmetric each pair of genomes is only processed once. That is,
        genome A and genome B will be processed with skani as (A,B) or (B,A), but not
        both since they will give the same result. Representative genomes are first
        sketched and then the query genomes searched against this DB. This is slower
        than using "dist" but keeps memory usage much more reasonable.
        """

        os.makedirs(output_dir, exist_ok=True)

        log_str = f'Calculating ANI between {len(query_paths):,} query and {len(ref_paths):,} reference genomes'
        log_str += f' ({preset}; min_af = {min_af}; s = {min_sketch_ani}):'
        self.log.info(log_str)

        try:
            # create file with path to reference genomes
            self.log.info(' - creating file with path to reference genomes')
            rep_path_file = os.path.join(output_dir, 'skani_ref_genome_paths.tsv')
            fout = open(rep_path_file, 'w')
            for gp in ref_paths.values():
                fout.write(f'{gp}\n')
            fout.close()

            # create file with path to query genomes
            self.log.info(' - creating file with path to query genomes')
            query_path_file = os.path.join(output_dir, 'skani_query_genome_paths.tsv')
            fout = open(query_path_file, 'w')
            for gid, gp in query_paths.items():
                if gid in ref_paths and ref_paths[gid] != gp:
                    self.log.error(f'Genome ID {gid} has different genomic FASTA files for query and reference.')
                    sys.exit(1)
                fout.write(f'{gp}\n')
            fout.close()
            
            # sketch the reference genomes
            self.log.info(' - sketching reference genomes')
            sketch_dir = os.path.join(output_dir, 'skani_ref_sketches')
            cmd = ['skani', 'sketch',
                    '-t', str(self.cpus), 
                    '-l', rep_path_file,
                    '-o', sketch_dir]

            if preset is not None:
                cmd += [preset]

            proc = subprocess.Popen(cmd,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        encoding='utf-8')
            _stdout, stderr = proc.communicate()
            if proc.returncode != 0:  # skani returned an error code
                print(stderr)
                raise ANIError(f"skani sketch exited with code {proc.returncode}.")
            
            # search query genomes against reference genome sketches
            self.log.info(' - calculating ANI and AF between reference and query genomes')
            results_file = os.path.join(output_dir, 'skani_query_vs_ref.tsv')
            cmd = ['skani', 'search',
                    '-t', str(self.cpus),
                    '-d', sketch_dir,
                    '--ql', query_path_file,
                    '-o', results_file,
                    '--min-af', str(min_af),
                    '-s', str(min_sketch_ani),
                    '--short-header']
            returncode, stderr = self._run_with_query_progress(cmd, len(query_paths))
            if returncode != 0:  # skani returned an error code
                print(stderr)
                raise ANIError(f"skani search exited with code {returncode}.")
        except Exception as e:
            print(e)
            raise

        # get map from genome path to genome ID
        gp_to_gid_map = {}
        for gid, gp in query_paths.items():
            gp_to_gid_map[gp] = gid
        for gid, gp in ref_paths.items():
            gp_to_gid_map[gp] = gid

        # parse skani results
        self.log.info(' - parsing skani results')
        ani_af = defaultdict(dict)
        with open(results_file) as f:
            header = f.readline().strip().split('\t')

            ref_file_idx = header.index('Ref_file')
            query_file_idx = header.index('Query_file')
            ani_idx = header.index('ANI')
            af_r_idx = header.index('Align_fraction_ref')
            af_q_idx = header.index('Align_fraction_query')

            for line in f:
                tokens = line.strip().split('\t')

                rid = gp_to_gid_map[tokens[ref_file_idx]]
                qid = gp_to_gid_map[tokens[query_file_idx]]
                ani = float(tokens[ani_idx])
                af_r = float(tokens[af_r_idx])
                af_q = float(tokens[af_q_idx])

                ani_af[qid][rid] = (ani, af_r, af_q)

        # remove sketch directory if requested
        if remove_sketches:
            shutil.rmtree(sketch_dir)

        return ani_af
