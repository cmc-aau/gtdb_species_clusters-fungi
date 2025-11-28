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
import csv
import logging

from gtdblib.util.shell.gtdbshutil import check_file_exists, make_sure_path_exists

from gsc_fungi.qc_genomes import QcGenomes, QcCriteria
from gsc_fungi.select_sp_rep import SelectSpeciesRepresentative
from gsc_fungi.species_clusters import SpeciesClusters

from gsc_fungi.taxa_by_rank import TaxaByRank
from gsc_fungi.plots.plot_genome_stats import PlotGenomeStats


if sys.version_info >= (3, 11):
    import tomllib

    def parse_toml_file(toml_file: str):
        with open(toml_file, 'rb') as f:
            return tomllib.load(f)
else:
    import toml as tomllib

    def parse_toml_file(toml_file: str):
        with open(toml_file, 'r') as f:
            return tomllib.load(f)


class OptionsParser():
    """Validate and execute command-line interface."""

    def __init__(self):
        """Initialization"""

        self.log = logging.getLogger('rich')

    def args_abs_path(self, params, argument: str) -> str:
        """Get absolute path for argument."""

        if argument in params['data_files']:
            return os.path.abspath(os.path.join(params['root_dir'], params['data_files'][argument]))
        elif argument in params['nomenclature_files']:
            return os.path.abspath(os.path.join(params['root_dir'], params['nomenclature_files'][argument]))
        elif argument in params['ledgers']:
            return os.path.abspath(os.path.join(params['root_dir'], params['ledgers'][argument]))
        elif argument in params['output_dirs']:
            return os.path.abspath(os.path.join(params['root_dir'], params['output_dirs'][argument]))
        else:
            self.log.error(f"Unknown argument: {argument}")
            sys.exit(-1)

    def qc_genomes(self, args) -> None:
        """Quality check genomes."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        check_file_exists(self.args_abs_path(params, 'genome_metadata_file'))
        check_file_exists(self.args_abs_path(params, 'busco_file'))
        check_file_exists(self.args_abs_path(params, 'ncbi_assembly_summary_genbank_file'))
        check_file_exists(self.args_abs_path(params, 'ncbi_taxonomy_file'))

        out_dir = self.args_abs_path(params, args.subparser_name)
        make_sure_path_exists(out_dir)
        self.log.info(f"Output directory: {out_dir}")

        qc_criteria = QcCriteria(
            args.min_comp,
            args.max_cont,
            args.min_quality,
            args.max_contigs,
            args.min_N50,
            args.max_ambiguous_perc)

        p = QcGenomes(out_dir)
        p.run(self.args_abs_path(params, 'genome_metadata_file'),
                self.args_abs_path(params, 'busco_file'),
                self.args_abs_path(params, 'ncbi_assembly_summary_genbank_file'),
                self.args_abs_path(params, 'ncbi_taxonomy_file'),
                qc_criteria)

    def cdn_select_sp_reps(self, args) -> None:
        """Select species representative genomes."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        qc_dir = self.args_abs_path(params, 'cdn_qc_genomes')
        qc_pass_file = os.path.join(qc_dir, 'qc_passed.tsv')

        check_file_exists(self.args_abs_path(params, 'ncbi_assembly_summary_genbank_file'))
        check_file_exists(self.args_abs_path(params, 'ncbi_taxonomy_file'))
        check_file_exists(self.args_abs_path(params, 'cur_genome_path_file'))

        out_dir = self.args_abs_path(params, args.subparser_name)
        make_sure_path_exists(out_dir)
        self.log.info(f"Output directory: {out_dir}")

        p = SelectSpeciesRepresentative(out_dir)
        p.run(qc_pass_file,
              self.args_abs_path(params, 'ncbi_assembly_summary_genbank_file'),
              self.args_abs_path(params, 'ncbi_taxonomy_file'),
              self.args_abs_path(params, 'cur_genome_path_file'))

    def cdn_sp_clusters(self, args) -> None:
        """Create ANI-based species clusters."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        qc_dir = self.args_abs_path(params, 'cdn_qc_genomes')
        qc_pass_file = os.path.join(qc_dir, 'qc_passed.tsv')

        select_sp_reps_dir = self.args_abs_path(params, 'cdn_select_sp_reps')
        sp_rep_file = os.path.join(select_sp_reps_dir, 'sp_reps.tsv')

        check_file_exists(self.args_abs_path(params, 'ncbi_taxonomy_file'))
        check_file_exists(self.args_abs_path(params, 'cur_genome_path_file'))

        out_dir = self.args_abs_path(params, args.subparser_name)
        make_sure_path_exists(out_dir)
        self.log.info(f"Output directory: {out_dir}")

        p = SpeciesClusters(args.cpus, out_dir)
        p.run(qc_pass_file,
                sp_rep_file,
                self.args_abs_path(params, 'ncbi_taxonomy_file'),
                self.args_abs_path(params, 'cur_genome_path_file'))

    def merge_sp_naming_priority(self, args) -> None:
        """Determine name with priority for merged species."""

        check_file_exists(args.input_params_file)
        params = parse_toml_file(args.input_params_file)

        check_file_exists(self.args_abs_path(params, 'mycobank_file'))
        check_file_exists(self.args_abs_path(params, 'mycobank_sanctioned_names_persoon'))
        check_file_exists(self.args_abs_path(params, 'mycobank_sanctioned_names_ef'))

        # TBD


    def taxa_by_rank(self, args) -> None:
        """Get number of genomes and taxa at each rank."""

        check_file_exists(args.ncbi_taxonomy_file)
        check_file_exists(args.assembly_summary_file)
        os.makedirs(args.out_dir, exist_ok=True)

        p = TaxaByRank(args.out_dir)
        p.run(args.ncbi_taxonomy_file, args.assembly_summary_file)

    def plot_genome_stats(self, args) -> None:
        """Plot common genomic statistics."""

        check_file_exists(args.busco_file)
        check_file_exists(args.assembly_summary_file)
        os.makedirs(args.out_dir, exist_ok=True)

        p = PlotGenomeStats(args.out_dir)
        p.run(args.busco_file, args.assembly_summary_file)

    def run(self, args) -> None:
        """Parse user arguments and call the correct pipeline(s)"""

        if args.subparser_name == 'cdn_qc_genomes':
            self.qc_genomes(args)
        elif args.subparser_name == 'cdn_select_sp_reps':
            self.cdn_select_sp_reps(args)
        elif args.subparser_name == 'cdn_sp_clusters':
            self.cdn_sp_clusters(args)
        elif args.subparser_name == 'u_qc_genomes':
            self.qc_genomes(args)
        elif args.subparser_name == 'taxa_by_rank':
            self.taxa_by_rank(args)
        elif args.subparser_name == 'plot_genome_stats':
            self.plot_genome_stats(args)
        else:
            self.log.error(
                f'Unknown gtdb_species_clusters command: {args.subparser_name}\n')
            sys.exit()
