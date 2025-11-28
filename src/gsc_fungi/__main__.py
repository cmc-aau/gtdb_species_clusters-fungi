#!/usr/bin/env python3

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

__author__ = "Donovan Parks"
__copyright__ = "Copyright 2025"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

import os
import sys
import argparse
import logging
import ntpath
import time

from gtdblib.cli.custom_help_formatter import CustomHelpFormatter

from gsc_fungi import __version__
from gsc_fungi.main import OptionsParser, parse_toml_file
from gsc_fungi import defaults as Defaults


def print_help():
    """Help menu."""

    print('')
    print('                ...::: GTDB Species Clusters: Fungi v' + __version__ + ' :::...''')
    print('''\

    Create de novo species clusters:
     cdn_qc_genomes     -> Quality check genomes
     cdn_select_sp_reps -> Select species representative genomes
     cdn_sp_clusters    -> Create ANI-based species clusters

    Update species clusters:
      u_qc_genomes -> Quality check genomes

    Post-manual curation
      ...

    Inspect species clusters:
      ...

    Inspect genomes:
      ...

    Utilities:
      taxa_by_rank      -> Get number of genomes and taxa at each rank
      plot_genome_stats -> Plot gcommon genomic statistics

  Use: gsc_fungi <command> -h for command specific help.

  Feature requests or bug reports can be sent to Donovan Parks (donovan.parks@gmail.com)
    or posted on GitHub (https://github.com/cmc-aau/gtdb_species_clusters-fungi).
    ''')


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')

    # quality check genomes
    cdn_qc_genomes_parser = subparsers.add_parser('cdn_qc_genomes',
                                                formatter_class=CustomHelpFormatter,
                                                description='Quality check genomes.')
    cdn_qc_genomes_parser.add_argument('input_params_file', help='TOML file indicating input parameters')
    cdn_qc_genomes_parser.add_argument(
        '--min_comp', help='minimum completeness of genome [0, 100]', type=float, default=Defaults.QC_MIN_COMP)
    cdn_qc_genomes_parser.add_argument(
        '--max_cont', help='maximum contamination of genome [0, 100]', type=float, default=Defaults.QC_MAX_CONT)
    cdn_qc_genomes_parser.add_argument(
        '--min_quality', help='minimum genome quality (comp - 5*cont) [0, 100]', type=float, default=Defaults.QC_MIN_QUALITY)
    cdn_qc_genomes_parser.add_argument(
        '--max_contigs', help='maximum number of contigs', type=int, default=Defaults.QC_MAX_CONTIGS)
    cdn_qc_genomes_parser.add_argument(
        '--min_N50', help='minimum N50 of scaffolds', type=int, default=Defaults.QC_MIN_N50)
    cdn_qc_genomes_parser.add_argument(
        '--max_ambiguous_perc', help='maximum percentage of ambiguous bases', type=float, default=Defaults.QC_MAX_AMBIGUOUS_PERC)
    cdn_qc_genomes_parser.add_argument('--silent', help="suppress output", action='store_true')

    # quality check genomes
    cdn_select_sp_reps_parser = subparsers.add_parser('cdn_select_sp_reps',
                                                formatter_class=CustomHelpFormatter,
                                                description='Select species representative genomes.')
    cdn_select_sp_reps_parser.add_argument('input_params_file', help='TOML file indicating input parameters')
    cdn_select_sp_reps_parser.add_argument('--silent', help="suppress output", action='store_true')

    # quality check genomes
    cdn_sp_clusters_parser = subparsers.add_parser('cdn_sp_clusters',
                                                formatter_class=CustomHelpFormatter,
                                                description='Create ANI-based species clusters.')
    cdn_sp_clusters_parser.add_argument('input_params_file', help='TOML file indicating input parameters')
    cdn_sp_clusters_parser.add_argument('-c', '--cpus', help='number of CPUs', default=1)
    cdn_sp_clusters_parser.add_argument('--silent', help="suppress output", action='store_true')

    # quality check genomes
    u_qc_genomes_parser = subparsers.add_parser('u_qc_genomes',
                                                formatter_class=CustomHelpFormatter,
                                                description='Quality check genomes.')
    u_qc_genomes_parser.add_argument('input_params_file', help='TOML file indicating input parameters')
    u_qc_genomes_parser.add_argument(
        '--min_comp', help='minimum completeness of genome [0, 100]', type=float, default=Defaults.QC_MIN_COMP)
    u_qc_genomes_parser.add_argument(
        '--max_cont', help='maximum contamination of genome [0, 100]', type=float, default=Defaults.QC_MAX_CONT)
    u_qc_genomes_parser.add_argument(
        '--min_quality', help='minimum genome quality (comp - 5*cont) [0, 100]', type=float, default=Defaults.QC_MIN_QUALITY)
    u_qc_genomes_parser.add_argument(
        '--max_contigs', help='maximum number of contigs', type=int, default=Defaults.QC_MAX_CONTIGS)
    u_qc_genomes_parser.add_argument(
        '--min_N50', help='minimum N50 of scaffolds', type=int, default=Defaults.QC_MIN_N50)
    u_qc_genomes_parser.add_argument(
        '--max_ambiguous_perc', help='maximum percentage of ambiguous bases', type=float, default=Defaults.QC_MAX_AMBIGUOUS_PERC)
    u_qc_genomes_parser.add_argument('--silent', help="suppress output", action='store_true')

    # get number of genomes and taxa at each rank
    taxa_by_rank_parser = subparsers.add_parser('taxa_by_rank',
                                                add_help=False,
                                                formatter_class=CustomHelpFormatter,
                                                description='Get number of genomes and taxa at each rank.')
    req_parser = taxa_by_rank_parser.add_argument_group('Required')
    req_parser.add_argument('-i', '--ncbi_taxonomy_file', help='standardized NCBI taxonomy file with subranks', required=True)
    req_parser.add_argument('-a', '--assembly_summary_file', help='GenBank assembly summary file from NCBI', required=True)
    req_parser.add_argument('-o', '--out_dir', help='output directory', required=True)

    opt_parser = taxa_by_rank_parser.add_argument_group('Optional')
    opt_parser.add_argument('--silent', help="suppress output", action='store_true')
    opt_parser.add_argument('-h', '--help', action='help', help='show this help message and exit')

    # get number of genomes and taxa at each rank
    plot_genome_stats_parser = subparsers.add_parser('plot_genome_stats',
                                                add_help=False,  
                                                formatter_class=CustomHelpFormatter,
                                                description='Plot common genomic statistics.')
    req_parser = plot_genome_stats_parser.add_argument_group('Required')
    req_parser.add_argument('-b', '--busco_file', help='file with BUSCO genome quality estimates', required=True)
    req_parser.add_argument('-a', '--assembly_summary_file', help='GenBank assembly summary file from NCBI', required=True)
    req_parser.add_argument('-o', '--out_dir', help='output directory', required=True)

    opt_parser = plot_genome_stats_parser.add_argument_group('Optional')
    opt_parser.add_argument('--silent', help="suppress output", action='store_true')
    opt_parser.add_argument('-h', '--help', action='help', help='show this help message and exit')

    # get and check options
    args = None
    if (len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv == '--help'):
        print_help()
        sys.exit(0)
    else:
        args = parser.parse_args()

    # parse input parameters file
    params = None
    if hasattr(args, 'input_params_file'):
        params = parse_toml_file(args.input_params_file)

    # get output directory
    output_dir = None
    if params and args.subparser_name in params['output_dirs']:
        output_dir = os.path.abspath(os.path.join(params['root_dir'], params['output_dirs'][args.subparser_name]))
    elif hasattr(args, 'output_dir'):
        output_dir = args.output_dir

    # setup logging
    log = logging.getLogger("rich")
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        log_file = os.path.join(output_dir, 'gsc_fungi.log')
        file_logger = logging.FileHandler(log_file, 'a')
        file_logger.setFormatter(logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                                   datefmt="%Y-%m-%d %H:%M:%S"))
        log.addHandler(file_logger)

    log.info('GTDB Species Clusters: Fungi v%s' % __version__)
    log.info(ntpath.basename(sys.argv[0]) + ' ' + ' '.join(sys.argv[1:]))

    # run specified command
    try:
        start = time.time()
        OptionsParser().run(args)
        end = time.time()

        elapsed_time = (end - start) / 60
        log.info(f'Elapsed time (min): {elapsed_time:.1f}')
        log.info('Done.')
    except SystemExit:
        print("Controlled exit resulting from an unrecoverable error or warning.")
        raise
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise
