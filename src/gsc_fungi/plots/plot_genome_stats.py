
import os
import logging
from pathlib import PurePath
from collections import namedtuple
from dataclasses import dataclass
from typing import Dict, List

from numpy import (median as np_median,
                    mean as np_mean,
                    std as np_std,
                    percentile as np_percentile,
                    ones_like as np_ones_like)
                    
from matplotlib.ticker import FuncFormatter, MaxNLocator

from gsc_fungi.plots.abstract_plot import AbstractPlot


@dataclass
class GenomeMetadata:
    completeness: float
    contamination: float
    genome_quality: float
    genome_size: int
    contig_count: int
    n50_contigs: int
    ambiguous_bases: int


class GenomeStatsHistogram(AbstractPlot):
    """Create histogram of common genomic statistics."""

    def __init__(self, options, rows: int, cols: int):
        """Initialize."""
        
        AbstractPlot.__init__(self, options)
        
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        self.rows = rows
        self.cols = cols

    def plot(self, 
                plot_num: int, 
                data: list, 
                xlabel: str, 
                ylabel: str, 
                txt_position: str = None, 
                center_xticks: bool = False) -> None:
        """Create histogram for statistic."""
            
        self.axis = self.fig.add_subplot(self.rows, self.cols, plot_num)
        
        align = 'mid'
        if center_xticks:
            # not intuative, but setting the alignment to left
            # puts labels in the middle of each bar
            align = 'left' 
                        
        weights = np_ones_like(data)/float(len(data))
        num_bins = min(20, len(set(data))-1)
        counts, bins, patches = self.axis.hist(data, 
                                                bins=num_bins, 
                                                rwidth=0.9, 
                                                weights=weights, 
                                                color='#fdae6b',
                                                align=align)

        self.axis.set_xlabel(xlabel)
        self.axis.set_ylabel(ylabel)
        
        self.axis.xaxis.set_major_locator(MaxNLocator(integer=True))
        self.axis.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.1%}'.format(y)))
        
        # report summary statistics
        stat_txt = f'median = {np_median(data):.1f}\n'
        stat_txt += f'mean = {np_mean(data):.1f}\n'
        stat_txt += f'std = {np_std(data):.1f}'
        if txt_position == 'left':
            self.axis.text(0.05, 0.95, 
                            stat_txt, 
                            transform=self.axis.transAxes,
                            fontsize=self.options.tick_font_size,
                            verticalalignment='top')
        elif txt_position == 'right':
            self.axis.text(0.95, 0.95, 
                            stat_txt, 
                            transform=self.axis.transAxes,
                            fontsize=self.options.tick_font_size,
                            verticalalignment='top',
                            horizontalalignment='right')

        self.prettify(self.axis)
        for loc, spine in self.axis.spines.items():
            if loc in ['right', 'top']:
                spine.set_color('none')

        self.fig.tight_layout(pad=0.1, w_pad=1.0, h_pad=1.0)
        self.draw()
   

class PlotGenomeStats(object):
    """Plot of common genomic statistics."""

    def __init__(self, out_dir: str):
        """Initialization."""
        
        self.log = logging.getLogger('rich')
        self.out_dir = out_dir

    def read_metadata(self, busco_file: str, gb_assembly_file: str) -> Dict[str, GenomeMetadata]:
        """Read metadata for genomes from BUSCO TSV file."""

        # read genome metadata from NCBI assembly summary files
        ambig_bases = {}
        with open(gb_assembly_file) as f:
            f.readline()
            header = f.readline().strip().split('\t')

            genome_size_idx = header.index('genome_size')
            genome_size_ungapped_idx = header.index('genome_size_ungapped')

            for line in f:
                tokens = line.strip().split('\t')

                gid = tokens[0]

                cur_ambig_bases = int(tokens[genome_size_idx]) - int(tokens[genome_size_ungapped_idx])
                assert cur_ambig_bases >= 0
                ambig_bases[gid] = cur_ambig_bases

        # read genome metadata from BUSCO genome quality files
        metadata = {}
        with open(busco_file) as f:
            header = f.readline().strip().split('\t')

            genome_idx = header.index('gid')
            comp_idx = header.index('complete')
            cont_idx = header.index('multi_copy')
            genome_size_idx = header.index('genome_size')
            contig_count_idx = header.index('num_contigs')
            n50_idx = header.index('n50_contigs')
                
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = line_split[genome_idx]
                comp = float(line_split[comp_idx])
                cont = float(line_split[cont_idx])
                quality = comp-5*cont  
                genome_size = int(line_split[genome_size_idx])
                contig_count = int(line_split[contig_count_idx])
                n50 = int(line_split[n50_idx])
                
                metadata[gid] = GenomeMetadata(comp,
                                                cont,
                                                quality,
                                                genome_size,
                                                contig_count,
                                                n50,
                                                ambig_bases[gid])

        return metadata
        
    def get_stat_data(self, 
                        metadata: Dict[str, GenomeMetadata],
                        field: str, 
                        scale_factor: float, 
                        min_value: float,
                        max_value: float) -> List[float]:
        """Get statistic across desired genomes."""
        
        data = []
        for d in metadata.values():
            v = getattr(d, field)/scale_factor
            
            if (min_value and v < min_value) or (max_value and v > max_value):
                continue
            
            data.append(v)
                
        return data
 
    def run(self, busco_file: str, gb_assemfly_file: str) -> None:
        """Create plot of common genomic statistics."""

        # get genome metadata
        self.log.info('Reading genomi statistics from BUSCO and GenBank Assembly files.')
        genome_metadata = self.read_metadata(busco_file, gb_assemfly_file)
        
        # create plot for each genomic statistic
        options = AbstractPlot.Options(width=2, 
                                        height=2, 
                                        label_font_size=7, 
                                        tick_font_size=6, 
                                        dpi=600)
    
        table_data = []
        table_data.append(('', 
                            'Median', 
                            'Mean', 
                            'Std. deviation',
                            'Min.',
                            'Max.',
                            '5th percentile',
                            '95th percentile'))
        
        # plot genome completeness
        self.log.info(f'Creating genome completeness plot.')
        data = self.get_stat_data(genome_metadata, 'completness', 1, 0, 100)
        plot = GenomeStatsHistogram(options, 1, 1)
        plot.plot(1, data, 'Completeness (%)', f'Genomes ({len(data):,})')
        table_data.append(('Completeness', 
                            f'{np_median(data):.1f}',
                            f'{np_mean(data):.1f}',
                            f'{np_std(data):.1f}',
                            f'{min(data):.1f}',
                            f'{max(data):.1f}',
                            f'{np_percentile(data, 5):.1f}',
                            f'{np_percentile(data, 95):.1f}'))

        plot.save_plot(os.path.join(self.out_dir,'qc_completeness.png'), dpi=600)
        plot.save_plot(os.path.join(self.out_dir,'qc_completeness.svg'), dpi=600)

        # plot genome contamination
        self.log.info(f'Creating genome contamination plot.')
        data = self.get_stat_data(genome_metadata, 'contamination', 1, 0, 100)
        plot = GenomeStatsHistogram(options, 1, 1)
        plot.plot(1, data, 'Contamination (%)', f'Genomes ({len(data):,})')
        table_data.append(('Contamination', 
                            f'{np_median(data):.1f}',
                            f'{np_mean(data):.1f}',
                            f'{np_std(data):.1f}',
                            f'{min(data):.1f}',
                            f'{max(data):.1f}',
                            f'{np_percentile(data, 5):.1f}',
                            f'{np_percentile(data, 95):.1f}'))

        plot.save_plot(os.path.join(self.out_dir,'qc_contamination.png'), dpi=600)
        plot.save_plot(os.path.join(self.out_dir,'qc_contamination.svg'), dpi=600)

        # plot genome quality
        self.log.info(f'Creating genome quality plot.')
        data = self.get_stat_data(genome_metadata, 'genome_quality', 1, 0, 100)
        plot = GenomeStatsHistogram(options, 1, 1)
        plot.plot(1, data, 'Genome Quality', f'Genomes ({len(data):,})')
        table_data.append(('Genome Quality', 
                            f'{np_median(data):.1f}',
                            f'{np_mean(data):.1f}',
                            f'{np_std(data):.1f}',
                            f'{min(data):.1f}',
                            f'{max(data):.1f}',
                            f'{np_percentile(data, 5):.1f}',
                            f'{np_percentile(data, 95):.1f}'))

        plot.save_plot(os.path.join(self.out_dir,'qc_genome_quality.png'), dpi=600)
        plot.save_plot(os.path.join(self.out_dir,'qc_genome_quality.svg'), dpi=600)

        # plot genome contig count
        self.log.info(f'Creating number of contigs plot.')
        data = self.get_stat_data(genome_metadata, 'contig_count', 1, 0, None)
        plot = GenomeStatsHistogram(options, 1, 1)
        plot.plot(1, data, 'Contig Count', f'Genomes ({len(data):,})')
        table_data.append(('Contig count', 
                            f'{np_median(data):.1f}',
                            f'{np_mean(data):.1f}',
                            f'{np_std(data):.1f}',
                            f'{min(data):.1f}',
                            f'{max(data):.1f}',
                            f'{np_percentile(data, 5):.1f}',
                            f'{np_percentile(data, 95):.1f}'))

        plot.save_plot(os.path.join(self.out_dir,'qc_contig_count.png'), dpi=600)
        plot.save_plot(os.path.join(self.out_dir,'qc_contig_count.svg'), dpi=600)

        # plot genome N50
        self.log.info(f'Creating N50 plot.')
        data = self.get_stat_data(genome_metadata, 'n50_contigs', 1_000_000, 0, None)
        plot = GenomeStatsHistogram(options, 1, 1)
        plot.plot(1, data, 'N50 contigs (1e6)', f'Genomes ({len(data):,})')
        table_data.append(('N50 contigs (1e6)', 
                            f'{np_median(data):.1f}',
                            f'{np_mean(data):.1f}',
                            f'{np_std(data):.1f}',
                            f'{min(data):.1f}',
                            f'{max(data):.1f}',
                            f'{np_percentile(data, 5):.1f}',
                            f'{np_percentile(data, 95):.1f}'))

        plot.save_plot(os.path.join(self.out_dir,'qc_n50_contigs.png'), dpi=600)
        plot.save_plot(os.path.join(self.out_dir,'qc_n50_contigs.svg'), dpi=600)

        # plot ambiguous bases
        self.log.info(f'Creating ambiguous bases plot.')
        data = self.get_stat_data(genome_metadata, 'ambiguous_bases', 1_000_000, 0, None)
        plot = GenomeStatsHistogram(options, 1, 1)
        plot.plot(1, data, 'Ambiguous bases', f'Genomes ({len(data):,})')
        table_data.append(('Ambiguous_bases', 
                            f'{np_median(data):.1f}',
                            f'{np_mean(data):.1f}',
                            f'{np_std(data):.1f}',
                            f'{min(data):.1f}',
                            f'{max(data):.1f}',
                            f'{np_percentile(data, 5):.1f}',
                            f'{np_percentile(data, 95):.1f}'))

        plot.save_plot(os.path.join(self.out_dir,'qc_ambiguous_bases.png'), dpi=600)
        plot.save_plot(os.path.join(self.out_dir,'qc_ambiguous_bases.svg'), dpi=600)

        # write out table
        self.log.info('Creating QC table.')
        fout = open(os.path.join(self.out_dir, 'qc_table.tsv'), 'w')
        for row in table_data:
            row_str ='\t'.join(row)
            fout.write(f'{row_str}\n')
        fout.close()
