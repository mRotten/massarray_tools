# ma_tools.py

"""
Tools for working with MassARRAY xml result files, and some limited
methods for designing new SBE probes around existing SBE probes.
"""

# Want List:
#   - Fully-functional GUI (maybe not worth it).
#   - Data grapher: user-defined sample order.
#   - Define a new minimum snr, peak uncertainty for calling + assays.
#   - Calculate new well calls with peak uncertainties and snr.
#   - Method to average replicate signals.
#   - Error bars for graphs (when experiment has replicates).
#   - Normalized height - height of every peak relative to control
#       oligo peak heights.
#   - Ability to plot control signals (e.g. W1_SAP_QC).

from random import randrange
import itertools as itls
from Bio import SeqIO as sio
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tkinter
import tkinter.filedialog as fd
from lxml import etree
import os
import re


class Control_DNA:

    def __init__(self, path='C:/Users/25676/AppData/Local/Programs/' \
                 'Python/Python36/control_dna_library_v01.txt'):
        self.__read_csv(path)


    def __read_csv(self, path):
        self.data = pd.read_csv(path, sep="\t")


    def get_known_mutations(self, gene, codon):

        """
        Input gene and codon, returns c-dot (string variable).
        """

        a = self.data[self.data['gene'] == gene]
        b = a[a['codon_num'] == codon]
        if len(b) > 0:
            return b.cds_mutation.iloc[0]
        else:
            return "wild-type"


class MA4_Experiment:

    def __init__(self, xml_path):

        self.xml_path = xml_path

        self.sample_groups = None
        self.assays = pd.DataFrame()
        self.records = pd.DataFrame()
        self.spectra = pd.DataFrame()
        self.qc_recs = pd.DataFrame()

        self.flt_data = ['area',
                         'areaUncert',
                         'callScore',
                         'frequency',
                         'frequencyUncert',
                         'height', 'mass',
                         'maxShift',
                         'noiseStdDev',
                         'peakScore',
                         'resolution',
                         'snr']

        self.str_data = ['assayId',
                         'genotypeId',
                         'peakType',
                         'sampleId',
                         'wellPosition',
                         'description',
                         'entryOperator',
                         'direction',
                         'group_id',
                         'id', 'pcr1',
                         'pcr2',
                         'project_id',
                         'terminators',
                         'oligo_id',
                         'sequence',
                         'type',
                         'well-position']

        self.int_data = ['assayPk',
                         'calibration',
                         'callPk',
                         'rasters',
                         'samplePk',
                         'spectraPk',
                         'status',
                         'use-call',
                         'use-well',
                         'wellsPk',
                         'pk',
                         'plex_id',
                         'some-pk',
                         'well=pk']

        self.__parse_records()
        self.__parse_assays()
        self.__annotate_from_assay()


    ########### PUBLIC DEFS FOR GETTING STUFF ########################


    def get_sample_name(self, plate_loc):
        
        """
        From plate location (A01-H12), returns sample name.

        Plate location must be a 3-character string of (A-H)
        followed by ints 1-12, padded to two characters (e.g.
        01, 03, 12, etc.).
        """
        
        a = self.records[self.records['wellPosition'] == plate_loc]
        return list(a.sampleId)[0]


    def get_sample_names(self):
        """
        Returns a list of all unique sample names in the experiment.
        """
        return list(self.records.sampleId.unique())


    def get_sample_groups(self):
        """
        Returns a list of sample groups. Each sample group is a list comprising
        the sample names of all group members.
        """        
        return self.sample_groups


    def get_probe_names(self):
        """
        Returns a list of all unique probe names in the assay.
        """        
        return list(self.records.assayId.unique())


    def get_all_records(self, qc_records=False):
        
        """
        Returns the full dataset as a pandas dataframe. Equivalent to
        calling "get_subset" with no arguments.
        """
        if qc_records == True:
            return pd.concat([self.records, self.qc_recs])
        else:
            return self.records


    def export_data(self, include_spectra=False):
        """
        Exports assay data and annotated record data to a csv file in the
        same directory as the source data.
        Setting 'include_spectra' to 'True' will parse spectra and save
        them in a third csv file.
        """
        
        path = os.path.dirname(self.xml_path)
        self.records.to_csv(path + "/records.csv")
        self.assays.to_csv(path + "/assays.csv")
        print("All records and assays have been saved to:")
        print(path)

        if include_spectra:

            self.__parse_spectra()
            self.spectra.to_csv(path + "/spectra.csv")

            # Clearing spectra from memory.
            self.spectra = pd.DataFrame()

            print("All spectra have been saved to:")
            print(path)


    def get_subset(self, probe_name=None, sample_list=None, \
                   sample_group=None):
        """
        Will return the full dataset if no arguments are given.
        """
        
        a = self.records.copy()
        
        if type(sample_group) == int:
            sg_name = 'sg_' + str(sample_group)
            a = a[a[sg_name] == True]

        if probe_name:
            a = a[a['assayId'] == probe_name]

        if sample_list:
            a = a[a['sampleId'].isin(sample_list)]

        return a


    ########### PUBLIC DEFS FOR SETTING STUFF ########################


    def define_sample_groups(self, sample_groups):

        """
        Each sample group is a list list/tuple of sample
        names that belong to that sample group.
        Membership in one does not preclude membership in
        another.
        """

        for sg_num, sg in enumerate(sample_groups):
            self.__add_sg(sg_num, sg)
            self.sample_groups = sample_groups


    def rename_sample(self, new_name, plate_locations):
        """
        Given a list of plate locations, assigns "new_name" as the
        sample name.

        Plate locations must be a 3-character string comprising one
        letter (A-H) followed by an interger (1-12) padded to 2
        characters.
        """

        for i, row in self.records.iterrows():
            if row['wellPosition'] in plate_locations:
                self.records.loc[i,'sampleId'] = new_name


    ########### PRIVATE DEFS FOR PARSING FILES #######################


    def __parse_records(self):

        tav = etree.parse(self.xml_path).getroot().getchildren()[0]
        for element in tav:
            if element.tag == 'all-records':
                a = element
                break

        header = list(a[0].attrib.keys())
        b, qc_recs = [], []

        # regex to recognize & ignore a qc record:
        qc_oligo_name = r'\bW[1-4]_[PSE][CAX][RPT]_QC\b'

        for el in a:
            c = [self.__convert(x, el.attrib[x]) for x in header]
            if re.match(qc_oligo_name, el.attrib['assayId']):
                qc_recs.append(c)
            else:
                b.append(c)

        for i in range(len(header)):
            d = [x[i] for x in b]
            e = [x[i] for x in qc_recs]
            self.records[header[i]] = pd.Series(d)
            self.qc_recs[header[i]] = pd.Series(e)


    def __parse_assays(self):

        # Annotate self.records as assays are being parsed?
        #   - Not yet implemented, but may not be any faster.

        tav = etree.parse(self.xml_path).getroot().getchildren()[0]
        for element in tav:
            if element.tag == 'all-assays':
                a = element
                break

        assay_header = a[0].keys()
        peak_header = a[0][0].keys()
        data_lists = []

        for assay in a:
            b = [assay.attrib[x] for x in assay_header]

            for peak in assay:
                c = [self.__convert(x, peak.attrib[x]) for x in peak_header]
                data_lists.append( b + c )

        assay_header[assay_header.index('id')] = 'assay_id'
        peak_header[peak_header.index('id')] = 'analyte_call'

        column_names = assay_header + peak_header

        for i, c_name in enumerate(column_names):
            self.assays[c_name] = [x[i] for x in data_lists]


    def __parse_spectra(self):

        tav = etree.parse(self.xml_path).getroot().getchildren()[0]
        for element in tav:
            if element.tag == 'all-spectra':
                all_spec = element
                break

        for spec in all_spec:
            plate_pos = spec.attrib['well-position']
            if plate_pos != 'CALIBRATION':
                df1 = self.records[self.records['wellPosition'] == plate_pos]
                s_name, plex_id = df1.sampleId.iloc[0], df1.plex_id.iloc[0]
                data = [float(x) for x in spec.text.split(',')]
                spec_ser = pd.Series(data)
                self.spectra[s_name + '_well_' + plex_id] = spec_ser


    ########### PRIVATE DEFS FOR ANNOTATING DFs ######################


    def __annotate_from_assay(self):

        pos_neg = []
        call_contribs = []
        plex_ids = []

        # The following nested loop could be faster.

        for i, row1 in self.records.iterrows():

            a = self.assays[self.assays['assay_id'] == row1['assayId']]
            b, c = [], []

            for j, row2 in a.iterrows():

                if self.__is_same_mass(row1['mass'], row2['mass']):
                    b.append(row2['analyte_call'])
                    c.append(row2['oligo_id'])

            # Add to list of plex_ids.

            plex_ids.append(row2['plex_id'])

            # Add to list of pos/neg calls.

            if row1['genotypeId'] in b:
                pos_neg.append('pos')
            else:
                pos_neg.append('neg')

            # Add to list of call contributions.

            if c[0][:3] in ('PCR', 'SAP', 'MUT', 'EXT'):
                call_contribs.append("")
            elif c[0][:3] == "UEP":
                call_contribs.append("UEP")
            else:
                call_contribs.append(c[0])

        self.records['pos_neg'] = pd.Series(pos_neg, index=self.records.index)
        self.records['call_if_pos'] = pd.Series(call_contribs,
                                                 index=self.records.index)
        self.records['plex_id'] = pd.Series(plex_ids, index=self.records.index)


    def __add_sg(self, sg_num, sg):

        memberships = []
        name = 'sg_' + str(sg_num)

        for i, row1 in self.records.iterrows():
            if row1['sampleId'] in sg:
                memberships.append(1)
            else:
                memberships.append(0)
        self.records[name] = pd.Series(memberships, \
                                       index=self.records.index)


    ########### PRIVATE UTILITIES ####################################


    def __convert(self, data_label, data):

        if data_label in self.flt_data:
            return float(data)

        elif data_label in self.int_data:
            return int(data)

        elif data_label in self.str_data:
            return str(data)


    def __is_same_mass(self, m1, m2):
        
        return m1 >= (m2 - 0.5) and m1 <= (m2 + 0.5)


class MA_Data_Grapher:

    """
    How to use this class:
    1) Instantiate the class with no arguments.
    2) Get a one-probe, one-sample-group dataframe from MA_Experiment
        with MA_Experiment.get_subset().
    3) Pass the dataframe from step 2 to MA_Data_Grapher.new_dataset().
    4) Figure out what title your graph should have.
    5) Decide what parameter to graph (e.g. 'height', 'snr', etc.).
    6) Pass parameters from steps 4-5 to function MA_Data_Grapher.graph().
    7) Figure out where you want to save the image (target directory).
    8) Pass the target directory to "save_graph" (must include file
        name and extension.
    """


    def __init__(self):

        self.graph_data = None
        self.current_probe = None
        self.current_sample_group = None
        self.sample_list = None
        self.fig = None
        
        self.uep_color = 'orange'
        self.neg_color = 'y'
        self.pos_color = 'g'

                        # analyte count:  2     3       4       5
                        #----------------------------------------------
        self.chart_params = {'cols':     (1,    1,      2,      2),
                             'rows':     (2,    3,      2,      3),
                             'width':    (12,   12,     24,     24),
                             'height':   (18,   27,     18,     27),
                             'b_marg':   (0.2,  0.15,   0.2,    0.15)}

        self.font = {'weight' : 'bold',
                     'size'   : 24}


    ############ PUBLIC DEFs FOR SETTING STUFF #######################


    def new_dataset(self, data):
        """
        Arg 'data' should come from 'get_subset' function in "MA_Experiment".
        Assumes 'data' only contains info from one probe.
        """

        if self.fig != None:
            plt.close(self.fig)
        self.graph_data = data
        self.current_probe = None
        self.current_sample_group = None
        self.__map_colors()


    def set_colors(self, neg='y', pos='g', uep='orange'):
        """
        Determines what color graph bars will be. Positive calls will be
        'pos' (default is green), negative calls will be 'neg' (default
        is yellow) and UEP bars will be 'uep' (default is orange).
        """

        self.uep_color = uep
        self.neg_color = neg
        self.pos_color = pos


    ############ PUBLIC DEFs FOR GRAPHING ############################


    def graph(self, graph_title, y_axis='height'):
        """
        Given a graph title and y_axis ('height', 'snr', etc.), function
        will plot the data passed to 'new_dataset'.

        Graph title can be any string.
        """

        self.y_axis = y_axis

        plt.rc('font', **self.font)
        data, colormap = self.__format_for_graphing()

        c, r, h, w = self.__get_chart_dimensions(len(data.columns))
        self.fig, ax_obj = plt.subplots(nrows=r, ncols=c, sharey=True)

        plt.suptitle(graph_title)

        axes = self.__get_axis_list(c, r, ax_obj)

        for i, col in enumerate(data.columns):
            data[col].plot(kind='bar', ax=axes[i], figsize=(w, h), \
                           title=col, color=colormap[i])
            axes[i].set_ylabel("Peak " + y_axis.capitalize())

        self.__clean_up_plots(len(data.columns), ax_obj)

        bottom_margin = self.chart_params['b_marg'][len(data.columns)-2]
        
        plt.gcf().subplots_adjust(bottom=bottom_margin)


    def save_graph(self, target_dir):
        """
        Function will save the graph data to 'target_dir' directory.
        Arg 'target_dir' must include the file name and extension
        (e.g. '.png'). Extension .png is recommended.
        """

        plt.savefig(target_dir)
        plt.close(self.fig)


    ############ PRIVATE DEFs FOR FORMATTING PLOTS/DATA ##############


    def __clean_up_plots(self, analyte_count, axes):

        ax_off = ()

        if analyte_count == 3:
            xt_off = (axes[0], axes[1])
        elif analyte_count == 4:
            xt_off = (axes[0][0], axes[0][1])
        elif analyte_count == 5:
            xt_off = (axes[0][0], axes[0][1], axes[1][0])
            ax_off = (axes[2][1],)
        elif analyte_count == 2:
            xt_off = (axes[0],)

        for ax in xt_off:
            ax.set_xticklabels([])
        for ax in ax_off:
            ax.axis('off')


    def __format_for_graphing(self):

        col_names, colormap = [], []
        gr_data, color_df = pd.DataFrame(), pd.DataFrame()
        
        #df = self.graph_data.sort_values(by=['sampleId'])
        df = self.graph_data
        a = list(set(df['call_if_pos']))

        for el in a:
            if el[:3] == "WT-":
                col_names.insert(0, el)
            elif el != "UEP":
                col_names.append(el)
        col_names.insert(0, "UEP")

        for cn in col_names:

            b = df[df.call_if_pos == cn]

            yvals = list(b[self.y_axis])
            names = list(b['sampleId'])
            colors = list(b['col_map'])

            gr_data[cn] = pd.Series(yvals, index=names)
            color_df[cn + '_color'] = pd.Series(colors, index=names)

        for c in color_df.columns:

            colormap.append(list(color_df[c]))

        return gr_data, colormap


    def __map_colors(self):

        col_map = []

        for i, row in self.graph_data.iterrows():
            if row['peakType'] == 'P':
                col_map.append(self.uep_color)
            elif row['pos_neg'] == 'neg':
                col_map.append(self.neg_color)
            elif row['pos_neg'] == 'pos':
                col_map.append(self.pos_color)
            else:
                print("Color mapping found an error. (line 273)")
                return

        a = pd.Series(col_map, index=self.graph_data.index)
        self.graph_data['col_map'] = a


    ############ PRIVATE DEFs FOR GETTING INFO #######################


    def __get_axis_list(self, col_count, row_count, axes):

        c = col_count
        r = row_count
        a = axes

        if c == 1:
            al = [a[x] for x in range(r)]
        else:
            al = [a[x][y] for x in range(r) for y in range(c)]

        return al


    def __get_chart_dimensions(self, analyte_count):

        index = analyte_count - 2

        c = self.chart_params['cols'][index]
        r = self.chart_params['rows'][index]
        w = self.chart_params['width'][index]
        h = self.chart_params['height'][index]

        return c, r, h, w


class Masshole_Finder:

    def __init__(self, path, pad=30):

        """
        Var 'path' is the path for an csv-formatted MassARRAY
        assay definition file.

        Although functional MassARRAY assay definition files are
        Excel-formatted, this class does not yet support parsing
        Excel-formatted files.
        """

        self.probe_set = pd.DataFrame()
        self.gaps = {}
        self.massholes = []
        self.well_list = []
        self.pad = pad
        self.source_path = path
        self.term_masses = {'none': 0,
                            'C': 247.2,
                            'A': 271.2,
                            'G': 287.2,
                            'T': 327.2}
        self.__parse_ma_file(path)
        self.__calc_all_gaps()


    ############ PUBLIC DEFs FOR GETTING MASSHOLES ###################


    def get_all_uep_holes(self, m_min=4500, m_max=11000):

        """
        Returns a dictionary of available spectrum positions for
        each well. These indicate what masses are available for
        designing new SBE probes.
        """

        complete_massholes = {}

        for well in self.well_list:
            wn = int(well[1:])
            complete_massholes[well] = self.get_uep_holes(wn, m_min, m_max)

        return complete_massholes
    

    def get_uep_holes(self, well, m_min=4500, m_max=11000):

        """
        Returns a list of spectrum positions for a new probe in
        a given well.
        """

        well = "W" + str(int(well))

        massholes = []
        gaps = self.gaps[well].copy()

        m_min, m_max = int(m_min), int(m_max)

        result = m_min

        while result < m_max:
            result = self.__get_next_masshole(m_min, m_max, gaps)
            if result >= m_min and result < m_max:
                massholes.append(result)
                gaps = self.__remove_all_uep_analytes(result, gaps)
                m_min = result

        return massholes


    ############ PRIVATE DEFs FOR DOING STUFF ########################


    def __parse_ma_file(self, path):
        
        data = pd.read_csv(path, sep='\t')
        self.probe_set = data[data['UEP_MASS'] > 4500]


    def __calc_all_gaps(self):

        self.well_list = list(self.probe_set.WELL.unique())
    
        old_spectrum = list(range(4500, 11000))

        for w in self.well_list:
            a = self.probe_set[self.probe_set['WELL'] == w]
            b = list(a.UEP_MASS)

            for mass in b:
                new_spectrum = self.__remove_all_uep_analytes(mass, old_spectrum)
                old_spectrum = new_spectrum
                
            self.gaps[w] = new_spectrum


    def __remove_all_uep_analytes(self, uep_mass, gap_list):

        for t, m in self.term_masses.items():

            gap_list = self.__remove_mass(uep_mass + m, gap_list)

        return gap_list


    def __remove_mass(self, mass, gap_list):

        mass, new_gap_list = int(mass), []

        mass_swath = [mass + x for x in range(-self.pad, self.pad)]

        for mg in gap_list:
            if mg not in mass_swath:
                new_gap_list.append(mg)

        return new_gap_list


    def __get_next_masshole(self, m_min, m_max, gaps):

        g = set(gaps)
        test = m_min

        while test < m_max:
            ms = self.__get_mass_set(test)
            if ms.issubset(g):
                return test
            else:
                test += 1

        return m_max + 1


    ############ PRIVATE UTILITIES ###################################


    def __get_mass_set(self, uep_mass):
        a = uep_mass
        b = a + 247
        c = a + 271
        d = a + 287
        e = a + 327
        return set((a, b, c, d, e))


class Mass_Blaster:

    """
    Class for generating candidate appendages for a-la-carte SBE
    design.
    
    Given an assay description file (csv format), and output from
    Masshole_Finder.get_all_uep_holes(), this class will find nucleotide
    compositions for an appendage that will give the target mass.

    Useful for determining what sequence you can add to the 5' of a
    single-base extension oligo such that the oligo will have some
    target mass.

    Arguments:
        'tolerance' is the maximum acceptable mass difference
        'max_gc' is the maximum acceptable G/C content for an appendage (0-1).
        'max_gc_count' is the maximum acceptable number of G/C residues.
    """

    def __init__(self, tolerance=10, max_gc=0.4, max_gc_count=5):

        self.dna_masses = {'A':313.21,'C':289.18,'G':329.21,'T':304.2}
        self.tolerance = tolerance
        self.max_gc = max_gc
        self.max_gc_count = max_gc_count

        self.source_path = None
        self.result_path = None
        self.mass_lib = None
        self.base_probes = None

        self.mass_lib = self.get_mass_lib()


    ############ PUBLIC GETTERS ######################################
        

    def save_results(self, path, how_many=5):

        """
        Kwarg 'how_many' is the number of top candidates that will be
        reported for each probe.

        Arg "path" is the path where the results should be saved, must
        include the file name and extension (.txt is recommended).
        """

        self.result_path = path
        open(self.result_path, "w").close()
        self.__write_result_file_headers()
        self.__analyze_all_base_seqs(5)


    def get_mass_lib(self):
        
        its = self.__get_iterators(2,6)
        
        base_lib = self.__make_base_lib(its)
        comb_lib1 = self.__make_combination_lib(base_lib)
        comb_lib2 = self.__make_combination_lib(comb_lib1)
        mass_lib = self.__cull_gc_metrics(comb_lib2)

        return mass_lib


    ############ PUBLIC SETTERS ######################################
        

    def input_base_probe_sequences(self, path):

        """
        Var 'path' should be a path to a multi-fasta file, i.e.:
            >sequence_1_name
            >AGCTAGCTAGCT
            >sequence_2_name
            >ACGTACGTACGT
            ...
        Arg 'path' must include the file name.
        """
        self.source_path = path

        a = sio.parse(path, 'fasta')
        self.base_probes = []
        
        for el in a:
            name = str(el.name)
            seq = str(el.seq)
            mass = self.__get_ntlist_mass(list(seq))
            print(name, seq, mass)
            self.base_probes.append((mass, name, seq))


    def input_massholes(self, massholes):

        """
        Arg 'massholes' should be the dict output from
        Masshole_Finder.get_all_uep_holes().
        """
        
        self.massholes = []

        for w, ml in massholes.items():
            entries = [(m, w) for m in ml]
            self.massholes.extend(entries)


    ############ PRIVATE THINKERS ####################################


    def __analyze_all_base_seqs(self, how_many):

        for bp in self.base_probes:
            self.__analyze_probe(bp, how_many)


    def __analyze_probe(self, pr_rec, how_many):
    
        with open(self.result_path, 'a') as rfile:
            
            for m_rec in self.massholes:
                
                t_mass = m_rec[0] - pr_rec[0]
                seqs, masses = self.__get_optimal_nt_comp(t_mass, how_many)
                
                r_txt = pr_rec[1] + "\t" + pr_rec[2] + "\t" + str(pr_rec[0]) + \
                        "\t" + m_rec[1] + "\t" + str(m_rec[0]) + "\t" + \
                        str(t_mass) + "\t"
                
                for i in range(how_many):
                    index = how_many - 1 - i
                    
                    if len(seqs[index]) > 0:
                        entry = r_txt + seqs[index] + "\t" + str(masses[index]) + \
                                "\t" + str(pr_rec[0] + masses[index]) + "\n"
                        rfile.write(entry)


    ############ PRIVATE WRITERS #####################################


    def __write_result_file_headers(self):
            
        headers = ["probe_name",
                   "base_sequence",
                   "base_sequence_mass",
                   "well_number",
                   "target_mass", 
                   "appendage_target_mass",
                   "appendage_composition",
                   "appendage_mass",
                   "probe_mass"]
        
        header_line = "\t".join(headers) + "\n"
        
        with open(self.result_path, 'a') as rfile:
            rfile.write(header_line)


    ############ PRIVATE GETTERS #####################################


    def __get_optimal_nt_comp(self, mass, how_many):
        
        seqs = [''] * how_many
        masses = [''] * how_many

        if mass < 289.2:
            return seqs, masses
        
        for m, s in self.mass_lib.items():
            
            md = abs(mass - m)

            if md < self.tolerance:

                seqs.append(self.mass_lib[m])
                masses.append(m)
                seqs = seqs[1:]
                masses = masses[1:]
                min_mass_diff = md

        return seqs, masses


    def __get_ntlist_mass(self, nts):
        
        mass = 0

        for nt in nts:
            mass += self.dna_masses[nt]
            
        return mass - 61.962


    def __make_base_lib(self, iterator_list):

        ml = {}

        for it in iterator_list:
            for nt_seq in it:
                m = self.__get_ntlist_mass(nt_seq) + 61.962
                if m not in ml:
                    ml[m] = ''.join(nt_seq)

        return ml
    

    def __make_combination_lib(self, ml):        

        ml_temp = {}        

        for m1, s1 in ml.items():

            for m2, s2 in ml.items():

                m = m1 + m2
                s = s1 + s2

                if m not in ml:
                    ml_temp[m] = s

        ml.update(ml_temp)     

        return ml


    ############ PRIVATE THINKERS ####################################


    def __cull_gc_metrics(self, lib):

        new_lib = {}

        for m, seq in lib.items():
            a = seq.count("G") + seq.count("C") < self.max_gc_count
            b = a/len(seq) < self.max_gc
            if a and b:
                new_lib[m] = seq

        return new_lib


    def __get_iterators(self, minimum, maximum):
        
        its = []

        for i in range(minimum, maximum + 1):
            its.append(itls.product('ACGT', repeat=i))

        return its


################### DUMMY DATAFRAME FOR TESTING ######################


def get_test_dataframe(namelist=0):
    """
    Returns a small DataFrame for testing Pandas DataFrame operations.
    """
    a = (('tom','dick','harry','sally','laverne','shirley'),
         ('andria','sau','tyson','argelia','jarod','shona'),
         ('rasheeda','angla','angie','catrina','nico','sam')) 
    b = ('bench','squat','curl')
    c = ((150, 180, 165, 130, 120, 150), \
         (180, 95,  200, 80,  80,  95), \
         (40,  35,  45,  20,  35,  20))
    d = {}
    for i in range(3):
        d[b[i]] = c[i]
    return pd.DataFrame(d, index=a[namelist])
    



################### CONTROL CODE #####################################


def main_1(export=False):
    """
    Currently, this function when the script is run directly. Ideally,
    this will be replaced with a GUI.
    """
    root = tkinter.Tk()
    data_path = fd.askopenfilename(title="Where is the exported MassARRAY data file?")
    root.destroy()
    ma4expt = MA4_Experiment(data_path)

######################################################################
##
## Dealing with the NTC replicate in expt 012418:
##
##    a = ma4expt.get_sample_name("G01")
##    b = ma4expt.get_sample_name("C05")
##    
##    if a == b:
##        ma4expt.rename_sample("NTC_1", ('G01','G02','G03','G04'))
##        ma4expt.rename_sample("NTC_2", ('C05','C06','C07','C08'))
##
######################################################################


######################################################################
##
## Uncomment this (below) for samples Mut_50, _40, _30, etc.
##
##    sgs = [[],[]]
##    for sn in c:
##        if sn[:3] == "Mut":
##            sgs[0].append(sn)
##        elif sn[:3] == "NTC":
##            sgs[0].append(sn)
##            sgs[1].append(sn)
##        else:
##            sgs[1].append(sn)
##
##    for l in sgs:
##        l.sort()
##
##    groups = [x if len(x) > 1 else [] for x in sgs]
##
######################################################################

    # Assigning sample groups (specific for expt 012418).

    if export == True:

        ma4expt.export_data(include_spectra=True)

    elif export == False:
    
        c = ma4expt.get_sample_names()
        groups = [c,]
                
        ma4expt.define_sample_groups(groups)
        
        grapher = MA_Data_Grapher()
        probe_names = ma4expt.get_probe_names()
        sample_groups = ma4expt.get_sample_groups()

        group_names = ['con', 'clin']

        for i, probe_name in enumerate(probe_names):
            
            for j, sample_list in enumerate(sample_groups):

                d = ma4expt.get_subset(probe_name=probe_name, sample_list=sample_list)
                grapher.new_dataset(d)
                grapher.graph(probe_name)

                t_dir = os.path.dirname(data_path)

                res_fname = probe_name + "_" + group_names[j] + ".png"

                grapher.save_graph(t_dir + "/" + res_fname)
                break


def main_2():

    intro_text = \
               "This facility requires two files:\n" \
               "    1) The MassARRAY assay description file (.txt).\n" \
               "    2) A fasta-formatted file containing the allele-\n" \
               "        specific regions of your new SBE probes.\n"
    print(intro_text)
    a = input("Hit enter when you are ready.")
    
    root = tkinter.Tk()

    path1 = fd.askopenfilename(title="Where is the MassARRAY assay description file?")
    
    mhf = Masshole_Finder(path1)
    massholes = mhf.get_all_uep_holes()
    
    print("Your new UEP can have one of the following masses (by well):")
    for k, v in massholes.items():
        print(k)
        print(v)

    path2 = fd.askopenfilename(title="Where is the Fasta file with new base probes?")
    
    mb = Mass_Blaster()
    mb.input_massholes(massholes)
    mb.input_base_probe_sequences(path2)

    path3 = fd.asksaveasfilename(title="Where do you want to save the results?", \
                                 defaultextension=".txt")

    mb.save_results(path3)

    print("Results have been saved to:")
    print(path3)
    
    root.destroy()


def main_0():
    print("What do you want to do?")
    print("    1) Graph experiment results.")
    print("    2) Design a new SBE probe.")
    print("    3) Export data to CSV file.")
    choice = int(input("Enter your choice: "))
    if choice == 1:
        main_1()
    elif choice == 2:
        main_2()
    elif choice == 3:
        main_1(export=True)


if __name__ == '__main__':
    main_0()
