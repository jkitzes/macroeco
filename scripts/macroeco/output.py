#!/usr/bin/python

'''This module provides functions for outputting results of macroeco 
analyses'''


from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import logging
from macroeco.utils.form_func import output_form, add_field
import copy as cp
import os
import shutil

__author__ = "Mark Wilber"
__copyright__ = "Copyright 2012, Regents of the University of California"
__license__ = None
__version__ = "0.5"
__maintainer__ = "Mark Wilber"
__email__ = "mqw@berkeley.edu"
__status__ = "Development"

readme_info_plots =\
'''
FOLDER DESCRIPTION
-------------------

The folder {3} contains {0} files.  There are {1} {4} represented as png
files and {2} csv files which contain the data required to generate each
plot.  The csv files have identical names to the png files to which they
correspond. Each file name is a concatenation of the following strings:
analysis name, run name, data name, and {5}.  An additional identifier is
appended to the file name after {5} in order to make each file unique.  It is
either a species identifier or a number.

On the right hand side of each plot, you will see a string that begins
'Criteria for plot'.  The criteria are either a species name or string that 
looks like

'y': [('>=', 0.0), ('<', 150.0)], 'x': [('>=', 0.0), ('<', 50.0)]

This can be interpreted as follows.  The plot under consideration has 'y' values
greater than or equal to 0 and less than 150 and 'x' values greater than or
equal to 0 and less than 50.  Similarly a criteria string of the form

'year' : ('==' , 1998)

can be interpreted as the plot under consideration has 'year' values equal to
1998. The criteria is determined by how you decided to divide your plot for the
analysis. A criteria string of the form

'temperature' : ('==', 'cool')

can be interpreted as the plot under consideration has 'temperature' values
equal to 'cool'. '''

readme_info_summary=\
u"""
FOLDER DESCRIPTION
------------------

The folder {0} contains {1} txt file(s) and {1} csv file(s).  Each .txt file
contains a summary for the plot generate by the criteria at the header of the
file.  The criteria are either a species name or string that looks like

'y': [('>=', 0.0), ('<', 150.0)], 'x': [('>=', 0.0), ('<', 50.0)]

This can be interpreted as follows.  The plot under consideration has 'y' values
greater than or equal to 0 and less than 150 and 'x' values greater than or
equal to 0 and less than 50.  Similarly a criteria string of the form

'year' : ('==' , 1998)

can be interpreted as the plot under consideration has 'year' values equal to
1998. 

Each txt file has a corresponding csv plot with the AIC values in tabular form
for easy analysis.

Each summary file contains summary statistics for the observed data and each
distribution to which the observed data was compared. Each file name is a
concatenation of the following strings: analysis name, data name and
summary_table or AIC_table.  An additional identifier is appended to the file
name after summary_table in order to make each file unique.  It is either a
species identifier, a number, or both."""

readme_info_rarity =\
'''
FOLDER DESCRIPTION
------------------

The folder {0} contains {1} csv files.  Each file contains the
columns 'data_name', 'criteria', 'observed', and any number of columns with
distribution names.  These are the distributions to which the data was
compared.  The column data_name gives the name of the data being examined, the
column criteria describes the specifications that made the given plot, the
remaining columns describe the number of items that had a value below a
prespecified minimum.  The prespecified minimum can be found in the file name
immediately after '_<=_'.   Each file name is a concatenation of the following
strings: analysis name, data name and 'rarity_<=_' some minimum.
'''

readme_info_sar=\
'''
FOLDER DESCRIPTION
------------------

The folder {0} contains {1} png files and {2} csv files.  The png file(s) are 
log-log SAR plot(s) with area_fraction on the x-axis and species on the y-axis.
The names of the png file(s) are a concatenation of the following strings:
analysis name, run_name, data_name, and SAR_plot. A number is appended to the
end of the plot to ensure the filename is unique. The csv
files contain the data required to make the given plot(s). Each csv file
contains two columns, species and area_fraction.  area_fraction assigns the
base area a value of 1 and represents all other areas as a fraction of the base
area. The csv file name(s) are a concatenation of the following strings:
analysis_name, run_name, data_name, SAR_plot_, a unique number, and the SAR
name.

'''


class DistributionOutput(object):
    '''
    This formats and outputs analyses on distributions

    '''

    def __init__(self, out_dir):
        '''
        Parameters
        ----------
        out_dir : string
            Output directory of object
        '''

        self.out_dir = out_dir
        self.urns = 'Urns'
        self.balls = 'Balls'
        self.rad_x_axis = 'Rank'
        self.rad_y_axis = 'Abundance'
        self.cdf_x_axis = 'Abundance'
        self.cdf_y_axis = 'Cumulative Probability'
        self.variable = 'abundance'


    def write_summary_table(self, smry, criteria=None):
        '''
        Parameters
        ---------
        smry : tuple
            A tuple of length two in which the first object is a dictionary as
            returned by the function compare_summary within the 
            CompareDistribution class.  The second object is  dictionary with
            the keyword 'mins' that refers to the abundance that the 'tot_min'
            keyword is less than OR equal to. If mins = 1, tot_min describes the
            number of items with counts <= 1. 

        criteria : array-like object
            An array-like object in which contains either string or dicts that
            tell how each dataset was generated.  Describes the subsetting of
            an sad and the species ID of an ssad.

        Notes
        -----
        Writes out a formatted txt file to self.out_dir 

        '''
        # Make output folder
        folder_name = 'summary_statistics_' + self.out_dir
        make_directory(folder_name)
        
        tot_sad = len(smry['observed']['balls'])
        if criteria != None:
            assert len(criteria) == tot_sad, "len(criteria) must  equal" + \
                                   " number of data arrays under consideration"
        ob = smry['observed']

        count = 0
        for i in xrange(tot_sad):
            if criteria != None and np.all([type(crt) != dict for crt in
                                                                  criteria]):
                filename = os.path.join(folder_name, self.out_dir + \
                        '_summary_table_' + str(criteria[i]) + '.txt')
                filename_aic = os.path.join(folder_name, self.out_dir + \
                        '_AIC_table_' + str(criteria[i]))

                
            else:
                filename = os.path.join(folder_name, self.out_dir + 
                                    '_summary_table_' + str(i) + '.txt')
                filename_aic = os.path.join(folder_name, self.out_dir + 
                                    '_AIC_table_' + str(i))


            fout = open(filename, 'w')
            logging.info('Writing summary table %s' % filename)


            if criteria != None:
                fout.write('CRITERIA: ' + str(criteria[i]) + '\n\n')
            else:
                fout.write('CRITERIA: NONE ' + str(i) + '\n\n')

            # Getting rarity 
            ob_rare = {} 
            for mins in ob['tot_min'].iterkeys():
                ob_rare['<=' + str(mins)] = ob['tot_min'][mins][i]

            fout.write('EMPIRICAL VALUES:\n' + self.urns + ' = ' + 
                    str(ob['urns'][i]) + '\n' + self.balls + ' = ' + 
                    str(ob['balls'][i]) + '\nObserved Nmax= ' + 
                    str(ob['max'][i]) + '\nObserved Rarity = ' +
                    str(ob_rare) + '\n\n')


            # Also output AIC values in for each table. Could add other other
            # measures to this table as well. 
            # Might break this out later
            aic_vals = {}

            for kw in smry.iterkeys():
                if kw != 'observed':
                    dt= smry[kw]
                    # set relevant aic values for table output
                    aic_vals[kw]={'AIC_weights' : dt['aic_w'][i], 'Delta_AIC' :
                                   dt['aic_d'][i], 'Parameter_number' :
                                   dt['par_num'][i], 'Corrected_AIC' :
                                   dt['aic'][i]}
                    # Getting rarity
                    dt_rare = {}
                    for mins in dt['tot_min'].iterkeys():
                        dt_rare['<=' + str(mins)] = dt['tot_min'][mins][i]

                    fout.write('PREDICTED DISTRIBUTION : ' + kw + '\n' + 
                    self.urns + ' = ' + str(dt['urns'][i]) + '\n' + 
                    self.balls + ' = ' + str(dt['balls'][i]) + 
                    '\nAIC = ' + str(dt['aic'][i]) + '\nDelta_AIC = ' +
                    str(dt['aic_d'][i]) + '\nAIC_weight = ' +
                    str(dt['aic_w'][i]) + '\nNumber of Parameters = ' + 
                    str(dt['par_num'][i]) + '\nPredicted Nmax = ' + 
                    str(dt['max'][i]) + '\nPredicted Rarity = ' +
                    str(dt_rare) + '\n\n')
            fout.close()
            count += 1

            # Make and print AIC table
            dtype = [('Model', 'S30'), ('Parameter_number', np.float),
                        ('Corrected_AIC', np.float), ('AIC_weights', np.float),
                        ('Delta_AIC', np.float)]
            aic_array = np.empty(len(aic_vals), dtype=dtype)
            for j, model_name in enumerate(aic_vals.iterkeys()):
                aic_array['Model'][j] = model_name
                aic_array['Parameter_number'][j] =\
                                       aic_vals[model_name]['Parameter_number']
                aic_array['Corrected_AIC'][j] =\
                                       aic_vals[model_name]['Corrected_AIC']
                aic_array['AIC_weights'][j] =\
                                       aic_vals[model_name]['AIC_weights']
                aic_array['Delta_AIC'][j] =\
                                       aic_vals[model_name]['Delta_AIC']
            output_form(aic_array, filename_aic)

        fout = open(os.path.join(folder_name, 'README'), 'w')
        fout.write(readme_info_summary.format(folder_name, count))
        fout.close()

    

    def plot_rads(self, rads, criteria=None, species=None):
        '''
        Plotting the observed and predicted rank abundance distributions

        Parameters
        ----------
        rads : dict
            A dictionary that is returned from the function compare_rads in the
            CompareDistribution class.

        criteria : list of objects
            If not none, the objects in criteria will be printed a strings in
            the plots and file names.

        Notes
        -----
        Saves RAD plots to given out_dir.  Saves as many plots as there are
        observed distributions.

        '''
        folder_name = 'rank_abundance_plots_' + self.out_dir
        make_directory(folder_name)

        tot_sad = len(rads['observed'])
        recs = make_rec_from_dict(rads, tot_sad, species=species)

        if criteria != None:
            assert len(criteria) == tot_sad, "len(criteria) must  equal" + \
                                   " number of data arrays under consideration"
        count = 0
        for i, data in enumerate(recs):
            
            # Plot all columns of the rec array
            plot_rec_columns(data)
            plt.semilogy()
            plt.ylabel('Log ' + self.rad_y_axis)
            plt.xlabel(self.rad_x_axis)
            
            if criteria != None and np.all([type(crt) != dict for crt in
                                                                    criteria]):
                plt.title('Rank abundance distribution for ' + str(criteria[i]))
                filename = os.path.join(folder_name, self.out_dir +
                                    '_rank_abundance_plot_' + str(criteria[i]))

                logging.info('Saving figure and csv ' + filename)
                plt.savefig(filename)
                output_form(recs[i], filename)
                count += 2

            elif criteria != None and np.all([type(crt) == dict for crt in
                                                                    criteria]):
                plt.title('Rank abundance distribution')
                plt.figtext(.97, .5, 'Criteria for plot: ' + str(criteria[i]), 
                                rotation='vertical', size=8,
                                horizontalalignment='center',
                                verticalalignment='center')

                filename = os.path.join(folder_name,  self.out_dir + 
                                              '_rank_abundance_plot_' + str(i))
                logging.info('Saving figure ' + filename)
                plt.savefig(filename)
                output_form(recs[i], filename)
                count += 2

            else:
                plt.title('Rank abundance distribution: plot number ' + str(i))
                
                filename = os.path.join(folder_name, self.out_dir + 
                                              '_rank_abundance_plot_' + str(i))
                logging.info('Saving figure ' + filename)
                plt.savefig(filename)
                output_form(recs[i], filename)
                count += 2
            
            plt.clf()

        fout = open(os.path.join(folder_name, 'README'), 'w')
        fout.write(readme_info_plots.format(count, count /2, count/2,
            folder_name, 'rank abundance plots (RAD)', 'rank_abundance_plot'))
        fout.close()

    
    def plot_cdfs(self, cdfs, obs_sads, criteria=None, species=None):
        '''

        Plots observed vs predicted cdfs and returns a csv file with values
        used for plotting.


        Parameters
        ----------
        cdfs : dict
            A dictionary that is returned from the function compare_cdfs in the
            CompareDistribution class. 

        obs_sads : list
            A list of arrays.  The observed sad(s)

        '''
        # Make directory
        folder_name = 'cdf_plots_' + self.out_dir
        make_directory(folder_name)

        tot_sad = len(cdfs['observed'])
        recs = make_rec_from_dict(cdfs, tot_sad, add_rank=False)
        if criteria != None:
            assert len(criteria) == tot_sad, "len(criteria) must  equal" + \
                                   " number of data arrays under consideration"

        count = 0
        for i, data in enumerate(recs):
            
            names = data.dtype.names
            for nm in names:
                fig = plt.plot(np.sort(obs_sads[i]), np.sort(data[nm]), '-o')
            
            # Formatting
            fig[0].axes.xaxis.tick_bottom()
            fig[0].axes.yaxis.tick_left()
            ylim = list(plt.ylim())
            if ylim[0] == 0:
                ylim[0] = -.1
            plt.ylim((ylim[0], 1.1))            
            xlim = plt.xlim()
            plt.xlim((.9, xlim[1] + 10))
            plt.legend(names, loc='best')
            plt.semilogx()
            plt.ylabel(self.cdf_y_axis)
            plt.xlabel('Log ' + self.cdf_x_axis)
            
            # Add observed to cdf array
            if species != None:
                sorted_ab, sorted_spp = sort_rank_abund([obs_sads[i]],
                                                                [species[i]])
                n_rec = add_field(data, [(self.variable, np.float)])
                n_rec = add_field(n_rec, [('species', 'S40')])
                n_rec[self.variable] = sorted_ab[0]
                n_rec['species'] = sorted_spp[0]
            else:
                n_rec = add_field(data, [(self.variable, np.float)])
                n_rec[self.variable] = np.sort(obs_sads[i])

            if criteria != None and np.all([type(crt) != dict for crt in
                                                                    criteria]):
                plt.title('Cumulative density function for species ' + str(criteria[i]))

                filename = os.path.join(folder_name, self.out_dir +
                                               '_cdf_plot_' + str(criteria[i]))
                logging.info('Saving figure and csv ' + filename)
                plt.savefig(filename)
                output_form(n_rec, filename)
                count += 2

            elif criteria != None and np.all([type(crt) == dict for crt in
                                                                    criteria]):
                plt.title('Cumulative Density Function')
                plt.figtext(.97, .5, 'Criteria for plot: ' + str(criteria[i]), 
                                rotation='vertical', size=8,
                                horizontalalignment='center',
                                verticalalignment='center')

                filename = os.path.join(folder_name, self.out_dir + 
                                                   '_cdf_plot_' + str(i)) 
                logging.info('Saving figure ' + filename)
                plt.savefig(filename)
                output_form(n_rec, filename)
                count += 2

            else:
                plt.title('CDF: plot number ' + str(i))
                filename = os.path.join(folder_name, self.out_dir +
                                                        '_cdf_plot_' + str(i)) 
                logging.info('Saving figure and csv ' + filename)
                plt.savefig(filename)
                output_form(n_rec, filename)
                count += 2

            plt.clf()

        fout = open(os.path.join(folder_name, 'README'), 'w')
        fout.write(readme_info_plots.format(count, count/2, count/2,
                folder_name, 'cumulative density plots (cdf)', 'cdf_plot'))
        fout.close()

class SADOutput(DistributionOutput):
    '''
    Derived class for SAD output
    '''

    def __init__(self, out_dir):
        '''
        Parameters
        ----------
        out_dir : string
            Output directory of object

        '''
        self.out_dir = out_dir
        self.urns = 'Species'
        self.balls = 'Total Individuals'
        self.rad_x_axis = 'Rank'
        self.rad_y_axis = 'Abundance'
        self.cdf_x_axis = 'Abundance'
        self.cdf_y_axis = 'Cumulative Probability'
        self.variable = 'abundance'

class SSADOutput(DistributionOutput):
    '''
    Derived class for SSAD output
    '''

    def __init__(self, out_dir):
        '''
        Parameters
        ----------
        out_dir : string
            Output directory of object

        '''
        self.out_dir = out_dir
        self.urns = 'Cells'
        self.balls = 'Individuals'
        self.rad_x_axis = 'Rank'
        self.rad_y_axis = 'Abundance'
        self.cdf_x_axis = 'Abundance'
        self.cdf_y_axis = 'Cumulative Probability'
        self.variable = 'abundance'

class SAROutput(object):
    '''
    This object interacts with CompareSARCurves
    '''

    def __init__(self, out_dir):
        '''
        Parameters
        ----------
        out_dir : string
            Output directory of object
        '''
        self.out_dir = out_dir

    def plot_sars(self, sars, names=[]):
        '''
        Plots observed vs predicted sars

        Parameters
        ----------
        sars : list of dicts
            The output of CompareSARCurve method compare_curves

        names : list or strings
            If not None, names is a list of the same length as sars.  Gives the
            desired names for the plots.

        '''

        folder_name = 'sar_plots_' + self.out_dir
        make_directory(folder_name)

        if len(names) != 0:
            assert len(names) == len(sars); "Length of names must equal" + \
                                           "length of sars"
        count = 0
        for i, sar in enumerate(sars):
            filename = os.path.join(folder_name, self.out_dir + '_SAR_plot_' +
                                                                        str(i))
            legend = []
            for kw in sar.iterkeys():
                legend.append(kw)
                if kw == 'observed':
                    fig = plt.plot(sar[kw]['area'], sar[kw]['items'], '-o')
                else:
                    fig = plt.plot(sar[kw]['area'], sar[kw]['items'])

                # Change dtype names and output
                defnm = sar[kw].dtype.names
                sar[kw].dtype.names = ('species', 'area_fraction')
                output_form(sar[kw], filename + '_' + kw)
                sar[kw].dtype.names = defnm

            # Plot formatting 
            fig[0].axes.xaxis.tick_bottom()
            fig[0].axes.yaxis.tick_left()

            plt.loglog()
            plt.legend(tuple(legend), loc='best')
            plt.xlabel('log(Area Fraction)')
            plt.ylabel('log(Species Number)')
            if len(names) != 0:
                plt.title(names[i])
            else:
                plt.title('SAR plot %i' % (i))
            filename = os.path.join(folder_name, self.out_dir + '_SAR_plot_' +
                                                                        str(i))
            logging.info('Saving figure ' + filename)
            plt.savefig(filename)
            plt.clf()
            count += 1

        fout = open(os.path.join(folder_name, 'README'), 'w')
        fout.write(readme_info_sar.format(folder_name, count, count * len(sar))) 
        fout.close()


class ASEDOutput(object):
    '''
    Class outputs the average species energy distributions by interacting with
    CompareASED
    '''

    def __init__(self, out_dir):
        '''
        Parameters
        ----------
        out_dir : string
            Output directory of object
        '''
        self.out_dir = out_dir         

    def plot_reds(self, reds, criteria=None, species=None):
        '''
        Plotting the observed and predicted rank abundance distributions

        Parameters
        ----------
        reds : dict
            A dictionary that is returned from the function compare_reds in the
            CompareASED class.

        criteria : list of objects
            If not none, the objects in criteria will be printed a strings in
            the plots and file names.

        Notes
        -----
        Saves RAD plots to given out_dir.  Saves as many plots as there are
        observed distributions.

        '''
        folder_name = 'ased_rank_energy_plots_' + self.out_dir
        make_directory(folder_name)

        tot_sad = len(reds['observed'])
        recs = make_rec_from_dict(reds, tot_sad, species=species)

        if criteria != None:
            assert len(criteria) == tot_sad, "len(criteria) must  equal" + \
                                   " number of data arrays under consideration"
        count = 0
        for i, data in enumerate(recs):
            
            # Plot all columns of the rec array
            plot_rec_columns(data)
            plt.semilogy()
            plt.ylabel('Log Energy')
            plt.xlabel('Rank')
            
            if criteria != None and np.all([type(crt) != dict for crt in
                                                                    criteria]):
                plt.title('ASED rank energy distribution for ' + str(criteria[i]))
                filename = os.path.join(folder_name, self.out_dir +
                                    '_rank_abundance_plot_' + str(criteria[i]))

                logging.info('Saving figure and csv ' + filename)
                plt.savefig(filename)
                output_form(recs[i], filename)
                count += 2

            elif criteria != None and np.all([type(crt) == dict for crt in
                                                                    criteria]):
                plt.title('ASED rank energy distribution')
                plt.figtext(.97, .5, 'Criteria for plot: ' + str(criteria[i]), 
                                rotation='vertical', size=8,
                                horizontalalignment='center',
                                verticalalignment='center')

                filename = os.path.join(folder_name,  self.out_dir + 
                                              '_ased_rank_energy_plot_' + str(i))
                logging.info('Saving figure ' + filename)
                plt.savefig(filename)
                output_form(recs[i], filename)
                count += 2

            else:
                plt.title('ASED rank energy distribution: plot number ' + str(i))
                
                filename = os.path.join(folder_name, self.out_dir + 
                                              '_ased_rank_energy_plot_' + str(i))
                logging.info('Saving figure ' + filename)
                plt.savefig(filename)
                output_form(recs[i], filename)
                count += 2
            
            plt.clf()

        fout = open(os.path.join(folder_name, 'README'), 'w')
        fout.write(readme_info_plots.format(count, count /2, count/2,
            folder_name, 
            'average species energy distribution (ASED) rank' + 
            ' energy plots', 'ased_rank_energy_plot'))
        fout.close()
        


class IEDOutput(object):
    '''
    Class outputs individual energy distributions by interacting with
    CompareIED

    '''

    def __init__(self, out_dir):
        '''
        Parameters
        ----------
        out_dir : string
            Output directory of object
        '''
        self.out_dir = out_dir 

    def plot_reds(self, reds, criteria=None):
        '''
        Saves plot and csv file with predicted and empirical rank energy data

        Parameters
        ----------
        reds : tuple
            The output from the CompareIED.compare_rads method
        criteria : list or None
            A list of dicts with the criteria for divisions.  See Patch.sad

        Output
        ------
        This method outputs both a plot and a csv that compare observed and
        predicted individual rank energy curves for the entire community at the
        given subset.  

        '''
        folder_name = 'ied_rank_energy_plots_' + self.out_dir
        make_directory(folder_name)

        
        tot_reds = len(reds['observed'])
        recs = make_rec_from_dict(reds, tot_reds)
        if criteria != None:
            assert len(criteria) == tot_reds, "len(criteria) must  equal" + \
                                      " number of reds under consideration"
        count = 0                             
        for i, data in enumerate(recs):
            
            #Plot all data in a single rec array
            plot_rec_columns(data)

            # Make appropriate title for figure
            if criteria != None:
                plt.title('Rank Energy Distribution')
                plt.figtext(.97, .5, 'Criteria for plot: ' + str(criteria[i]), 
                                rotation='vertical', size=8,
                                horizontalalignment='center',
                                verticalalignment='center')
            else:
                plt.title('Rank Energy Distribution')
                plt.figtext(.97, .5, 'Plot number: ' + str(i), 
                                rotation='vertical', size=8,
                                horizontalalignment='center',
                                verticalalignment='center')

            plt.loglog()
            plt.ylabel('Log Energy')
            plt.xlabel('Log Rank')

            filename = os.path.join(folder_name, self.out_dir +
                                                  '_ied_rank_energy_' + str(i))
 
            logging.info('Saving figure ' + filename)
            plt.savefig(filename)
            plt.clf()
            output_form(recs[i], filename)
            count += 2

        fout = open(os.path.join(folder_name, 'README'), 'w')
        fout.write(readme_info_plots.format(count, count/2, count/2,
                   folder_name, 
                   'individual energy distribution (IED) rank energy plots',
                   'ied_rank_energy'))
        fout.close()

class SEDOutput(object):
    '''
    Class outputs species-level energy distributions by interacting with
    CompareSED

    '''

    def __init__(self, out_dir):
        '''
        Parameters
        ----------
        out_dir : string
            Output directory of object
        '''
        self.out_dir = out_dir  

    def plot_reds(self, reds, criteria=None):
        '''
        Saves plot and csv file with predicted and empirical rank energy data

        Parameters
        ----------
        reds : tuple
            The output from the CompareSED.compare_rads method
        criteria : list or None
            A list of dicts with the criteria for divisions.  See Patch.sad

        Output
        ------
        This method outputs both a plot and a csv that compare observed and
        predicted species-level rank energy curves.  

        '''
        folder_name = 'sed_rank_energy_plots_' + self.out_dir
        make_directory(folder_name)

        spp = reds[1]
        tot_reds = len(reds[0]['observed'])
        recs = make_rec_from_dict(reds[0], tot_reds)
        if criteria != None:
            assert len(criteria) == tot_reds, "len(criteria) must  equal" + \
                                      " number of reds under consideration"
        count = 0
        for i, data in enumerate(recs):

            plot_rec_columns(data)
            plt.semilogx()
            plt.ylabel('Energy')
            plt.xlabel('Log Rank')

            if spp != None:
                if criteria != None:
                    plt.title('Rank Energy Distribution for species ' +
                                                                   str(spp[i]))
                    plt.figtext(.97, .5, 'Criteria for plot: ' + 
                                str(criteria[i]), rotation='vertical', size=8,
                                horizontalalignment='center',
                                verticalalignment='center')
                else:
                    plt.title('Rank Energy Distribution for species ' + 
                                                                   str(spp[i]))

                filename = os.path.join(folder_name, self.out_dir + 
                              '_sed_rank_energy_' + str(spp[i]) + '_' + str(i))

                logging.info('Saving figure ' + filename)
                plt.savefig(filename)
                output_form(recs[i], filename)
                count += 2

            elif spp == None:
                if criteria != None:
                    plt.title('Criteria: ' + str(criteria[i]))
                else:
                    plt.title('Plot number ' + str(i))
                
                filename = os.path.join(folder_name, self.out_dir + 
                                                  '_sed_rank_energy_' + str(i))
                logging.info('Saving figure ' + filename)
                plt.savefig(filename)
                output_form(recs[i], filename)
            plt.clf()

        fout = open(os.path.join(folder_name, 'README'), 'w')
        fout.write(readme_info_plots.format(count, count/2, count/2,
                   folder_name, 
                   'species-level energy distribution (SED) rank energy plots',
                   'sed_rank_energy'))
        fout.close()

class OutputRarity(object):
    '''
    This object accepts output from the Compare.compare_rarity method to 
    output rarity

    '''

    def __init__(self, out_dir):
        '''

        Parameters
        ----------
        out_dir : string
            The output directory ID

        '''

        self.out_dir = out_dir

    def output_rarity(self, rarity, data_path, data, criteria=None):
        '''
        Outputs csv files containing rarity measures

        Parameters
        ----------
        rarity : a CompareRarity object
        
        data_path : str
            data_path string for identifying data in csv file

        data : list 
            A list of observed species abundance distributions

        criteria : dict or None
            The criteria for how the plot was split

        '''
        folder_name = 'rarity_values_' + self.out_dir
        make_directory(folder_name)

        keys = list(rarity.viewkeys())
        dtype = [(kw, np.int) for kw in keys]
        dtype.insert(0, ('criteria', 'S90')) # arbitrary length
        dtype.insert(0, ('data_name', 'S90')) # arbitrary length

        # Get a list of my minimums
        rare_list = []
        mins = list(rarity['observed'].viewkeys())
        for mn in mins:
            rarity_array = np.empty(len(data), dtype=dtype)
            rarity_array['criteria'] = criteria
            nm = os.path.split(data_path)[1].split('.')[0]
            rarity_array['data_name'] = np.repeat(nm, len(rarity_array))
            for kw in keys:
                rarity_array[kw] = rarity[kw][mn]
            rare_list.append(rarity_array)
        
        # Output results
        count = 0
        for i, rare in enumerate(rare_list):
            filename = os.path.join(folder_name, self.out_dir + '_rarity_<=_' +
                                                                  str(mins[i]))
            logging.info('Saving rarity data ' + filename)
            output_form(rare, filename)
            count += 1

        fout = open(os.path.join(folder_name, 'README'), 'w')
        fout.write(readme_info_rarity.format(folder_name, count))
        fout.close()

def make_rec_from_dict(dist_dict, num, species=None, dt=np.float, add_rank=True):
    '''
    Makes a structured/rec array from a dictionary

    Parameters
    ----------
    dist_dict : dict
        A dictionary with each keyword referencing a list of arrays

    num : int
        Number of rec_arrays to return in list

    species : None or list of iterables
        If not None, species should be a list of iterables that is the same
        length as the list of iterables in any keyword in dist_dict.

    Returns
    -------
    : structured array

    '''
    
    # Check that species has the appropriate length
    if species != None:
        species = cp.deepcopy(species)
        for val in dist_dict.itervalues():
            if len(species) != len(val):
                raise TypeError('Species must contain the same number of ' +
                                 'iterables as each value in dist_dict')
    # Sort Observed and species list
    if species != None:
        dist_dict['observed'], species = sort_rank_abund(dist_dict['observed'],
                                                                       species)
    recs = []
    names = list(dist_dict.viewkeys())
    dtype = zip(names, np.repeat(dt, len(names)))
    if species != None:
        dtype.insert(0, ('species', 'S40'))
    if add_rank:
        dtype.insert(0, ('rank', dt))
    for i in xrange(num):
        temp = np.empty(len(dist_dict[names[0]][i]), dtype=dtype)
        if species != None:
            temp['species'] = species[i]
        if add_rank:
            temp['rank'] = np.arange(1,len(temp) + 1)[::-1]
        for kw in dist_dict.iterkeys():
            temp[kw] = np.sort(dist_dict[kw][i])
        recs.append(temp)
    return recs

def sort_rank_abund(abund_list, spp_list):
    '''
    Sorts and returns two lists based on abundance
    
    Parameters
    ----------
    abund_list : list of arrays
    
    spp_list : list of arrays

    Returns
    -------
    :tuple
        sorted_abund, sorted_spp

    '''

    assert len(abund_list) == len(spp_list), 'Lengths of arguments not equal'
    assert np.all([len(a) == len(b) for a,b in zip(abund_list, spp_list)]),\
                        'Lengths of all corresponding iterables not equal'
    abund_list = [np.array(t) for t in abund_list]
    spp_list = [np.array(t) for t in spp_list]

    sorted_abund = []
    sorted_spp = []
    for i in xrange(len(abund_list)):
        temp = np.array(zip(abund_list[i], spp_list[i]), dtype=[('a',
                               abund_list[i].dtype), ('s', spp_list[i].dtype)])
        temp_sorted = np.sort(temp, order='a')
        sorted_abund.append(temp_sorted['a'])
        sorted_spp.append(temp_sorted['s'])

    return sorted_abund, sorted_spp

def plot_rec_columns(rec_array):
    '''
    Function plots the columns in a rec array.
    '''

    # Available plotting symbols
    plot_symbols = ['+', 's', 'd', '*', 'x', '8', 'H', '1', 'p', '2', '3',
                                                        '4', '|', 4, 5, 6, 7]
    names = rec_array.dtype.names
    legend = []

    # If there are more arrays than symbols just change colors of lines
    if len(names) > len(plot_symbols):
        for nm in names:
            if nm != 'species' and nm != 'rank':
                if nm == 'observed':
                    fig = plt.plot(np.arange(1, len(rec_array) + 1),
                                        np.sort(rec_array[nm])[::-1], '-o',
                                        color='black')
                    legend.append(nm)
                else:
                    fig = plt.plot(np.arange(1, len(rec_array) + 1),
                                        np.sort(rec_array[nm])[::-1], '-o')
                    legend.append(nm)

    # Else, use different symbols/markers for each line
    elif len(names) <= len(plot_symbols):

        # Counter is 0
        cnt = 0
        for nm in names:
            if nm != 'species' and nm != 'rank':
                if nm == 'observed':
                    
                    fig = plt.plot(np.arange(1, len(rec_array) + 1),
                                        np.sort(rec_array[nm])[::-1], '-o',
                                        color='black')
                    legend.append(nm)
                else:
                    fig = plt.plot(np.arange(1, len(rec_array) + 1),
                                        np.sort(rec_array[nm])[::-1], '-' +
                                        str(plot_symbols[cnt]))
                    legend.append(nm)
                cnt += 1
    # Include ticks only on bottom and left 
    fig[0].axes.xaxis.tick_bottom()
    fig[0].axes.yaxis.tick_left()
    
    plt.legend(tuple(legend), loc='best')

def make_directory(folder_name):
    '''Makes a directory named folder_name.  If the directory exists it
    is overwritten

    folder_name - Name of the directory
    '''

    try:
        os.mkdir(folder_name)
    except OSError:
        shutil.rmtree(folder_name)
        os.mkdir(folder_name)

            


    













           
            
            
                        

            


        



"""
def sad_cdf_obs_pred_plot(sad, outputID, params={}, interactive=False):
    '''Makes a plot of observed vs predicted cdf values
    based on the given sad

    Parameters
    ----------
    sad : ndarray
        SAD
    outputID : string
        The name of the saved plot
    params : dict
        If empty uses default values, if not incorporates given params
        into params into dict
    interactive : bool
        If True, figure is shown in interactive window.  
        If false, no figure is shown.

    Notes
    -----
    Saves plot to current working directory. Default distribution to graph
    is METE. Need top specify in parameters.xml if another distribution is
    desired.

    '''

    #Allowing the user to set parameters if they pass some in
    prm = {'clr_jit':'grey', 'clr_sct':'black', 'ln_clr':'grey', 'jit_scl':0.007,\
          'ylbl':'Observed cdf', 'xlbl':'Predicted cdf',\
          'title': 'Observed vs. predicted values of SAD cdf', 'distr':'mete'}
    if len(params) != 0:
        inter = set(params.viewkeys()).intersection(set(prm.viewkeys()))
        if len(inter) != 0:
            module_logger.debug('Setting parameters ' + str(inter))
            for key in inter:
                if key is 'jit_scl':
                    prm[key] = eval(params[key])
                else:
                    prm[key] = params[key]
          
    cdf = sad_analysis.get_obs_vs_pred_cdf(sad, prm['distr'])
    jitt_x, jitt_y = jitter(cdf['pred'], cdf['observed'], jitter_scale=prm['jit_scl'])
    plt.plot([0] + list(cdf['observed']),[0] +  list(cdf['observed']),\
                                            color=prm['ln_clr'], linestyle='--')
        
    plt.scatter(jitt_x, jitt_y, color=prm['clr_jit']) 
    plt.scatter(cdf['pred'], cdf['observed'], color=prm['clr_sct'])
    '''for s in xrange(len(cdf)):
        plt.text(cdf['pred'][s], cdf['observed'][s] + 0.01, "n=" + str(cdf['n'][s]),\
                                fontsize=7)'''
    plt.ylabel(prm['ylbl'])
    plt.xlabel(prm['xlbl'])
    plt.title(prm['title'])
    plt.xlim(-0.1,1.1)
    plt.ylim(-0.1,1.1)
    plt.legend(('Slope 1', 'All species cum.\nprob. (jittered)',\
                    'Cum. prob.\nfor given n'),loc='best', prop={'size':12})
    plt.savefig(outputID + '_cdf')
    form_func.output_form(cdf, outputID + '_cdf.csv')
    module_logger.debug('Plot and csv saved to cwd')

    if interactive:
        plt.show()
    else:
        plt.clf()


def jitter(x, y, jitter_scale=0.001):
    '''Function takes in x and y values and jitters
    the identical points (x,y) in the x direction
    
    Parameters
    ----------
    x : ndarray
      x points
    y : ndarray
      y points
    jitter_scale : float
        The distance a point is jittered in the x direction

    Returns
    -------
    : ndarrays
        returns jittered x and y as separate ndarrays



    '''

    #TODO: Add option to jitter in y direction
    assert len(x) == len(y), "x and y must be the same length"
    jitt_x = np.copy(x)
    jitt_y = np.copy(y)
    xy_tuple = []
    for i in xrange(len(x)):
        xy_tuple.append((x[i], y[i]))

    xy_tuple = np.array(xy_tuple)
    unq_xy = np.unique(xy_tuple)
    for xy in unq_xy:
        indices = np.where(xy_tuple == xy)[0]
        if len(indices > 1):
            for ind in xrange(len(indices)):
                if ind > 0:
                    jitt_x[indices[ind]] += (jitter_scale * ind)
    return jitt_x, jitt_y

"""
