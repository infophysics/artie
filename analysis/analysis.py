"""
Set of analysis functions for ARTIE sim/data
"""
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import uproot
import sys
import os
import logging

# set up logger
logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("artie_analysis.log"),
        logging.StreamHandler(sys.stdout)
    ]
)

simulation_variables = [
    'gen_energy',
    'arrival_time',
    'arrival_e',
    'num_elastic',
    'num_inelastic',
    'num_ncapture',
    'num_fission',
    'num_scatter',
    'num_scatout',
    'z_scatter',
    'first_gas',
    'max_dphi',
    'max_dp',
    'max_de',
]

class ArtieData:
    """
    Class for storing ARTIE data.
    """
    def __init__(self,
        input_file: str,
        variables:  list=simulation_variables,
    ):
        logging.info(f"Attempting to open file {input_file}.")
        if not os.path.isfile(input_file):
            logging.error(f"Specified input directory '{input_file}' not in path.")
            raise ValueError(f"Specified input directory '{input_file}' not in path.")
        logging.info(f"Loading data from input file: {input_file}.")
        try:
            self.input_file = uproot.open(input_file)['artie']
        except Exception as e:
            logging.error(f"Error loading file {input_file}: {e}")
            raise ValueError(f"Error loading file {input_file}: {e}")
        self.variables = {}
        # check that all simulation variables are present
        for item in variables:
            if item not in self.input_file.keys():
                logging.error(f"Variable {item} not present in input file: {input_file}.")
                raise ValueError(f"Variable {item} not present in input file: {input_file}.")
            else:
                logging.info(f"Adding variable '{item}' to variables.")
                self.variables[item] = self.input_file[item].array(library="np")
        self.num_events = len(self.variables['gen_energy'])
        logging.info(f"Loaded {self.num_events} events in {len(self.variables)} variables.")
        if 'arrival_time' in self.variables.keys():
            self.target_mask = np.where(self.variables['arrival_time'] > 0)
            self.num_reach_target = len(self.variables['arrival_time'][self.target_mask])
            logging.info(f"Out of {self.num_events} events, {self.num_reach_target} reached the target.")
        

    def __len__(self):
        """Return the number of events"""
        return self.num_events

class ArtieAnalysis:
    """
    Class for handling the analysis of ARTIE data.  The user
    must pass in a ArtieData dictionary, for example one might do:

        >>> vacuum_data = ArtieData(...)
        >>> argon_data  = ArtieData(...)
        >>> data_dict = {
                'vacuum_data':  vacuum_data,
                'argon_data':   argon_data
            }
        >>> analysis = ArtieAnalysis(data_dict)

    """
    def __init__(self,
        data_dict
    ):
        """
        Essentially we just copy the dictionary.  The first argument
        to essentially all functions below must be the name
        of the dataset(s) you want to use.
        """
        self.artie_data = data_dict

    def check_for_dataset(self,
        dataset:    str,
    ):
        """Check if 'dataset' is in the data dictionary"""
        if dataset not in self.artie_data.keys():
            logging.error(f"Specified dataset '{dataset}' not in analysis set: {self.artie_data.keys()}.")
            raise ValueError(f"Specified dataset '{dataset}' not in analysis set: {self.artie_data.keys()}.")

    def convert_tof_to_energy(self,
        tof,
    ):
        """
        Convert a time of flight value to energy using the formula:
            E = m(gamma - 1),           where
            gamma = 1/sqrt(1 - beta^2), and
            beta = v/c = L/ct
        """
        neutron_mass    = 939.5654133   # MeV 
        speed_c         = .299792458    # m / ns
        artie_length    = 70.0          # m
        beta_sq         = np.power((artie_length / (tof * speed_c)), 2.0)
        gamma           = 1.0 / np.sqrt(1.0 - beta_sq)
        return (gamma - 1.0) * neutron_mass    

    def smear_time(self,
        tof:  float,
        resolution: float=0.125
    ):
        """
        Smear a time value(s) with a uniform random number
        """
        smear = np.random.uniform(
            -resolution/2.0, 
            resolution/2.0, 
            size=tof.size
        ) 
        return tof + smear

    def simulate_tof_data_from_mc(self,
        dataset:    str,
        resolution: float=125   # ns
    ):
        """
        Simulate the effect of the detector by smearing the tof.
        The input resolution should be in [ns]
        """
        self.check_for_dataset(dataset)
        logging.info(f"Generating smeared tof and energy.")
        tof = self.artie_data[dataset].variables['arrival_time']
        # now, throw away any tof's = 0
        tof = tof[(tof > 0)]
        sim_tof = self.smear_time(tof,resolution)
        sim_energy  = self.convert_tof_to_energy(sim_tof)
        logging.info(f"Simulated {len(sim_energy)} events by smearing tof.")
        return sim_tof, sim_energy
    
    def plot_energy_spectrum(self,
        dataset:    str,
        plot_type:  str='target',
        num_bins:   int=200,
        energy_min: float=40e-3, # MeV
        energy_max: float=70e-3, # MeV
        max_events: int=1000000,
        title:      str='ARTIE: Energy Spectrum',
        save:       str='',
        show:       bool=True,
    ):
        """
        Plot the energy spectrum for all events, and events which reach the target.
        The dataset argument can be a single string: 'vacuum', or a list: ['vacuum', 'argon'].
        The plot_type argument can be either 'target', 'all', or 'both'.
        """
        fig, axs = plt.subplots()
        if plot_type not in ['target','all','both']:
            logging.warning(f"Specified plot type: {plot_type} not allowed, using 'target'.")
            plot_type = 'target'
        if type(dataset) == str:
            self.check_for_dataset(dataset)
            energy_all = self.artie_data[dataset].variables['gen_energy'][:max_events]
            energy_target = energy_all[self.artie_data[dataset].target_mask]
            if plot_type == 'all' or plot_type == 'both':
                axs.hist(energy_all, bins=num_bins, histtype='step', label=dataset+'-produced', density=True)
            if plot_type == 'target' or plot_type == 'both':
                axs.hist(energy_target, bins=num_bins, histtype='step', label=dataset+'-detected', density=True)
        else:
            for item in dataset:
                self.check_for_dataset(item)
                energy_all = self.artie_data[item].variables['gen_energy'][:max_events]
                energy_target = energy_all[self.artie_data[item].target_mask]
                if plot_type == 'all' or plot_type == 'both':
                    axs.hist(energy_all, bins=num_bins, histtype='step', label=item+'-produced', density=True)
                if plot_type == 'target' or plot_type == 'both':
                    axs.hist(energy_target, bins=num_bins, histtype='step', label=item+'-detected', density=True)
        axs.set_xlabel(f"Energy (MeV) [{energy_min},{energy_max}]")
        axs.set_ylabel("Neutrons")
        axs.set_title(title)
        plt.legend()
        plt.tight_layout()
        if save != '':
            plt.savefig('plots/' + save + '.png')
        if show:
            plt.show()
    
    def plot_transmission_coeff(self,
        dataset:    str,
        num_bins:   int=200,
        energy_min: float=40e-3, # MeV
        energy_max: float=70e-3, # MeV
        max_events: int=1000000,
        title:      str='ARTIE: Transmission Coefficient vs. Energy',
        save:       str='',
        show:       bool=True,
    ):
        """
        Plot the transmission coefficient as a function of energy
        """
        fig, axs = plt.subplots()
        if type(dataset) == str:
            self.check_for_dataset(dataset)
            energy_all = self.artie_data[dataset].variables['gen_energy'][:max_events]
            energy_target = energy_all[self.artie_data[dataset].target_mask]
            hall, edges = np.histogram(energy_all, bins=num_bins,range=[energy_min, energy_max])
            hpss, edges = np.histogram(energy_target, bins=num_bins,range=[energy_min, energy_max])
            cbins = (edges[1:] + edges[:-1])/2.0
            trans_ratio = np.zeros(cbins.size)
            trans_mask = (hall > 0)
            # find trans ratio
            trans_ratio[trans_mask] = hpss[trans_mask] / hall[trans_mask]
            axs.plot(cbins, trans_ratio, linestyle='--', label=dataset+'-transmission ratio')
        else:
            for item in dataset:
                self.check_for_dataset(item)
                energy_all = self.artie_data[item].variables['gen_energy'][:max_events]
                energy_target = energy_all[self.artie_data[item].target_mask]
                hall, edges = np.histogram(energy_all, bins=num_bins,range=[energy_min, energy_max])
                hpss, edges = np.histogram(energy_target, bins=num_bins,range=[energy_min, energy_max])
                cbins = (edges[1:] + edges[:-1])/2.0
                trans_ratio = np.zeros(cbins.size)
                trans_mask = (hall > 0)
                # find trans ratio
                trans_ratio[trans_mask] = hpss[trans_mask] / hall[trans_mask]
                axs.plot(cbins, trans_ratio, linestyle='--', label=item+'-transmission ratio')
        axs.set_xlabel(f"Energy (MeV) [{energy_min},{energy_max}]")
        axs.set_ylabel(f"Transmission Ratio")
        axs.set_title(title)
        plt.legend()
        plt.tight_layout()
        if save != '':
            plt.savefig('plots/' + save + '.png')
        if show:
            plt.show()
    
    def plot_tof(self,
        dataset:    str,
        plot_einstein:  bool=True,
        num_bins:   int=200,
        max_events: int=1000000,
        title:      str='ARTIE: Time-of-flight vs. Energy',
        save:       str='',
        show:       bool=True,
    ):
        fig, axs = plt.subplots()
        if type(dataset) == str:
            self.check_for_dataset(dataset)
            tof_all = self.artie_data[dataset].variables['arrival_time'][:max_events]
            energy_all = self.artie_data[dataset].variables['gen_energy'][:max_events] 
            tof_target = tof_all[self.artie_data[dataset].target_mask]
            energy_target = energy_all[self.artie_data[dataset].target_mask]
            # generate fake tof to compare
            tof_compare = np.linspace(min(tof_target),max(tof_target), num_bins)
            tof_indices = np.digitize(tof_target, tof_compare)
            energy_target_bins = [[] for ii in range(num_bins)]
            for ii in range(len(energy_target)):
                energy_target_bins[tof_indices[ii]-1].append(energy_target[ii])
            energy_target_means = [np.mean(e_bins) for e_bins in energy_target_bins]
            energy_target_stds  = [np.std(e_bins)/np.sqrt(len(e_bins)) for e_bins in energy_target_bins]
            axs.errorbar(
                energy_target_means, 
                tof_compare, 
                yerr=energy_target_stds, 
                linestyle='--',  
                label=dataset+'-geant4', 
                capsize=2.0
            )
        else:
            for item in dataset:
                self.check_for_dataset(item)
                tof_all = self.artie_data[item].variables['arrival_time'][:max_events]
                energy_all = self.artie_data[item].variables['gen_energy'][:max_events] 
                tof_target = tof_all[self.artie_data[item].target_mask]
                energy_target = energy_all[self.artie_data[item].target_mask]
                # generate fake tof to compare
                tof_compare = np.linspace(min(tof_target),max(tof_target), num_bins)
                tof_indices = np.digitize(tof_target, tof_compare)
                energy_target_bins = [[] for ii in range(num_bins)]
                for ii in range(len(energy_target)):
                    energy_target_bins[tof_indices[ii]-1].append(energy_target[ii])
                energy_target_means = [np.mean(e_bins) for e_bins in energy_target_bins]
                energy_target_stds  = [np.std(e_bins)/np.sqrt(len(e_bins)) for e_bins in energy_target_bins]
                axs.errorbar(
                    energy_target_means, 
                    tof_compare, 
                    yerr=energy_target_stds, 
                    linestyle='--', 
                    label=item+'-geant4', 
                    capsize=2.0
                )
        if plot_einstein:
            energy_compare = self.convert_tof_to_energy(tof_compare)
            axs.plot(energy_compare, tof_compare, linestyle='--', color='r', label='Einstein')
        axs.set_xlabel(f"Energy (MeV)")
        axs.set_ylabel(f"Time-of-flight (ns)")
        axs.set_title(title)
        plt.legend()
        plt.tight_layout()
        if save != '':
            plt.savefig('plots/' + save + '.png')
        if show:
            plt.show()

    def plot_tof_comparison(self,
        dataset:    str,
        num_bins:   int=200,
        max_events: int=1000000,
        title:      str='ARTIE: Energy Comparison',
        save:       str='',
        show:       bool=True,
    ):
        fig, axs = plt.subplots()
        if type(dataset) == str:
            self.check_for_dataset(dataset)
            tof_all = self.artie_data[dataset].variables['arrival_time'][:max_events]
            energy_all = self.artie_data[dataset].variables['gen_energy'][:max_events] 
            tof_target = tof_all[self.artie_data[dataset].target_mask]
            energy_target = energy_all[self.artie_data[dataset].target_mask]
            # compute energy from tof
            energy_target_tof = self.convert_tof_to_energy(tof_target)
            # generate fake tof to compare
            tof_compare = np.linspace(min(tof_target),max(tof_target))
            energy_compare = self.convert_tof_to_energy(tof_compare)
            # compare the truth
            delta = energy_target_tof - energy_target
            # create plots
            axs.hist(delta, bins=num_bins, label=dataset)
        else:
            for item in dataset:
                self.check_for_dataset(item)
                tof_all = self.artie_data[item].variables['arrival_time'][:max_events]
                energy_all = self.artie_data[item].variables['gen_energy'][:max_events] 
                tof_target = tof_all[self.artie_data[item].target_mask]
                energy_target = energy_all[self.artie_data[item].target_mask]
                # compute energy from tof
                energy_target_tof = self.convert_tof_to_energy(tof_target)
                # generate fake tof to compare
                tof_compare = np.linspace(min(tof_target),max(tof_target))
                energy_compare = self.convert_tof_to_energy(tof_compare)
                # compare the truth
                delta = energy_target_tof - energy_target
                # create plots
                axs.hist(delta, bins=num_bins, label=item)
        axs.set_xlabel(rf"Energy (measured - actual) (MeV) ($\mu$,$\sigma$): ({np.mean(delta):.2e}, {np.std(delta):.2e})")
        axs.set_title(title)
        plt.legend()
        plt.tight_layout()
        if save != '':
            plt.savefig('plots/' + save + '.png')
        if show:
            plt.show()

    def plot_xsec(self,
        max_events: int=1000000,
    ):
        # collect time and energy information
        tof_all = self.variables['arrival_time'][:max_events]
        energy_all = self.variables['gen_energy'][:max_events] 
        tof_target = tof_all[self.target_mask]
        energy_target = energy_all[self.target_mask]

if __name__ == "__main__":

    vacuum_data = ArtieData(
        input_file="../build/vacuum.root",
        variables=simulation_variables,
    )
    argon_data = ArtieData(
        input_file="../build/argon.root",
        variables=simulation_variables,
    )
    data_dict = {
        'vacuum':   vacuum_data,
        'argon':    argon_data,
    }

    analysis = ArtieAnalysis(data_dict)

    analysis.plot_energy_spectrum(
        dataset=['vacuum','argon'],
        plot_type='target',
        save='example_energy_spectrum'
    )
    analysis.plot_transmission_coeff(
        dataset=['vacuum','argon'],
        save='example_transmission_coeff'
    )
    analysis.plot_tof(
        dataset=['vacuum','argon'],
        plot_einstein=False,
        save='example_tof'
    )
    analysis.plot_tof_comparison(
        dataset=['vacuum','argon'],
        save='example_tof_comparison'
    )