import numpy as np
import tqdm
import pandas as pd
import plotting
import bayesModels as bayes
import concurrent.futures as concurrent
import collections as col
import pymol_visualisation as pym
import numerics

class PolymerTypeSims(plotting.PolymerTypeSims, bayes.PolymerTypeSims):
    '''Stores all simulation runs of one polymer type.
    '''
    def __init__(self, project, polymer_sims, ic50=None, with_pymol=True):
        self.project = project
        self.sims = polymer_sims
        self.register(self.sims)
        self.ic50 = ic50
        self._distance_probability = None
        self._energy_distance_distribution = None
        self.name = self.sims[0].meta['poly_name']
        self.weights = self.sims[0].weights
        if with_pymol:
            self.pymol = pym.PymolVisPolyType(self, protein_path=self.project.protein_path)
        super(PolymerTypeSims, self).__init__()

    def distance_probability():
        doc = "The distance_probability property."
        def fget(self):
            if self._energy_distance_distribution == None:
                results = self._calc_distance_distribution(self.sims)
                self.distance_probability = results[['distance', 'density']]
                self.energy_distance_distribution = results[['distance', 'energy']]
            return self._distance_probability
        def fset(self, value):
            self._distance_probability = value
        def fdel(self):
            self._distance_probability = None
        return locals()
    distance_probability = property(**distance_probability())


    def energy_distance_distribution():
        doc = "The energy_distance_distribution property."
        def fget(self):
            if self._energy_distance_distribution == None:
                results = self._calc_distance_distribution(self.sims)
                self.distance_probability = results[['distance', 'density']]
                self.energy_distance_distribution = results[['distance', 'energy']]
            return self._energy_distance_distribution
        def fset(self, value):
            self._energy_distance_distribution = value
        def fdel(self):
            self._energy_distance_distribution = None
        return locals()
    energy_distance_distribution = property(**energy_distance_distribution())

    def _calc_distance_distribution(self, sims, thread_num=2):
        dist_v = np.zeros(2000, dtype=[('count', np.int), ('energy', np.float)])
        with concurrent.ThreadPoolExecutor(max_workers=thread_num) as executor:
            for job in tqdm.tqdm(sims, desc='Calculating distance distribution of Simulation'):
                job_count = job.distance_frequency['frequency']
                dist_v['count'][:len(job_count)] += job_count
                dist_v['energy'][:len(job_count)] += job.energy_distance_distribution['energy']
        # Shortening the array by trailing zeros
        reduced_hist = numerics.discard_tailing_zeros(dist_v['count'])
        reduced_len = reduced_hist.shape[0]
        
        # Mean Energy
        energy_mean = dist_v['energy'][:reduced_len] / len(self.sims)
        # set elements to zero, who are devided by zero
        energy_mean[np.isnan(energy_mean)] = 0.0

        # Calculating norm
        n = np.sum(reduced_hist).astype(np.float)
        hist_normed = reduced_hist/n
        # results
        results_dtype = [('distance', np.float), ('density', np.float),
                         ('energy', np.float)]
        results = np.empty(reduced_len, results_dtype)
        results['distance'] = np.arange(reduced_len, dtype=np.float)/10
        results['density'] = hist_normed
        results['energy'] = energy_mean
        return results

    def cumulative_binding_probability(self, distance):
        return numerics.cdf(self.distance_probability['density'], distance)

    def cumulative_mean_energy(self, distance):
        mean_energy = (numerics.cdf(self.distance_probability['density'] * self.energy_distance_distribution['energy'], distance)
                                   / self.cumulative_binding_probability(distance))
        return mean_energy

    def distance_cumulative_probability(self):
        '''Cumulative distance function: 
        Return a 2xN Array that denotes the cumulative density within a 
        given radius around the active site.
        '''
        cum_distance = np.empty_like(self.distance_probability)
        cum_distance['distance'] = self.distance_probability['distance'] + 0.1
        cum_distance['density'] = numerics.cumulative_bins(self.distance_probability['density'])
        return cum_distance

    def register(self, children):
        for child in children:
            child.poly_type = self

    def data_frame(self):
        sub_frames = np.empty(len(self.sims), dtype=object)
        for i, job in enumerate(self.sims):
            sub_frames[i] = job.to_dataFrame()
            job._parse.close()
        data_frame = pd.concat(sub_frames)
        col_index = pd.MultiIndex.from_product([[self.name], data_frame.columns], 
                                                names=['PolyType', 'Results'])
        data_frame.set_axis(1, col_index)
        return data_frame

    def __repr__(self):
        out  = []
        out += ['Polymer Type %s.' % self.name]
        out += ['Containing %d Jobs.' % len(self.sims)]
        return '\n'.join(out)

    def __str__(self):
        out  = []
        out += ['Polymer Type %s.' % self.name]
        out += ['Containing %d Jobs.' % len(self.sims)]
        out += ['Composition:']
        for name,weight in self.weights.items():
            monomer_charge = self.monomer_charge_for_all_beads(name)
            out += ['  %s: %0.2f (charge: %d)' % (
                name, weight, monomer_charge)]
        out += ['Average Charge: %0.2f (monomer), %0.2f (sequence)' % self.charge_average()]
        if self.ic50:
            out += ['Experimental inhibition: %0.2f (inverse IC50)' % self.ic50]
        return '\n'.join(out)

    def charge_average(self):
        '''Return the average charges of the polymer type.
        "monomer" gives the average charge per monomer.
        "sequence" returns the average charge from the generates sequences
        '''
        ChargeAverage = col.namedtuple('ChargeAverage', 'monomer, sequence')

        per_monomer = sum((self.monomer_charge_for_all_beads(name) * weight
                           for name,weight in self.weights.items())) / sum(self.weights.values())
        sequence = sum((sim.charge for sim in self.sims))/len(self.sims)

        return ChargeAverage(monomer=per_monomer, sequence=sequence)

    def monomer_charge_for_all_beads(self, name):
        monomer_charge = 0
        for beads in self.project.parameters['Atoms'][name].keys():
            monomer_charge += self.project.parameters['Atoms'][name][beads]['charge']
        return monomer_charge