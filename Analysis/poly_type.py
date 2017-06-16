import numpy as np
import tqdm
import pandas as pd
import plotting
import bayesModels as bayes
import concurrent.futures as concurrent
import numerics

class PolymerTypeSims(plotting.PolymerTypeSims, bayes.PolymerTypeSims):
    '''Stores all simulation runs of one polymer type.
    '''
    def __init__(self, project, polymer_sims, ic50=None):
        self.project = project
        self.sims = polymer_sims
        self.register(self.sims)
        self.ic50 = ic50
        self._distance_probability = None
        self._energy_distance_distribution = None
        self.polymer_type = self.sims[0].meta['poly_name']
        self.weights = self.sims[0].weights
        super(PolymerTypeSims, self).__init__()

    def distance_probability():
        doc = "The distance_probability property."
        def fget(self):
            if self._energy_distance_distribution == None:
                results = self._calc_distance_distribution(self.sims)
                self._distance_probability = results[['distance', 'density']]
                self._energy_distance_distribution = results[['distance', 'energy']]
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
                self._distance_probability = results[['distance', 'density']]
                self._energy_distance_distribution = results[['distance', 'energy']]
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
        col_index = pd.MultiIndex.from_product([[self.polymer_type], data_frame.columns], 
                                                names=['PolyType', 'Results'])
        data_frame.set_axis(1, col_index)
        return data_frame

    def __repr__(self):
        out  = []
        out += ['Polymer Type %s.' % self.polymer_type]
        out += ['Containing %d Jobs.' % len(self.sims)]
        return '\n'.join(out)
