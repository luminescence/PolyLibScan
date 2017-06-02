from matplotlib.artist import setp
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import numba as nb
import pymc as mc
import numerics as num_

class Project(object):

    def scatter_plot(self, with_errors=False, with_labels=False, confidence_interval=0.95, ax=None, save_path=None, min_dist_to_ac=5):
        '''create a scatter plot with the probability
        of binding (x-axis) and the mean strength of binding 
        at the active site.

        input:
            ax: pyplot axis
            with_errors: bool
            confidence_interval: float
            min_dist_to_ac: float
            save_path: [string]

        output:
            none
        '''
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            plt.title(self.jobs[0].meta['protein'], size=20)
            
        if self.experimental_data is not None:
            polymer_list = list(set(self.endstate_matrix.columns.levels[0]) & set(self.experimental_data.index))
            experimental = self.experimental_data[polymer_list]
            experimental.sort_values(inplace=True)
        else:
            polymer_list = self.endstate_matrix.columns.levels[0]
            
        full_distance_matrix = self.endstate_matrix.swaplevel(0,1, axis=1)['Distance'][polymer_list]
        full_energy_matrix = self.endstate_matrix.swaplevel(0,1, axis=1)['Energy'][polymer_list]
        
        distance_matrix = full_distance_matrix[full_distance_matrix<min_dist_to_ac]
        energy_matrix = full_energy_matrix[full_distance_matrix<min_dist_to_ac]
            
        if with_errors:
            df1 = num_.distance_with_error(distance_matrix)
            df2 = num_.energy_with_error(energy_matrix, 0.95)
            results = pd.concat([df1, df2], axis=1)
            error = results['dist_max_error'].max()
        else:
            results = pd.DataFrame(index=distance_matrix.columns)
            results['energy_mean'] = energy_matrix.mean()
            results['dist_mean'] = distance_matrix.count()/float(full_distance_matrix.shape[0])
            error = 0.0

        if self.experimental_data is not None:
            results['color_by_inhibition'] = experimental.apply(lambda x:'b' if x>0 else 'r')
            results.plot(kind='scatter', x='dist_mean', y='energy_mean',
                         ax=ax, c=results.dropna()['color_by_inhibition'], s=100)
            if with_errors:
                ax.errorbar(results['dist_mean'] ,results['energy_mean'],
                            xerr=[results['dist_min_error'], results['dist_max_error']], 
                            yerr=[results['energy_min_error'], results['energy_max_error']], 
                            capsize=6, fmt=' ', color='grey', zorder=-1)
            inhib_leg = mpatches.Patch(color='red', label='not inhibiting')
            non_inhib_leg = mpatches.Patch(color='blue', label='inhibiting')
            ax.legend(handles=[inhib_leg, non_inhib_leg], fontsize=20, loc=1)
        else:
            results.plot(kind='scatter', x='dist_mean', y='energy_mean', 
                                  ax=ax, s=100)
        if with_labels:
            self._annotate(ax, results, 'dist_mean', 'energy_mean')
        ax.set_ylabel('energy', size=20)
        ax.set_xlabel('binding probability within %.2fA to active site' % round(min_dist_to_ac,1), size=20)
        ax.set_xlim([-0.2*results['dist_mean'].max(), 1.2*results['dist_mean'].max()+error])
        if save_path:
            plt.savefig(save_path)

    def _annotate(self, ax, df, x_name, y_name):
        for key, val in df.iterrows():
            ax.annotate(key, (val[x_name], val[y_name]),
                xytext=(5,-10), textcoords='offset points',
                family='sans-serif', fontsize=12, color='darkslategrey')

    def histogramm(self, min_dist_to_ac=5, ax=None, save_path=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            plt.title(self.jobs[0].meta['protein'], size=20)
        if self.experimental_data is not None:
            near_active_site = self.endstate_matrix.stack(0).loc[self.endstate_matrix.stack(0)['Distance']<min_dist_to_ac, :].unstack()['Energy']
            a = near_active_site.loc[:, self.experimental_data[near_active_site.columns].isnull()].mean()
            weights1 = np.ones_like(a)/len(a)
            a.hist(ax=ax, weights=weights1, label='not inhibiting', color='red')
            b = near_active_site.loc[:, self.experimental_data[near_active_site.columns].isnull()==False].mean()
            weights2 = np.ones_like(b)/len(b)
            b.hist(ax=ax, weights=weights2, label='inhibiting', alpha=0.7, color='blue')
            ax.legend(loc='upper right', fontsize=22)
        else:
            a = near_active_site.mean()
            a.hist(ax=ax)
        ax.set_xlabel('Energy', size=25)
        ax.set_ylabel('# of copolymer-types', size=25)
        if save_path:
            plt.savefig(save_path)

    def multibox_plot(self, experimental_data=None, ax=None, save_path=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            plt.title(self.jobs[0].meta['protein'], size=20)
        if experimental_data is None:
            experimental = self.experimental_data
        else:
             experimental = experimental_data[self.jobs[0].meta['protein']]
        energy_matrix = self.endstate_matrix.stack(0)['Energy'].unstack()
        col_substitutes = {columns: columns.replace('[x]', '') for columns in energy_matrix.columns}
        energy_matrix = energy_matrix.rename(columns=col_substitutes)
        shared_columns = list(set(energy_matrix.columns) & set(experimental.index))
        p = energy_matrix.loc[:, shared_columns].boxplot(ax=ax, return_type='dict')
        ax.set_xticks(rotation=90, size=18)
        setp(p['whiskers'], linewidth=2)
        setp([item for i,item in enumerate(p['whiskers']) if experimental[shared_columns].isnull()[i/2]], color='red')
        setp([item for i,item in enumerate(p['boxes']) if experimental[shared_columns].isnull()[i]], color='red')
        setp(p['boxes'], linewidth=4)
        setp(p['caps'], linewidth=4)
        ax.set_ylabel('Energy', size=20)
        ax.set_xlabel('Copolymer', size=20)
        inhib_leg = mpatches.Patch(color='red', label='Not inhibiting')
        non_inhib_leg = mpatches.Patch(color='blue', label='inhibiting')
        ax.legend(handles=[inhib_leg, non_inhib_leg], fontsize=22, loc=1)
        # fig.tight_layout()
        if save_path:
            plt.savefig(save_path)

    def plot_distance_density(self, cumulative=False, max_distance_range=None, ax=None, save_path=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            plt.title(self.jobs[0].meta['protein'], size=20)
        for name, poly_type in self.polymer_types.items():
            poly_type.distance_histogram(max_distance_range=max_distance_range, cumulative=cumulative, ax=ax)
        ax.set_xlabel('Distance $\AA$', size=20)
        ax.set_ylabel('Probability Density', size=20)
        ax.legend()
        if save_path:
            plt.savefig(save_path)

    def plot_results(self, min_dist_to_ac=5, save_path=None):
        fig = plt.figure(figsize=(18,12))
        ax = plt.subplot2grid((2,2),(0,0))
        self.scatter_plot(ax=ax, min_dist_to_ac=min_dist_to_ac)
        ax2 = plt.subplot2grid((2,2),(0,1))
        self.plot_bars(ax=ax2, distance_cutoff=min_dist_to_ac)
        ax3 = plt.subplot2grid((2,2),(1,0), colspan=2)
        self.plot_distance_density(ax=ax3, max_distance_range=min_dist_to_ac*10)
        if save_path:
            plt.savefig(save_path)

    def plot_bars(self, distance_cutoff=5, ax=None, save_path=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            plt.title(self.jobs[0].meta['protein'], size=20)
        width = 0.35
        index = np.arange(len(self.polymer_types.keys()))
        inhib_polys = filter(lambda x:self.polymer_types[x].ic50>=0,self.polymer_types.keys())
        non_inhib_polys = filter(lambda x:np.isnan(self.polymer_types[x].ic50),self.polymer_types.keys())
        poly_by_inhib = sorted(inhib_polys, key=lambda x:self.polymer_types[x].ic50, reverse=True) + non_inhib_polys
        ic_50s = [self.polymer_types[poly_name].ic50 for poly_name in poly_by_inhib]
        # experimental Data
        experimental = ax.bar(index, ic_50s, width)
        binding = [self.polymer_types[poly_name].cumulative_binding_probability(distance_cutoff) 
                    for poly_name in poly_by_inhib]
        # simulation data
        simulation = ax.bar(index + width, binding, width)
        ax.set_xticks(index + width / 2)
        ax.set_xticklabels(poly_by_inhib)
        ax.set_xlabel('Polymer Types', size=20)
        ax.set_ylabel('Inhibition/Binding', size=20)
        ax.legend([experimental, simulation], [r'$ic50^{-1}$ (Experiment)', 'Binding (Simulation)'])
        if save_path:
            plt.savefig(save_path)

class PolymerTypeSims(object):

    def scatter_plot(self, ax=None, save_path=None, with_error=False, min_dist_to_ac=5):
        '''create a scatter plot with the probability
        of binding (x-axis) and the mean strength of binding 
        at the active site.

        input:
            ax: pyplot axis
            save_path: [string]

        output:
            none
        '''
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            plt.title('%s - %s' % (self.sims[0].meta['protein'], self.polymer_type), size=20)
    
        digested_results = pd.DataFrame()
        if with_error:
            bayes_results_energy = {}
            for sim in self.sims:
                bayes_results_energy[sim.Id] = sim.bayes_results(model_name='energy')
            digested_results['Energy'] = pd.Series({key: val[1] for key,val in bayes_results_energy.items()})
            digested_results['EnergyErrMin'] = pd.Series({key: val[0] for key,val in bayes_results_energy.items()})
            digested_results['EnergyErrMax'] = pd.Series({key: val[2] for key,val in bayes_results_energy.items()})

            bayes_results_distance = {}
            for sim in self.sims:
                bayes_results_distance[sim.Id] = sim.bayes_results(model_name='distance')
            digested_results['BindingProbability'] = pd.Series({key: val[1] for key,val in bayes_results_distance.items()})
            digested_results['BindingProbabilityErrMin'] = pd.Series({key: val[0] for key,val in bayes_results_distance.items()})
            digested_results['BindingProbabilityErrMax'] = pd.Series({key: val[2] for key,val in bayes_results_distance.items()})
        else:
            data = self.data_frame()
            number_of_polymer_sims = data.shape[0]/len(data.index.get_level_values(0).unique())
            near_active_site = data.stack(0).loc[data.stack(0)['Distance']<min_dist_to_ac, :].unstack(0)
            digested_results['Energy'] = near_active_site['Energy'].mean()
            digested_results['BindingProbability'] = near_active_site['Energy'].count()/number_of_polymer_sims

        if self.project.parameters:
            digested_results['charge'] = pd.Series({sim.Id: sim.charge for sim in self.sims})
            digested_results.plot(kind='scatter', x='BindingProbability', y='Energy',
                                  ax=ax, c=digested_results['charge'], s=100, cmap='gist_rainbow')
        else:
            digested_results.plot(kind='scatter', x='BindingProbability', y='Energy', 
                                  ax=ax, s=100)
        if with_error:
            ax.errorbar(digested_results['BindingProbability'] ,digested_results['Energy'],
                            xerr=[digested_results['BindingProbabilityErrMin'], digested_results['BindingProbabilityErrMax']], 
                            yerr=[digested_results['EnergyErrMin'], digested_results['EnergyErrMax']], 
                            fmt=' ', color='grey', zorder=-1)
        ax.set_ylabel('Energy', size=25)
        ax.set_xlabel('Binding Probability within %.2fA to Active Site' % round(min_dist_to_ac,1), size=20)
        ax.set_xlim([0.0-0.2*digested_results['BindingProbability'].max(), 1.2*digested_results['BindingProbability'].max()])
        if save_path:
            plt.savefig(save_path)

    def distance_histogram(self, max_distance_range=None, cumulative=False, ax=None, save_path=None):
        # discretised bin container - hopefully nothing higher than
        # 200 Angstrom comes along
        
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            plt.title('%s - %s' % (self.sims[0].meta['protein'], self.polymer_type), size=20)
        if self.ic50:
            if self.ic50 >= 0:
                line_style = '-'
                tag = 'inhibiting'
            else:
                line_style = ':'
                tag = 'non-inhibiting'
        else:
            line_style = '-'
            tag = ' '
        if not max_distance_range:
            max_distance_range = self.distance_probability['distance'].shape[0]
        if cumulative:
            density = num_.cumulative_bins(self.distance_probability['density'])
        else:
            density = self.distance_probability['density']
        ax.plot(self.distance_probability['distance'][:max_distance_range], 
                density[:max_distance_range], 
                alpha=0.6, linestyle=line_style, lw=3,
                label='%s - %s' % (self.polymer_type, tag) )

        ax.set_xlabel('Distance [$\AA$]')
        ax.set_ylabel('Probability')
        if save_path:
            plt.savefig(save_path)



class Run(object):

    def plot_energies(self, ax=None, save_path=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            ax.set_title('%s - %s - Job Id: %d - Run Id: %d' % (self.meta()['protein'], 
                self.meta()['poly_name'],
                self.job.Id, self.Id), size=20)
        else:
            ax.set_title('Energies of Run %d' % self.Id , size=25)

        total_ene_data = self.total_energy()
        sim_time_in_ps = total_ene_data[:,1]/1000 * self.job.lmp_parameters['timestep']
        ax.plot(sim_time_in_ps, self.total_energy()[:,0], label='Total Energy')
        ax.plot(sim_time_in_ps, self.potential_energy()[:,0], label='Potential Energy')
        ax.plot(sim_time_in_ps, self.kinetic_energy()[:,0], label='Kinetic Energy')
        ax.set_ylabel('Energy', size=20)
        ax.set_xlabel('Time [ps]', size=20)
        legend1_info = ax.get_legend_handles_labels()

        ax2 = ax.twinx()
        ax2.plot(sim_time_in_ps, self.binding_energy()[:,0], color='m', label='Binding Energy')
        ax2.set_ylabel('Binding Energy', size=20)
        ax2.tick_params('y', color='m')
        legend2_info = ax2.get_legend_handles_labels()

        ax2.legend(legend1_info[0]+legend2_info[0], legend1_info[1]+legend2_info[1],)
        if save_path:
            plt.savefig(save_path)

    def plot_temperature(self, ax=None, save_path=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            ax.set_title('%s - %s - Job Id: %d - Run Id: %d' % (self.meta()['protein'], 
                self.meta()['poly_name'],
                self.job.Id, self.Id), size=20)
        else:
            ax.set_title('Temperature of Run %d' % self.Id , size=25)
        data = self.temperature()
        time_series = data[:,1] * self.job.lmp_parameters['timestep']
        ax.plot(time_series, data[:, 0], label='Temperature')
        ax.set_ylabel('Temperture [K]', size=20)
        ax.set_xlabel('Time [ps]', size=20)
        legend1_info = ax.get_legend_handles_labels()

        if save_path:
            plt.savefig(save_path)