from matplotlib.artist import setp
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import numba as nb
import pymc as mc
import numerics as num_
try:
    import adjustText as at
    adjust_text_enabled = True
except ImportError:
    adjust_text_enabled = False

class Project(object):

    def scatter_plot(self, subset=None, with_errors=False, with_labels=False, with_crossvalidation=False, 
                           confidence_interval=0.95, ax=None, save_path=None, min_dist_to_ac=10, property_='charge',
                           ignore_experiment=False, label_only_misclassified=False, error_method='bootstrap'):
        '''create a scatter plot with the probability
        of binding (x-axis) and the mean strength of binding 
        at the active site.

        input:
            ax: pyplot axis
            with_errors: bool
            with_labels: bool
            with_crossvalidation: bool
            confidence_interval: float
            ax: pyplot subplot object
            save_path: [string]
            min_dist_to_ac: float

        output:
            none
        '''
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            ax.set_title('%s (PDB: %s)' % (self.jobs[0].meta['protein_name'], 
                                           self.jobs[0].meta['protein'].upper()), size=20)
        
        results = self._scatter_data(subset=subset, 
                                     with_errors=with_errors, 
                                     with_labels=with_labels, 
                                     label_only_misclassified=label_only_misclassified,
                                     with_crossvalidation=with_crossvalidation, 
                                     property_=property_,
                                     confidence_interval=confidence_interval, 
                                     min_dist_to_ac=min_dist_to_ac,
                                     ignore_experiment=ignore_experiment,
                                     error_method=error_method)

        if with_errors:
            error = results['dist_max_error'].max()
        else:
            error = 0.0

        if self.experimental_data is not None:
            c_min = results['property'].min()
            c_max = results['property'].max()
            cm = plt.cm.get_cmap('bwr_r')
            non_inhib = results[results.color_by_inhibition=='r']

            legend_items = []
            if property_:
                ax.scatter(non_inhib['dist_mean'], non_inhib['energy_mean'], edgecolor='k', c=non_inhib['property'],
                       vmin=c_min, vmax=c_max, s=100, marker='o', cmap=cm)
                inhibitors = results[results.color_by_inhibition=='b']
                plot = ax.scatter(inhibitors['dist_mean'], inhibitors['energy_mean'], edgecolor='k', c=inhibitors['property'],
                           vmin=c_min, vmax=c_max, s=100, marker='s', cmap=cm)
                ax.figure.colorbar(plot)
                cbar = ax.figure.axes[1]
                cbar.axes.get_yaxis().labelpad = 20
                cbar.axes.set_ylabel(property_, rotation=270, fontsize=25)
                legend_items.append(mlines.Line2D([], [], markeredgecolor='k', color='w', marker='s', linestyle='None', markersize=10,
                                             label='inhibiting'))
                legend_items.append(mlines.Line2D([], [], markeredgecolor='k', color='w', marker='o', linestyle='None', markersize=10,
                                             label='not inhibiting'))
            else:
                results.plot(kind='scatter', x='dist_mean', y='energy_mean', alpha=0.7, 
                         ax=ax, c=results.dropna()['color_by_inhibition'], s=100) 
                legend_items = [mpatches.Patch(color='red', label='not inhibiting')] 
                legend_items.append(mpatches.Patch(color='blue', label='inhibiting')) 

            if with_errors:
                ax.errorbar(results['dist_mean'] ,results['energy_mean'],
                            xerr=[results['dist_min_error'], results['dist_max_error']], 
                            yerr=[results['energy_min_error'], results['energy_max_error']], 
                            capsize=6, fmt=' ', color='grey', zorder=-1)
            if with_crossvalidation:
                classification = results['color_by_inhibition'].apply(lambda x:x=='b')
                roc_auc_score = self._roc_auc(classification, results['probabilities'])
                kappa = self._kappa(classification, results['model_predictions'])
                # plotting black x on false predictions
                if np.any(~results['true_predictions']):
                    results[~results['true_predictions']].plot(kind='scatter', x='dist_mean', 
                                                       y='energy_mean', ax=ax, marker='x', c='lightgreen', s=35, linewidth=2)
                legend_items.append(mlines.Line2D([], [], markeredgecolor='lightgreen', color='w', marker='x', 
                                                         linewidth=3, linestyle='None', 
                                                         markeredgewidth=3, markersize=10, label='misclassified'))
                legend_items.append(mpatches.Patch(color='white', label='{:<9s} {:>4.2f}'.format('ROC-AUC:', roc_auc_score) ))
                legend_items.append(mpatches.Patch(color='white', label='{:<9s} {:>4.2f}'.format('kappa:', kappa) ))

            ax.legend(handles=legend_items, loc='best', prop={'size': 20, 'family': 'monospace'})
        else:
            results.plot(kind='scatter', x='dist_mean', y='energy_mean', alpha=0.7, 
                         ax=ax, s=100)
        if with_labels:
            self._annotate(ax, results, 'dist_mean', 'energy_mean', only_misclassified=label_only_misclassified)
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.set_ylabel('Energy', size=25)
        ax.set_xlabel(r'Binding probability within $%d\AA$ to interface' % round(min_dist_to_ac,1), size=25)
        ax.set_xlim([-0.2*results['dist_mean'].max(), 1.2*results['dist_mean'].max()+error])
        if save_path:
            plt.savefig(save_path)

    def _annotate(self, ax, df, x_name, y_name, only_misclassified=False):
        if only_misclassified and 'true_predictions' in df.columns:
            data = df[~df['true_predictions']]
        else:
            data =  df
        texts = [ax.text(val[x_name], val[y_name], key, fontsize=14) 
                                    for key, val in data.iterrows()]
        if adjust_text_enabled:
            at.adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black'), 
                                  expand_points=(1.2, 1.75),
                                  force_points=0.5)

    def histogramm(self, min_dist_to_ac=5, ax=None, save_path=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            ax.set_title('%s (PDB: %s)' % (self.jobs[0].meta['protein_name'], 
                                           self.jobs[0].meta['protein']), size=20)
        near_active_site = self.endstate_matrix.stack(0).loc[self.endstate_matrix.stack(0)['Distance'] < min_dist_to_ac, :].unstack()['Energy']
        if self.experimental_data is not None:
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
            ax.set_title('%s (PDB: %s)' % (self.jobs[0].meta['protein_name'], 
                                           self.jobs[0].meta['protein']), size=20)
        if experimental_data is None:
            experimental = self.experimental_data
        else:
             experimental = experimental_data[self.jobs[0].meta['protein']]
        energy_matrix = self.endstate_matrix.stack(0)['Energy'].unstack()
        col_substitutes = {columns: columns.replace('[x]', '') for columns in energy_matrix.columns}
        energy_matrix = energy_matrix.rename(columns=col_substitutes)
        shared_columns = list(set(energy_matrix.columns) & set(experimental.index))
        p = energy_matrix.loc[:, shared_columns].boxplot(ax=ax, return_type='dict')
        # ax.set_xticks(rotation=90, size=18)
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
            ax.set_title('%s (PDB: %s)' % (self.jobs[0].meta['protein_name'], 
                                           self.jobs[0].meta['protein']), size=20)
        for name, poly_type in self.polymer_types.items():
            poly_type.distance_histogram(max_distance_range=max_distance_range, cumulative=cumulative, ax=ax)
        ax.set_xlabel('Distance $\AA$', size=20)
        ax.set_ylabel('Probability Density', size=20)
        ax.legend(loc=2, ncol=2)
        if save_path:
            plt.savefig(save_path)

    def plot_results(self, min_dist_to_ac=5, save_path=None):
        fig = plt.figure(figsize=(18,12))
        ax = plt.subplot2grid((2,2),(0,0))
        self.scatter_plot(ax=ax, min_dist_to_ac=min_dist_to_ac)
        ax2 = plt.subplot2grid((2,2),(0,1))
        self.plot_bars(ax=ax2, distance_cutoff=min_dist_to_ac)
        ax3 = plt.subplot2grid((2,2),(1,0), colspan=2)
        self.plot_distance_density(ax=ax3, max_distance_range=min_dist_to_ac)
        if save_path:
            plt.savefig(save_path)

    def plot_bars(self, distance_cutoff=10, ax=None, save_path=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            ax.set_title('%s (PDB: %s)' % (self.jobs[0].meta['protein_name'], 
                                           self.jobs[0].meta['protein']), size=20)
        
        model_inhibition = self._inhibition_series(distance_cutoff)
        ic50_data = self._ic50_series()
        
        combined = pd.DataFrame(data={'model': model_inhibition, 'ic50': ic50_data})
        combined.sort(columns=['ic50'], ascending=False, inplace=True)
        
        colors = self._bar_colors(combined)
        # set unknowns (marked with -1) to NaN
        combined.loc[combined['ic50']==-1, 'ic50'] = np.nan
        
        combined.plot(kind='bar', ax=ax, color=colors)
        ax.set_xlabel('Polymer Types', size=20)
        ax.set_ylabel('Inhibition/Binding', size=20)
        if save_path:
            plt.savefig(save_path)

    def _bar_colors(self, data):
        ic50_inhibitor_color = 'steelblue'
        inhibitor_color = 'g'
        non_inhibitors = 'r'
        unknowns = 'k'
        color_list = []
        for name,d in data.iterrows():
            if d.ic50 > 0:
                color_list.append(ic50_inhibitor_color)
                color_list.append(inhibitor_color)
            elif np.isnan(d.ic50):
                color_list.append(non_inhibitors)
                color_list.append(non_inhibitors)
            else:
                color_list.append(unknowns)
                color_list.append(unknowns)
        return color_list

    def plot_experimental_model_comparison(self, distance, ax=None, save_path=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            ax.set_title('%s (PDB: %s)' % (self.jobs[0].meta['protein_name'], 
                                           self.jobs[0].meta['protein']), size=20)
        
        model_inhibition = self._inhibition_series(distance)
        ic50_data = self._ic50_series()
        
        combined = pd.DataFrame(data={'model': model_inhibition, 'ic50': ic50_data})
        combined.plot(kind='scatter', x='ic50', y='model')
        if save_path:
            plt.savefig(save_path)

    def _ic50_series(self):
        ic_50s = {poly_name: self.polymer_types[poly_name].ic50 for poly_name in self.polymer_types}
        return pd.Series(ic_50s).sort(inplace=False, ascending=False)

    def _inhibition_series(self, distance_cutoff):
        data = {name: polytype.cumulative_binding_probability(distance_cutoff) 
                    for name,polytype in self.polymer_types.iteritems()}
        return pd.Series(data)

    def plot_inhibition(self, distance=10, ax=None, save_path=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            ax.set_title('%s (PDB: %s)' % (self.jobs[0].meta['protein_name'], 
                                           self.jobs[0].meta['protein']), size=20)
        data = self._inhibition_series(distance_cutoff=distance)
        data.plot(kind='bar', ax=ax)
        if save_path:
            plt.savefig(save_path)

    def plot_ic50(self, ax=None, save_path=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(18,12))
            ax.set_title('%s (PDB: %s)' % (self.jobs[0].meta['protein_name'], 
                                           self.jobs[0].meta['protein']), size=20)
        
        data = self._ic50_series()
        data_with_experimental_value = data[data!=-1]
        if len(data_with_experimental_value) == 0:
            raise ValueError('No PolymerTypes with known inhibition in project.')
        data_with_experimental_value.plot(kind='bar', ax=ax)
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
            ax.set_title('Poly Type: %s | Protein %s (PDB: %s)' % (
                self.name,
                self.sims[0].meta['protein_name'], 
                self.sims[0].meta['protein'].upper()
            ), 
            size=20)
    
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
            ax.set_title('%s (PDB: %s)' % (self.sims[0].meta['protein_name'], self.sims[0].meta['protein']), size=20)
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
        else:
            max_distance_range = int(max_distance_range*10)
        if cumulative:
            density = num_.cumulative_bins(self.distance_probability['density'])
        else:
            density = self.distance_probability['density']

        ax.plot(self.distance_probability['distance'][:max_distance_range], 
                density[:max_distance_range], 
                alpha=0.6, linestyle=line_style, lw=3,
                label='%s - %s' % (self.name, tag) )

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
