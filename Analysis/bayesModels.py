import pymc as mc
import pandas as pd
import numpy as np
from sklearn import linear_model, metrics
from sklearn.preprocessing import StandardScaler

from PolyLibScan.helpers.jupyter_compatibility import agnostic_tqdm


class Project(object):

    def bayes_model_setup(self, distance=5):
        for p_type in self.polymer_types.values():
            p_type.model_setup(distance)

    def bayes_sampling(self, sampling_rate=1000, burn_in=1000):
        self.mcmc = {}
        for p_type in agnostic_tqdm(self.polymer_types.values(), desc='Sampling Polymer Types:'):
            p_type.mcmc = {}
            p_type.bayes_sample_distance(show_progress=False)
            p_type.bayes_sample_energy(show_progress=False)

    @staticmethod
    def _cross_validate(input_data, classification):
        # Normalize the data
        X = StandardScaler().fit_transform(input_data)
        data = pd.DataFrame(index=input_data.index, columns=input_data.columns, data=X)

        logReg = linear_model.LogisticRegression(max_iter=500, tol=1e-6)
        loo_predictions = pd.Series(index=data.index, name='predictions')
        loo_score = pd.Series(index=data.index, dtype='bool', name='score')
        loo_probabilities = pd.Series(index=data.index, name='probabilities')

        for idx in data.index:
            # leave one out
            logReg.fit(data.drop([idx]), classification.drop([idx]))
            loo_predictions[idx] = logReg.predict([data.loc[idx]])
            loo_score[idx] = (loo_predictions[idx] == classification.loc[idx])
            # check at what index the true-class resides
            true_index = np.argwhere(logReg.classes_)[0][0]
            loo_probabilities[idx] = logReg.predict_proba([data.loc[idx]])[0][true_index]
        return loo_score, loo_predictions, loo_probabilities
    
    @staticmethod
    def _roc_auc(classification, probabilities):
        return metrics.roc_auc_score(classification, probabilities)

    @staticmethod
    def _kappa(true_classification, model_classification):
        return metrics.cohen_kappa_score(true_classification, model_classification)

    @staticmethod
    def _matthews(true_classification, model_classification):
        return metrics.matthews_corrcoef(true_classification, model_classification)


class PolymerTypeSims(object):

    def model_setup(self, distance):
        self.models = {}
        for sim in self.sims:
            sim.models = {}
        self.models['distance'] = self.distance_model_setup(distance)
        self.models['energy'] = self.energy_model_setup(distance)

    def create_bayes_energy_model(self):

        # Priors for polymer-Type
        pt_mu = mc.Normal('pt_mu', mu=0.0, tau=10.0)
        pt_sigma = mc.Uniform('pt_sigma', lower=0.0, upper=20.0)
        
        return {'pt_mu': pt_mu, 'pt_sigma': pt_sigma}

    def energy_model_setup(self, distance):
        model = {}
        model['distance'] = distance
        model = self.create_bayes_energy_model()
        model['p_theta'] = []
        model['p_sigma'] = []
        model['obs'] = []
        for sim in self.sims:
            sim.energy_model_setup(model, distance)
            model['p_theta'].append(sim.models['energy']['p_theta'])
            model['p_sigma'].append(sim.models['energy']['p_sigma'])
            model['obs'].append(sim.models['energy']['obs'])
        return model

    def distance_model_setup(self, distance):
        model = {}
        model['distance'] = distance
        model = self.create_bayes_distance_model()
        model['ptheta'] = []
        model['binding_obs'] = []
        for sim in self.sims:
            sim.distance_model_setup(model, distance)
            model['ptheta'].append(sim.models['distance']['p_theta'])
            model['binding_obs'].append(sim.models['distance']['binding_obs'])
        return model

    def create_bayes_distance_model(self):

        pt_mu = mc.Beta('pt_mu', alpha=2, beta=2)
        K = 6
        return {'mu': pt_mu, 'K': K}

    def bayes_results(self, model_name='distance', confidence_interval=0.95): 
        M = self.mcmc[model_name]
        confidence_interval = mc.utils.hpd(M.trace('pt_mu')[:], 1-confidence_interval)
        mean = M.trace('pt_mu')[:].mean()
        error = np.abs((confidence_interval- mean))
        return np.array([error[0], mean, error[1]])

    def bayes_sample_distance(self, sampling_rate=5000, burn_in=1000, show_progress=True):
        self.mcmc['distance'] = mc.MCMC(self.models['distance'])
        self.mcmc['distance'].sample(sampling_rate, burn_in, 10, True, save_interval=10, progress_bar=show_progress)

    def bayes_sample_energy(self, sampling_rate=5000, burn_in=1000, show_progress=True):
        self.mcmc['energy'] = mc.MCMC(self.models['energy'])
        self.mcmc['energy'].sample(sampling_rate, burn_in, 10, True, save_interval=10, progress_bar=show_progress)

class Job(object):

    def distance_model_setup(self, p_type_priors, max_distance):
        '''Create priors for distance model.
        '''
        distance_data = self.distances().apply(binding, dist=max_distance).values
        self.models['distance'] = self.create_bayes_distance_model(p_type_priors, distance_data)

    def energy_model_setup(self, p_type_priors, max_distance):
        energy_data = self.energies()[self.distances()<max_distance].dropna().values
        self.models['energy'] = self.create_bayes_energy_model(p_type_priors, energy_data)

    def create_bayes_energy_model(self, p_type_priors, data):
        pt_mu = p_type_priors['pt_mu']
        pt_sigma = p_type_priors['pt_sigma']

        p_theta = mc.Normal('p_theta_%i'%(self.Id-1), mu=pt_mu, tau=pt_sigma)
        p_sigma = mc.Lognormal('p_sigma_%i'%(self.Id-1), mu=3.0, tau=0.2)

        obs = mc.Normal('obs_%i'%(self.Id-1), mu=p_theta, tau=p_sigma, 
                            value=data, 
                            observed=True)
        return {'p_theta': p_theta, 'p_sigma': p_sigma, 'obs': obs}

    def create_bayes_distance_model(self, p_type_priors, data):
        # uninformed Prior
        pt_mu = p_type_priors['mu']
        K = p_type_priors['K']

        p_theta = mc.Beta('p_theta_%i'%(self.Id-1), alpha=pt_mu*K, beta=(1-pt_mu)*K)
        
        binding_obs = mc.Bernoulli('binding_obs_%i'%(self.Id-1), 
                                            p=p_theta, 
                                            value=data, 
                                            observed=True)

        return {'p_theta': p_theta, 'binding_obs': binding_obs}


    def bayes_results(self, model_name='distance', confidence_interval=0.95): 
        M = self.poly_type.mcmc[model_name]
        confidence_interval = mc.utils.hpd(M.trace('p_theta_%i'%(self.Id-1))[:], 1-confidence_interval)
        mean = M.trace('p_theta_%i'%(self.Id-1))[:].mean()
        error = np.abs((confidence_interval- mean))
        return np.array([error[0], mean, error[1]])
    

def binding(value, dist=15.0):
    if value <= dist:
        return 1
    else:
        return 0
