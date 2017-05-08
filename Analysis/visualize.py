import matplotlib.pyplot as plt
import numpy as np
import PolyLibScan.Tools.lmp_helpers
from lightning import Lightning
from IPython.display import display, HTML

class Visualize(object):

    def __init__(self, parent):
        self.parent = parent
        self.lgn = Lightning(ipython=True, local=True)

    def trajectory(self, run):
        scene = self._generate_scene(run)
        display(HTML(scene))

    def _generate_scene(self, run):
        points = self._get_components(run)
        return self._setup_scene(points, run)

    def _setup_scene(self, points, run):
        # devides the points into groups of polymer
        # protein and active site particles
        group = self._identify_groups(points, run)
        # ball size
        size = size = np.zeros(points.shape, np.int)
        size[:] = 4
        # opacity
        alpha = np.array([0.4 if g==2 else 1.0 for g in group])
        
        graph = self.lgn.scatter3(points['x'], points['y'], points['z'], 
                             group=group, size=size, alpha=alpha).get_html()
        return graph

    def _get_components(self, run):
        '''returns the coordinates of protein, polymer, interaction particles
        and active site in the start and end pose.
        '''
        points_end = run.coordinates()['end']
        points_start = run.coordinates()['start']
        poly_mask = np.in1d(points_start['atom_type'], run.job._particle_ids['polymer'])
        polymer_points_start = points_start[poly_mask]
        complete_points = np.append(points_end, polymer_points_start)
        return complete_points
    
    def _identify_groups(self, points, run):
        group = np.zeros(points.shape, np.int)
        for i,at in enumerate(points['atom_type']):
            if at in run.job._particle_ids['polymer']:
                group[i] = 1
            else:
                group[i] = 2
        group[run.job.active_site['xyz']+1] = 3
        return group
