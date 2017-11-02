import tempfile as temp
# import moviepy.editor as mpy
import time
import tqdm
import os
from xmlrpclib import ServerProxy
import pathlib2 as pl

def render_contact_picture(out_path, pdb_name, pymol_handle, gif=False):
    '''Creates a series of slightly rotated pictures with pymol and 
    creates a gif.

    input:
        out_path:   output path [string]
        
    '''
    pymol_handle.set("depth_cue", 0)
    pymol_handle.set("ray_opaque_background", "off")
    if gif:
        file_list = []
        folder = temp.mkdtemp()
        for i in tqdm.tqdm(range(120)):
            filename = '%s/pymolRotation%03d.png'% (folder, i)
            pymol_handle.do("cmd.png('%s', width=800, ray=1)" % filename )
            time.sleep(0.5)
            pymol_handle.do("cmd.rotate(\"y\", 3, \"%s\")" % pdb_name)
            time.sleep(0.5)
            file_list.append(filename)
        clip = mpy.ImageSequenceClip(file_list, fps=30)
        clip.write_gif(out_path)
        # cleaning up
        [os.remove(file_) for file_ in file_list]
        os.rmdir(folder)
    else:
        pymol_handle.do("cmd.png('%s', width=800, ray=1)" % out_path)


def connect(ip='132.252.170.144', port=9123):
    pymol = ServerProxy(uri='http://%s:%d/RPC2' % (ip, port) )
    return pymol

def load_pose(lig_path, receptor_path, pose_id, pymol_handle, interaction_res=None, clear_before=True):
    lig = pl.Path(lig_path)
    rec = pl.Path(receptor_path)
    pymol = pymol_handle
    if clear_before:
        pymol.delete('all')
    pymol.load(receptor_path)
    rec_name = rec.name.split('.')[0]
    lig_name = lig.name.split('.')[0]
    pymol.show_as('cartoon', rec_name)
    # ligand
    pymol.load(lig.as_posix())
    pymol.split_states(lig_name, pose_id, pose_id)
    pymol.delete(lig_name)
    pose_name = lig_name+ '_%04d' % pose_id
    pymol.show_as('sticks', pose_name)
    # interaction partner
    if interaction_res:
        pymol.show('sticks', '(resi %d)'%interaction_res)
        pymol.color('red', '(resi %d)'%interaction_res)
    pymol.zoom(rec_name)

def pymol_pic(path, pymol, size=(800,600), dpi=300):
    pymol.png(path, size[0], size[1], dpi, 1)
    return path
