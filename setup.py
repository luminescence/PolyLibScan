from setuptools import setup

setup(name='PolyLibScan',
      version='0.1',
      description='A toolkit to model polymer-protein interaction.',
      url='https://github.com/BioinformaticsBiophysicsUDE/PolyLibScan',
      author='Ludwig Ohl, Niklas Toetsch',
      author_email='Ludwig.Ohl@uni-due.de, Niklas.Toetsch@uni-due.de',
      license='MIT',
      packages=['PolyLibScan'],
      install_requires=[
      	'numba==0.30.1+0.g8c1033f.dirty',
		'mock==2.0.0',
		'pandas==0.19.2',
		'numpy==1.14.2',
		'pathlib2==2.3.2',
		'Jinja2==2.10',
		'xarray==0.10.2',
		'matplotlib==2.0.0',
		'tqdm==4.11.2',
		'scipy==0.18.1',
		'tables==3.3.0',
		'MDAnalysis==0.18.0',
		'adjustText==0.7.3',
		'futures==3.2.0',
		'pymc==2.3.6',
		'scikit_learn==0.19.2',
		'PyYAML==3.13']
      zip_safe=False)