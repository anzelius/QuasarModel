from distutils.core import setup
setup(
  name = 'QuasarModel',
  packages = ['QuasarModel'],
  version = '0.1',
  license='MIT',
  description = 'Program that models quasars from fits files as gaussian functions.',
  author = 'Tuss Anzelius & Ludvig Rodung',
  author_email = 'tuss.anzelius@gmail.com, ludvig@rodung.se',
  url = 'https://github.com/anzelius/QuasarModel',
  download_url = 'https://github.com/anzelius/QuasarModel/archive/refs/tags/v_01.tar.gz', 
  keywords = ['VLBI', 'QUASAR', 'LEASTSQUARE', 'MODEL'],
  install_requires=[           
          'astropy',
          'matplotlib',
          'numpy',
          'Pillow',
          'PySimpleGUI',
          'scikit-image',
          'keyboard',
          'scipy',
          'pandas'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3.9',
  ],
)
