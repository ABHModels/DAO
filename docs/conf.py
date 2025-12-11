# docs/conf.py

import os
import sys
sys.path.insert(0, os.path.abspath('..')) 

project = 'DAO'
copyright = '2025, Yimin Huang'
author = 'Yimin Huang'
version = '0.1'


extensions = [
    'sphinx.ext.autodoc',      
    'sphinx.ext.napoleon',      #
    'sphinx.ext.mathjax',       #
    'sphinx.ext.viewcode',      # 
    'myst_parser',   
    'sphinxcontrib.bibtex'         
]
bibtex_bibfiles = ['ref.bib']
bibtex_reference_style = 'author_year'
bibtex_default_style = 'unsrt'

myst_enable_extensions = [
    "dollarmath", 
    "amsmath",
]

import sphinx_rtd_theme
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_theme_options = {
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',

    
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}
