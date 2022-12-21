# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'Cytopath'
copyright = '2022, Revant Gupta'
author = 'Revant Gupta'

release = '0.1.9'
version = '0.1.9'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
github_repo = "cytopath"
github_nb_repo = "cytopath-notebooks"

# -- Options for EPUB output
epub_show_urls = 'footnote'

# Width
html_static_path = ['_static']

def setup(app):
    app.add_css_file('custom.css')
