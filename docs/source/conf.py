# Configuration file for the Sphinx documentation builder.
# See https://www.sphinx-doc.org/en/master/usage/configuration.html

project = "HeteroTyper Pipeline"
copyright = "2026, avicanLAB"
author = "avicanLAB"
release = "1.0"

# -- General configuration ---------------------------------------------

extensions = [
    "myst_parser",
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "substitution",
]

myst_heading_anchors = 3

source_suffix = {
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output ---------------------------------------------

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]

html_theme_options = {
    "github_url": "https://github.com/avicanlab/HeteroTyper",
    "show_toc_level": 2,
    "navigation_with_keys": True,
    "navbar_align": "left",
    "use_edit_page_button": True,
    "footer_start": ["copyright"],
    "footer_end": ["theme-version"],
    "footer_center": ["sphinx-version"]
}

html_context = {
    "github_user": "avicanlab",
    "github_repo": "HeteroTyper",
    "github_version": "main",
    "doc_path": "docs/source",
}

html_title = "HeteroTyper Pipeline"
