from setuptools import setup, find_packages

setup(
    name='astro-utils',
    use_scm_version=True,  # Enables versioning based on Git
    setup_requires=['setuptools_scm'],  # Ensures setuptools_scm is available for building
    packages=find_packages(),
    install_requires=[
        'ipywidgets',
        'IPython',
        'numpy',
        'matplotlib',
        'astropy'
    ],
)
