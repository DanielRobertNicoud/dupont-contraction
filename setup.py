import setuptools

long_description = """
Please find the package documentation in the original Github repository_. We 
recommend the use of a LaTeX rendered for Github.

.. _repository: https://github.com/DanielRobertNicoud/dupont-contraction
"""
    
setuptools.setup(
    name='dupont-contraction',
    version='0.1.1',
    author='Daniel Robert-Nicoud',
    author_email='daniel.robertnicoud@gmail.com',
    description="A package for Sullivan and Dupont forms, and the Dupont " \
        "contraction",
    long_description=long_description,
    url='https://github.com/DanielRobertNicoud/dupont-contraction',
    install_requires=[
        'numpy'
    ],
    packages=setuptools.find_packages(),
    namespace_packages=['dupontcontraction']
)