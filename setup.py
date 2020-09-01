import setuptools

long_description = """
Please find the package documentation in the original 
`Github repository<http://www.python.org/>`_.
"""
    
setuptools.setup(
    name='dupont-contraction',
    version='0.0.0',
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