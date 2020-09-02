import setuptools

long_description = """
Please find the package documentation in the original Github repository_. We 
recommend the use of a LaTeX viewer for Github.

+---------+-------------------------------------------+
| Version | Comments                                  |
+=========+===========================================+
| 1.0.0   | First productive version.                 |
+---------+-------------------------------------------+
| 1.0.1   | No code changes, version to generate DOI. |
+---------+-------------------------------------------+
| 1.0.2   | Minor changes/bug fixes:                            |
|         |                                           |
|         | - Fixed sign error in a_infinity_product. |
|         | - Improved README.md                      |
+---------+-------------------------------------------+

.. _repository: https://github.com/DanielRobertNicoud/dupont-contraction
"""
    
setuptools.setup(
    name='dupont-contraction',
    version='1.0.2',
    author='Daniel Robert-Nicoud',
    author_email='daniel.robertnicoud@gmail.com',
    description="A package for computations using Sullivan and Dupont forms," \
        " and the Dupont contraction.",
    long_description=long_description,
    url='https://github.com/DanielRobertNicoud/dupont-contraction',
    install_requires=[
        'numpy'
    ],
    packages=setuptools.find_packages(),
    namespace_packages=['dupontcontraction']
)