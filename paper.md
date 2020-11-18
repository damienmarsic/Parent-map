---
title: 'Parent-map: analysis of parental contributions to evolved or engineered protein or DNA sequences'
tags:
  - Python
  - Directed evolution
  - AAV
  - Sequence analysis
authors:
  - name: Damien Marsic
    orcid: 0000-0003-0847-8095
    affiliation: 1
affiliations:
 - name: Porton Biologics, Suzhou Industrial Park, Jiangsu, China
   index: 1
date: 18 November 2020
bibliography: paper.bib
---

# Summary

Parent-map analyzes protein or DNA sequences which are derived from one or multiple parent sequences, and shows parental contributions as well as differences from relevant parents. Originally developed to analyze capsid protein sequences obtained by directed evolution, parent-map can be used in any case where variant sequences are to be compared to parent sequences from which they are derived. Parent-map detects sequence shuffling as well as substitutions, insertions and deletions, and displays results in user-friendly formats. Parent-map is an open-source, platform-independent Python 3 script, available as a Bioconda package as well as a Windows program.
Source code: <https://github.com/damienmarsic/Parent-map>
Python package: <https://pypi.org/project/parent-map/>
Bioconda recipe and package: <http://bioconda.github.io/recipes/parent-map/README.html>
Windows installer: <https://sourceforge.net/projects/parent-map/>
Documentation: <https://parent-map.readthedocs.io/>

# Statement of need

Adeno-associated virus (AAV) capsid directed evolution projects typically generate multiple enriched variant sequences after 2 to 5 rounds of selection starting from complex capsid libraries. For libraries developed from a single parental serotype, through random peptide insertion at a specific position or surface loop diversification in well-defined variable regions for example, a single multiple alignment of all enriched variant sequences against the parent sequence conveniently shows how each variant differs from the parent. However, when more than one parental sequence is involved, such as when different libraries are mixed together, or when a library design involves DNA shuffling from several parents, such alignments can quickly become illegible, particularly when the complete capsid gene is sequenced. In such cases, in the absence of appropriate software tools, each variant needs to be separately aligned against all possible parents, a time-consuming and cumbersome process. An added difficulty in the case of shuffled libraries is that, because of high sequence homology between parents, multiple regions will share sequence identities with more than one parent, complicating attempts at comprehensively defining the variant sequences in terms of parental contributions. To date, SALANTO [@herrmann_robust_2019] seems to be the only relevant publicly available software. However, it only applies to shuffled libraries, and its user-friendliness is limited as it requires the user to perform a multiple sequence alignment beforehand, and to further process the data manually after analysis. The software described in this article, parent-map, provides a user-friendly and comprehensive solution. It can be used with sequences derived from any type of library, or even with naturally-occurring mutants or rationally engineered variants. It is not limited to protein sequences. It only requires one file containing the variant sequences to be analyzed, and one file containing parental sequences, without any prior manipulation. It generates a set of five files covering most end-users’ needs, in directly usable formats. Finally, although it was developed to address a need in the field of AAV capsid directed evolution, parent-map can be used whenever protein or DNA sequences, whether originating from natural evolution, directed evolution or rational design, are to be compared with one or more possible parental sequences.

Parent-map was written under Python 3.7 as both a command-line interface (CLI) and a graphical user interface (GUI) application, by allowing parser modules [argparse](https://docs.python.org/3/library/argparse.html) and [Gooey](https://github.com/chriskiehl/Gooey) to coexist within a single file (the GUI will start if no argument is present, while any argument will cause parent-map to start in CLI mode). A parent-map Python package was created and uploaded to the Python Package Index (PyPI) according to [packaging instructions](https://packaging.python.org/tutorials/packaging-projects/). A parent-map Bioconda [@gruning_bioconda_2018] recipe based on the PyPI package was written and submitted according to [instructions](https://bioconda.github.io/contributor/index.html). A stand-alone Windows executable and its installation program were created using respectively [PyInstaller](https://www.pyinstaller.org/) and [Inno Setup](https://jrsoftware.org/isinfo.php). The documentation was written using [Sphinx](https://www.sphinx-doc.org/en/master/).

Parent-map is a multi-platform Python script that generates a set of five output files from two input files. Input file names and options can be entered as arguments at launch time, resulting in parent-map running in CLI mode, or within the GUI, which starts if parent-map is launched without arguments. This flexibility allows parent-map to be deployed in a variety of settings, as a simple desktop application or even as a bioinformatics pipeline component. The first input file contains the variant sequences, typically the most frequent or the most enriched sequences obtained at the completion of a directed evolution experiment. The other input file is a set of potential parental sequences to the variant sequences. The most useful files generated by parent-map, particularly in the case of variants derived from DNA shuffling, are parental contribution maps (file names ending in –par.txt and –par.html, the latter being a colorized version of the former). Instead of all possible combinations, the simplest map that can accurately describe the variant is shown, using as few parents and as few fragments as possible. Other output files include a statistics file summarizing the variant sequences main features, a sequence definition file comprehensively defining each variant in terms of its parents, and an alignment file showing how variants differ from their common parent.

Parent-map can be tested using the provided [variant](https://github.com/damienmarsic/Parent-map/blob/master/example_variants.fasta) and [parent](https://github.com/damienmarsic/Parent-map/blob/master/example_parents.fasta) sample files, based on available literature describing evolved and rationally designed AAV capsid variants. Variants AAV-DJ [@grimm_vitro_2008], AAV2.5T [@excoffon_directed_2009], NP84 [@paulk_bioengineered_2018] and OLIG001 [@powell_characterization_2016] are derived from shuffled DNA libraries. Variants AAV-F [@hanlon_selection_2019], AAV-PHP.B [@deverman_cre-dependent_2016], 7m8 [@dalkara_vivo-directed_2013] and rAAV2-retro [@tervo_designer_2016] are derived from peptide insertion libraries. Variants SCH2, SCH9 [@ojala_vivo_2018], LI-A and LI-C [@marsic_vector_2014] are derived from more complex rationally designed libraries. Variants AAV2i8 [@asokan_reengineering_2010] and AAV2-sept-Y-F [@petrs-silva_novel_2011] were rationally designed.

A comprehensive description of parent-map is provided in the [documentation](https://parent-map.rtfd.io).

# Acknowledgements

We thank Yan Chen and Oleksandr Kondratov for testing parent-map and providing valuable feedback.

# References
