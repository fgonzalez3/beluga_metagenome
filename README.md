# Beluga Metagenome
This repo contains documentation for Beluga whale metagenome assembly. 

Initial metagenomics QC and assembly can be found [here](https://github.com/fgonzalez3/beluga_metagenome/blob/main/workflow.md). Genome annotation can be found [here](https://github.com/fgonzalez3/beluga_metagenome/blob/main/Genome_annotation.md). Downstream R analysis can be found [here](https://github.com/fgonzalez3/beluga_metagenome/blob/main/abundance.md). 

There was one paper ([Bai et al 2021](https://www.frontiersin.org/articles/10.3389/fmicb.2021.769012/full)) that looked at Beluga feces using shotgun metagenomics. I analyzed that data in a similar fashion to our data above, found [here](https://github.com/fgonzalez3/beluga_metagenome/blob/main/china_metagenome.md). 

The initial SeqCenter data showed hits to Morbillivirus, which is something we're interested in for this project. However, this was shotgun DNA data so that shouldn't have happened. I wrote some quick code to sort this out [here](https://github.com/fgonzalez3/beluga_metagenome/blob/main/Morbillivirus.md). 

This pipeline broadly followed that laid out [here](https://carpentries-lab.github.io/metagenomics-analysis/). 
