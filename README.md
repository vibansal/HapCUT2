
### HapCUT2 

HapCUT2 is a haplotype assembly algorithm that uses a likelihood model and is efficient for long read and Hi-C sequence technologies. The paper describing this method and a comparison to other tools for haplotype assembly has been published recently in Genome Research: http://genome.cshlp.org/content/early/2016/12/09/gr.213462.116.abstract 

A tarball of the source code (hapcut2.tar.gz) can be downloaded here. For the latest source code, check out: https://github.com/pjedge/hapcut2


### FragmentCut 

FragmentCut is a likelihood based method for the detection of long DNA fragments from dilution pool sequencing experiments (e.g. fosmid pooling, CGI LFR). FragmentCut identifies haplotype fragments that can subsequently be assembled into haplotype blocks using HapCUT2 or other tools. It is currently unpublished.


### CODE for generating haplotype fragment files from 10X barcoded BAM files 

See README file in the '10X' folder for instructions 


