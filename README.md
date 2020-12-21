# degeneracy-decomposition
Decomposes heterogeneous reads caused by indel

This script was designed for a very specific problem. We were doing a CRISPR knockout of a gene which had an indel close to the gRNA cut site. This was the only appropriate gRNA cut site without off-target effects. We couldn't design primers to fit between the indel and gRNA site, and the reverse primers were too far to get a quality read from that direction.

A normal CRISPR tool (TIDE) was confused by the WT having degeneracy and threw a fit.

This was a lower-priority gene, and we didn't have a very large throughput. So these quick analyses were done by hand.
