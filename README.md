# ISO-scSLAM

Snakemake workflow for event-level alternative-splicing analysis from completed
long-read SLAM results. The public tree contains reusable rules, report
builders, environment definitions, and tests.

External ISO-SLAM inputs and completed result trees are not included. Create a
local configuration with external read-only source results and repository-local
output, log, and temporary roots.

~~~bash
make lint
make smoke
~~~
