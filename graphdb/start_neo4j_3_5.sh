#!/bin/bash
# Current issues:  port map not working, ulimit warning persisting
cd /srv/rbd/covid19/bisystem
singularity instance start --fakeroot --writable-tmpfs --net --bind neo4jdata_3_5:/var/lib/neo4j/data,neo4jrun_3_5:/var/lib/neo4j/run,neo4jlogs_3_5:/var/lib/neo4j/logs,neo4jplugins_3_5:/var/lib/neo4j/plugins,neo4jcerts_3_5:/var/lib/neo4j/certificates  neo4j_3_5.sif ne4j2
singularity exec --fakeroot instance://ne4j2 neo4j start
